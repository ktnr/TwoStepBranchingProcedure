#include "BranchAndBound.h"

#include <algorithm>
#include <cassert>
#include <execution>
#include <iostream>
#include <mutex>
#include <numeric>
#include <string>

#include "ltalloc/ltalloc.h"
#include "taskflow/taskflow/taskflow.hpp"

namespace tsbp
{

// Use std::min<T> and std::max<T> explicitly to avoid definition clash: https://stackoverflow.com/a/30924806/5587903.
////#include <psapi.h>
////#include <windows.h>

void Packing2D::AddItem(Rectangle&& rectangle, int itemId)
{
    Items.emplace_back(std::move(rectangle));

    if (itemId == -1)
    {
        return;
    }

    PlacedItems.set((size_t)itemId);
}

bool LeftmostActiveOnly::Node::AllItemsPlaced(const std::vector<Rectangle>& items) const
{
    return items.size() == PlacedItems().count();
}

#pragma region Branch and bound logic

void LeftmostActiveOnly::Preprocess()
{
    size_t domainReducedItemIndex = FindMaximumItemDx();

    this->domainReducedItemIndex = domainReducedItemIndex;

    DomainReduction(domainReducedItemIndex);
}

size_t LeftmostActiveOnly::FindMinimumItemDx()
{
    int minDx = this->items[0].Dx;

    size_t domainReducedItemIndex = 0;
    for (size_t i = 1; i < this->items.size(); i++)
    {
        const Rectangle& item = this->items[i];

        if (item.Dx < minDx)
        {
            minDx = item.Dx;
            domainReducedItemIndex = i;
        }
    }

    return domainReducedItemIndex;
}

size_t LeftmostActiveOnly::FindMaximumItemDx()
{
    int maxDx = this->items[0].Dx;

    size_t domainReducedItemIndex = 0;
    for (size_t i = 1; i < this->items.size(); i++)
    {
        const Rectangle& item = this->items[i];

        if (item.Dx > maxDx)
        {
            maxDx = item.Dx;
            domainReducedItemIndex = i;
        }
    }

    return domainReducedItemIndex;
}

void LeftmostActiveOnly::DomainReduction(size_t domainReducedItemIndex)
{
    if (!this->fixedItemCoordinatesX.empty())
    {
        DomainReductionY(domainReducedItemIndex);
    }
    else
    {
        DomainReductionXY(domainReducedItemIndex);
    }
}

void LeftmostActiveOnly::DomainReductionY(const size_t domainReducedItemIndex)
{
    itemSpecificPlacementPointsY.reserve(this->items.size());

    int reducedDomainThresholdX = std::floor((this->container.Dx - this->items[domainReducedItemIndex].Dx) / 2.0);
    int reducedDomainThresholdY = std::floor((this->container.Dy - this->items[domainReducedItemIndex].Dy) / 2.0);

    for (size_t i = 0; i < this->items.size(); i++)
    {
        const Rectangle& item = this->items[i];
        assert(item.InternId == i);

        if (i == domainReducedItemIndex)
        {
            auto& bitsetY = this->itemSpecificPlacementPointsY.emplace_back(boost::dynamic_bitset<>(this->container.Dy));
            bitsetY.reset();

            for (size_t j = 0; j <= reducedDomainThresholdY; j++)
            {
                bitsetY.set(j);
            }

            this->reducedItemDomainX = reducedDomainThresholdX;
            this->reducedItemDomainY = reducedDomainThresholdY;
            this->reducedItemFeasiblePlacementPoints = (reducedDomainThresholdX + 1) * (reducedDomainThresholdY + 1);

            continue;
        }

        // TODO: reduce by item.Dy to exclude overlap with the container. Enables simplification of no-overlap checks.
        auto& bitsetY = this->itemSpecificPlacementPointsY.emplace_back(boost::dynamic_bitset<>(this->container.Dy));
        bitsetY.set();
    }
}

void LeftmostActiveOnly::DomainReductionXY(const size_t domainReducedItemIndex)
{
    itemSpecificPlacementPointsX.reserve(this->items.size());
    itemSpecificPlacementPointsY.reserve(this->items.size());

    int reducedDomainThresholdX = std::floor((this->container.Dx - this->items[domainReducedItemIndex].Dx) / 2.0);
    int reducedDomainThresholdY = std::floor((this->container.Dy - this->items[domainReducedItemIndex].Dy) / 2.0);

    // TODO.Performance: implement further domain reduction according to Soh (2008). In particular, deactivate placement points
    // when the $lr$ clause can be set to false (the second part of domain reduction).
    for (size_t i = 0; i < this->items.size(); i++)
    {
        const Rectangle& item = this->items[i];
        assert(item.InternId == i); // Caution: assert does not throw with NDEBUG flag.

        if (i == domainReducedItemIndex)
        {
            auto& bitsetX = this->itemSpecificPlacementPointsX.emplace_back(boost::dynamic_bitset<>(this->container.Dx));
            auto& bitsetY = this->itemSpecificPlacementPointsY.emplace_back(boost::dynamic_bitset<>(this->container.Dy));
            bitsetX.reset();
            bitsetY.reset();

            for (size_t j = 0; j <= reducedDomainThresholdX; j++)
            {
                bitsetX.set(j);
            }

            for (size_t j = 0; j <= reducedDomainThresholdY; j++)
            {
                bitsetY.set(j);
            }

            this->reducedItemDomainX = reducedDomainThresholdX;
            this->reducedItemDomainY = reducedDomainThresholdY;
            this->reducedItemFeasiblePlacementPoints = (reducedDomainThresholdX + 1) * (reducedDomainThresholdY + 1);

            continue;
        }

        // TODO: reduce by item.Dx/Dy to exclude overlap with the container. Enables simplification of no-overlap checks.
        auto& bitsetX = this->itemSpecificPlacementPointsX.emplace_back(boost::dynamic_bitset<>(this->container.Dx));
        auto& bitsetY = this->itemSpecificPlacementPointsY.emplace_back(boost::dynamic_bitset<>(this->container.Dy));
        bitsetX.set();
        bitsetY.set();
    }
}

SearchStatus LeftmostActiveOnly::Solve()
{
    auto start = std::chrono::steady_clock::now();

    int numberOfThreads = this->parameters.LMAOThreads;
    SearchStatus searchStatus = SearchStatus::None;

    std::vector<int> itemOrder(items.size());
    std::iota(itemOrder.begin(), itemOrder.end(), 0);
    std::sort(itemOrder.begin(), itemOrder.end(), CompareLexicographicDxDy(&items));

    this->branchingOrder = itemOrder;

    if (numberOfThreads == 1)
    {
        searchStatus = SolveSequential();
        this->statistics.NodeCount = this->tree.m_vertices.size();
    }
    else if (numberOfThreads >= 1)
    {
        searchStatus = SolveParallelTaskflow();
        // Slightly inferior performance compared to taskflow with this memory allocator. Cannot limit threads.
        ////searchStatus = SolveParallelNative();
    }
    else
    {
        std::cout << "Invalid number of threads parameter, solve in sequential mode.\n";
        searchStatus = SolveSequential();
    }

    auto end = std::chrono::steady_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    this->statistics.TimeOverall = time;
    ////this->statistics.MemoryUsage = MemoryUsedByCurrentProcess();

    return searchStatus;
}

SearchStatus LeftmostActiveOnly::SolveParallelNative()
{
    std::atomic<std::chrono::steady_clock::time_point> begin = std::chrono::steady_clock::now();
    std::atomic<std::chrono::steady_clock::time_point> endLast = begin.load();

    std::atomic<bool> feasibleSolutionCancellationToken = false;

    std::vector<std::unique_ptr<LeftmostActiveOnly>> subtrees;
    subtrees.reserve(items.size());
    ////for (size_t i = 0; i < 1; i++)
    for (size_t itemId: this->branchingOrder)
    {
        std::unique_ptr<LeftmostActiveOnly> lmao = std::make_unique<LeftmostActiveOnly>(this->parameters);
        lmao->AddContainer(container);
        lmao->AddItems(items, itemId, this->fixedItemCoordinatesX);
        lmao->SetBranchingOrder(this->branchingOrder);
        lmao->AddCancellationToken(&feasibleSolutionCancellationToken); // Optional: comment out for sequential mode.
        ////lmao->AddFixedItem(fixedItem);

        subtrees.push_back(std::move(lmao));
    }

    SearchStatus status = SearchStatus::None;

    std::atomic<size_t> exploredSubtrees(0);
    std::atomic<int64_t> exploredNodes(0);
    std::atomic<bool> isFeasible(false);
    Packing2D solution;

    std::mutex mutex;

    std::for_each(
        std::execution::par_unseq,
        subtrees.begin(),
        subtrees.end(),
        [&](std::unique_ptr<LeftmostActiveOnly>& lmao)
        {
            SearchStatus localStatus = lmao->SolveSequential();
            if (localStatus == SearchStatus::Abort)
            {
                status = localStatus;
            }

            if (localStatus == SearchStatus::Feasible)
            {
                mutex.lock();
                solution = lmao->GetSolution();
                mutex.unlock();

                isFeasible.store(true);
            }

            exploredNodes += lmao->GetTreeSize();
            exploredSubtrees++;

            std::cout << "Explored subtrees = " + std::to_string(exploredSubtrees) + " / " + std::to_string(items.size()) + "\n";

            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            endLast = end;
            std::cout << "Elapsed time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin.load()).count() << "s\n";

            lmao.reset();
        });

    if (status == SearchStatus::Abort)
    {
        std::cout << "Search was aborted.\n";
        return status;
    }

    if (isFeasible.load())
    {
        this->solution = solution;
        status = SearchStatus::Feasible;
        std::cout << "Feasible solution found.\n";
    }
    else
    {
        status = SearchStatus::Infeasible;
        std::cout << "Problem is infeasible.\n";
    }

    std::chrono::steady_clock::time_point endParallelFor = std::chrono::steady_clock::now();

    std::cout << "Elapsed time until after for_each() = " << std::chrono::duration_cast<std::chrono::seconds>(endParallelFor - begin.load()).count() << "s\n";
    std::cout << "Termination delay = " << std::chrono::duration_cast<std::chrono::seconds>(endParallelFor - endLast.load()).count() << "s\n";

    std::cout << "Total explored nodes = " + std::to_string(exploredNodes.load()) + "\n";

    return status;
}

SearchStatus LeftmostActiveOnly::SolveParallelTaskflow()
{
    std::atomic<std::chrono::steady_clock::time_point> begin = std::chrono::steady_clock::now();
    std::atomic<std::chrono::steady_clock::time_point> endLast = begin.load();

    size_t numberOfThreads = this->parameters.LMAOThreads;
    tf::Executor executor(numberOfThreads);
    tf::Taskflow taskflow;

    SearchStatus status = SearchStatus::None;

    std::atomic<size_t> exploredSubtrees(0);
    std::atomic<size_t> exploredNodes(0);
    std::atomic<bool> feasibleSolutionCancellationToken = false;
    Packing2D solution;

    std::mutex mutex;

    for (size_t itemId: this->branchingOrder)
    {
        tf::Task itemSpecificBranch = taskflow.emplace(
            [&, itemId]()
            {
                    // TODO.Performance: write copy assignment constructor for Node to run the preprocess only once. It has its copy assignment constructor disabled because it has a unique_ptr member.
                    ////LeftmostActiveOnly lmaoSubtree = *this;
                    LeftmostActiveOnly lmao(this->parameters);
                    lmao.AddContainer(container);
                    lmao.AddItems(items, itemId, this->fixedItemCoordinatesX);
                    lmao.SetBranchingOrder(this->branchingOrder);
                    lmao.AddCancellationToken(&feasibleSolutionCancellationToken);

                    SearchStatus localStatus = lmao.SolveSequential();
                    if (localStatus == SearchStatus::Abort)
                    {
                        status = localStatus;
                    }

                    if (localStatus == SearchStatus::Feasible)
                    {
                        mutex.lock();
                        solution = lmao.GetSolution();
                        mutex.unlock();
                    }

                    exploredNodes += lmao.GetTreeSize();
                    exploredSubtrees++;

                    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
                    endLast = end;

                    ////if (this->parameters.LMAOEnableLogging)
                    if (false)
                    {
                        std::cout << "LMAO \t Explored subtrees = " + std::to_string(exploredSubtrees) + " / " + std::to_string(items.size()) + "\n";

                        std::cout << "Elapsed time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin.load()).count() << "s\n";
                    } });
    }

    // TODO.Performance: Consider prioritizing the reduced domain item and smaller items. It takes the longest to compute and should not come late in the process, s.t. the other threads cannot be used as no computations are left.
    // TODO.Performance: Execution with the taskflow library takes longer in sequential and parallel mode.
    // In sequential mode, the issue seems to be that a task/thread is blocked unitl after it has finished computation *and* until its memory fully deallocated.
    // This also has an impact in parallel mode.
    // Look into async tasks https://taskflow.github.io/taskflow/AsyncTasking.html and parallel algorithms.

    tf::Future<void> fu = executor.run(taskflow);
    fu.wait(); // block until the execution completes

    std::chrono::steady_clock::time_point endTaskflow = std::chrono::steady_clock::now();

    if (status == SearchStatus::Abort)
    {
        std::cout << "Search was aborted.\n";
        return status;
    }

    if (feasibleSolutionCancellationToken.load())
    {
        status = SearchStatus::Feasible;
        this->solution = solution;
        std::cout << "Feasible solution found.\n";
    }
    else
    {
        status = SearchStatus::Infeasible;
        std::cout << "Problem is infeasible.\n";
    }

    ////if (this->parameters.LMAOEnableLogging)
    if (false)
    {
        std::cout << "LMAO \t Total explored nodes = " + std::to_string(exploredNodes.load()) + "\n";
    }

    this->statistics.NodeCount = exploredNodes.load();
    this->statistics.TimeMemoryDeallocation = std::chrono::duration_cast<std::chrono::milliseconds>(endTaskflow - endLast.load()).count();

    return status;
}

void LeftmostActiveOnly::UpdateSolverStatistics(SolverStatistics& const statistics) const
{
    statistics.CallCountLMAO += 1;
    statistics.NodeCountLMAO += this->statistics.NodeCount;
    statistics.TimeLMAO += this->statistics.TimeOverall;
    statistics.TimeMemoryDeallocation += this->statistics.TimeMemoryDeallocation;
    statistics.MemoryUsage = std::max<double>(statistics.MemoryUsage, this->statistics.MemoryUsage);
}

Packing2D LeftmostActiveOnly::GetSolution() const
{
    return this->solution;
}

size_t LeftmostActiveOnly::InitializeSearchTree()
{
    size_t rootNodeId = boost::add_vertex(this->tree);
    this->activeNodes.push_back(rootNodeId);

    Node& rootNode = this->tree[rootNodeId];
    rootNode.ItemIdToPlace = -1;
    rootNode.NodeStatus = BaseNode::Status::Root;

    std::unique_ptr<Packing2D> layout2D = std::make_unique<Packing2D>();
    layout2D->Items.reserve(items.size());
    layout2D->PlacedItems = boost::dynamic_bitset<>(items.size());
    layout2D->RemainingItemAreaToPlace = std::accumulate(items.begin(), items.end(), 0.0, [](double volume, const Rectangle& i)
                                                         { return volume + i.Area; });
    layout2D->PlacedAreaVector = std::vector<size_t>(container.Dy + 1, 0);
    layout2D->DeactivatedAreaVector = std::vector<size_t>(container.Dy + 1, 0);
    layout2D->ReducedItemInfeasiblePlacementPoints = 0;

    rootNode.Packing = std::move(layout2D);

    if (this->fixedItemAtOrigin.has_value())
    {
        // Restrict placement of fixedItem
        this->itemSpecificPlacementPointsX[this->fixedItemAtOrigin.value()].reset();
        this->itemSpecificPlacementPointsY[this->fixedItemAtOrigin.value()].reset();
        this->itemSpecificPlacementPointsX[this->fixedItemAtOrigin.value()].set(0);
        this->itemSpecificPlacementPointsY[this->fixedItemAtOrigin.value()].set(0);

        ////this->tree.m_vertices.reserve(250000000);
    }

    return rootNodeId;
}

bool LeftmostActiveOnly::IsCancelled() const
{
    return this->cancellationToken != nullptr && this->cancellationToken->load();
}

void LeftmostActiveOnly::SignalCancelling()
{
    if (this->cancellationToken != nullptr)
    {
        this->cancellationToken->store(true);
    }
}

void LeftmostActiveOnly::AddItems(const std::vector<Rectangle>& items, const std::optional<int>& fixedItemAtOrigin, const std::vector<int>& itemCoordinatesX)
{
    if (this->container.Dx == 0)
    {
        throw std::exception("Container not assigned.");
    }

    this->tree.m_vertices.reserve(items.size() + 1);
    this->items = items;

    if (fixedItemAtOrigin.has_value())
    {
        this->fixedItemAtOrigin = fixedItemAtOrigin;
    }

    // The intern IDs of the items must correpond to the indices of items.
    for (size_t i = 0; i < itemCoordinatesX.size(); i++)
    {
        const Rectangle& item = items[i];

        if (i != item.InternId)
        {
            throw std::exception("Internal item IDs do not match.");
        }

        auto& bitsetX = this->itemSpecificPlacementPointsX.emplace_back(boost::dynamic_bitset<>(this->container.Dx));
        bitsetX.reset();

        int coordinate = itemCoordinatesX[i];
        if (coordinate < 0)
        {
            throw std::exception("Negative fixed coordinate.");
        }

        bitsetX.set(coordinate);
    }
}

void LeftmostActiveOnly::AddContainer(const Bin& container)
{
    this->container = container;
}

void LeftmostActiveOnly::AddCancellationToken(std::atomic<bool>* cancellationToken)
{
    this->cancellationToken = cancellationToken;
}

SearchStatus LeftmostActiveOnly::SolveSequential()
{
    Preprocess();

    size_t nodeId = InitializeSearchTree();

    size_t previousNodeId = nodeId;

    // TODO.Performance: parallelize tree search, e.g. a separate sub-tree for each node starting at depth 1
    this->searchStatus = SearchStatus::InProgress;
    while (this->searchStatus == SearchStatus::InProgress && !IsCancelled())
    {
        previousNodeId = nodeId;
        std::optional<size_t> newNodeId = Branch(nodeId);

        if (!newNodeId.has_value())
        {
            // No new node could be created: an infeasible sequence has been found.
            std::optional<size_t> backtrackedNodeId = Backtrack(nodeId);

            if (!backtrackedNodeId.has_value())
            {
                this->searchStatus = SearchStatus::Infeasible;
                break;
            }

            nodeId = backtrackedNodeId.value();
            continue;
        }

        Node& newNode = this->tree[newNodeId.value()];

        switch (newNode.NodeType)
        {
            case BaseNode::Type::UsePlacement:
                EvaluateLeaf(newNode, newNodeId.value());
                break;
            case BaseNode::Type::DeactivatePlacement:
                DeactivatePlacement(newNode, newNodeId.value());
                break;
            default:
                throw std::exception("Invalid node type.");
        }

        if (newNode.AllItemsPlaced(this->items))
        {
            this->searchStatus = SearchStatus::Feasible;
            this->feasibleLeafNodeId = newNodeId.value(); // newNode.Id;
            break;
        }

        // SelectNextNode(nodeStatus, previousNodeId) return nodeId to branch on
        switch (newNode.NodeStatus)
        {
            case BaseNode::Status::InfeasiblePlacement:
                // Prune node. Backtrack a single edge and select next most promising node to branch on.
                nodeId = previousNodeId;
                break;
            case BaseNode::Status::FeasiblePlacement:
                // Fallthrough.
            case BaseNode::Status::DeactivatedPlacement:
                // Fallthrough.
            case BaseNode::Status::InfeasibleSequence:
                // Continue search at current node (DFS).
                nodeId = newNodeId.value();
                break;
            default:
                throw std::exception("Invalid node status.");
        }
    }

    size_t bestNodeId = this->maxPlacedItemsInfeasibleNodeId;
    if (IsCancelled())
    {
        this->searchStatus = SearchStatus::Cancelled;
        bestNodeId = 0;
    }

    if (this->searchStatus == SearchStatus::Feasible)
    {
        bestNodeId = this->feasibleLeafNodeId;

        SignalCancelling();
    }

    std::cout << "Search status: " << (int)this->searchStatus << "\n";

    if (bestNodeId == -1)
    {
        throw std::exception("Search status is feasible but feasible leaf node is -1.");
    }

    if (fixedItemAtOrigin.has_value())
    {
        std::cout << "Fixed item (intern id): " + std::to_string(this->fixedItemAtOrigin.value()) + ", explored nodes = " + std::to_string(GetTreeSize()) + "\n";
    }
    else
    {
        std::cout << "Explored nodes = " + std::to_string(GetTreeSize()) + "\n";
    }

    switch (this->searchStatus)
    {
        case SearchStatus::Feasible:
        case SearchStatus::Infeasible:
            break;
        case SearchStatus::Cancelled:
            break;
        case SearchStatus::Abort:
        default:
            throw std::exception("Search status after search is neither feasible nor infeasible.");
    }

    if (this->searchStatus == SearchStatus::Feasible)
    {
        const auto& feasibleNode = this->tree[feasibleLeafNodeId];
        this->solution = *feasibleNode.Packing;
    }

    return this->searchStatus;
}

LeftmostActiveOnly::SearchTree::vertex_descriptor LeftmostActiveOnly::AddLeafNode(Node&& node)
{
    SearchTree::vertex_descriptor newNodeId = boost::add_vertex(this->tree);
    ////assert(node.Id == newNodeId);

    Node& newNode = this->tree[newNodeId];
    size_t remainingItemsToPlace = this->items.size() - node.PlacedItems().count();

    // TODO.Performance: this is helpful with directionalS and bidirectionalS. With bidirectionalS, edges are additionally (?) stored as members of tree.m_edges.
    this->tree.m_vertices[newNodeId].m_out_edges.reserve(remainingItemsToPlace + 1);

    newNode = std::move(node);

    return newNodeId;
}

std::optional<size_t> LeftmostActiveOnly::Branch(size_t nodeId)
{
    Node& const currentNode = this->tree[nodeId];
    const boost::dynamic_bitset<>& placedItems = currentNode.PlacedItems();

    if (currentNode.NodeStatus == BaseNode::Status::InfeasibleSequence)
    {
        return std::nullopt;
    }

    // Determine number of outgoing edges.
    size_t numberOfOutgoingEdges = boost::out_degree(nodeId, this->tree);
    size_t numberOfRemainingItems = items.size() - placedItems.count();

    ////bool isFixedSequenceDeactivation = numberOfOutgoingEdges == 1 && IsFixedSequence(currentNode, -1);
    ////if (numberOfOutgoingEdges == numberOfRemainingItems || isFixedSequenceDeactivation)
    if (numberOfOutgoingEdges == numberOfRemainingItems)
    {
        // TODO.Logic: currently the DeactivatePlacement node is the last node that is generated. Rework to allow it to be explored in any sequence.

        this->nodesToDeactivate.emplace(nodeId);

        const Packing2D& packing = *currentNode.Packing;

        int placementX = packing.ActiveX;
        int placementY = packing.ActiveY;

        if (currentNode.NodeStatus == Node::Status::Root)
        {
            return std::nullopt;
        }

        if (placementY == 0 && placementX == packing.MaxX)
        {
            currentNode.NodeStatus == BaseNode::Status::InfeasibleSequence;
            return std::nullopt;
        }

        /*
        typename SearchTree::in_edge_iterator inEdgeIterator;
        typename SearchTree::in_edge_iterator inEdgeEndIterator;
        typename SearchTree::edge_descriptor inEdge;

        for (boost::tie(inEdgeIterator, inEdgeEndIterator) = boost::in_edges(this->tree, nodeId); inEdgeIterator != inEdgeEndIterator; ++inEdgeIterator)
        {
            inEdge = *inEdgeIterator;
            typename SearchTree::vertex_descriptor targetId = boost::source(inEdge, this->tree);
            const Node& targetNode = this->tree[targetId];
        }
        */

        typename SearchTree::vertex_descriptor leftNodeId = boost::num_vertices(this->tree);
        AddLeafNode(Node{-1, BaseNode::Type::DeactivatePlacement, BaseNode::Status::NotEvaluated, std::make_unique<Packing2D>(*currentNode.Packing)});
        boost::add_edge(nodeId, leftNodeId, this->tree);

        // TODO.Safety: Caution, possible dangling reference to currentNode because the underlying graph might resize on AddLeafNode().
        // TODO.Performance: Consider using std::unique_ptr<Node> to reduce copy effort.

        this->activeNodes.push_back(leftNodeId);

        return leftNodeId;
    }

    typename SearchTree::out_edge_iterator outEdgeIterator;
    typename SearchTree::out_edge_iterator outEdgeEndIterator;
    typename SearchTree::edge_descriptor outEdge;

    boost::dynamic_bitset<> attemptedItemIds(items.size());

    for (boost::tie(outEdgeIterator, outEdgeEndIterator) = boost::out_edges(nodeId, this->tree); outEdgeIterator != outEdgeEndIterator; ++outEdgeIterator)
    {
        outEdge = *outEdgeIterator;
        typename SearchTree::vertex_descriptor targetId = boost::target(outEdge, this->tree);
        const Node& targetNode = this->tree[targetId];

        if (targetNode.ItemIdToPlace >= 0)
        {
            attemptedItemIds.set(targetNode.ItemIdToPlace);
        }
    }

    std::optional<size_t> newNodeId;
    // TODO.Logic: different branching logic, e.g. largest item or deliberately disable placement.
    ////for (size_t i = 0; i < items.size(); i++)
    for (int i: this->branchingOrder)
    {
        ////if (placedItems.contains(i) || attemptedItemIds[i])
        if (placedItems[i] || attemptedItemIds[i])
        {
            continue;
        }

        /*
        if (!IsFixedSequence(currentNode, i))
        {
            continue;
        }
        */

        typename SearchTree::vertex_descriptor leftNodeId = boost::num_vertices(this->tree);
        AddLeafNode(Node{(int)i, BaseNode::Type::UsePlacement, BaseNode::Status::NotEvaluated, std::make_unique<Packing2D>(*currentNode.Packing)});
        boost::add_edge(nodeId, leftNodeId, this->tree);

        // TODO.Safety: Caution, possible dangling reference to currentNode because the underlying graph might resize on AddLeafNode().

        this->activeNodes.push_back(leftNodeId);
        newNodeId = leftNodeId;

        break;
    }

    return newNodeId;
}

bool LeftmostActiveOnly::IsFixedSequence(const Node& node, size_t nextItemId)
{
    const auto& packing = node.Packing;
    size_t currentItemId = node.ItemIdToPlace;

    size_t numberOfPlacedItems = node.PlacedItems().count();

    std::vector<size_t> fixedSequence = {0, 10, 5, 9, 14, 11, 2, 4, 7, 3, 6, 16, 12, 1, 15, 17, 13};

    const auto& placedItems = node.Packing->Items;
    for (size_t i = 0; i < numberOfPlacedItems; i++)
    {
        const Rectangle& placedItem = placedItems[i];

        size_t itemIdToPlace = fixedSequence[i];

        if (placedItem.InternId != itemIdToPlace)
        {
            return false;
        }
    }

    if (nextItemId == -1)
    {
        return true;
    }

    size_t nextItemIdToPlace = fixedSequence[numberOfPlacedItems];
    if (nextItemIdToPlace != nextItemId)
    {
        return false;
    }

    return true;
}

std::optional<size_t> LeftmostActiveOnly::Backtrack(size_t node)
{
    // Check here instead of IsCancelled() to avoid excessive number of memory checks.
    ////CheckMemoryUsage();

    // Remove deactivated leaves.
    if (!this->nodesToDeactivate.empty())
    {
        // Delete ordersToRemove from orders according to https://en.wikipedia.org/wiki/Erase%E2%80%93remove_idiom
        this->activeNodes.erase(
            std::remove_if(
                std::begin(this->activeNodes),
                std::end(this->activeNodes),
                [&](size_t nodeToRemove)
                {
                    return this->nodesToDeactivate.contains(nodeToRemove);
                }),
            std::end(this->activeNodes));
    }

    // Clear memory from nodes that will never be used again (fathomed/pruned? What's the terminology here?).
    for (int nodeIdToDeactivate: this->nodesToDeactivate)
    {
        Node& nodeToDeactivate = this->tree[nodeIdToDeactivate];

        if (nodeIdToDeactivate == this->maxPlacedItemsInfeasibleNodeId)
        {
            continue;
        }

        nodeToDeactivate.Packing.reset();
    }

    // TODO.Performance: is this the most efficient way to clear the set?
    ////this->nodesToDeactivate = std::unordered_set<size_t>{ this->maxPlacedItemsInfeasibleNodeId };
    this->nodesToDeactivate.clear();
    this->nodesToDeactivate.emplace(this->maxPlacedItemsInfeasibleNodeId);

    if (this->activeNodes.size() == 0)
    {
        return std::nullopt;
    }

    if ((int)this->tree.m_vertices.size() - (int)this->exploredNodeDisplayCount > (int)this->exploredNodeDisplayInterval)
    {
        this->exploredNodeDisplayCount = this->tree.m_vertices.size();
        size_t outDegree = boost::out_degree(0, this->tree);

        std::cout << "LMAO \t Explored / root outedges / active nodes: " << this->tree.m_vertices.size() << " / " << outDegree << " / " << this->activeNodes.size() << "\n";
    }

    // Backtrack
    return this->activeNodes.back();

    /*
    typename SearchTree::out_edge_iterator outEdgeIterator;
    typename SearchTree::out_edge_iterator outEdgeEndIterator;
    typename SearchTree::edge_descriptor outEdge;
    for (size_t activeNodeId: this->activeNodes)
    {
        const Node& activeNode = this->tree[activeNodeId];

        // Check some condition on the active node to see if it is a candidate to continue the search.

        return activeNodeId;

        for (boost::tie(outEdgeIterator, outEdgeEndIterator) = boost::out_edges(activeNodeId, this->tree); outEdgeIterator != outEdgeEndIterator; ++outEdgeIterator)
        {
            outEdge = *outEdgeIterator;
            typename SearchTree::vertex_descriptor targetId = boost::target(outEdge, this->tree);

            ////Node& targetNode = this->tree[targetId];

            return targetId;
        }
    }
    */

    /*
    typename SearchTree::out_edge_iterator inEdgeIterator;
    typename SearchTree::out_edge_iterator inEdgeEndIterator;
    typename SearchTree::edge_descriptor inEdge;

    typename SearchTree::out_edge_iterator outEdgeIterator;
    typename SearchTree::out_edge_iterator outEdgeEndIterator;
    typename SearchTree::edge_descriptor outEdge;

    size_t currentNode = node;

    while (boost::in_degree(currentNode, this->tree) > 0)
    {
        assert(boost::in_degree(currentNode, this->tree) == 1); // There must only be on in-edge.
        inEdge = *inEdgeIterator;
        typename SearchTree::vertex_descriptor sourceId = boost::source(inEdge, this->tree);

        Node& sourceNode = this->tree[sourceId];
        auto inEdges = boost::in_edges(currentNode, this->tree);
        size_t outDegree = boost::out_degree(sourceId, this->tree);

        size_t placedItems = sourceNode.PlacedItems.size();
        size_t remainingItems = this->items.size() - placedItems;

        if (outDegree < remainingItems - placedItems)
        {
            currentNode = sourceId;
            continue;
        }

        // Logic ...
    }
    */
}

void LeftmostActiveOnly::EvaluateLeaf(Node& const node, size_t nodeId)
{
    ////Node& node = this->tree[nodeId];

    const Rectangle& item = this->items[node.ItemIdToPlace];

    bool isFeasiblePlacement = PlaceBottomLeftPlacement(node, item);

    if (isFeasiblePlacement)
    {
        node.NodeStatus = BaseNode::Status::FeasiblePlacement;
        node.Packing->PlacedItems.set((size_t)node.ItemIdToPlace);
        node.Packing->RemainingItemAreaToPlace -= item.Area;

        const Node& bestNode = this->tree[this->maxPlacedItemsInfeasibleNodeId];

        if (node.PlacedItems().count() > bestNode.PlacedItems().count())
        {
            this->maxPlacedItemsInfeasibleNodeId = nodeId;
        }

        if (this->domainReducedItemIndex == item.InternId)
        {
            node.Packing->IsReducedItemPlaced = true;
        }

        if (this->fixedItemAtOrigin.has_value() && item.InternId == this->fixedItemAtOrigin.value())
        {
            node.Packing->IsFixedItemPlaced = true;
        }
    }
    else
    {
        node.NodeStatus = BaseNode::Status::InfeasiblePlacement;

        nodesToDeactivate.emplace(nodeId);
    }

    if (node.Packing->RemainingItemAreaToPlace + node.Packing->GetDeactivatedArea() > this->container.Area
        || (!node.Packing->IsReducedItemPlaced && node.Packing->ReducedItemInfeasiblePlacementPoints >= this->reducedItemFeasiblePlacementPoints)
        || this->fixedItemAtOrigin.has_value() && !node.Packing->IsFixedItemPlaced)
    {
        node.NodeStatus = BaseNode::Status::InfeasibleSequence;

        nodesToDeactivate.emplace(nodeId);
    }
}

#pragma endregion

#pragma region Packing logic

void LeftmostActiveOnly::DeactivatePlacement(Node& const node, size_t nodeId)
{
    ////Node& node = this->tree[nodeId];
    node.NodeStatus = BaseNode::Status::DeactivatedPlacement;

    Packing2D& packing = *node.Packing;

    int placementX = packing.ActiveX;
    int placementY = packing.ActiveY;

    if (placementY == 0 && placementX == packing.MaxX)
    {
        this->nodesToDeactivate.emplace(nodeId);

        node.NodeStatus = BaseNode::Status::InfeasibleSequence;
        return;
    }

    // Determine the closest item above the current placement point. The east edge of the item is one of two bounds on the maximum x coordinate of the dummy item.
    size_t aboveMinY = packing.PlacedAreaVector.size() - 1; // TODO.Performance: - minDy
    for (size_t y = placementY; y < packing.PlacedAreaVector.size() - 1; y++)
    {
        size_t x = packing.PlacedAreaVector[y];
        ////assert(packing.Layout.DeactivatedAreaVector[y] < placementX);

        if (x > placementX)
        {
            aboveMinY = y;
            break;
        }
    }

    // Determine the other bound for the maximum x coordinate of the dummy item, i.e. the closest item below the current one.
    size_t previousX = placementY == 0 ? this->container.Dx : packing.PlacedAreaVector[placementY - 1];
    size_t xAtMinY = previousX;

    if (aboveMinY < packing.PlacedAreaVector.size() - 1)
    {
        xAtMinY = packing.PlacedAreaVector[aboveMinY];
    }

    // The maximum x coordinate of the dummy item.
    size_t disabledDx = std::min<int>(xAtMinY, previousX);

    // Determine maximum item.Dy to disable the maximum area possible.
    int minItemDy = container.Dy;
    bool infeasibleFixedSequence = false;
    const auto& placedItems = packing.PlacedItems;

    for (size_t i = 0; i < this->items.size(); i++)
    {
        const Rectangle& item = this->items[i];

        if (placedItems[i])
        {
            continue;
        }

        if (this->areItemCoordinatesFixedX
            ////&& !placedItems[i]
            && this->fixedItemCoordinatesX[i] < packing.MinX)
        {
            // Must be compared against MinX because after placing the current item, a new placement point might be generated with x < newActiveX, e.g. inverted L-shape.
            // Can be improved by accounting for item.Dy but requires more complicated checks.
            this->nodesToDeactivate.emplace(nodeId);
            node.NodeStatus = BaseNode::Status::InfeasibleSequence;
            return;
        }

        minItemDy = std::min<int>(minItemDy, item.Dy);
    }

    size_t maxY = placementY + 1;
    if (minItemDy > (aboveMinY)-placementY)
    {
        maxY = aboveMinY;
    }

    // Deactivate complete area between the two neighboring items.
    // Complete area between items (or the top of the bin) can only be deactivated if the distance to the item above (aboveMinY - placementY) is smaller than the smallest remaining item.
    // Otherwise, it cannot because feasible solutions might be excluded. Only one horizontal bar (disabledDx) at placementY can be deactivated.
    double oldCoveredArea = 0.0;
    double updatedCoveredArea = 0.0;
    size_t oldInfeasiblePlacementPoints = 0;
    size_t updatedInfeasiblePlacementPoints = 0;
    for (size_t y = placementY; y < maxY; y++)
    ////for (size_t y = placementY; y < aboveMinY; y++)
    {
        size_t placedAreaHorizontalBarAtY = packing.PlacedAreaVector[y];
        size_t deactivatedHorizontalBarAtY = packing.DeactivatedAreaVector[y];

        size_t coveredHorizontalBarAtY = std::max<int>(placedAreaHorizontalBarAtY, deactivatedHorizontalBarAtY);

        oldCoveredArea += coveredHorizontalBarAtY;
        updatedCoveredArea += disabledDx;

        if (y <= this->reducedItemDomainY)
        {
            oldInfeasiblePlacementPoints += std::min<int>(coveredHorizontalBarAtY, this->reducedItemDomainX + 1);
            updatedInfeasiblePlacementPoints += std::min<int>(disabledDx, this->reducedItemDomainX + 1);
        }

        packing.DeactivatedAreaVector[y] = disabledDx;
    }

    packing.DeactivatedArea += (updatedCoveredArea - oldCoveredArea);
    packing.ReducedItemInfeasiblePlacementPoints += (updatedInfeasiblePlacementPoints - oldInfeasiblePlacementPoints);

    if (packing.RemainingItemAreaToPlace + packing.GetDeactivatedArea() > this->container.Area
        || (!packing.IsReducedItemPlaced && packing.ReducedItemInfeasiblePlacementPoints >= this->reducedItemFeasiblePlacementPoints)
        || fixedItemAtOrigin.has_value() && !packing.IsFixedItemPlaced)
    {
        nodesToDeactivate.emplace(nodeId);

        node.NodeStatus = BaseNode::Status::InfeasibleSequence;
        return;
    }

    // Determine new active coordinates.
    auto [newActiveX, newActiveY] = FindNewBottomLeft(packing);

    packing.ActiveX = newActiveX;
    packing.ActiveY = newActiveY;

    /*
    // Snippet for managing deactivated area with boost::geometry. Worse performance than the current solution.
    // Deprecated.Logic: account for case when activeX = binDx
    if (packing.ActiveY == 0)
    {
        // TODO.Logic: add correct deactivated area.
        return node.NodeStatus;
    }

    // TODO.Logic: must this not be before if (packing.ActiveY) ?
    packing.ActiveX = packing.BackupX;
    packing.ActiveY = packing.BackupY;

    // Add dummy item of height 1.
    const Rectangle& previouslyPlacedItem = packing.Layout.Items.back();

    // Update deactivated area. Use https://www.boost.org/doc/libs/1_76_0/libs/geometry/doc/html/geometry/reference/algorithms/union_/union__3.html
    using boost::assign::tuple_list_of;
    using boost::make_tuple;
    using bg::append;

    ////typedef bg::model::polygon<bg::model::d2::point_xy<double> > polygon;
    ////typedef bg::model::polygon<boost::tuple<int, int> > polygon;

    // Deactivate the area to the left.
    Polygon newItemDummy;
    append(newItemDummy, tuple_list_of(previouslyPlacedItem.X, previouslyPlacedItem.Y)(previouslyPlacedItem.X, previouslyPlacedItem.Y + previouslyPlacedItem.Dy + 1)(previouslyPlacedItem.X + previouslyPlacedItem.Dx, previouslyPlacedItem.Y + previouslyPlacedItem.Dy + 1)(previouslyPlacedItem.X + previouslyPlacedItem.Dx, previouslyPlacedItem.Y)(previouslyPlacedItem.X, previouslyPlacedItem.Y));

    std::vector<Polygon> output;
    bg::union_(packing.Layout.DeactivatedArea, newItemDummy, output);

    ////for (const auto& p : output)
    ////{
    ////    std::cout << "Covered area = " << bg::area(p) << "\n";
    ////}

    // For MiM, output.size() > 1 can occur.
    packing.Layout.DeactivatedArea = output.front();

    return node.NodeStatus;
    */
}

bool LeftmostActiveOnly::PlaceBottomLeftPlacement(Node& node, const Rectangle& item)
{
    Packing2D& packing = *node.Packing;

    if (IsPlacementFeasible(packing, item))
    {
        packing.AddItem(Rectangle{packing.ActiveX, packing.ActiveY, item.Dx, item.Dy, item.InternId, item.ExternId}, node.ItemIdToPlace);

        int placedX = packing.ActiveX;
        int placedY = packing.ActiveY;

        // Add dummy item for area left of current placement.
        AddSuccessfulPlacementDummy(packing);

        size_t xAboveNewY = std::max<size_t>(packing.PlacedAreaVector[placedY + item.Dy], packing.DeactivatedAreaVector[placedY + item.Dy]);

        if (placedY + item.Dy == this->container.Dy || xAboveNewY > placedX + item.Dx)
        {
            auto [newActiveX, newActiveY] = FindNewBottomLeft(packing);

            packing.ActiveX = newActiveX;
            packing.ActiveY = newActiveY;
        }
        else
        {
            int newActiveX = xAboveNewY;

            packing.ActiveX = newActiveX;
            packing.ActiveY += item.Dy;
        }

        packing.MaxX = std::max<int>(packing.MaxX, placedX + item.Dx);
        if (placedY == 0)
        {
            assert(placedX + item.Dx > packing.AbscissaMaxX);
            packing.AbscissaMaxX = placedX + item.Dx;
        }

        return true;
    }

    // Placement was infeasible: Prune node, deactivate placement i, activate placement ii.
    return false;
}

bool LeftmostActiveOnly::IsPlacementFeasible(Packing2D& packing, const Rectangle& itemToPlace)
{
    if (packing.ActiveX + itemToPlace.Dx > this->container.Dx || packing.ActiveY + itemToPlace.Dy > this->container.Dy)
    {
        return false;
    }

    // TODO.Logic: For compatibility with outer B&B scheme, check if placement point is in set of MiM. Does this exclude solutions? Cannot be used in conjunction with domain reduction according to Soh (2010).
    if (!(this->itemSpecificPlacementPointsX[itemToPlace.InternId][packing.ActiveX] && this->itemSpecificPlacementPointsY[itemToPlace.InternId][packing.ActiveY]))
    {
        return false;
    }

    if (packing.AbscissaMaxX == packing.MaxX && packing.ActiveX == packing.AbscissaMaxX)
    {
        return true;
    }

    size_t minY = packing.ActiveY;
    size_t maxY = packing.ActiveY + itemToPlace.Dy;

    for (size_t y = minY; y < maxY; y++)
    {
        size_t x = packing.PlacedAreaVector[y];

        if (x > packing.ActiveX)
        {
            return false;
        }
    }

    return true;

    /*
    // Checking intersections/crossing with boost geometry is about 18% slower than the manual overlap check below.
    typedef bg::model::point<int, 2, bg::cs::cartesian> Point2D;
    typedef bg::model::linestring<Point2D> Linestring2D;

    // Checking for crossing only with the west edge will not always produce expected true/false values. Need to add north and east edge as well.
    Linestring2D newRectangleWestEdge{ { packing.ActiveX, packing.ActiveY }, { packing.ActiveX, packing.ActiveY + item.Dy }, { packing.ActiveX + item.Dx, packing.ActiveY + item.Dy }, { packing.ActiveX + item.Dx, packing.ActiveY } };

    const bool isCrossing = boost::geometry::crosses(newRectangleWestEdge, packing.Layout.DeactivatedArea);
    if (isCrossing)
    {
        return false;
    }
    */

    /*
    // Naive overlap check. Worse performance than current solution.
    for (const Rectangle& rectangle : packing.Layout.Items)
    {
        // overlapDx > 0 && overlapDy > 0
        if (std::min(rectangle.X + rectangle.Dx, packing.ActiveX + item.Dx) - std::max(rectangle.X, packing.ActiveX) > 0
            && std::min(rectangle.Y + rectangle.Dy, packing.ActiveY + item.Dy) - std::max(rectangle.Y, packing.ActiveY) > 0)
        {
            return false;
        }
    }
    */
}

std::tuple<size_t, size_t> LeftmostActiveOnly::FindNewBottomLeft(Packing2D& packing)
{
    // TODO.Test: what if abcissaMaxX = container.Dx. Test with item.Dx = container.Dx, container.Dy > item.Dy and (0, item.Dy) deactivated.

    size_t abcissaMaxX = std::max<int>(packing.PlacedAreaVector[0], packing.DeactivatedAreaVector[0]);
    size_t minX = abcissaMaxX;
    size_t minY = 0;
    size_t previousX = minX;

    size_t overallMinX = abcissaMaxX;

    for (size_t y = 0; y < container.Dy; y++)
    {
        size_t xBarPlaced = packing.PlacedAreaVector[y];
        size_t deactivatedX = packing.DeactivatedAreaVector[y];

        size_t currentMinX = std::max<int>(xBarPlaced, deactivatedX);
        overallMinX = std::min<int>(overallMinX, currentMinX);

        if (currentMinX < minX
            && xBarPlaced < previousX
            && deactivatedX < previousX)
        {
            minX = currentMinX;
            minY = y;
        }

        previousX = xBarPlaced;
    }

    packing.MinX = overallMinX;

    return std::make_tuple(minX, minY);
}

void LeftmostActiveOnly::AddSuccessfulPlacementDummy(Packing2D& packing)
{
    const Rectangle& newlyPlacedItem = packing.Items.back();

    double oldCoveredArea = 0.0;
    double updatedCoveredArea = 0.0;
    size_t oldInfeasiblePlacementPoints = 0;
    size_t updatedInfeasiblePlacementPoints = 0;
    for (size_t y = newlyPlacedItem.Y; y < newlyPlacedItem.Y + newlyPlacedItem.Dy; y++)
    {
        size_t placedAreaHorizontalBarAtY = packing.PlacedAreaVector[y];
        size_t deactivatedHorizontalBarAtY = packing.DeactivatedAreaVector[y];

        size_t coveredHorizontalBarAtY = std::max<int>(placedAreaHorizontalBarAtY, deactivatedHorizontalBarAtY);

        oldCoveredArea += coveredHorizontalBarAtY;

        size_t newHorizontalBarAtY = (size_t)newlyPlacedItem.X + (size_t)newlyPlacedItem.Dx;

        if (y <= this->reducedItemDomainY)
        {
            // Count the infeasible placement points of the domain reduced item. This information can be used to prematurely abort when there are
            // feasible placement points for left for the reduced item to be placed, i.e. when infeasiblePLacementPoints == reducedDomainArea.

            oldInfeasiblePlacementPoints += std::min<int>(coveredHorizontalBarAtY, this->reducedItemDomainX + 1);
            updatedInfeasiblePlacementPoints += std::min<int>(newHorizontalBarAtY, this->reducedItemDomainX + 1);
        }

        packing.PlacedAreaVector[y] = newHorizontalBarAtY;

        updatedCoveredArea += newHorizontalBarAtY;
    }

    assert(updatedCoveredArea >= oldCoveredArea);
    packing.DeactivatedArea += (updatedCoveredArea - oldCoveredArea);
    packing.ReducedItemInfeasiblePlacementPoints += (updatedInfeasiblePlacementPoints - oldInfeasiblePlacementPoints);

    /*
    // Snippet for managing deactivated area with boost::geometry. Worse performance than the current solution.
    // Update deactivated area. Use https://www.boost.org/doc/libs/1_76_0/libs/geometry/doc/html/geometry/reference/algorithms/union_/union__3.html
    using boost::assign::tuple_list_of;
    using boost::make_tuple;
    using bg::append;

    ////typedef boost::gg::model::polygon<bg::model::d2::point_xy<double> > polygon;
    ////typedef bg::model::polygon<boost::tuple<int, int> > polygon;

    // Deactivate the area to the left.
    Polygon newItemDummy;
    append(newItemDummy, tuple_list_of(0, newlyPlacedItem.Y)(0, newlyPlacedItem.Y + newlyPlacedItem.Dy)(newlyPlacedItem.X + newlyPlacedItem.Dx, newlyPlacedItem.Y + newlyPlacedItem.Dy)(newlyPlacedItem.X + newlyPlacedItem.Dx, newlyPlacedItem.Y)(0, newlyPlacedItem.Y));

    std::vector<Polygon> output;
    bg::union_(packing.Layout.DeactivatedArea, newItemDummy, output);

    for (const auto& item: packing.Layout.Items)
    {
        std::cout << item.X << " + " << item.Dx << ", " << item.Y << " + " << item.Dy << "\n";
    }
    ////for (const auto& p: output)
    ////{
    ////    std::cout <<  "Covered area = " << bg::area(p) << "\n";
    ////}

    // For MiM, output.size() > 1 can occur.
    packing.Layout.DeactivatedArea = output.front();
    */
}

#pragma endregion

// TODO.Logic: before DeactivatedAreaVector is now separated into PlacedAreaVector and DeactivatedAreaVector. Has not been verified in TwoStepBranchingProcedure.

#pragma region Two - step branching procedure

bool TwoStepBranchingProcedure::Node::AllItemsPlaced(const std::vector<Rectangle>& items) const
{
    return items.size() == PlacedItems().count();
}

void TwoStepBranchingProcedure::AddItems(const std::vector<Rectangle>& items, int fixedItemAtOrigin)
{
    if (this->container.Dx == 0)
    {
        throw std::exception("Container not assigned.");
    }

    this->tree.m_vertices.reserve(items.size() + 1);
    this->items = items;

    this->fixedItemAtOrigin = fixedItemAtOrigin;

    // The intern IDs of the items must correpond to the indices of items.
    for (size_t i = 0; i < items.size(); i++)
    {
        const Rectangle& item = items[i];

        if (i != item.InternId)
        {
            throw std::exception("Internal item IDs do not match.");
        }
    }
}

void TwoStepBranchingProcedure::AddContainer(const Bin& container)
{
    this->container = container;
}

void TwoStepBranchingProcedure::AddCancellationToken(std::atomic<bool>* cancellationToken)
{
    this->cancellationToken = cancellationToken;
}

bool TwoStepBranchingProcedure::IsCancelled() const
{
    // ((this->cancellationToken == nullptr) || (this->cancellationToken != nullptr && !(this->cancellationToken->load())))
    return this->cancellationToken != nullptr && this->cancellationToken->load();
}

void TwoStepBranchingProcedure::SignalCancelling()
{
    if (this->cancellationToken != nullptr)
    {
        this->cancellationToken->store(true);
    }
}

/*
bool TwoStepBranchingProcedure::IsFixedSequence(const Node& node, size_t nextItemId)
{
    const auto& packing = node.Packing;
    size_t currentItemId = node.ItemIdToPlace;

    size_t numberOfPlacedItems = node.PlacedItems().count();

    std::vector<size_t> fixedSequence = {1, 5, 12, 4, 10, 14, 6, 8, 3, 7, 0, 11, 13, 15, 16, 2, 9};
    ////std::vector<size_t> fixedSequence = {1, 5, 12, 4, 10, 14, 6, 8, 7, 0, 11, 13, 15, 16, 2, 9};

    const auto& placedItems = node.Packing->ItemSequence;
    for (size_t i = 0; i < numberOfPlacedItems; i++)
    {
        int placedItemId = placedItems[i];

        size_t itemIdToPlace = fixedSequence[i];

        if (placedItemId != itemIdToPlace)
        {
            return false;
        }
    }

    if (nextItemId == -1)
    {
        return true;
    }

    size_t nextItemIdToPlace = fixedSequence[numberOfPlacedItems];
    if (nextItemIdToPlace != nextItemId)
    {
        return false;
    }

    return true;
}
*/

SearchStatus TwoStepBranchingProcedure::Solve()
{
    auto start = std::chrono::steady_clock::now();

    SearchStatus searchStatus = SearchStatus::None;
    int numberOfThreads = this->parameters.TSBPThreads;

    std::vector<int> itemOrder(items.size());
    std::iota(itemOrder.begin(), itemOrder.end(), 0);
    std::sort(itemOrder.begin(), itemOrder.end(), CompareLexicographicDxDy(&items));

    this->branchingOrder = itemOrder;

    if (numberOfThreads == 1)
    {
        searchStatus = SolveSequential();
        this->statistics.NodeCountTSBP = this->tree.m_vertices.size();
    }
    else if (numberOfThreads >= 1)
    {
        searchStatus = SolveParallelTaskflow();
        ////searchStatus = SolveParallelNative();
    }
    else
    {
        std::cout << "Invalid number of threads parameter, solve in sequential mode.\n";
        searchStatus = SolveSequential();
    }

    auto end = std::chrono::steady_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    this->statistics.TimeOverall = time;
    this->statistics.TimeTSBP = time - this->statistics.TimeLMAO;
    ////this->statistics.MemoryUsage = MemoryUsedByCurrentProcess();

    return searchStatus;
}

SearchStatus TwoStepBranchingProcedure::SolveParallelTaskflow()
{
    std::atomic<std::chrono::steady_clock::time_point> begin = std::chrono::steady_clock::now();
    std::atomic<std::chrono::steady_clock::time_point> endLast = begin.load();

    size_t numberOfThreads = this->parameters.TSBPThreads;

    tf::Executor executor(numberOfThreads);
    tf::Taskflow taskflow;

    std::atomic<size_t> exploredSubtrees(0);
    std::atomic<int64_t> exploredNodes(0);
    std::atomic<int64_t> callCountLMAO(0);
    std::atomic<int64_t> exploredLMAO(0);

    std::atomic<bool> feasibleSolutionCancellationToken = false;

    std::mutex mutex;
    Packing2D solution;

    SearchStatus status = SearchStatus::None;

    for (int itemId: this->branchingOrder)
    {
        tf::Task itemSpecificBranch = taskflow.emplace(
            [&, itemId]()
            {
                    // TODO.Performance: write copy assignment constructor for Node to run the preprocess only once. It has its copy assignment constructor disabled because it has a unique_ptr member.
                    ////LeftmostActiveOnly lmao = *this;
                    TwoStepBranchingProcedure tsbp(this->parameters);
                    tsbp.AddContainer(container);
                    tsbp.AddItems(items, itemId);
                    tsbp.SetBranchingOrder(this->branchingOrder);
                    tsbp.AddCancellationToken(&feasibleSolutionCancellationToken);

                    SearchStatus localStatus = tsbp.SolveSequential();
                    if (localStatus == SearchStatus::Abort)
                    {
                        status = localStatus;
                    }

                    if (localStatus == SearchStatus::Feasible)
                    {
                        mutex.lock();
                        solution = tsbp.GetSolution();
                        mutex.unlock();
                    }

                    exploredNodes += tsbp.GetTreeSize();
                    callCountLMAO += tsbp.GetStatistics().CallCountLMAO;
                    exploredLMAO += tsbp.GetStatistics().NodeCountLMAO;
                    exploredSubtrees++;

                    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
                    endLast = end;

                    ////if (this->parameters.TSBPEnableLogging)
                    if (false)
                    {
                        std::cout << "Explored subtrees = " + std::to_string(exploredSubtrees) + " / " + std::to_string(items.size()) + "\n";

                        std::cout << "Elapsed time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin.load()).count() << "s\n";
                    } });
    }

    // TODO.Performance: Consider prioritizing the reduced domain item and smaller items. It takes the longest to compute and should not come late in the process, s.t. the other threads cannot be used as no computations are left.
    // TODO.Performance: Execution with the taskflow library takes longer in sequential and parallel mode.
    // In sequential mode, the issue seems to be that a task/thread is blocked unitl after it has finished computation *and* until its memory fully deallocated.
    // This also has an impact in parallel mode.
    // Look into async tasks https://taskflow.github.io/taskflow/AsyncTasking.html and parallel algorithms.

    tf::Future<void> fu = executor.run(taskflow);
    fu.wait(); // block until the execution completes

    std::chrono::steady_clock::time_point endTaskflow = std::chrono::steady_clock::now();

    if (status == SearchStatus::Abort)
    {
        std::cout << "Search was aborted.\n";
        return status;
    }

    if (feasibleSolutionCancellationToken.load())
    {
        status = SearchStatus::Feasible;
        this->solution = solution;
        std::cout << "Feasible solution found.\n";
    }
    else
    {
        status = SearchStatus::Infeasible;
        std::cout << "Problem is infeasible.\n";
    }

    ////if (this->parameters.TSBPEnableLogging)
    if (false)
    {
        std::cout << "TSBP \t Total explored nodes = " + std::to_string(exploredNodes.load()) + "\n";
    }

    this->statistics.NodeCountTSBP = exploredNodes.load();
    this->statistics.NodeCountLMAO = exploredLMAO.load();
    this->statistics.TimeMemoryDeallocation = std::chrono::duration_cast<std::chrono::milliseconds>(endTaskflow - endLast.load()).count();

    // TODO.Logic: retrieve return value of different tasks.

    return status;
}

SearchStatus TwoStepBranchingProcedure::SolveParallelNative()
{
    std::atomic<std::chrono::steady_clock::time_point> begin = std::chrono::steady_clock::now();
    std::atomic<std::chrono::steady_clock::time_point> endLast = begin.load();

    std::atomic<bool> feasibleSolutionCancellationToken = false;

    std::vector<std::unique_ptr<TwoStepBranchingProcedure>> subtrees;
    subtrees.reserve(items.size());
    ////for (size_t i = 0; i < 1; i++)
    for (size_t itemId: this->branchingOrder)
    {
        std::unique_ptr<TwoStepBranchingProcedure> tsbp = std::make_unique<TwoStepBranchingProcedure>(this->parameters);
        tsbp->AddContainer(container);
        tsbp->AddItems(items, itemId);
        tsbp->SetBranchingOrder(this->branchingOrder);
        tsbp->AddCancellationToken(&feasibleSolutionCancellationToken); // Optional: comment out for sequential mode.
        ////tsbp->AddFixedItem(fixedItem);

        subtrees.push_back(std::move(tsbp));
    }

    SearchStatus status = SearchStatus::None;

    std::atomic<size_t> exploredSubtrees(0);
    std::atomic<int64_t> exploredNodes(0);
    std::atomic<int64_t> callCountLMAO(0);
    std::atomic<int64_t> exploredLMAO(0);
    std::atomic<bool> isFeasible(false);
    Packing2D solution;

    std::mutex mutex;

    std::for_each(
        std::execution::par_unseq,
        subtrees.begin(),
        subtrees.end(),
        [&](std::unique_ptr<TwoStepBranchingProcedure>& tsbp)
        {
            SearchStatus localStatus = tsbp->SolveSequential();
            if (localStatus == SearchStatus::Abort)
            {
                status = localStatus;
            }

            if (localStatus == SearchStatus::Feasible)
            {
                mutex.lock();
                solution = tsbp->GetSolution();
                mutex.unlock();

                isFeasible.store(true);
            }

            exploredNodes += tsbp->GetTreeSize();
            callCountLMAO += tsbp->GetStatistics().CallCountLMAO;
            exploredLMAO += tsbp->GetStatistics().NodeCountLMAO;
            exploredSubtrees++;

            std::cout << "Explored subtrees = " + std::to_string(exploredSubtrees) + " / " + std::to_string(items.size()) + "\n";

            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            endLast = end;
            std::cout << "Elapsed time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin.load()).count() << "s\n";

            tsbp.reset();
        });

    if (status == SearchStatus::Abort)
    {
        std::cout << "Search was aborted.\n";
        return status;
    }

    if (isFeasible.load())
    {
        this->solution = solution;
        status = SearchStatus::Feasible;
        std::cout << "Feasible solution found.\n";
    }
    else
    {
        status = SearchStatus::Infeasible;
        std::cout << "Problem is infeasible.\n";
    }

    std::chrono::steady_clock::time_point endParallelFor = std::chrono::steady_clock::now();

    std::cout << "Elapsed time until after for_each() = " << std::chrono::duration_cast<std::chrono::seconds>(endParallelFor - begin.load()).count() << "s\n";
    std::cout << "Termination delay = " << std::chrono::duration_cast<std::chrono::seconds>(endParallelFor - endLast.load()).count() << "s\n";

    std::cout << "Total explored nodes = " + std::to_string(exploredNodes.load()) + "\n";

    this->statistics.NodeCountTSBP = exploredNodes.load();
    this->statistics.NodeCountLMAO = exploredLMAO.load();
    this->statistics.TimeMemoryDeallocation = std::chrono::duration_cast<std::chrono::milliseconds>(endParallelFor - endLast.load()).count();

    return status;
}

void TwoStepBranchingProcedure::UpdateSolverStatistics(SolverStatistics& const statistics) const
{
    statistics.CallCountLMAO += this->statistics.CallCountLMAO;
    statistics.NodeCountLMAO += this->statistics.NodeCountLMAO;
    statistics.TimeLMAO += this->statistics.TimeLMAO;

    statistics.NodeCountTSBP += this->statistics.NodeCountTSBP;
    statistics.TimeTSBP += this->statistics.TimeTSBP;

    statistics.TimeMemoryDeallocation += this->statistics.TimeMemoryDeallocation;

    statistics.MemoryUsage = std::max<double>(statistics.MemoryUsage, this->statistics.MemoryUsage);
}

Packing2D TwoStepBranchingProcedure::GetSolution() const
{
    return this->solution;
}

void TwoStepBranchingProcedure::Preprocess()
{
    size_t domainReducedItemIndex = FindMaximumItemDx();

    this->domainReducedItemIndex = domainReducedItemIndex;

    DomainReductionX(domainReducedItemIndex);
}

size_t TwoStepBranchingProcedure::FindMaximumItemDx()
{
    int maxDx = this->items[0].Dx;

    size_t domainReducedItemIndex = 0;
    for (size_t i = 1; i < this->items.size(); i++)
    {
        const Rectangle& item = this->items[i];

        if (item.Dx > maxDx)
        {
            maxDx = item.Dx;
            domainReducedItemIndex = i;
        }
    }

    return domainReducedItemIndex;
}

void TwoStepBranchingProcedure::DomainReductionX(size_t domainReducedItemIndex)
{
    itemSpecificPlacementPointsX.reserve(this->items.size());

    int reducedDomainThreshold = std::floor((this->container.Dx - this->items[domainReducedItemIndex].Dx) / 2.0);

    for (size_t i = 0; i < this->items.size(); i++)
    {
        const Rectangle& item = this->items[i];
        if (item.InternId != i)
        {
            throw std::exception("Item IDs do not match in TSBP.");
        }

        if (i == domainReducedItemIndex)
        {
            auto& bitsetX = this->itemSpecificPlacementPointsX.emplace_back(boost::dynamic_bitset<>(this->container.Dx));
            bitsetX.reset();

            for (size_t j = 0; j <= reducedDomainThreshold; j++)
            {
                bitsetX.set(j);
            }

            this->reducedItemDomainX = reducedDomainThreshold;

            continue;
        }

        // TODO: reduce by item.Dy to exclude overlap with the container. Enables simplification of no-overlap checks.
        auto& bitsetX = this->itemSpecificPlacementPointsX.emplace_back(boost::dynamic_bitset<>(this->container.Dx));
        bitsetX.set();
    }
}

size_t TwoStepBranchingProcedure::InitializeSearchTree()
{
    size_t rootNodeId = boost::add_vertex(this->tree);
    this->activeNodes.push_back(rootNodeId);

    Node& rootNode = this->tree[rootNodeId];
    rootNode.ItemIdToPlace = -1;
    rootNode.NodeStatus = BaseNode::Status::Root;

    std::unique_ptr<PackingRelaxed2D> layout2D = std::make_unique<PackingRelaxed2D>();
    layout2D->ItemCoordinatesX.reserve(items.size());
    layout2D->PlacedItems = boost::dynamic_bitset<>(items.size());
    layout2D->PlaceableItemsAtCurrentX = ~layout2D->PlacedItems;
    ////layout2D.PlacedItems = std::unordered_set<size_t>();
    ////layout2D.PlacedItems.reserve(this->items.size());
    layout2D->RemainingItemAreaToPlace = std::accumulate(
        items.begin(),
        items.end(),
        0.0,
        [](double volume, const Rectangle& i)
        { return volume + i.Area; });
    layout2D->PlacedAreaVector = std::vector<size_t>(container.Dx + 1, 0);
    layout2D->ItemCoordinatesX = std::vector<int>(items.size(), -1);
    layout2D->MinimumRemainingItemDy = container.Dy;
    layout2D->MinimumRemainingItemDyAtCurrentX = container.Dy;

    rootNode.Packing = std::move(layout2D);

    if (this->fixedItemAtOrigin.has_value())
    {
        // Restrict placement of fixedItem
        this->itemSpecificPlacementPointsX[this->fixedItemAtOrigin.value()].reset();
        this->itemSpecificPlacementPointsX[this->fixedItemAtOrigin.value()].set(0);

        ////this->tree.m_vertices.reserve(250000000);
    }

    return rootNodeId;
}

SearchStatus TwoStepBranchingProcedure::SolveSequential()
{
    Preprocess();

    size_t nodeId = InitializeSearchTree();

    size_t previousNodeId = nodeId;

    // TODO.Performance: parallelize tree search, e.g. a separate sub-tree for each node starting at depth 1
    this->searchStatus = SearchStatus::InProgress;
    while (this->searchStatus == SearchStatus::InProgress && !IsCancelled())
    {
        previousNodeId = nodeId;
        std::optional<size_t> newNodeId = Branch(nodeId);

        if (!newNodeId.has_value())
        {
            // No new node could be created: an infeasible sequence has been found.
            std::optional<size_t> backtrackedNodeId = Backtrack(nodeId);

            if (!backtrackedNodeId.has_value())
            {
                this->searchStatus = SearchStatus::Infeasible;
                break;
            }

            nodeId = backtrackedNodeId.value();
            continue;
        }

        Node& newNode = this->tree[newNodeId.value()];

        switch (newNode.NodeType)
        {
            case BaseNode::Type::UsePlacement:
                EvaluateLeaf(newNode, newNodeId.value());
                break;
            case BaseNode::Type::DeactivatePlacement:
                DeactivatePlacement(newNode, newNodeId.value());
                break;
            default:
                throw std::exception("Invalid node type.");
        }

        if (newNode.AllItemsPlaced(this->items))
        {
            bool isFeasible = SolveSubproblem(newNode, newNodeId.value());

            if (isFeasible)
            {
                break;
            }
        }

        // SelectNextNode(nodeStatus, previousNodeId) return nodeId to branch on
        switch (newNode.NodeStatus)
        {
            case BaseNode::Status::InfeasiblePlacement:
                // Prune node. Backtrack a single edge and select next most promising node to branch on.
                nodeId = previousNodeId;
                break;
            case BaseNode::Status::FeasiblePlacement:
                // Fallthrough.
            case BaseNode::Status::DeactivatedPlacement:
                // Fallthrough.
            case BaseNode::Status::InfeasibleSequence:
                // Continue search at current node (DFS).
                nodeId = newNodeId.value();
                break;
            default:
                throw std::exception("Invalid node status.");
        }
    }

    size_t bestNodeId = this->maxPlacedItemsInfeasibleNodeId;
    if (this->cancellationToken != nullptr && this->cancellationToken->load())
    {
        this->searchStatus = SearchStatus::Cancelled;
        bestNodeId = 0;
    }

    if (this->searchStatus == SearchStatus::Feasible)
    {
        bestNodeId = this->feasibleLeafNodeId;

        SignalCancelling();
    }

    std::cout << "Search status: " << (int)this->searchStatus << "\n";

    if (bestNodeId == -1)
    {
        throw std::exception("Search status is feasible but feasible leaf node is -1.");
    }

    if (fixedItemAtOrigin.has_value())
    {
        std::cout << "Fixed item (intern id): " + std::to_string(this->fixedItemAtOrigin.value()) + ", explored nodes = " + std::to_string(GetTreeSize()) + "\n";
    }
    else
    {
        std::cout << "Explored nodes = " + std::to_string(GetTreeSize()) + "\n";
    }

    switch (this->searchStatus)
    {
        case SearchStatus::Feasible:
        case SearchStatus::Infeasible:
            break;
        case SearchStatus::Cancelled:
            break;
        case SearchStatus::Abort:
        default:
            throw std::exception("Search status after search is neither feasible nor infeasible.");
    }

    return this->searchStatus;
}

void TwoStepBranchingProcedure::DeactivatePlacement(Node& const node, size_t nodeId)
{
    ////Node& node = this->tree[nodeId];
    node.NodeStatus = BaseNode::Status::DeactivatedPlacement;

    PackingRelaxed2D& packing = *node.Packing;

    // Fill with dummy item.
    double oldArea = 0.0;
    double newArea = 0.0;

    size_t lastX = packing.MinX;
    size_t lastY = packing.PlacedAreaVector[packing.MinX];
    for (size_t x = packing.MinX; x < this->container.Dx; x++)
    {
        size_t currentY = packing.PlacedAreaVector[x];
        if (currentY < lastY)
        {
            assert(x > packing.MinX);
            // TODO.Performance: check against packing.MinimumRemainingItemDy. If packing.MinimumRemainingItemDy > container.Dy - currentY, dont break out of the for loop and keep iterating to find the new x. Also, dont contine, must update deactivated area!
            // This way, many child nodes will not be enumerated when no placement at the new x is possible anyways.
            // There should be no benefit in updating the lower bound packing.MinimumRemainingItemDyAtCurrentX here because it differs from packing.MinimumRemainingItemDy only in that
            // the former excludes items with IDs less than items currently placed at x. However, by construction, no item can be placed at x, yet.
            lastX = x;
            lastY = currentY;
            break;
        }

        oldArea += currentY;
        newArea += this->container.Dy;
        packing.PlacedAreaVector[x] = this->container.Dy;
    }

    // TODO: consider updating PlaceableItemsAtCurrentX if item.Dy < container.Dy - currentY to reduce number of child nodes created.
    packing.DeactivatedArea += (newArea - oldArea);
    packing.MinX = lastX;
    packing.MinimumRemainingItemDyAtCurrentX = packing.MinimumRemainingItemDy;
    packing.PlaceableItemsAtCurrentX = ~packing.PlacedItems;

    if (packing.RemainingItemAreaToPlace + packing.GetDeactivatedArea() > this->container.Area
        || (!packing.IsReducedItemPlaced && packing.MinX > this->reducedItemDomainX)
        || fixedItemAtOrigin.has_value() && !packing.IsFixedItemPlaced)
    {
        node.NodeStatus = BaseNode::Status::InfeasibleSequence;

        nodesToDeactivate.emplace(nodeId);
    }
}

std::optional<size_t> TwoStepBranchingProcedure::Branch(size_t nodeId)
{
    Node& const currentNode = this->tree[nodeId];
    const boost::dynamic_bitset<>& placedItems = currentNode.PlacedItems();
    ////const std::unordered_set<size_t>& placedItems = currentNode.PlacedItems();

    if (currentNode.NodeStatus == BaseNode::Status::InfeasibleSequence)
    {
        return std::nullopt;
    }

    const PackingRelaxed2D& packing = *currentNode.Packing;

    // Determine number of outgoing edges.
    size_t numberOfOutgoingEdges = boost::out_degree(nodeId, this->tree);
    size_t numberOfRemainingItems = packing.PlaceableItemsAtCurrentX.count();
    ////size_t numberOfRemainingItems = items.size() - placedItems.count();

    // TODO.Logic: currently the DeactivatePlacement node is the last node that is generated. Rework to allow it to be explored in any sequence.
    ////bool isFixedSequenceDeactivation = numberOfOutgoingEdges == 1 && IsFixedSequence(currentNode, -1);
    ////if (numberOfOutgoingEdges == numberOfRemainingItems || isFixedSequenceDeactivation)
    if (numberOfOutgoingEdges == numberOfRemainingItems)
    {
        // AddDeactivatedPlacementNode()
        ////currentNode.Packing.ActiveX = currentNode.Packing.BackupX;
        ////currentNode.Packing.ActiveY = currentNode.Packing.BackupY;

        this->nodesToDeactivate.emplace(nodeId);

        if (currentNode.NodeStatus == BaseNode::Status::Root)
        {
            this->nodesToDeactivate.emplace(0);
            return std::nullopt;
        }

        if (currentNode.NodeType == BaseNode::Type::DeactivatePlacement)
        {
            return std::nullopt;
        }

        typename SearchTree::vertex_descriptor leftNodeId = boost::num_vertices(this->tree);
        AddLeafNode(Node{-1, BaseNode::Type::DeactivatePlacement, BaseNode::Status::NotEvaluated, std::make_unique<PackingRelaxed2D>(packing)});
        boost::add_edge(nodeId, leftNodeId, this->tree);

        // Caution, possible dangling reference to currentNode because the underlying graph might resize on AddLeafNode().
        // TODO.Performance: Consider using std::unique_ptr<Node> to reduce copy effort.

        this->activeNodes.push_back(leftNodeId);

        return leftNodeId;
    }

    typename SearchTree::out_edge_iterator outEdgeIterator;
    typename SearchTree::out_edge_iterator outEdgeEndIterator;
    typename SearchTree::edge_descriptor outEdge;

    boost::dynamic_bitset<> attemptedItemIds(items.size());

    for (boost::tie(outEdgeIterator, outEdgeEndIterator) = boost::out_edges(nodeId, this->tree); outEdgeIterator != outEdgeEndIterator; ++outEdgeIterator)
    {
        outEdge = *outEdgeIterator;
        typename SearchTree::vertex_descriptor targetId = boost::target(outEdge, this->tree);
        const Node& targetNode = this->tree[targetId];

        if (targetNode.ItemIdToPlace >= 0)
        {
            attemptedItemIds.set(targetNode.ItemIdToPlace);
        }
    }

    std::optional<size_t> newNodeId;
    // TODO.Logic: different branching logic, e.g. largest item or deliberately disable placement.
    ////for (size_t i = 0; i < items.size(); i++)
    for (int i: this->branchingOrder)
    {
        ////if (placedItems.contains(i) || attemptedItemIds[i])
        if (placedItems[i] || attemptedItemIds[i] || !packing.PlaceableItemsAtCurrentX[i])
        {
            continue;
        }

        /*
        if (!IsFixedSequence(currentNode, i))
        {
            continue;
        }
        */

        typename SearchTree::vertex_descriptor leftNodeId = boost::num_vertices(this->tree);
        AddLeafNode(Node{(int)i, BaseNode::Type::UsePlacement, BaseNode::Status::NotEvaluated, std::make_unique<PackingRelaxed2D>(packing)});
        boost::add_edge(nodeId, leftNodeId, this->tree);

        // Caution, possible dangling reference to currentNode because the underlying graph might resize on AddLeafNode().

        this->activeNodes.push_back(leftNodeId);
        newNodeId = leftNodeId;

        break;
    }

    return newNodeId;
}

std::optional<size_t> TwoStepBranchingProcedure::Backtrack(size_t currentNode)
{
    // Check here instead of IsCancelled() to avoid excessive number of memory checks.
    // According to https://stackoverflow.com/a/64166/5587903.
    ////CheckMemoryUsage();

    // Remove deactivated leaves.
    if (!this->nodesToDeactivate.empty())
    {
        // Delete ordersToRemove from orders according to https://en.wikipedia.org/wiki/Erase%E2%80%93remove_idiom
        this->activeNodes.erase(
            std::remove_if(
                std::begin(this->activeNodes),
                std::end(this->activeNodes),
                [&](size_t nodeToRemove)
                {
                    return this->nodesToDeactivate.contains(nodeToRemove);
                }),
            std::end(this->activeNodes));
    }

    // Clear memory from nodes that will never be used again (fathomed/pruned? What's the terminology here?).
    for (int nodeIdToDeactivate: this->nodesToDeactivate)
    {
        Node& nodeToDeactivate = this->tree[nodeIdToDeactivate];

        if (nodeIdToDeactivate == this->maxPlacedItemsInfeasibleNodeId)
        {
            continue;
        }

        nodeToDeactivate.Packing.reset();
    }

    // TODO.Performance: is this the most efficient way to clear the set?
    ////this->nodesToDeactivate = std::unordered_set<size_t>{ this->maxPlacedItemsInfeasibleNodeId };
    this->nodesToDeactivate.clear();
    this->nodesToDeactivate.emplace(this->maxPlacedItemsInfeasibleNodeId);

    if (this->activeNodes.size() == 0)
    {
        return std::nullopt;
    }

    if ((int)this->tree.m_vertices.size() - (int)this->exploredNodeDisplayCount > (int)this->exploredNodeDisplayInterval)
    {
        this->exploredNodeDisplayCount = this->tree.m_vertices.size();
        size_t outDegree = boost::out_degree(0, this->tree);

        std::cout << "TSBP \t Explored / root outedges / active nodes: " << this->tree.m_vertices.size() << " / " << outDegree << " / " << this->activeNodes.size() << "\n";
    }

    // Backtrack
    return this->activeNodes.back();
}

void TwoStepBranchingProcedure::EvaluateLeaf(Node& const node, size_t nodeId)
{
    ////Node& node = this->tree[nodeId];
    PackingRelaxed2D& packing = *node.Packing;

    const Rectangle& item = this->items[node.ItemIdToPlace];

    bool isFeasiblePlacement = PlaceAtCurrentX(node, item);

    if (isFeasiblePlacement)
    {
        node.NodeStatus = BaseNode::Status::FeasiblePlacement;
        packing.PlacedItems.set((size_t)node.ItemIdToPlace);
        ////packing.Layout.PlacedItems.emplace((size_t)node.ItemIdToPlace);
        packing.RemainingItemAreaToPlace -= item.Area;

        ////packing.ItemSequence.push_back(item.InternId);

        const Node& bestNode = this->tree[maxPlacedItemsInfeasibleNodeId];

        if (node.PlacedItems().count() > bestNode.PlacedItems().count())
        {
            this->maxPlacedItemsInfeasibleNodeId = nodeId;
        }

        if (this->domainReducedItemIndex == item.InternId)
        {
            packing.IsReducedItemPlaced = true;
        }

        if (this->fixedItemAtOrigin.has_value() && item.InternId == this->fixedItemAtOrigin.value())
        {
            node.Packing->IsFixedItemPlaced = true;
        }
    }
    else
    {
        node.NodeStatus = BaseNode::Status::InfeasiblePlacement;

        // Is reset unnecessary because search node is pruned?
        packing.MaximumRemainingItemDx = std::max<int>(packing.MaximumRemainingItemDx, item.Dx);

        nodesToDeactivate.emplace(nodeId);
    }

    if (packing.RemainingItemAreaToPlace + packing.GetDeactivatedArea() > this->container.Area
        || (packing.MinX + packing.MaximumRemainingItemDx > this->container.Dx)
        || (!packing.IsReducedItemPlaced && packing.MinX > this->reducedItemDomainX)
        || this->fixedItemAtOrigin.has_value() && !node.Packing->IsFixedItemPlaced)
    {
        node.NodeStatus = BaseNode::Status::InfeasibleSequence;

        nodesToDeactivate.emplace(nodeId);
    }
}

bool TwoStepBranchingProcedure::IsPlacementFeasible(PackingRelaxed2D& packing, const Rectangle& itemToPlace)
{
    // No need to check feasibility in y-direction for x > packing.MinX, because, by construction, packing.PlacedAreaVector (K(x)) is non-increasing.
    if (packing.MinX + itemToPlace.Dx > this->container.Dx || packing.PlacedAreaVector[packing.MinX] + itemToPlace.Dy > this->container.Dy)
    {
        return false;
    }

    if (!this->itemSpecificPlacementPointsX[itemToPlace.InternId][packing.MinX])
    {
        return false;
    }

    int maximumRemainingItemDx = 0;
    int minimumRemainingItemDy = this->container.Dy;
    int minimumRemainingItemDyAtCurrentX = this->container.Dy;
    for (size_t i = 0; i < items.size(); i++)
    {
        const Rectangle& item = items[i];
        // Caution, requires i == items[i].InternId
        if (packing.PlacedItems.test(item.InternId))
        {
            if (packing.ItemCoordinatesX[item.InternId] == packing.MinX
                && itemToPlace.InternId < item.InternId)
            {
                return false;
            }
        }
        else
        {
            if (item.InternId != itemToPlace.InternId)
            {
                minimumRemainingItemDy = std::min<int>(minimumRemainingItemDy, item.Dy);
                maximumRemainingItemDx = std::max<int>(maximumRemainingItemDx, item.Dx);
            }

            if (itemToPlace.InternId < item.InternId)
            {
                minimumRemainingItemDyAtCurrentX = std::min<int>(minimumRemainingItemDyAtCurrentX, item.Dy);
            }
        }
    }

    ////packing.MinimumRemainingItemDy = 0;
    ////packing.MinimumRemainingItemDyAtCurrentX = 0;
    packing.MinimumRemainingItemDy = minimumRemainingItemDy;
    packing.MinimumRemainingItemDyAtCurrentX = minimumRemainingItemDyAtCurrentX;
    packing.MaximumRemainingItemDx = maximumRemainingItemDx;

    return true;
}

bool TwoStepBranchingProcedure::SolveSubproblem(Node& newNode, const size_t newNodeId)
{
    ////LeftmostActiveOnly leftmostActiveOnlyAlgorithm(this->items, this->container, newNode.Packing->ItemCoordinatesX, this->parameters);
    ////SearchStatus status = leftmostActiveOnlyAlgorithm.Solve();

    LeftmostActiveOnly leftmostActiveOnlyAlgorithm(this->parameters);
    leftmostActiveOnlyAlgorithm.AddContainer(container);
    leftmostActiveOnlyAlgorithm.AddItems(items, std::nullopt, newNode.Packing->ItemCoordinatesX);

    SearchStatus status = leftmostActiveOnlyAlgorithm.Solve();

    switch (status)
    {
        case SearchStatus::Feasible:
            this->searchStatus = SearchStatus::Feasible;
            this->feasibleLeafNodeId = newNodeId;
            break;
        case SearchStatus::Infeasible:
            this->nodesToDeactivate.emplace(newNodeId);
            newNode.NodeStatus = BaseNode::Status::InfeasibleSequence;
            break;
        case SearchStatus::Cancelled:
        case SearchStatus::Abort:
        default:
            throw std::exception("Lower level branch and bound search returned neither feasible nor infeasible. The search cannot continue.");
    }

    if (this->searchStatus == SearchStatus::Feasible)
    {
        leftmostActiveOnlyAlgorithm.UpdateSolverStatistics(this->statistics);
        this->solution = leftmostActiveOnlyAlgorithm.GetSolution();
        return true;
    }

    return false;
}

size_t TwoStepBranchingProcedure::UpdatePlacedAreaAndDetermineLoss(PackingRelaxed2D& packing, const Rectangle& item)
{
    // This lower bound can be raised by determining the remaining minY at newX at a cost of O(n) or by solving a subset sum problem; cf. SS_1 in [CJCM2008] section 5.2.
    // Instead, what we do here is to see if there exists any item at current x that can be placed, which was already determined in IsPlacementFeasible().
    // If none can be placed, no SSP must be solved. The space between the item and the container can be marked as lost by deactivating the area.
    // By solving SSPs, infeasibility by detected lost space might be detected earlier during search and save a lot of child nodes. When SSPs are solved, there should be a
    // condition for when to solve them to not unnecessarily expend time. For example, it would make sense to only do so if the used area plus remaining area of items to place is
    // close to the container area.
    // TODO.Performance: subset sum can be solved faster than classically in O(nC), where n is the number of items and C is the (remaining) capacity/height at at current x,
    // namely in O(n max Z), \tilde{O}(sqrt(n)C), \tilde{O}(n + C); cp. https://arxiv.org/abs/1610.04712.
    // Further references:
    // - https://en.wikipedia.org/wiki/Subset_sum_problem#Pseudo-polynomial_time_dynamic_programming_solutions
    // - https://en.wikipedia.org/wiki/Subset_sum_problem#Polynomial_time_approximate_algorithm
    // - For the \tilde{O}(sqrt(n)C) algorithm and a good comparison of algorithms, see https://arxiv.org/pdf/1507.02318v1.pdf and table 1.1 therein. A Java implementation is available: https://github.com/shtratos/subsetsum.
    // - For the \tilde{O}(n + C) algorithm, see https://arxiv.org/abs/1610.04712 (the error term means that is is not guaranteed to be exact?).
    //   It also states that "Using number-theoretic arguments, certain dense cases of SubsetSum are solvable in near-linear time, e.g., if t << n^2 then there is an \tilde{O}(n) algorithm [11]",
    //   and t << n^2 might be satisfied in our case.
    // - The complexities stated in the more recent paper https://arxiv.org/abs/1610.04712 use a tilde above big O, which means that logarithmic factors are ignored, see https://www.johndcook.com/blog/2019/01/17/big-o-tilde-notation/.
    //   Here, we deal mostly with small instances of SSPs, so ignored logarithmic factors should be accounted for as it might have a disproportionate impact for small n or a considerable tradeoff in space complexity.
    // Probably just go for the classical O(nC).
    int deactivatedDy = packing.PlacedAreaVector[packing.MinX] + item.Dy > container.Dy - packing.MinimumRemainingItemDyAtCurrentX
                            ? container.Dy - packing.PlacedAreaVector[packing.MinX]
                            : item.Dy;

    size_t newX = packing.MinX + item.Dx;
    size_t lastY = packing.PlacedAreaVector[packing.MinX];
    for (size_t x = packing.MinX; x < packing.MinX + item.Dx; x++)
    {
        size_t currentY = packing.PlacedAreaVector[x];

        if (currentY < lastY)
        {
            // K(x) is non-increasing, so the first occurence where currentY < lastY is the next placement position.
            // When currentY < lastY corresponds to the values that the function \roh(Z) in [CJCM2008] returns: "the set of coordinates Z_1 at which the value H_x moves".

            // Mus be reset because the next placement point x > packing.MinX.
            // This lower bound (MinimumRemainingItemDy) can be raised by determining the remaining minY at newX at a cost of O(n) or by solving a subset sum problem.
            // The latter would yield SS_1 of [CJCM2008] if it was performed for [packing.MinX, container.Dx] instead of only [packing.MinX, packing.MinX + item.Dx].
            packing.MinimumRemainingItemDyAtCurrentX = packing.MinimumRemainingItemDy;

            if (packing.PlacedAreaVector[x] + item.Dy > container.Dy - packing.MinimumRemainingItemDyAtCurrentX)
            {
                lastY = currentY;
                deactivatedDy = container.Dy - packing.PlacedAreaVector[x];
            }
            else
            {
                // We set lastY to zero to be able to complete the loop but not re-enter this branch, in order to correctly update DeactivatedArea and PlacedAreaVector.
                newX = x;
                lastY = 0;
                deactivatedDy = item.Dy;
            }
        }

        packing.DeactivatedArea += deactivatedDy;
        packing.PlacedAreaVector[x] += deactivatedDy;
    }

    return newX;
}

bool TwoStepBranchingProcedure::PlaceAtCurrentX(Node& node, const Rectangle& item)
{
    PackingRelaxed2D& packing = *node.Packing;

    if (IsPlacementFeasible(packing, item))
    {
        packing.ItemCoordinatesX[node.ItemIdToPlace] = packing.MinX;

        size_t newCandidateX = UpdatePlacedAreaAndDetermineLoss(packing, item);
        ////size_t newX = UpdatePlacedAreaAndDetermineLossSS1(packing, item);

        packing.PlaceableItemsAtCurrentX.reset(0, (size_t)item.InternId + 1);

        if (packing.PlacedAreaVector[packing.MinX] == container.Dy)
        {
            packing.MinX = newCandidateX;

            // Must be reset because newX > packingMinX.
            // This lower bound can be raised by determining the remaining minY at newX at a cost of O(n) or by solving a subset sum problem in O(nC).
            packing.MinimumRemainingItemDyAtCurrentX = packing.MinimumRemainingItemDy;
            packing.PlaceableItemsAtCurrentX = ~packing.PlacedItems;
            packing.PlaceableItemsAtCurrentX.reset((size_t)item.InternId);
        }

        return true;
    }

    // Placement was infeasible.
    return false;
}

typename TwoStepBranchingProcedure::SearchTree::vertex_descriptor TwoStepBranchingProcedure::AddLeafNode(Node&& node)
{
    SearchTree::vertex_descriptor newNodeId = boost::add_vertex(this->tree);

    Node& newNode = this->tree[newNodeId];
    size_t remainingItemsToPlace = this->items.size() - node.PlacedItems().count();

    // TODO.Performance: this is helpful with directionalS and bidirectionalS. With bidirectionalS, edges are additionally (?) stored as members of tree.m_edges.
    this->tree.m_vertices[newNodeId].m_out_edges.reserve(remainingItemsToPlace + 1);

    newNode = std::move(node);

    return newNodeId;
}

#pragma endregion

/*
double MemoryUsedByCurrentProcess()
{
    PROCESS_MEMORY_COUNTERS_EX pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));

    SIZE_T physicalMemoryUsedByCurrentProcess = pmc.WorkingSetSize;

    return (double)physicalMemoryUsedByCurrentProcess;
}

void CheckMemoryUsage()
{
    // According to https://stackoverflow.com/a/64166/5587903.
    MEMORYSTATUSEX MemInfo;
    MemInfo.dwLength = sizeof(MEMORYSTATUSEX);
    GlobalMemoryStatusEx(&MemInfo);
    double ramUsage = (double)(MemInfo.ullTotalPhys - MemInfo.ullAvailPhys) / (double)MemInfo.ullTotalPhys;
    if (ramUsage > 0.95)
    {
        std::cout << "Memory at " + std::to_string(ramUsage) + "%: abort.\n";
        throw std::exception("System out of memory!");
    }
}
*/

}