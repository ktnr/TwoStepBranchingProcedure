#pragma once

#include "Model.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/exterior_property.hpp>
#include <boost/graph/graph_traits.hpp>

#include <boost/property_map/property_map.hpp>

#include <boost/dynamic_bitset/dynamic_bitset.hpp>

namespace tsbp
{

////static void CheckMemoryUsage();
////static double MemoryUsedByCurrentProcess();

enum class SearchStatus
{
    None = 0,
    InProgress,
    Feasible,
    Infeasible,
    Abort,
    Cancelled
};

struct SolverStatistics
{
  public:
    double TimeOverall = 0.0;
    double TimeLMAO = 0.0;
    double TimeTSBP = 0.0;
    double TimeMemoryDeallocation = 0.0;
    double MemoryUsage = 0.0;

    int64_t NodeCountLMAO = 0;
    int64_t NodeCountTSBP = 0;
    int64_t CallCountLMAO = 0;

  public:
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(SolverStatistics,
                                   TimeOverall,
                                   TimeLMAO,
                                   TimeTSBP,
                                   TimeMemoryDeallocation,
                                   MemoryUsage,
                                   NodeCountLMAO,
                                   NodeCountTSBP,
                                   CallCountLMAO)
};

struct Packing2D
{
  public:
    int ActiveX = 0;
    int ActiveY = 0;

    size_t MinX = 0;
    int MaxX = 0;
    int AbscissaMaxX = 0;

    boost::dynamic_bitset<> PlacedItems;
    std::vector<Rectangle> Items;

    // Worse performance than the solution with std::vector<size_t>.
    std::vector<size_t> PlacedAreaVector;
    std::vector<size_t> DeactivatedAreaVector;

    double DeactivatedArea = 0.0;
    double RemainingItemAreaToPlace;

    size_t ReducedItemInfeasiblePlacementPoints = std::numeric_limits<size_t>::max();

    bool IsReducedItemPlaced = false;
    bool IsFixedItemPlaced = false;

    void Initialize(const std::vector<Rectangle>& items, const Bin& bin);
    void AddItem(Rectangle&& rectangle, int itemId);

    double GetDeactivatedArea() const { return DeactivatedArea; }

  public:
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(Packing2D,
                                   Items,
                                   PlacedAreaVector,
                                   DeactivatedAreaVector,
                                   DeactivatedArea)
};

struct BaseNode
{
  public:
    enum class Type
    {
        UsePlacement = 0,
        DeactivatePlacement
    };

    enum class Status
    {
        None = 0,
        Root,
        NotEvaluated,
        FeasiblePlacement,
        InfeasiblePlacement,
        DeactivatedPlacement,
        InfeasibleSequence
    };

    // TODO.Performance: move member into packing class to reduce the memory footprint after unique_ptr.reset() when nodes have been deactivated.
    ////size_t Id;
    ////size_t Depth;
    int ItemIdToPlace;
    Type NodeType;
    Status NodeStatus;
    ////size_t PredecessorId;
};

class IBranchAndBoundSolver
{
  public:
    virtual ~IBranchAndBoundSolver() = default;
    
    virtual SearchStatus Solve() = 0;
    virtual SearchStatus GetSearchStatus() const = 0;

    virtual void UpdateSolverStatistics(SolverStatistics& statistics) const = 0;
    virtual Packing2D GetSolution() const = 0;
};

/// The lower level B&B procedure from Clautiaux, F., Carlier, J., & Moukrim, A. (2007). A new exact method for the two-dimensional orthogonal packing problem.
/// European Journal of Operational Research, 183(3), 1196-1211.
class LeftmostActiveOnly : public IBranchAndBoundSolver
{
    struct Node : BaseNode
    {
      public:
        std::unique_ptr<Packing2D> Packing;

        const boost::dynamic_bitset<>& PlacedItems() const { return Packing->PlacedItems; }

        bool AllItemsPlaced(const std::vector<Rectangle>& items) const;
    };

    struct Statistics
    {
      public:
        double TimeOverall = 0.0;
        double TimeMemoryDeallocation = 0.0;
        double MemoryUsage = 0.0;

        int NodeCount = 0;
    };

    // bidirectionalS is needed to iterate in-edges.
    ////using SearchTree = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, Node>;
    using SearchTree = boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, Node>;

  public:
    SearchStatus Solve() override;
    SearchStatus SolveSequential();
    SearchStatus SolveParallelTaskflow();
    SearchStatus SolveParallelNative();

    SearchStatus GetSearchStatus() const override { return this->searchStatus; };
    void UpdateSolverStatistics(SolverStatistics& statistics) const override;

    const Statistics& GetStatistics() const { return this->statistics; };
    Packing2D GetSolution() const override;

    LeftmostActiveOnly() = default;
    LeftmostActiveOnly(const InputParameters& inputParameters) : parameters(inputParameters) {}

    // Assigns new InternIds to items in the order in which they appear in items.
    LeftmostActiveOnly(const std::vector<Rectangle>& items, const Bin& container, const InputParameters& inputParameters)
    : container(container), parameters(inputParameters)
    {
        this->tree.m_vertices.reserve(items.size() + 1);
        this->items.reserve(items.size());

        // The intern IDs of the items must correpond to the indices of items.
        for (size_t i = 0; i < items.size(); i++)
        {
            const Rectangle& item = items[i];

            // InternId is the index of items, externId is the internId of the passed items.
            const Rectangle& addedItem = this->items.emplace_back(0, 0, item.Dx, item.Dy, i, item.ExternId);
        }
    }

    /// Requires that the each item.InternId corresponds to its index in items and is consistent with itemCoordinatesX.
    void AddItems(const std::vector<Rectangle>& items, const std::optional<int>& fixedItemAtOrigin, const std::vector<int>& itemCoordinatesX);
    void AddContainer(const Bin& container);
    void AddCancellationToken(std::atomic<bool>* cancellationToken);

    size_t GetTreeSize() const { return this->tree.m_vertices.size(); }
    const std::optional<int>& GetFixedItem() const { return this->fixedItemAtOrigin; }
    void SetBranchingOrder(const std::vector<int>& branchingOrder) { this->branchingOrder = branchingOrder; }

  private:
    int exploredNodeDisplayInterval = 100000000;
    int feasibleLeafNodeId = -1;

    SearchStatus searchStatus = SearchStatus::None;

    size_t domainReducedItemIndex = 0;
    size_t reducedItemDomainX = std::numeric_limits<size_t>::max();
    size_t reducedItemDomainY = std::numeric_limits<size_t>::max();
    size_t reducedItemFeasiblePlacementPoints = std::numeric_limits<size_t>::max();

    size_t maxPlacedItemsInfeasibleNodeId = 0;
    size_t exploredNodeDisplayCount = 0;

    bool areItemCoordinatesFixedX = false;
    std::vector<int> fixedItemCoordinatesX;
    std::vector<int> branchingOrder;

    // Should this be https://en.cppreference.com/w/cpp/atomic/atomic_flag?
    std::atomic<bool>* cancellationToken = nullptr;

    std::optional<int> fixedItemAtOrigin = std::nullopt;

    SearchTree tree;
    Bin container;
    InputParameters parameters;
    Statistics statistics;

    Packing2D solution;

    std::vector<Rectangle> items;
    std::vector<size_t> activeNodes;
    std::unordered_set<size_t> nodesToDeactivate;

    std::vector<boost::dynamic_bitset<>> itemSpecificPlacementPointsX;
    std::vector<boost::dynamic_bitset<>> itemSpecificPlacementPointsY;

    size_t InitializeSearchTree();

    size_t FindMinimumItemDx();
    size_t FindMaximumItemDx();

    /// Reduced rotational and mirror symmetries by domain reduction according to Soh, T., Inoue, K., Tamura, N., Banbara, M., & Nabeshima, H. (2010). A SAT-based method for solving the two-dimensional strip packing problem.
    /// Fundamenta Informaticae, 102(3-4), 467-487.
    void DomainReduction(size_t domainReducedItemIndex);

    void DomainReductionXY(const size_t domainReducedItemIndex);
    void DomainReductionY(const size_t domainReducedItemIndex);

    typename SearchTree::vertex_descriptor AddLeafNode(Node&& node);

    std::optional<size_t> Branch(size_t nodeId);
    std::optional<size_t> Backtrack(size_t currentNode);

    void EvaluateLeaf(Node& node, size_t nodeId);

    void DeactivatePlacement(Node& node, size_t nodeId);

    bool PlaceBottomLeftPlacement(Node& node, const Rectangle& item);

    std::tuple<size_t, size_t> FindNewBottomLeft(Packing2D& packing);
    void AddSuccessfulPlacementDummy(Packing2D& packing);

    bool IsFixedSequence(const Node& node, size_t itemId);

    bool IsCancelled() const;
    void SignalCancelling();

    /// These methods should actually be private, but they are used in tests.
    public:
    void Preprocess();
    bool IsPlacementFeasible(Packing2D& packing, const Rectangle& itemToPlace);
};

struct PackingRelaxed2D
{
  public:
    size_t MinX = 0;
    int MinimumRemainingItemDyAtCurrentX = std::numeric_limits<int>::max();
    int MinimumRemainingItemDy = std::numeric_limits<int>::max();
    int MaximumRemainingItemDx = 0;
    boost::dynamic_bitset<> PlacedItems;
    boost::dynamic_bitset<> PlaceableItemsAtCurrentX;
    std::vector<int> ItemCoordinatesX;
    ////std::vector<int> ItemSequence;

    /// Corresponds to K(x) in Clautiaux, Carlier, Moukrim (2007).
    std::vector<size_t> PlacedAreaVector;
    double DeactivatedArea = 0.0;

    double RemainingItemAreaToPlace;

    bool IsReducedItemPlaced = false;
    bool IsFixedItemPlaced = false;

    void Initialize(const std::vector<Rectangle>& items, const Bin& container);
    ////double GetDeactivatedArea() const { return boost::geometry::area(DeactivatedArea); }
    double GetDeactivatedArea() const { return DeactivatedArea; }
};

class TwoStepBranchingProcedure : public IBranchAndBoundSolver
{
    struct Node : BaseNode
    {
      public:
        ////size_t Id;
        ////size_t Depth;
        std::unique_ptr<PackingRelaxed2D> Packing;

        const boost::dynamic_bitset<>& PlacedItems() const { return Packing->PlacedItems; }

        bool AllItemsPlaced(const std::vector<Rectangle>& items) const;
    };

    // bidirectionalS is needed to iterate in-edges (?).
    ////using SearchTree = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, Node>;
    using SearchTree = boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, Node>;

  public:
    SearchStatus Solve() override;
    SearchStatus SolveSequential();
    SearchStatus SolveParallelTaskflow();
    SearchStatus SolveParallelNative();

    SearchStatus GetSearchStatus() const override { return this->searchStatus; };
    void UpdateSolverStatistics(SolverStatistics& statistics) const override;

    const SolverStatistics& GetStatistics() const { return this->statistics; };
    Packing2D GetSolution() const override;

    TwoStepBranchingProcedure() = default;
    TwoStepBranchingProcedure(const InputParameters& inputParameters) : parameters(inputParameters) {}
    TwoStepBranchingProcedure(const std::vector<Rectangle>& items, const Bin& container, const InputParameters& inputParameters)
    : container(container), parameters(inputParameters)
    {
        this->tree.m_vertices.reserve(items.size() + 1);
        this->items.reserve(items.size());

        // The intern IDs of the items must correpond to the indices of items.
        for (size_t i = 0; i < items.size(); i++)
        {
            const Rectangle& item = items[i];

            // InternId is the index of items, externId is the internId of the passed items.
            const Rectangle& addedItem = this->items.emplace_back(0, 0, item.Dx, item.Dy, i, item.ExternId);
        }
    }

    /// Requires that the each item.InternId corresponds to its index in item.
    void AddItems(const std::vector<Rectangle>& items, int fixedItemAtOrigin);
    void AddContainer(const Bin& container);
    void AddCancellationToken(std::atomic<bool>* cancellationToken);

    size_t GetTreeSize() const { return this->tree.m_vertices.size(); }
    const std::optional<int>& GetFixedItem() const { return this->fixedItemAtOrigin; }
    void SetBranchingOrder(const std::vector<int>& branchingOrder) { this->branchingOrder = branchingOrder; }

  private:
    int exploredNodeDisplayInterval = 100000000;
    int feasibleLeafNodeId = -1;
    size_t maxPlacedItemsInfeasibleNodeId = 0;
    SearchStatus searchStatus = SearchStatus::None;
    size_t exploredNodeDisplayCount = 0;

    size_t domainReducedItemIndex = 0;
    size_t reducedItemDomainX = std::numeric_limits<size_t>::max();

    std::atomic<bool>* cancellationToken = nullptr;
    std::optional<int> fixedItemAtOrigin = std::nullopt;
    std::vector<int> branchingOrder;

    SearchTree tree;
    Bin container;
    InputParameters parameters;
    SolverStatistics statistics;

    Packing2D solution;

    std::vector<Rectangle> items;
    std::vector<size_t> activeNodes;
    std::unordered_set<size_t> nodesToDeactivate;

    std::vector<boost::dynamic_bitset<>> itemSpecificPlacementPointsX;

    void Preprocess();

    size_t InitializeSearchTree();

    size_t FindMaximumItemDx();

    /// Reduced rotational and mirror symmetries by domain reduction according to Soh, T., Inoue, K., Tamura, N., Banbara, M., & Nabeshima, H. (2010). A SAT-based method for solving the two-dimensional strip packing problem.
    /// Fundamenta Informaticae, 102(3-4), 467-487.
    void DomainReductionX(size_t domainReducedItemIndex);

    bool PlaceAtCurrentX(Node& node, const Rectangle& item);
    size_t UpdatePlacedAreaAndDetermineLoss(PackingRelaxed2D& packing, const Rectangle& item);

    typename SearchTree::vertex_descriptor AddLeafNode(Node&& node);

    std::optional<size_t> Branch(size_t nodeId);
    std::optional<size_t> Backtrack(size_t currentNode);

    void EvaluateLeaf(Node& node, size_t nodeId);

    void DeactivatePlacement(Node& node, size_t nodeId);

    bool IsPlacementFeasible(PackingRelaxed2D& packing, const Rectangle& itemToPlace);

    bool SolveSubproblem(Node& newNode, const size_t newNodeId);

    bool IsCancelled() const;
    void SignalCancelling();

    bool IsFixedSequence(const Node& node, size_t itemId);
};

}