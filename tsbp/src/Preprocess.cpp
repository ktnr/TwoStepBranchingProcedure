#include "Preprocess.h"

#include "Serialization.h"

#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <iostream>

namespace tsbp
{

void PreprocessPacking2D::Run()
{
    std::vector<Rectangle> newItems = this->items;
    Bin newContainer = this->container;

    RemoveLargeItems(newItems, newContainer);
    FilterFrameConfigurations(newItems, newContainer);
    ////EnlargeItemsByOrthogonalPacking(newItems, newContainer);

    // Reconstruct to correctly initialize area and volume.
    newContainer = Bin(0, 0, newContainer.Dx, newContainer.Dy, -1, -1);

    std::cout << "Preprocess removed " << this->items.size() - newItems.size() << " items "
              << "and reduced container area by " << this->container.Area - newContainer.Area << ".\n";

    PreprocessedContainer = newContainer;
    ProcessedItems = newItems;
}

void PreprocessPacking2D::RemoveLargeItems(std::vector<Rectangle>& items, Bin& container)
{
    bool removedItem = true;
    while (removedItem)
    {
        removedItem = false;
        std::unordered_set<int> itemsToRemove;
        for (size_t i = 0; i < items.size(); i++)
        {
            const Rectangle& item = items[i];
            if (item.Dx == container.Dx)
            {
                container.Dy -= item.Dy;

                itemsToRemove.emplace(item.InternId);
                removedItem = true;

                break;
            }

            if (item.Dy == container.Dy)
            {
                container.Dx -= item.Dx;

                itemsToRemove.emplace(item.InternId);
                removedItem = true;

                break;
            }
        }

        items.erase(
            std::remove_if(
                items.begin(),
                items.end(),
                [&itemsToRemove](Rectangle& item)
                { return itemsToRemove.contains(item.InternId); }),
            items.end());
    }
}

void PreprocessPacking2D::FilterFrameConfigurations(std::vector<Rectangle>& items, Bin& container)
{
    for (size_t i = 0; i < items.size(); i++)
    {
        const Rectangle& itemI = items[i];
        for (size_t j = 0; j < items.size(); j++)
        {
            if (i == j)
            {
                continue;
            }

            const Rectangle& itemJ = items[j];
            for (size_t k = 0; k < items.size(); k++)
            {
                if (k == i || k == j)
                {
                    continue;
                }

                const Rectangle& itemK = items[k];
                for (size_t l = 0; l < items.size(); l++)
                {
                    if (l == i || l == j || l == k)
                    {
                        continue;
                    }

                    const Rectangle& itemL = items[l];

                    if (itemI.Dx + itemL.Dx == itemJ.Dx + itemK.Dx
                        && itemI.Dx + itemL.Dx == container.Dx
                        && itemI.Dy + itemK.Dy == itemJ.Dy + itemL.Dy
                        && itemI.Dy + itemK.Dy == container.Dy)
                    {
                        container.Dx = container.Dx - itemI.Dx - itemJ.Dx;
                        container.Dy = container.Dy - itemK.Dy - itemL.Dy;

                        std::unordered_set<int> itemsToRemove{itemI.InternId, itemJ.InternId, itemK.InternId, itemL.InternId};
                        items.erase(
                            std::remove_if(
                                items.begin(),
                                items.end(),
                                [&itemsToRemove](Rectangle& item)
                                { return itemsToRemove.contains(item.InternId); }),
                            items.end());

                        return;
                    }
                }
            }
        }
    }
}

/*
void PreprocessPacking2D::EnlargeItemsByOrthogonalPacking(std::vector<Rectangle>& items, Bin& container)
{
    std::unordered_map<int, std::vector<Rectangle>> feasibleThresholdToItemCount;
    int largestItemThreshold = -1;
    int largestItemCount = 0;

    for (size_t p = 1; p <= std::floor<int>(container.Dy / 2); p++)
    {
        std::vector<Rectangle> largeItemSubset;
        std::vector<Rectangle> smallItemSubset;
        std::vector<Rectangle> itemSubset;

        int largeItemDx = 0;
        int largestSmallItemDx = 0;

        for (size_t i = 0; i < items.size(); i++)
        {
            Rectangle& item = items[i];
            if (item.Dy >= container.Dy - p)
            {
                largeItemDx += item.Dx;

                itemSubset.emplace_back(item);
                largeItemSubset.emplace_back(item);
                continue;
            }

            if (item.Dy <= p)
            {
                largestSmallItemDx = std::max<int>(largestSmallItemDx, item.Dx);

                itemSubset.emplace_back(item);
                smallItemSubset.emplace_back(item);
                continue;
            }
        }

        if (largeItemSubset.empty() || smallItemSubset.empty() || itemSubset.empty() || largestSmallItemDx > largeItemDx)
        {
            continue;
        }

        Bin newContainer(largeItemDx, container.Dy, container.Dz, container.WeightLimit);

        ////BranchAndBound::OrthogonalPackingSolver2D exactAlgorithm(itemSubset, newContainer, this->parameters);
        ////auto status = exactAlgorithm.SolveWithoutPreprocessing();

        std::cout << "Preprocess CP call with " << itemSubset.size() << "/" << items.size() << "items\n";

        PackingCP2D packingCP;
        auto cpStatus = packingCP.Solve(newContainer, itemSubset);

        using namespace operations_research::sat;

        if (cpStatus.status() == CpSolverStatus::FEASIBLE || cpStatus.status() == CpSolverStatus::OPTIMAL)
        {
            if (itemSubset.size() > largestItemCount)
            {
                largestItemCount = itemSubset.size();
                largestItemThreshold = p;

                feasibleThresholdToItemCount.emplace(p, itemSubset);
            }
        }
    }

    if (largestItemThreshold == -1)
    {
        return;
    }

    auto maximalItemSet = feasibleThresholdToItemCount[largestItemThreshold];
    std::unordered_set<int> itemsToRemove;
    itemsToRemove.reserve(maximalItemSet.size());

    int reducedContainerDx = 0;
    for (auto& item :maximalItemSet)
    {
        itemsToRemove.emplace(item.InternId);
        reducedContainerDx += item.Dx;
    }

    container.Dx -= reducedContainerDx;

    items.erase(
    std::remove_if(
        items.begin(),
        items.end(),
        [&itemsToRemove](Rectangle& item)
        { return itemsToRemove.contains(item.InternId); }),
    items.end());
}
*/
}