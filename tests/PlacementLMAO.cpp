#include "doctest.h"

#include "tsbp/Model.h"
#include "tsbp/BranchAndBound.h"

#include <vector>

using namespace tsbp;

TEST_CASE("Test feasible item placement: inverted L")
{
    Bin bin(0, 0, 4, 4, 0, 0);
    std::vector<tsbp::Rectangle> items;

    items.emplace_back(tsbp::Rectangle(0, 0, 1, 3, 0, 0));
    items.emplace_back(tsbp::Rectangle(0, 3, 2, 1, 1, 1));
    
    Packing2D packing;
    packing.Initialize(items, bin);

    for (auto& item : items)
    {
        int itemId = item.InternId;

        packing.AddItem(std::move(item), itemId); // Caution.
        packing.PlacedItems.set(itemId);
    }

    // Build inverted L according to the dimensions and placements of "items".
    packing.PlacedAreaVector[0] = 1;
    packing.PlacedAreaVector[1] = 1;
    packing.PlacedAreaVector[2] = 1;
    packing.PlacedAreaVector[3] = 2;

    tsbp::Rectangle itemToPlace(0, 0, 1, 1, 2, 2);
    std::vector<tsbp::Rectangle> allItems = items;
    allItems.emplace_back(itemToPlace);

    InputParameters inputParameters;
    
    SUBCASE("Place at (0, 0)")
    {
        packing.ActiveX = 0;
        packing.ActiveY = 0;
        packing.MaxX = 2;
        packing.AbscissaMaxX = 1;

        LeftmostActiveOnly lmao(allItems, bin, inputParameters);
        lmao.Preprocess();
        bool isFeasible = lmao.IsPlacementFeasible(packing, itemToPlace);

        CHECK(!isFeasible);
    }
    
    SUBCASE("Place at (1, 0)")
    {
        packing.ActiveX = 1;
        packing.ActiveY = 0;
        packing.MaxX = 2;
        packing.AbscissaMaxX = 1;

        LeftmostActiveOnly lmao(allItems, bin, inputParameters);
        lmao.Preprocess();
        bool isFeasible = lmao.IsPlacementFeasible(packing, itemToPlace);

        CHECK(isFeasible);
    }
    
    SUBCASE("Place at (1, 2)")
    {
        packing.ActiveX = 1;
        packing.ActiveY = 2;
        packing.MaxX = 2;
        packing.AbscissaMaxX = 1;

        LeftmostActiveOnly lmao(allItems, bin, inputParameters);
        lmao.Preprocess();
        bool isFeasible = lmao.IsPlacementFeasible(packing, itemToPlace);

        CHECK(isFeasible);
    }
    
    SUBCASE("Place at (1, 3)")
    {
        packing.ActiveX = 1;
        packing.ActiveY = 3;
        packing.MaxX = 2;
        packing.AbscissaMaxX = 1;

        LeftmostActiveOnly lmao(allItems, bin, inputParameters);
        lmao.Preprocess();
        bool isFeasible = lmao.IsPlacementFeasible(packing, itemToPlace);

        CHECK(!isFeasible);
    }
    
    SUBCASE("Place at (2, 3)")
    {
        packing.ActiveX = 2;
        packing.ActiveY = 3;
        packing.MaxX = 2;
        packing.AbscissaMaxX = 1;

        LeftmostActiveOnly lmao(allItems, bin, inputParameters);
        lmao.Preprocess();
        bool isFeasible = lmao.IsPlacementFeasible(packing, itemToPlace);

        CHECK(isFeasible);
    }

    SUBCASE("Place at (3, 0)")
    {
        packing.ActiveX = 3;
        packing.ActiveY = 0;
        packing.MaxX = 2;
        packing.AbscissaMaxX = 1;

        LeftmostActiveOnly lmao(allItems, bin, inputParameters);
        lmao.Preprocess();
        bool isFeasible = lmao.IsPlacementFeasible(packing, itemToPlace);

        CHECK(isFeasible);
    }

    SUBCASE("Place at (4, 0)")
    {
        packing.ActiveX = 4;
        packing.ActiveY = 0;
        packing.MaxX = 2;
        packing.AbscissaMaxX = 1;

        LeftmostActiveOnly lmao(allItems, bin, inputParameters);
        lmao.Preprocess();
        bool isFeasible = lmao.IsPlacementFeasible(packing, itemToPlace);

        CHECK(!isFeasible);
    }

    SUBCASE("Place at (3, 3)")
    {
        packing.ActiveX = 3;
        packing.ActiveY = 3;
        packing.MaxX = 2;
        packing.AbscissaMaxX = 1;

        LeftmostActiveOnly lmao(allItems, bin, inputParameters);
        lmao.Preprocess();
        bool isFeasible = lmao.IsPlacementFeasible(packing, itemToPlace);

        CHECK(isFeasible);
    }

    SUBCASE("Place at (3, 4)")
    {
        packing.ActiveX = 3;
        packing.ActiveY = 4;
        packing.MaxX = 2;
        packing.AbscissaMaxX = 1;

        LeftmostActiveOnly lmao(allItems, bin, inputParameters);
        lmao.Preprocess();
        bool isFeasible = lmao.IsPlacementFeasible(packing, itemToPlace);

        CHECK(!isFeasible);
    }
    
    SUBCASE("Place at (4, 4)")
    {
        packing.ActiveX = 4;
        packing.ActiveY = 4;
        packing.MaxX = 2;
        packing.AbscissaMaxX = 1;

        LeftmostActiveOnly lmao(allItems, bin, inputParameters);
        lmao.Preprocess();
        bool isFeasible = lmao.IsPlacementFeasible(packing, itemToPlace);

        CHECK(!isFeasible);
    }
}