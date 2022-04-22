#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "tsbp/Model.h"
#include "tsbp/Preprocess.h"
#include "tsbp/BranchAndBound.h"

#include <vector>

using namespace tsbp;

////#define private public
// Use instead: https://stackoverflow.com/questions/27778908/define-private-to-public-in-c

TEST_CASE("Test frame configuration") 
{    
    using namespace tsbp;

    InputParameters inputParameters;
    Bin bin(0, 0, 4, 4, 0, 0);
    std::vector<tsbp::Rectangle> items;
    items.reserve(4);

    items.emplace_back(tsbp::Rectangle(0, 0, 1, 3, 0, 0));
    items.emplace_back(tsbp::Rectangle(0, 0, 3, 1, 1, 1));
    items.emplace_back(tsbp::Rectangle(0, 0, 1, 3, 2, 2));
    items.emplace_back(tsbp::Rectangle(0, 0, 3, 1, 3, 3));

    tsbp::PreprocessPacking2D preprocess(items, bin, inputParameters);
    preprocess.Run();

    CHECK(preprocess.ProcessedItems.size() == 0);
    CHECK(preprocess.PreprocessedContainer.Dx == 2);
    CHECK(preprocess.PreprocessedContainer.Dy == 2);
}