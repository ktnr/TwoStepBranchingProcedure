#pragma once

#include "Model.h"

namespace tsbp
{
class PreprocessPacking2D
{
  public:
    std::vector<Rectangle> ProcessedItems;
    Bin PreprocessedContainer;

    PreprocessPacking2D(const std::vector<Rectangle>& items, const Bin& container, const InputParameters& inputParameters)
    : items(items),
      container(container),
      parameters(inputParameters) {}

    void Run();

  private:
    std::vector<Rectangle> items;
    Bin container;
    InputParameters parameters;

    void RemoveLargeItems(std::vector<Rectangle>& items, Bin& container);
    void FilterFrameConfigurations(std::vector<Rectangle>& items, Bin& container);

    /// Enlarge item sizes by finding feasible packings within the container similar to the v_1^* procedure as described in section 3.2. in CJCM07.
    ////void EnlargeItemsByOrthogonalPacking(std::vector<Rectangle>& items, Bin& container);
};
}