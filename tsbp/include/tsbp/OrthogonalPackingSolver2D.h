#pragma once

#include "Model.h"
#include "BranchAndBound.h"

#include <memory>

namespace tsbp
{

class OrthogonalPackingSolver2D
{
public:
    OrthogonalPackingSolver2D(const std::vector<Rectangle>& items, const Bin& container, const InputParameters& inputParameters)
        :
        items(items),
        container(container),
        inputParameters(inputParameters) {};

    SearchStatus Solve();

private:
    ////ContainerLoadingInstance instance;
    std::vector<Rectangle> items;
    Bin container;
    InputParameters inputParameters;
    SolverStatistics solutionStatistics;
    Packing2D solution;

    // Must be unique_ptr<T>. Otherwise a use of deleted function error will be thrown (trying to reference deleted function operator=(const Node&)).
    std::unique_ptr<IBranchAndBoundSolver> branchAndBoundSearch;

    void Log(const Packing2D& solution, const SolverStatistics& statistics) const;
    void WriteFiles(SearchStatus status);
};

}
