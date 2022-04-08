#include "OrthogonalPackingSolver2D.h"

#include "Preprocess.h"
#include "Serialization.h"

#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>

namespace tsbp
{

SearchStatus OrthogonalPackingSolver2D::Solve()
{
    PreprocessPacking2D preprocess(this->items, this->container, this->inputParameters);
    preprocess.Run();

    std::vector<Rectangle> processedItems = preprocess.ProcessedItems;
    Bin processedContainer = preprocess.PreprocessedContainer;

    switch (this->inputParameters.SolutionAlgorithm)
    {
    case SolutionAlgorithm::LeftmostActiveOnlyAlgorithm:
        this->branchAndBoundSearch = std::make_unique<LeftmostActiveOnly>(processedItems, processedContainer, this->inputParameters);
        break;
    case SolutionAlgorithm::TwoStepBranchingProcedure:
        this->branchAndBoundSearch = std::make_unique<TwoStepBranchingProcedure>(processedItems, processedContainer, this->inputParameters);
        break;
    default:
        break;
    }

    auto start = std::chrono::steady_clock::now();

    SearchStatus status = this->branchAndBoundSearch->Solve();

    this->branchAndBoundSearch->UpdateSolverStatistics(this->solutionStatistics);
    this->solution = this->branchAndBoundSearch->GetSolution();

    auto end = std::chrono::steady_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    this->solutionStatistics.TimeOverall = time;

    Log(solution, this->solutionStatistics);

    WriteFiles(status);

    return status;
}

void OrthogonalPackingSolver2D::WriteFiles(SearchStatus status)
{
    const std::string statisticsFileName = "SolutionStatistics";
    Serializer::WriteToJson<SolverStatistics>(this->solutionStatistics, std::string(this->inputParameters.OutputPath + statisticsFileName));

    if (status == SearchStatus::Feasible)
    {
        const std::string solutionFileName = "Solution";
        Serializer::WriteToJson<Packing2D>(this->solution, this->inputParameters.OutputPath + solutionFileName);
    }
}

void OrthogonalPackingSolver2D::Log(const Packing2D& solution, const SolverStatistics& statistics) const
{
    if (!(this->inputParameters.LMAOEnableLogging || this->inputParameters.TSBPEnableLogging))
    {
        return;
    }

    for (const Rectangle& rectangle : solution.Items)
    {
        std::cout << "Id = " << rectangle.InternId << ":(" << rectangle.X << ", " << rectangle.Y << "), (" << rectangle.X + rectangle.Dx << ", " << rectangle.Y + rectangle.Dy << ")\n";
    }

    std::cout << "     \t call count \t # nodes \t time (ms)\n";
    std::cout << "LMAO \t " 
        << std::setw(10) << statistics.CallCountLMAO << " \t "
        << std::setw(10) << statistics.NodeCountLMAO << " \t "
        << std::setw(10) << statistics.TimeLMAO << "\n";

    if (this->inputParameters.SolutionAlgorithm != SolutionAlgorithm::LeftmostActiveOnlyAlgorithm)
    {
        std::cout << "TSBP \t "
            << std::setw(10) << 1 << " \t "
            << std::setw(10) << statistics.NodeCountTSBP << " \t "
            << std::setw(10) << statistics.TimeTSBP << "\n";
    }

    std::cout << "Total\t "
        << std::setw(26) << (statistics.NodeCountLMAO + statistics.NodeCountTSBP) << " \t "
        << std::setw(10) << statistics.TimeOverall << "\n";
}

}