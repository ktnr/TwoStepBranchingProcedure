#include "tsbp/Model.h"
#include "tsbp/Serialization.h"
#include "tsbp/OrthogonalPackingSolver2D.h"

#include "CLI11/CLI11.hpp"

#include <chrono>
#include <iomanip>
#include <iostream>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace tsbp;

void RunTwoDimensionalOrthogonalPackingSolver(const ProblemInstance& problemInstance, InputParameters parameters)
{
    OrthogonalPackingSolver2D exactAlgorithm(problemInstance.Items, problemInstance.Container, parameters);
    exactAlgorithm.Solve();
}

void RunSolver(std::string& inputFilePath, std::string& filename, std::string& parameterFile, std::string& outdir, bool enableTimeSuffix)
{
    InputParameters inputParameters;

    if (parameterFile != "")
    {
        inputParameters = Serializer::ReadFromJson<InputParameters>(parameterFile);
    }
    
    ProblemInstance problemInstance = Serializer::ReadFromJson<ProblemInstance>(inputFilePath + filename);

    // https://stackoverflow.com/questions/16357999/current-date-and-time-as-string/16358111
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d_%H-%M-%S");
    std::string dateTimeString = oss.str();

    std::string outputPath;

    if (enableTimeSuffix)
    {
        outputPath = outdir + dateTimeString + "/";
    }

    // https://stackoverflow.com/a/37524002
    if (!std::filesystem::is_directory(outdir) || !std::filesystem::exists(outdir))
    {
        std::filesystem::create_directory(outdir);
    }

    if (!std::filesystem::is_directory(outputPath) || !std::filesystem::exists(outputPath))
    {
        std::filesystem::create_directory(outputPath);
    }

    inputParameters.OutputPath = outputPath;

    RunTwoDimensionalOrthogonalPackingSolver(problemInstance, inputParameters);
}

int main(int argc, char** argv)
{
    // For example, call with: -i "data/input/CJCM08/" -f "E00N10.json" -o "data/output/" -p "data/2D-OPP/Converted/Parameters.json"
    CLI::App app;

    std::string inputFilePath = "default";
    std::string filename = "default";
    std::string outdir = "default";
    std::string parameterFile = "";
    bool enableTimeSuffix = true;

    app.add_option("-i,--inputdir", inputFilePath, "The directory where the input file -f resides")->required();
    app.add_option("-f,--file", filename, "The input file name")->required();
    app.add_option("-o,--outdir", outdir, "The output directory")->required();
    app.add_option("-p,--param", parameterFile, "The .json parameter full file path");
    app.add_option("-t,--timeSuffix", enableTimeSuffix, "If time should be appended to the output path (1=true, 0=false)");

    CLI11_PARSE(app, argc, argv);

    std::string inputFilePathDelimiter = inputFilePath.substr(inputFilePath.size() - 1, inputFilePath.size());
    std::string inputFileSuffix = filename.substr(filename.size() - 5, filename.size());
    std::string outdirDelimiter = outdir.substr(outdir.size() - 1, outdir.size());

    if (inputFilePathDelimiter != "/" && inputFilePathDelimiter != "\\")
    {
        throw CLI::ConversionError("-i directory path delimiter does neither math '/' nor '\\'");
    }

    if (inputFileSuffix != ".json")
    {
        throw CLI::ConversionError("-f does not have an .json suffix");
    }

    if (outdirDelimiter != "/" && outdirDelimiter != "\\")
    {
        throw CLI::ConversionError("-o directory path delimiter does neither math '/' nor '\\'");
    }
    try
    {
        RunSolver(inputFilePath, filename, parameterFile, outdir, enableTimeSuffix);
        return EXIT_SUCCESS;
    }
    catch (std::exception& e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return 0;
}