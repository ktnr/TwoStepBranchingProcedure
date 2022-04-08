#pragma once

#include "nlohmann/json.hpp"

#include <vector>

namespace tsbp
{

struct Rectangle
{
  public:
    int X = 0;
    int Y = 0;

    int Dx = 0;
    int Dy = 0;

    int InternId = -1;
    int ExternId = -1;

    double Area = 0.0;

    Rectangle() = default;
    Rectangle(const Rectangle& rectangle) = default;
    Rectangle(int x, int y, int dx, int dy, int internId, int externId)
    : X(x),
      Y(y),
      Dx(dx),
      Dy(dy),
      InternId(internId),
      ExternId(externId) 
      {
          Area = dx * dy;
      }

  public:
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(Rectangle,
                                   X,
                                   Y,
                                   Dx,
                                   Dy,
                                   Area,
                                   ExternId)
};

struct Bin : public Rectangle
{
    Bin() = default;
    Bin(const Bin& bin) = default;
    Bin(int x, int y, int dx, int dy, int internId, int externId)
    : Rectangle(x, y, dx, dy, internId, externId) {}

  public:
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(Bin,
                                   Dx,
                                   Dy,
                                   Area)
};

struct CompareLexicographicDxDy
{
  public:
    CompareLexicographicDxDy(std::vector<Rectangle>* itemsToPlace) : items(itemsToPlace) {}
    bool operator()(int i, int j) const
    {
        const Rectangle& itemA = (*this->items)[i];
        const Rectangle& itemB = (*this->items)[j];

        return std::tie(itemA.Dx, itemA.Dy) > std::tie(itemB.Dx, itemB.Dy);
    }

  private:
    std::vector<Rectangle>* items;
};

class ProblemInstance
{
  public:
    std::string Name;
    Bin Container;
    std::vector<Rectangle> Items;

    ProblemInstance() = default;

    ProblemInstance(std::string& name, Bin& container, std::vector<Rectangle>& items)
    : Name(name),
      Container(container),
      Items(items)
      {};

  public:
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(ProblemInstance,
                                   Name,
                                   Container,
                                   Items)
};

enum class SolutionAlgorithm
{
    None = 0,
    LeftmostActiveOnlyAlgorithm,
    TwoStepBranchingProcedure
};

struct InputParameters
{
  public:
    SolutionAlgorithm SolutionAlgorithm = SolutionAlgorithm::TwoStepBranchingProcedure;

    std::string OutputPath = "";

    // LMAO - B&B
    bool LMAOEnableLogging = false;
    int LMAOTimeLimit = 43200;
    int LMAOThreads = 1;

    std::string LMAOLogFileName = "LMAO";

    // TSBP - B&B
    bool TSBPEnableLogging = false;
    int TSBPTimeLimit = 43200;
    int TSBPThreads = 1;

    std::string TSBPLogFileName = "TSBP";

  public:
    // What we would really want is to serialize missing members to their default values. This was the behavior once, but changed.
    // Monitor https://github.com/nlohmann/json/pull/2819, which will provide the macro: NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(InputParameters,
                                   SolutionAlgorithm,
                                   // LMAO
                                   LMAOEnableLogging,
                                   LMAOTimeLimit,
                                   LMAOThreads,
                                   LMAOLogFileName,
                                   // TSBP
                                   TSBPEnableLogging,
                                   TSBPTimeLimit,
                                   TSBPThreads,
                                   TSBPLogFileName)
};

}