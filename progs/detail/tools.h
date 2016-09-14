#pragma once
#include "TGraph.h"

namespace ant
{
namespace progs
{

namespace tools
{


std::size_t FillGraph(TGraph* graph, const double x, const double y)
{
    graph->SetPoint(graph->GetN(),x,y);
    return graph->GetN();
}



}
}
}
