#include "TH2DrawTool.h"
#include "TH2Poly.h"

using namespace std;
using namespace ant;

TH2DrawTool::TH2DrawTool(TH2Poly *target)
{
    hist = target;
}

void TH2DrawTool::Draw(const matrixstack::Vector &vector)
{
    matrixstack::Vector t = Transform(vector);

    x.push_back(t(0));
    y.push_back(t(1));

}

void TH2DrawTool::Draw(const TH2DrawTool::point_list &points)
{
    x.reserve(x.size()+points.size());
    y.reserve(y.size()+points.size());

    point_list::const_iterator i;
    for(i=points.begin(); i!= points.end(); ++i) {
        Draw( *i );
    }
}

Int_t TH2DrawTool::FinishShape()
{
    Int_t bin = hist->AddBin(x.size(),&x[0],&y[0]);
    ResetShape();
    return bin;
}

void TH2DrawTool::ResetShape()
{
    x.clear();
    y.clear();
}
