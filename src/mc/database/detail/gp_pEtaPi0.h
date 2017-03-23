#include "ProductionDataBase.h"

using namespace ant;
using namespace ant::mc::data;



// In EPT-range: fit

// fit to cb/cbElsa data: f(x) = p0 + p1 * x + p1 * x^2
//  Minimizer is Linear
//  Chi2                      =      1.62506
//  NDf                       =           17
//  p0                        =     -33.1171   +/-   18.2706
//  p1                        =      46.9949   +/-   25.0073
//  p2                        =     -14.8619   +/-   8.5167
// IN GeV!!!!!!!!!!!!!!!!
class gp_pEtaPi0_Data
{
private:
    static constexpr double startFit =   1335.0;
    static constexpr double   endFit =   1600.0;
    static constexpr double      GeV =   1/1000.0;

    static constexpr double  p0 = - 33.1171;
    static constexpr double  p1 =   46.9949;
    static constexpr double  p2 = - 14.8619;

    static double fit(const double e)
    {
        return  p0
              + p1 * e
              + p2 * e * e;
    }
    static std::function<double(const double)> interpolated ;


public:
    static double Get(const double e)
    {
        const auto sigma = (e < startFit)? interpolated(e) : fit(e*GeV);
        return (e < endFit)? sigma : fit(endFit*GeV);
    }
};

std::function<double(const double)> gp_pEtaPi0_Data::interpolated =
        ProductionDataBase::MakeInterPolator(
{

                { 931,   0.0},  // threshold
                // cb/cbElsa data
                {937,   0.0010 },
                {949,   0.0020 },
                {961,   0.0060 },
                {972,   0.012 },
                {983,   0.025 },
                {995,   0.045 },
                {1006,   0.071 },
                {1017,   0.107 },
                {1028,   0.144 },
                {1039,   0.192 },
                {1050,   0.243 },
                {1061,   0.33 },
                {1072,   0.401 },
                {1082,   0.491 },
                {1093,   0.573 },
                {1103,   0.664 },
                {1114,   0.773 },
                {1124,   0.885 },
                {1134,   1.061 },
                {1144,   1.149 },
                {1154,   1.24 },
                {1164,   1.473 },
                {1174,   1.548 },
                {1184,   1.594 },
                {1193,   1.731 },
                {1203,   1.904 },
                {1212,   2.105 },
                {1221,   2.076 },
                {1230,   2.158 },
                {1239,   2.265 },
                {1248,   2.289 },
                {1257,   2.406 },
                {1266,   2.519 },
                {1275,   2.659 },
                {1283,   2.718 },
                {1291,   2.666 },
                {1300,   2.683 },
                {1308,   2.904 },
                {1316,   3.035 },
                {1322,   2.941 },
                {startFit,fit(startFit*GeV)}
                //                {1335,   3.127 },
                //                {1343,   3.025 }
            });





static ProductionDataBase::XSections_t::value_type gp_pEtaPi0 =
{ ParticleTypeTreeDatabase::Channel::gp_pEtaPi0,
  gp_pEtaPi0_Data::Get};

