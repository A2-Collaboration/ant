#include "IMPlots.h"
#include "utils/combinatorics.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;

IMPlots::IMPlots(const std::string& name, PhysOptPtr opts): Physics(name, opts),
//  cb("CB",HistFac),
//  taps("TAPS", HistFac),
  all("All", HistFac)
{
}

template <typename iter>
TLorentzVector sumlv(iter start, iter end) {
    TLorentzVector s;
    while(start != end ) {
        s += **(start);
        ++start;
    }
    return s;
}


void IMPlots::ProcessEvent(const data::Event& event)
{
    const auto& photons = event.Reconstructed().Particles().Get(ParticleTypeDatabase::Photon);

    for(unsigned n = all.MinNGamma(); n<all.MaxNGamma(); ++n) {
        for( auto comb = utils::makeCombination(photons,n); !comb.Done(); ++comb) {
            const TLorentzVector sum = sumlv(comb.begin(), comb.end());
            all.Fill(n,sum.M());
        }
    }
}

canvas& operator<<(canvas& c, const IMPlots::hist_set& s) {
    for(auto& h : s.m) {
        c << h;
    }
    c << endr;
    return c;
}

void IMPlots::ShowResult()
{
    canvas c("IMPlots");
    c << all << endc;
}

void IMPlots::hist_set::Fill(unsigned ngamma, double mm)
{
    m.at(ngamma-2)->Fill(mm);
}

IMPlots::hist_set::hist_set(const std::string& pref, SmartHistFactory& hf, std::size_t n):
    m(n-2)
{
    const BinSettings im(1600);

    for(size_t i=0;i <n-2; ++i) {
        m[i] = hf.makeTH1D(pref + " " + to_string(i+2) + " #gamma IM", "IM [MeV]", "", im, pref+"_"+to_string(i+2));
    }
}

AUTO_REGISTER_PHYSICS(IMPlots)
