#include "IMPlots.h"
#include "utils/combinatorics.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;

IMPlots::IMPlots(PhysOptPtr opts): Physics("IMPlots", opts),
  cb("CB",HistFac),
  taps("TAPS", HistFac)
{
}

void IMPlots::ProcessEvent(const data::Event& event)
{
    const auto& photons = event.Reconstructed().Particles().Get(ParticleTypeDatabase::Photon);

    for(unsigned n=2;n<=10;++n) {
    for( auto comb = utils::makeCombination(photons,n); !comb.Done(); ++comb) {
        TLorentzVector....
    }
    }
}

void IMPlots::hist_set::Fill(unsigned ngamma, double mm)
{
    m.at(ngamma)->Fill(mm);
}

IMPlots::hist_set::hist_set(const std::string& pref, SmartHistFactory& hf, std::size_t n):
    m(n-2)
{
    const BinSettings im(1600);

    for(size_t i=0;i <n-2; ++i) {
        m[i] = hf.makeTH1D(pref + " " + to_string(i+2) + " #gamma IM", "IM [MeV]", "", im, pref+"_"+to_string(i+2));
    }
}
