#include "JustPi0.h"

#include "utils/particle_tools.h"

#include "TH1D.h"

#include <memory>
#include <cassert>

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::analysis::data;
using namespace std;



JustPi0::JustPi0(const string& name, PhysOptPtr opts) :
    Physics(name, opts)
{
    for(unsigned mult=1;mult<=3;mult++) {
        multiPi0.emplace_back(HistFac, mult);
    }
}

void JustPi0::ProcessEvent(const Event& event)
{
    const auto& data = event.Reconstructed;
    for(auto& m : multiPi0)
        m.ProcessData(data);
}

void JustPi0::ShowResult()
{
    for(auto& m : multiPi0)
        m.ShowResult();
}




JustPi0::MultiPi0::MultiPi0(SmartHistFactory& histFac, unsigned nPi0) :
    multiplicity(nPi0),
    h_missingmass(promptrandom),
    IM_2g(promptrandom)
{
    std::string multiplicity_str = std_ext::formatter() << multiplicity << "Pi0";
    SmartHistFactory HistFac(multiplicity_str, histFac, multiplicity_str);

    promptrandom.AddPromptRange({-2.5,2.5});
    promptrandom.AddRandomRange({-50,-5});
    promptrandom.AddRandomRange({  5,50});

    steps = HistFac.makeTH1D("Steps","","#",BinSettings(10),"steps");

    h_missingmass.MakeHistograms(HistFac, "h_missingmass","Missing Mass",BinSettings(400,400, 1400),"MM / MeV","#");
    IM_2g.MakeHistograms(HistFac, "IM_2g","Invariant Mass 2#gamma",BinSettings(500,0,700),"IM / MeV","#");
}

void JustPi0::MultiPi0::ProcessData(const Event::Data& data)
{
    const auto nPhotons_expected = multiplicity*2;

    steps->Fill("Seen",1);

    // cut on energy sum and number of candidates

    if(data.Trigger.CBEnergySum <= 550)
        return;
    steps->Fill("CBESum>550MeV",1);

    const auto& cands = data.Candidates;
    const auto nCandidates = cands.size();
    const auto nCandidates_expected = nPhotons_expected+1;
    if(nCandidates != nCandidates_expected)
        return;
    std::string nCandidates_cutstr = std_ext::formatter() << "nCandidates==" << nCandidates_expected;
    steps->Fill(nCandidates_cutstr.c_str(),1);


    // use any candidate as proton, and do the analysis (ignore ParticleID stuff)

    for(auto i_proton=cands.begin();i_proton!=cands.end();i_proton++) {

        const auto proton = std::make_shared<Particle>(ParticleTypeDatabase::Proton, *i_proton);
        std::vector<ParticlePtr> photons;
        for(auto i_photon=cands.begin();i_photon!=cands.end();i_photon++) {
            if(i_photon == i_proton)
                continue;
            photons.emplace_back(make_shared<Particle>(ParticleTypeDatabase::Photon, *i_photon));
        }
        assert(photons.size() == nPhotons_expected);


        TLorentzVector photon_sum(0,0,0,0);
        for(const auto& p : photons) {
            photon_sum += *p;
        }

        for(const TaggerHit& taggerhit : data.TaggerHits) {
            promptrandom.SetTaggerHit(taggerhit.Time);
            if(promptrandom.State() == PromptRandom::Case::Outside)
                continue;

            const TLorentzVector beam_target = taggerhit.GetPhotonBeam() + TLorentzVector(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
            const TLorentzVector v_mm = beam_target - photon_sum;
            const double mm = v_mm.M();
            h_missingmass.Fill(mm);

            if(mm>850 && mm<1000) {
                steps->Fill("MM",1.0);
                utils::ParticleTools::FillIMCombinations([this] (double x) {IM_2g.Fill(x);},  2, photons);

            }

        }


    }
}

void JustPi0::MultiPi0::ShowResult()
{
    canvas(std_ext::formatter() << "JustPi0: " << multiplicity << "Pi0")
            << steps
            << h_missingmass.subtracted
            << IM_2g.subtracted
            << endc;
}

AUTO_REGISTER_PHYSICS(JustPi0)
