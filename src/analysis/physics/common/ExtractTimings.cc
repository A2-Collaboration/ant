#include "ExtractTimings.h"
#include "base/std_ext/math.h"

#include "TTree.h"

using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;


ExtractTimings::ExtractTimings(const string& name, PhysOptPtr opts):
    Physics(name, opts)
{
    const BinSettings timebins(opts->Get<int>("TimeBins",800),-40,40);

    EPT_CB    = HistFac.makeTH1D("EPT to CB Time",    "Time [ns]", "", timebins, "ept_cb");
    EPT_TAPS  = HistFac.makeTH1D("EPT to TAPS Time",  "Time [ns]", "", timebins, "ept_TAPS");
    CB_TAPS   = HistFac.makeTH1D("CB to TAPS Time",   "Time [ns]", "", timebins, "cb_taps");
    CB_CB     = HistFac.makeTH1D("CB to CB Time",     "Time [ns]", "", timebins, "cb_cb");
    TAPS_TAPS = HistFac.makeTH1D("TAPS to TAPS Time", "Time [ns]", "", timebins, "taps_taps");
    CB        = HistFac.makeTH1D("CB Time",           "Time [ns]", "", timebins, "cb");
    TAPS      = HistFac.makeTH1D("TAPS Time",         "Time [ns]", "", timebins, "taps");
    EPT       = HistFac.makeTH1D("EPT Time",          "Time [ns]", "", timebins, "ept");
}

template <typename T>
T inc(T x) {
    return ++x;
}

void ExtractTimings::ProcessEvent(const data::Event& event)
{
    const auto& photons = event.Reconstructed().Particles().Get(ParticleTypeDatabase::Photon);

    // EPT-CB and EPT-TAPS
    for(const auto& th : event.Reconstructed().TaggerHits()) {

        EPT->Fill(th->Time());

        for(const auto& p : photons) {

            if(p->Candidate()) {
                if(p->Candidate()->Detector() & Detector_t::Type_t::CB) {
                    EPT_CB->Fill(th->Time() - p->Candidate()->Time());
                } else if(p->Candidate()->Detector() & Detector_t::Type_t::TAPS) {
                    EPT_TAPS->Fill(th->Time() - p->Candidate()->Time());
                }

            }
        }
    }


    for(auto i = photons.cbegin(); i != photons.cend(); ++i) {

        const auto& pi = *i;

        if(pi->Candidate()) {

            if(pi->Candidate()->Detector() & Detector_t::Type_t::CB) {
                CB->Fill(pi->Candidate()->Time());
            } else if(pi->Candidate()->Detector() & Detector_t::Type_t::TAPS) {
                TAPS->Fill(pi->Candidate()->Time());
            }

            for(auto j = inc(i); j != photons.cend(); ++j) {

                const auto& pj = *j;

                // CB-TAPS
                if(pj->Candidate()) {
                    if((pi->Candidate()->Detector() & Detector_t::Type_t::CB && pj->Candidate()->Detector() & Detector_t::Type_t::TAPS)
                       || (pi->Candidate()->Detector() & Detector_t::Type_t::TAPS && pj->Candidate()->Detector() & Detector_t::Type_t::CB) ) {

                        double td = pi->Candidate()->Time() - pj->Candidate()->Time();

                        if(pi->Candidate()->Detector() & Detector_t::Type_t::TAPS) {
                            td *= -1.0;
                        }

                        CB_TAPS->Fill(td);
                    }
                }

                if(pi->Candidate()->Detector() & Detector_t::Type_t::CB && pj->Candidate()->Detector() & Detector_t::Type_t::CB) {
                    CB_CB->Fill(pi->Candidate()->Time() - pj->Candidate()->Time());
                }

                if(pi->Candidate()->Detector() & Detector_t::Type_t::TAPS && pj->Candidate()->Detector() & Detector_t::Type_t::TAPS) {
                    TAPS_TAPS->Fill(pi->Candidate()->Time() - pj->Candidate()->Time());
                }

            }
        }
    }
}

void ExtractTimings::ShowResult()
{

    canvas("Extract Timings")
            << EPT_CB
            << EPT_TAPS
            << CB_TAPS
            << CB_CB
            << TAPS_TAPS
            << CB
            << TAPS
            << EPT
            << endc;
}

AUTO_REGISTER_PHYSICS(ExtractTimings)
