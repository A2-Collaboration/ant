#include "FindProton.h"

#include "utils/ParticleTools.h"
#include "utils/Matcher.h"
#include "utils/ProtonPermutation.h"
#include "expconfig/ExpConfig.h"
#include "TTree.h"
#include <list>
#include <memory>
#include "base/Logger.h"

#include "analysis/utils/uncertainties/FitterSergey.h"

using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;


void FindProton::branches_t::SetBranchtes(TTree* tree)
{
    tree->Branch("chi2dof",     &chi2dof);
    tree->Branch("probability", &probability);
    tree->Branch("copl_angle",  &copl_angle);
    tree->Branch("angle_p_cp",  &angle_p_cp);
    tree->Branch("matched_p",   &matched_p);
    tree->Branch("isBest",      &isBest);
    tree->Branch("fitstatus",   &fitstatus);
    tree->Branch("mangle",      &mangle);
    tree->Branch("p_theta",     &p_theta);
    tree->Branch("p_phi",       &p_phi);
    tree->Branch("p_PSA_r",     &p_PSA_r);
    tree->Branch("p_PSA_a",     &p_PSA_a);
    tree->Branch("p_veto",      &p_veto);
    tree->Branch("p_detector",  &p_detector);
    tree->Branch("TagW",        &TagW);
    tree->Branch("TagE",        &TagE);
    tree->Branch("TagTime",     &TagTime);
    tree->Branch("TagCh",       &TagCh);
}

void FindProton::branches_t::Reset()
{
    chi2dof     = std_ext::NaN;
    probability = std_ext::NaN;
    copl_angle  = std_ext::NaN;
    angle_p_cp  = std_ext::NaN;
    mangle      = std_ext::NaN;
    matched_p   = -1;
    isBest      = -1;
    fitstatus   = -1;

    p_theta     = std_ext::NaN;
    p_phi       = std_ext::NaN;
    p_PSA_r     = std_ext::NaN;
    p_PSA_a     = std_ext::NaN;
    p_veto      = std_ext::NaN;
    p_detector  = -1;
    TagW        = std_ext::NaN;
    TagE        = std_ext::NaN;
    TagTime     = std_ext::NaN;
    TagCh       = -1;

}



FindProton::FindProton(const string& name, OptionsPtr opts):
    Physics(name,opts),
    nPhotons(opts->Get<unsigned>("nPhotons", 3)),
    fitter(make_shared<utils::UncertaintyModels::FitterSergey>())
{
    tree = HistFac.makeTTree("tree");
    tree_branches.SetBranchtes(tree);

    steps = HistFac.makeTH1D("Steps","","",BinSettings(3), "steps");


    promptrandom.AddPromptRange({-5,5});
    promptrandom.AddRandomRange({-20, -10});
    promptrandom.AddRandomRange({ 10,  20});
}

FindProton::~FindProton()
{

}

void FindProton::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event);

    steps->Fill("Total", 1.0);

    const auto& cands = event.Reconstructed().Candidates;

    if(cands.size() != nPhotons+1)
        return;

    steps->Fill("nCands", 1.0);

    TCandidatePtr matchedProton = nullptr;

    if(event.MCTrue().ParticleTree) {

        auto mctrue_particles = utils::ParticleTypeList::Make(event.MCTrue().ParticleTree);

        const auto mcparticles = mctrue_particles.GetAll();

        const auto ptree  = event.MCTrue().ParticleTree;

        auto true_proton  = utils::ParticleTools::FindParticle (ParticleTypeDatabase::Proton, ptree, 1);
        auto true_photons = utils::ParticleTools::FindParticles(ParticleTypeDatabase::Photon, ptree);

        if(true_proton && true_photons.size() == 3) {

            const auto matched  = utils::match1to1(mcparticles,
                                                   cands.get_ptr_list(),
                                                   [] (const TParticlePtr& p1, const TCandidatePtr& p2) {
                return p1->Angle(*p2);
            }, {0.0, degree_to_radian(15.0)});

            if(matched.size() == mcparticles.size()) {

                matchedProton = utils::FindMatched(matched, true_proton);
            }
        }
    }




    for(const auto taggerhit : event.Reconstructed().TaggerHits) {

        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));

        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        shared_ptr<branches_t> best_combo = nullptr;

        list<shared_ptr<branches_t>> b;

        for(utils::ProtonPermutation perm(event.Reconstructed().Candidates.get_ptr_list(), matchedProton); perm.Good(); perm.Next()) {

            auto branches = make_shared<branches_t>();

            b.emplace_back(branches);

            branches->TagW = promptrandom.FillWeight();
            branches->TagE    = taggerhit.PhotonEnergy;
            branches->TagCh   = taggerhit.Channel;
            branches->TagTime = taggerhit.Time;
            branches->p_detector  = 0;

            if(perm.Proton()->Candidate) {

                if(perm.Proton()->Candidate->Detector & Detector_t::Type_t::TAPS) {
                    branches->p_detector = 2;
                    const auto& cluster = perm.Proton()->Candidate->FindCaloCluster();
                    if(cluster) {
                        branches->p_PSA_a = std_ext::radian_to_degree(cluster->GetPSAAngle());
                        branches->p_PSA_r = cluster->GetPSARadius();
                    }
                } else if(perm.Proton()->Candidate->Detector & Detector_t::Type_t::CB) {
                    branches->p_detector = 1;
                }

            }

            if(matchedProton) {
                branches->matched_p = (perm.isTrueProton()) ? 1 : 0;
            } else {
                branches->matched_p = -2; // = No mc true info
            }


            const auto fitres = fitter.DoFit(taggerhit.PhotonEnergy, perm.Proton(), perm.Photons());

            branches->chi2dof     = fitres.ChiSquare / fitres.NDoF;
            branches->probability = fitres.Probability;
            branches->fitstatus   = fitres.Status == APLCON::Result_Status_t::Success ? 1 : 0;

            branches->p_theta = perm.Proton()->Theta();
            branches->p_phi   = perm.Proton()->Phi();
            branches->p_veto  = perm.Proton()->Candidate->VetoEnergy;

            if(!best_combo || branches->probability > best_combo->probability)
                    best_combo = branches;

        }

        for(auto& br : b) {

            br->isBest = (br==best_combo) ? 1 : 0;

            tree_branches.Reset();

            tree_branches = *br;
            VLOG(6) << "isBest=" << tree_branches.isBest;
            tree->Fill();

        }
    }




}


AUTO_REGISTER_PHYSICS(FindProton)
