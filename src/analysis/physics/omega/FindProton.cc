#include "FindProton.h"

#include "utils/particle_tools.h"
#include "utils/matcher.h"
#include "utils/ProtonPermutation.h"
#include "expconfig/ExpConfig.h"
#include "TTree.h"

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
    tree->Branch("mangle",   &mangle);
}

void FindProton::branches_t::Reset()
{
    chi2dof     = std_ext::NaN;
    probability = std_ext::NaN;
    copl_angle  = std_ext::NaN;
    angle_p_cp  = std_ext::NaN;
    double mangle      = std_ext::NaN;
    matched_p   = -1;
    isBest      = -1;
    fitstatus   = -1;
}



FindProton::FindProton(const string& name, OptionsPtr opts):
    Physics(name,opts),
    fitter("FindProton", 3)
{
    tree = HistFac.makeTTree("tree");
    branches.SetBranchtes(tree);

    const auto setup = ant::ExpConfig::Setup::GetLastFound();

    if(!setup) {
        throw std::runtime_error("No Setup found");
    }

    fitter.LoadSigmaData(setup->GetPhysicsFilesDirectory()+"/FitterSigmas.root");
}

FindProton::~FindProton()
{

}

void FindProton::ProcessEvent(const TEvent& event, Physics::manager_t&)
{

    const auto cands = event.Reconstructed->Candidates;

    if(cands.size() != 4)
        return;

    TCandidatePtr matchedProton = nullptr;

    if(event.MCTrue  && event.MCTrue->ParticleTree) {

        const auto mcparticles = event.MCTrue->Particles.GetAll();

        const auto ptree  = event.MCTrue->ParticleTree;

        auto true_proton  = utils::ParticleTools::FindParticle (ParticleTypeDatabase::Proton, ptree, 1);
        auto true_photons = utils::ParticleTools::FindParticles(ParticleTypeDatabase::Photon, ptree);

        if(true_proton && true_photons.size() == 3) {

            const auto matched  = utils::match1to1(mcparticles,
                                                   cands,
                                                   [] (const TParticlePtr& p1, const TCandidatePtr& p2) {
                return p1->Angle(*p2);
            }, {0.0, degree_to_radian(15.0)});

            if(matched.size() == mcparticles.size()) {

                matchedProton = utils::FindMatched(matched, true_proton);
            }
        }
    }

    branches.Reset();
    branches_t best_combo;
    best_combo.probability = -100;

    for(const auto taggerhit : event.Reconstructed->TaggerHits) {
        for(utils::ProtonPermutation perm(event.Reconstructed->Candidates, matchedProton); perm.Good(); perm.Next()) {

            if(matchedProton) {
                branches.matched_p = (perm.isTrueProton()) ? 1 : 0;
            } else {
                branches.matched_p = -2; // = No mc true info
            }


            fitter.SetEgammaBeam(taggerhit.PhotonEnergy);
            fitter.SetProton(perm.Proton());
            fitter.SetPhotons(perm.Photons());

            const auto fitres = fitter.DoFit();

            branches.chi2dof     = fitres.ChiSquare / fitres.NDoF;
            branches.probability = fitres.Probability;
            branches.fitstatus   = fitres.Status == APLCON::Result_Status_t::Success ? 1 : 0;

            if(branches.probability > best_combo.probability)
                best_combo = branches;

            tree->Fill();
        }
    }

    if(best_combo.probability != -100) {
        best_combo.isBest = 1;
        branches = best_combo;
        tree->Fill();
    }

}


AUTO_REGISTER_PHYSICS(FindProton)
