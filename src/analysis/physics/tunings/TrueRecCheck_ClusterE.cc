#include "TrueRecCheck_ClusterE.h"

#include "base/Logger.h"
#include "expconfig/ExpConfig.h"
#include "utils/ParticleTools.h"

// use some namespaces
using namespace std;
using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::physics;

TrueRecCheck_ClusterE::TrueRecCheck_ClusterE(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    CBThetaWindow(degree_to_radian(50.0), degree_to_radian(180.0-50.0)),
    TAPSThetaWindow(degree_to_radian(3.0), degree_to_radian(18.0)),
    CBHemisphereGap({degree_to_radian(interval<double>::CenterWidth(0.0,40.0)),degree_to_radian(interval<double>::CenterWidth(180.0,40.0))})
{
    string partname[] = {"p","ep","em","g"};
    string parttitle[] = {"p","e^{+}","e^{-}","#gamma"};
    string detname[] = {"CB","TAPS"};
    for(int i=0; i<nrParticles; i++){
        for(int j=0; j<2; j++){
            h_Ek_TrueRecvsRec[j][i] = HistFac.makeTH2D(Form("E_{k}^{true} - E_{k}^{rec} vs E_{k}^{rec} for %s in %s",parttitle[i].c_str(),detname[j].c_str()),"E_{k}^{rec}","E_{k}^{true} - E_{k}^{rec}",BinSettings(250,0,1000),BinSettings(100,-50.,50.),Form("h_Ek_TrueRecvsRec_%s_%s",partname[i].c_str(),detname[j].c_str()),true);
            h_PairedOpAngle[j][i] = HistFac.makeTH1D(Form("Opening angle for %s in %s",parttitle[i].c_str(),detname[j].c_str()),"#theta_{OA}","",BinSettings(90,0.,90.),Form("h_PairedOpAngle_%s_%s",partname[i].c_str(),detname[j].c_str()),true);

        }
    }

}

void TrueRecCheck_ClusterE::ProcessEvent(const TEvent& event, manager_t&)
{

    //-- Make sure that also a file with true particle information is provided
    if(!event.MCTrue().ParticleTree){
        throw Exception("Must also provide a file with true particle information as input");
    }
    auto mcparticleslist = utils::ParticleTypeList::Make(event.MCTrue().ParticleTree);

    //-- Check if the data is single particle MC or not
    bool singleMC = false;
    string decay = utils::ParticleTools::GetDecayString(event.MCTrue().ParticleTree);
    if (decay.find("ParticleGun") != std::string::npos) {
        singleMC = true;
    }

    //-- Single MC
    if(singleMC){
        if(mcparticleslist.GetAll().size() != 1){
            throw Exception("Shouldn't be more than 1 particle in single particle MC?!");
        }

        auto& p_true = mcparticleslist.GetAll().front();

        //-- Find out which type of particle it is
        int ptype = -1;
        if(p_true->Type() == ParticleTypeDatabase::Proton) ptype = en_p;
        else if(p_true->Type() == ParticleTypeDatabase::ePlus) ptype = en_ep;
        else if(p_true->Type() == ParticleTypeDatabase::eMinus) ptype = en_em;
        else if(p_true->Type() == ParticleTypeDatabase::Photon) ptype = en_g;
        else {LOG(INFO)<<"The simulated particle was neither a proton, lepton or photon.";  return;}

        //-- Fetch the reconstructed candidate, and only use 1-cluster events
        if(event.Reconstructed().Candidates.size() != 1)
            return;
        auto& cand = event.Reconstructed().Candidates.front();

        //-- Check if it is in CB or TAPS
        int dtype = -1;
        TClusterPtr clu = cand.FindCaloCluster();
        if(clu->DetectorType == Detector_t::Type_t::CB){
            dtype = en_cb;
            //-- Only care about candidate who are not too close to the holes and hemispheregap
            if(CBThetaWindow.Contains(clu->Position.Theta()) && !CBHemisphereGap.Contains(clu->Position.Phi())){
                h_Ek_TrueRecvsRec[dtype][ptype]->Fill(cand.CaloEnergy,p_true->Ek()-cand.CaloEnergy);
                h_PairedOpAngle[dtype][ptype]->Fill((p_true->Angle(clu->Position))*radtodeg);
            }
        }
        else {
            dtype = en_ta;
            //-- Only care about candidate who are not too close to the beamline or the edges of TAPS
            if(TAPSThetaWindow.Contains(clu->Position.Theta())){
                h_Ek_TrueRecvsRec[dtype][ptype]->Fill(cand.CaloEnergy,p_true->Ek()-cand.CaloEnergy);
                h_PairedOpAngle[dtype][ptype]->Fill((p_true->Angle(clu->Position))*radtodeg);
            }
        }
    }

    //-- MC from decay chains (should be the only other option, right?)
    else {

        //-- Match the true particles to the reconstructed candidates
        const auto mcparticles = mcparticleslist.GetAll();
        const auto& candidates = event.Reconstructed().Candidates;
        const auto matched  = utils::match1to1(mcparticles, candidates.get_ptr_list(),
                                               [] (const TParticlePtr& p1, const TCandidatePtr& p2) {return p1->Angle(*p2);},
                                               {0.0, std_ext::degree_to_radian(15.0)});

        //-- Fill their corresponding histograms for each true particle
        for(auto& p_true: mcparticles){

            //-- Find the particle type
            int ptype = -1;
            if(p_true->Type() == ParticleTypeDatabase::Proton) ptype = en_p;
            else if(p_true->Type() == ParticleTypeDatabase::ePlus) ptype = en_ep;
            else if(p_true->Type() == ParticleTypeDatabase::eMinus) ptype = en_em;
            else if(p_true->Type() == ParticleTypeDatabase::Photon) ptype = en_g;
            else {LOG(INFO)<<"The simulated particle was neither a proton, lepton or photon.";  return;}

            //-- Get the reconstructed candidate it matched with
            TCandidatePtr matchcand = utils::FindMatched(matched,p_true);
            if(matchcand){
                //-- Check if it is in CB or TAPS
                int dtype = -1;
                TClusterPtr clu = matchcand->FindCaloCluster();
                if(clu->DetectorType == Detector_t::Type_t::CB){
                    dtype = en_cb;
                    //-- Only care about candidate who are not too close to the holes and hemispheregap
                    if(CBThetaWindow.Contains(clu->Position.Theta()) && !CBHemisphereGap.Contains(clu->Position.Phi())){
                        h_Ek_TrueRecvsRec[dtype][ptype]->Fill(matchcand->CaloEnergy,p_true->Ek()-matchcand->CaloEnergy);
                        h_PairedOpAngle[dtype][ptype]->Fill((p_true->Angle(clu->Position))*radtodeg);
                    }
                }
                else {
                    dtype = en_ta;
                    //-- Only care about candidate who are not too close to the beamline or the edges of TAPS
                    if(TAPSThetaWindow.Contains(clu->Position.Theta())){
                        h_Ek_TrueRecvsRec[dtype][ptype]->Fill(matchcand->CaloEnergy,p_true->Ek()-matchcand->CaloEnergy);
                        h_PairedOpAngle[dtype][ptype]->Fill((p_true->Angle(clu->Position))*radtodeg);
                    }
                }
            }
        }
    }
}


AUTO_REGISTER_PHYSICS(TrueRecCheck_ClusterE)


