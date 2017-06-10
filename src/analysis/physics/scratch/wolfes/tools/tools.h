#pragma once

#include "tree/TEventData.h"
#include "tree/TCandidate.h"

#include "TH1D.h"
#include "TLorentzVector.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/EPT.h"
#include "expconfig/detectors/Tagger.h"

#include "analysis/physics/Plotter.h"

#include "TTree.h"
#include "tree/TSimpleParticle.h"


namespace ant {
namespace analysis {
namespace physics {

template<class WrapTree>
class TreePlotterBase_t: public Plotter{

protected:
    TTree* t = nullptr;
    WrapTree tree;
    unsigned nchannels;
    // Plotter interface
public:
    TreePlotterBase_t(const std::string& name, const WrapTFileInput& input, OptionsPtr opts):
        Plotter(name,input,opts),
        nchannels(ExpConfig::Setup::GetDetector<TaggerDetector_t>()->GetNChannels())
    {
        if(!input.GetObject(WrapTree::treeAccessName(),t))
            throw Exception("Input TTree not found");

        if(!tree.Matches(t))
            throw std::runtime_error("Tree branches don't match");
        tree.LinkBranches(t);
    }

    virtual long long GetNumEntries() const override {return t->GetEntries();}

    virtual ~TreePlotterBase_t(){}
};

template<class WrapTree>
class DetectionEffciencyBase_t: public TreePlotterBase_t<WrapTree>{

protected:
    TH1D* seenMC       = nullptr;
    TH1D* efficiencies = nullptr;

public:
    DetectionEffciencyBase_t(const std::string& name, const WrapTFileInput& input,
                             OptionsPtr opts):
        TreePlotterBase_t<WrapTree>(name,input,opts)
    {

        if(!input.GetObject("singlePi0/seenMC",seenMC))
            throw  std::runtime_error("Input TH1D for seenMC not found");

        efficiencies = this->HistFac.makeTH1D("reconstruction efficiency",
                                            "channel","eff [%]",
                                            BinSettings(this->nchannels),
                                            "eff",true);
        efficiencies->SetBit(TH1D::kIsAverage);
    }

    virtual ~DetectionEffciencyBase_t(){}
};

struct tools
{
    struct protonSelection_t
    {
        TParticlePtr      Proton;
        TParticleList     Photons;
        TCandidatePtr     CandP;
        TCandidatePtrList CandsG;
        LorentzVec     PhotonSum;
        LorentzVec     Proton_MM;
        LorentzVec     PhotonBeam;
        double         Copl_pg;
        double         Angle_pMM;
        double         Tagg_E;
        double         ProtonVetoE;
        double         PhotonVetoE;
        protonSelection_t(const TParticlePtr& proton, const TParticleList& photons,
                          const TCandidatePtr& candP, const TCandidatePtrList& candsG,
                          const LorentzVec& photonSum,
                          const LorentzVec& protonMM,
                          const LorentzVec& phtonBeam,
                          double copl_pg,
                          double angle_pMM,
                          double tagg_E,
                          double protonVetoE,
                          double photonVetoE):
            Proton(proton), Photons(photons),
            CandP(candP), CandsG(candsG),
            PhotonSum(photonSum),
            Proton_MM(protonMM),
            PhotonBeam(phtonBeam),
            Copl_pg(copl_pg),
            Angle_pMM(angle_pMM),
            Tagg_E(tagg_E),
            ProtonVetoE(protonVetoE),
            PhotonVetoE(photonVetoE){}
// neeed?
//        protonSelection_t() = default;

    };

    static std::vector<protonSelection_t> makeProtonSelections(
            const TCandidateList& candidates,
            const LorentzVec& photonBeam,
            double taggE,
            const IntervalD& imMMprotonCut = {-std_ext::inf,std_ext::inf}
            );

    template<class T>
    static bool cutOn(const std::string& fillName, const interval<T>& range, const T& val, TH1D* steps)
    {
        auto pass = range.Contains(val);
        if (pass)
        {
            std::string s = std_ext::formatter() << fillName << ": " << range;
            steps->Fill(s.c_str(),1);
        }
        return !pass;
    }

    static std::vector<TLorentzVector> MakeTLorenz(const TParticleList& particles)
    {
        std::vector<TLorentzVector> lg(particles.size());
        std::transform(particles.begin(),particles.end(),lg.begin(),
                       [](const TParticlePtr& ph){return TLorentzVector(*ph);});
        return lg;
    }


    static TCandidatePtrList getNeutral(const TEventData& data, const double threshold);

    static double getChargedClusterE(const TClusterList& clusters);
    static double getCandidateVetoE(const TCandidateList& cands);


    //strict stolen from oli:
    static double getCorrVetoEnergy(const TCandidate& photon, const TCandidate& proton);

    static double getPhotonVetoEnergy(const protonSelection_t& sel, const bool strict = false);


    static double getCBVetoEnergy(const TCandidatePtr& cand)
    {
        if (! (cand->Detector & Detector_t::Type_t::CB))
            return 0;
        return cand->VetoEnergy;
    }
    static double getCBProtonVeto(const protonSelection_t& sel)
    {
        return getCBVetoEnergy(sel.CandP);
    }
    static std::vector<double> getCBPhotonVeto(const protonSelection_t& sel)
    {
        std::vector<double> vetos(sel.CandsG.size());
        std::transform(sel.CandsG.begin(),sel.CandsG.end(),
                       vetos.begin(),
                       [](const TCandidatePtr& cptr) {return getCBVetoEnergy(cptr);});
        return vetos;
    }

    static std::vector<std::string> tokenize_cuts(const std::string& path)
    {
        std::stringstream ss(path);
        std::string s;
        std::vector<std::string> retval;

        while (getline(ss, s, '/')) {
         retval.push_back(s);
        }
        return retval;
    }

};


}
}
}
