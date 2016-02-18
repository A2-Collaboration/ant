#include "Fitter.h"

#include "expconfig/ExpConfig.h"

#include "base/std_ext/memory.h"
#include "base/std_ext/math.h"
#include "base/std_ext/system.h"
#include "base/ParticleType.h"
#include "base/WrapTFile.h"
#include "base/Logger.h"
#include "base/Paths.h"

#include "APLCON.hpp" // external project

#include "TTree.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TF1.h"

#include <cassert>

#include "utils/particle_tools.h"

using namespace std;
using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::utils;


Fitter::Fitter(const string& fittername)
{
    APLCON::Fit_Settings_t settings = APLCON::Fit_Settings_t::Default;
    settings.MaxIterations = 100;
    aplcon = make_unique<APLCON>(fittername, settings);
}

APLCON::Result_t Fitter::DoFit() {

    const auto res = aplcon->DoFit();

    result_chi2ndof    = res.ChiSquare / res.NDoF;
    result_iterations  = res.NIterations;
    result_status      = static_cast<int>(res.Status);
    result_probability = res.Probability;

    return res;
}

void Fitter::SetupBranches(TTree* tree, string branch_prefix)
{
    tree->Branch((branch_prefix+"_chi2dof").c_str(),     &result_chi2ndof);
    tree->Branch((branch_prefix+"_iterations").c_str(),  &result_iterations);
    tree->Branch((branch_prefix+"_status").c_str(),      &result_status);
    tree->Branch((branch_prefix+"_probability").c_str(), &result_probability);
}

double Fitter::EnergyResolution(const TParticlePtr& p) const
{
    assert(p->Candidate!=nullptr);

    if(p->Candidate->Detector & Detector_t::Type_t::CB) {

        return 1.07134e-02 * p->Ek();

    } else if(p->Candidate->Detector & Detector_t::Type_t::TAPS) {

        return 3.5E-2 * p->Ek();

    } else {
        LOG(WARNING) << "Photon not in (CB,TAPS?)";
    }


    return 0.0;
}

double Fitter::ThetaResolution(const TParticlePtr& p) const
{
    assert(p->Candidate!=nullptr);

    if(p->Candidate->Detector & Detector_t::Type_t::CB) {

        return cb_sigma_theta.GetSigma(p->Candidate->FindCaloCluster()->CentralElement, p->Ek());

    } if(p->Candidate->Detector & Detector_t::Type_t::TAPS) {

        return taps_sigma_theta.GetSigma(p->Candidate->FindCaloCluster()->CentralElement, p->Ek());

    } else {
        LOG(WARNING) << "Photon not in (CB,TAPS?)";
    }

    return 0.0;
}

double Fitter::PhiResolution(const TParticlePtr& p) const
{
    assert(p->Candidate!=nullptr);

    if(p->Candidate->Detector & Detector_t::Type_t::CB) {

        return cb_sigma_phi.GetSigma(p->Candidate->FindCaloCluster()->CentralElement, p->Ek());

    } else if(p->Candidate->Detector & Detector_t::Type_t::TAPS) {

        return taps_sigma_phi.GetSigma(p->Candidate->FindCaloCluster()->CentralElement, p->Ek());

    } else {
        LOG(WARNING) << "Photon not in (CB,TAPS?)";
    }

    return 0.0;
}

TLorentzVector Fitter::GetVector(const std::vector<double>& EkThetaPhi, const double m)
{
    const mev_t E = EkThetaPhi[0] + m;
    const mev_t p = sqrt( sqr(E) - sqr(m) );

    /// \bug This might be inefficient...

    TVector3 pv(1,0,0);

    pv.SetMagThetaPhi(p, EkThetaPhi[1], EkThetaPhi[2]);

    return TLorentzVector(pv,E);
}

double Fitter::fct_TaggerEGausSigma(double)
{
    return  3.0/sqrt(12.0);
}

void Fitter::FitParticle::Var_t::SetupBranches(TTree* tree, const string& prefix)
{
    tree->Branch(prefix.c_str(), addressof(Value));
    tree->Branch((prefix+"_pull").c_str(), addressof(Pull));
    tree->Branch((prefix+"_sigma").c_str(), addressof(Sigma));
}

void Fitter::FitParticle::SetupBranches(TTree* tree, const string& prefix)
{
    Ek.SetupBranches(tree, prefix+"_"+Name+"_Ek");
    Theta.SetupBranches(tree, prefix+"_"+Name+"_Theta");
    Phi.SetupBranches(tree, prefix+"_"+Name+"_Phi");
}

Fitter::angular_sigma::angular_sigma()
{
}

Fitter::angular_sigma::~angular_sigma()
{}

double Fitter::angular_sigma::GetSigma(const unsigned element, const double E) const
{
    const double vp0 = p0->GetBinContent(int(element));
    const double vp1 = p1->GetBinContent(int(element));
    const double vp2 = p2->GetBinContent(int(element));

    return f(E, vp0, vp1, vp2);

}

double Fitter::angular_sigma::f(const double x, const double p0, const double p1, const double p2) noexcept
{
    //return p0*exp(p1*x)+p2;
    return exp(p0 + p1*x) + p2;
}

double Fitter::angular_sigma::f_root(const double* x, const double* p) noexcept
{
    return f(x[0], p[0], p[1], p[2]);
}

TF1* Fitter::angular_sigma::GetTF1(const std::string& name)
{
    auto f = new TF1(name.c_str(), f_root, 0, 1600, 3);
    f->SetParName(0, "#alpha");
    f->SetParName(1, "#beta");
    f->SetParName(2, "Offset");
    return f;
}




void Fitter::angular_sigma::Load(WrapTFile& f, const string& prefix, const int bins)
{
    p0 = LoadHist(f, prefix+"_p0", bins);
    p1 = LoadHist(f, prefix+"_p1", bins);
    p2 = LoadHist(f, prefix+"_p2", bins);
}

Fitter::angular_sigma::Hist Fitter::angular_sigma::LoadHist(WrapTFile& f, const string& name, const int bins)
{
    auto h = f.GetSharedHist<TH1D>(name);

    if(!h)
        throw std::runtime_error("Could not load histogram " + name);

    if(h->GetNbinsX() != bins)
        throw std::runtime_error(formatter() << "TH1D " << name << " has wrong number of bins: " << h->GetNbinsX()  << ", expected: " << bins);

    return h;
}

void Fitter::LoadSigmaData(const string& filename)
{
    const auto& setup = ant::ExpConfig::Setup::GetLastFound();
    if(!setup)
        throw std::runtime_error("No Setup found!");

    unique_ptr<WrapTFileInput> f;
    if(!std_ext::system::testopen(filename)) {
        LOG(WARNING) << "Could not read sigmas from '" << filename << "', using default.";
        f = std_ext::make_unique<WrapTFileInput>(string(ANT_PATH_DATABASE)+"/default/physics_files/FitterSigmas.root");
    }
    else {
        f = std_ext::make_unique<WrapTFileInput>(filename);
    }


    const auto& CB = setup->GetDetector(Detector_t::Type_t::CB);
    if(!CB)
        throw std::runtime_error("No CB detector defined in setup");

    const auto nCB = CB->GetNChannels();

    const auto& TAPS = setup->GetDetector(Detector_t::Type_t::TAPS);
    if(!TAPS)
        throw std::runtime_error("No TAPS detector defined in setup");

    const auto nTAPS = TAPS->GetNChannels();

    cb_sigma_theta.Load(*f, "CB_sigma_Theta", int(nCB));
    cb_sigma_phi.Load(  *f, "CB_sigma_Phi",   int(nCB));

    taps_sigma_theta.Load(*f, "TAPS_sigma_Theta", int(nTAPS));
    taps_sigma_phi.Load(  *f, "TAPS_sigma_Phi",   int(nTAPS));
}

KinFitter::KinFitter(const std::string& name, unsigned numGammas) :
    Fitter(name)
{

    Photons.reserve(numGammas);
    for(unsigned i=0; i<numGammas;++i) {
        Photons.emplace_back(FitParticle("Photon"+to_string(i)));
    }

    aplcon->LinkVariable(Beam.Name,
                         Beam.Adresses() ,
                         Beam.Adresses_Sigma() );

    aplcon->LinkVariable(Proton.Name,
                         Proton.Addresses(),
                         Proton.Addresses_Sigma(),
                         Proton.Addresses_Pulls());

    vector<string> namesLInv      = { Beam.Name, Proton.Name };

    for ( auto& photon: Photons)
    {
        aplcon->LinkVariable(photon.Name,
                             photon.Addresses(),
                             photon.Addresses_Sigma(),
                             photon.Addresses_Pulls());

        namesLInv.push_back(photon.Name);
    }

    auto LorentzInvariance = [] (const vector<vector<double>>& values)
    {
        // beam    TLorentzVector(0.0, 0.0, PhotonEnergy(), PhotonEnergy());
        // target  TLorentzVector(0.0, 0.0, 0.0, ParticleTypeDatabase::Proton.Mass())
        const auto& Ebeam  = values[0][0];
        const auto& proton = values[1];

        //  Beam-LV:
        const TLorentzVector beam(0, 0, Ebeam, Ebeam);
        const TLorentzVector tg(0,0,0,ParticleTypeDatabase::Proton.Mass());

        TLorentzVector constraint = tg + beam;

        constraint -= GetVector(proton, ParticleTypeDatabase::Proton.Mass());

        const auto s = values.size();
        for ( unsigned photon = 0 ; photon < s-2 ; ++ photon)
            constraint -= GetVector(values[photon + 2], ParticleTypeDatabase::Photon.Mass());

        return vector<double>(
               { constraint.X(),
                 constraint.Y(),
                 constraint.Z(),
                 constraint.E()} );
    };

    aplcon->AddConstraint("LInv",{namesLInv},LorentzInvariance);

}

void KinFitter::SetEgammaBeam(const double &ebeam)
{
    Beam.E     = ebeam;
    Beam.Sigma = fct_TaggerEGausSigma(ebeam);
}

void KinFitter::SetProton(const TParticlePtr& proton)
{
    Proton.Ek.Value         = proton->Ek();
    Proton.Ek.Sigma    = 0.0; // unmeasured

    Proton.Theta.Value      = proton->Theta();
    Proton.Phi.Value        = proton->Phi();

    assert(proton->Candidate!=nullptr);

    if(proton->Candidate->Detector & Detector_t::Type_t::CB) {
        Proton.Theta.Sigma = degree_to_radian(5.43);
        Proton.Phi.Sigma   = degree_to_radian(5.31);
    } else if(proton->Candidate->Detector & Detector_t::Type_t::TAPS) {
        Proton.Theta.Sigma = degree_to_radian(2.89);
        Proton.Phi.Sigma   = degree_to_radian(4.45);
    } else {
        LOG(WARNING) << "Proton is not in (CB,TAPS)?";
    }

}

TParticlePtr KinFitter::GetFittedProton() const
{
    return make_shared<TParticle>(ParticleTypeDatabase::Proton, Proton.Ek.Value,
                                  Proton.Theta.Value, Proton.Phi.Value);
}

TParticleList KinFitter::GetFittedPhotons() const
{
    TParticleList photons;
    for(const auto& photon : Photons) {
        photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon,
                                                    photon.Ek.Value, photon.Theta.Value, photon.Phi.Value));
    }
    return photons;
}

double KinFitter::GetFittedBeamE() const
{
    return Beam.E;
}

void KinFitter::SetPhotons(const std::vector<TParticlePtr> &photons_data)
{
    if(Photons.size() != photons_data.size())
        throw Exception("Given number of photons does not match configured fitter");

    for ( unsigned i = 0 ; i < Photons.size() ; ++ i) {
        auto& photon = Photons.at(i);
        auto& data   = photons_data.at(i);

        photon.Ek.Value  = data->Ek();
        photon.Ek.Sigma = EnergyResolution(data);

        photon.Theta.Value  = data->Theta();
        photon.Theta.Sigma = ThetaResolution(data);

        photon.Phi.Value  = data->Phi();
        photon.Phi.Sigma = PhiResolution(data);
    }
}

void KinFitter::SetupBranches(TTree* tree, string branch_prefix)
{
    if(branch_prefix.empty())
        branch_prefix = aplcon->GetName();

    Proton.SetupBranches(tree, branch_prefix);
    for(auto& p : Photons) {
        p.SetupBranches(tree, branch_prefix);
    }

    Fitter::SetupBranches(tree, branch_prefix);
}

KinFitter::PhotonBeamVector::PhotonBeamVector(const string& name):
    Name(name)
{

}


TreeFitter::TreeFitter(const string& name, ParticleTypeTree ptree) :
    Fitter(name)
{
    tree = Tree<Node_t>::MakeNode(ptree);
    WalkTree(ptree, tree);

    vector<unsigned> p;
    for(unsigned i=0;i<tree_leaves.size();i++)
        p.push_back(i);

    struct unique_perm_t {
        vector<unsigned> Perm;
        tree_t Tree;
    };

    vector<unique_perm_t> unique_permutations;

    do {
        for(unsigned i=0;i<tree_leaves.size();i++)
            tree_leaves[i]->Get().leave_sum = 1 << p[i];

        CalcSum(tree);

        tree->Sort([] (const Node_t& a, const Node_t& b) {
            bool a_less_b = a.TypeTree->Get() < b.TypeTree->Get();
            bool b_less_a = b.TypeTree->Get() < a.TypeTree->Get();
            if(!a_less_b && !b_less_a) {
                return a.leave_sum < b.leave_sum;
            }
            return a_less_b;
        });

        bool found = false;
        for(const unique_perm_t& u : unique_permutations) {
            if(IsEqual(u.Tree, tree)) {
                found = true;
                break;
            }
        }
        if(!found) {
            unique_perm_t u;
            u.Perm = p;
            auto t  = Tree<Node_t>::MakeNode(tree->Get().TypeTree);
            CopyTree(tree, t);
            u.Tree = t;
            unique_permutations.emplace_back(move(u));
        }
    }
    while(std::next_permutation(p.begin(), p.end()));

    cout << "Unique Permutations:" << unique_permutations.size() << endl;
    for(const unique_perm_t& u : unique_permutations) {
        cout << u.Perm << endl;
    }

    cout << "Decay: " << ParticleTools::GetDecayString(ptree) << endl;
}

void TreeFitter::WalkTree(ParticleTypeTree ptree, tree_t tree)
{
    for(auto daughter : ptree->Daughters()) {
        auto& d = tree->CreateDaughter(daughter);
        if(daughter->Daughters().empty())
            tree_leaves.push_back(d);
        else
            WalkTree(daughter, d);
    }
}

void TreeFitter::CalcSum(tree_t tree)
{
    if(!tree->Daughters().empty())
        tree->Get().leave_sum = 0;
    for(auto daughter : tree->Daughters()) {
        CalcSum(daughter);
        tree->Get().leave_sum += daughter->Get().leave_sum;
    }
}

void TreeFitter::CopyTree(tree_t source, tree_t dest)
{
    dest->Get().leave_sum = source->Get().leave_sum;
    for(auto daughter : source->Daughters()) {
        CopyTree(daughter, dest->CreateDaughter(daughter->Get().TypeTree));
    }
}

bool TreeFitter::IsEqual(tree_t a, tree_t b)
{
    assert(a->Daughters().size() == b->Daughters().size());
    assert(a->Get().TypeTree->Get() == b->Get().TypeTree->Get());

    // ignore leaves
    if(a->Daughters().empty())
        return true;

    if(a->Get().leave_sum != b->Get().leave_sum)
        return false;

    auto it_a = a->Daughters().begin();
    auto it_b = b->Daughters().begin();

    for(; it_a != a->Daughters().end() && it_b != b->Daughters().end(); ++it_a, ++it_b) {
        if(!IsEqual(*it_a, *it_b))
            return false;
    }

    return true;
}
