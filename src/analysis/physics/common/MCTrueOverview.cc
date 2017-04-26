#include "MCTrueOverview.h"

#include "utils/ParticleTools.h"
#include "base/Logger.h"
#include "base/std_ext/string.h"

using namespace std;
using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::physics;


MCTrueOverview::MCTrueOverview(const std::string& name, OptionsPtr opts):
    Physics(name,opts)
{
}

void MCTrueOverview::ProcessEvent(const TEvent& event, manager_t&)
{
    auto& ptree = event.MCTrue().ParticleTree;
    if(!ptree)
        return;
    ParticleTypeTreeDatabase::Channel channel;
    ParticleTypeTree typetree;
    if(!utils::ParticleTools::TryFindParticleTypeTree(ptree, channel, typetree)) {
        LOG_N_TIMES(100, WARNING) << "Cannot find " << utils::ParticleTools::GetDecayString(ptree, false) << " in database (max 100x printed)";
        return;
    }

    // search for channel specific histograms
    auto it_perChannel = channels.find(channel);
    if(it_perChannel == channels.end()) {
        auto it = channels.emplace(make_pair(channel, perChannel_t(HistFac, typetree)));
        it_perChannel = it.first;
    }
    auto& perChannel = it_perChannel->second;
    perChannel.Fill(ptree, event.Reconstructed().Candidates);

}

void MCTrueOverview::ShowResult()
{
    for(auto& it_ch : channels) {
        auto decaystring = utils::ParticleTools::GetDecayString(ParticleTypeTreeDatabase::Get(it_ch.first), false);
        canvas c(GetName()+" "+decaystring);
        it_ch.second.Show(c);
        c << endc;
    }
}

inline int countLeaves(const ParticleTypeTree& t) {
    int n =0;
    t->Map_nodes([&n] (const ParticleTypeTree& node) { if(node->IsLeaf()) ++n;});
    return n;
}

MCTrueOverview::perChannel_t::perChannel_t(const HistogramFactory& histFac, const ParticleTypeTree& typetree):
    mult(histFac, countLeaves(typetree))
{
    HistogramFactory HistFac(utils::ParticleTools::GetDecayString(typetree, false),
                             histFac,
                             utils::ParticleTools::GetDecayString(typetree, true));

    histtree = typetree->DeepCopy<histnode_t>([HistFac] (const ParticleTypeTree& t) {
        // only create histograms at this node if a leaf daughter exists
        // indicated by non-empty leafTypes
        std::vector<histnode_t::typeptr_t> leafTypes;
        for(auto d  : t->Daughters()) {
            if(d->IsLeaf()) {
                leafTypes.push_back(addressof(d->Get()));
            }
        }
        auto histFacPtr = leafTypes.empty() ?
                              nullptr : // only create histogramfactory if really needed, prevents dummy directories
                              std_ext::make_unique<const HistogramFactory>(t->Get().Name(), HistFac, t->Get().PrintName());
        return histnode_t(move(histFacPtr), leafTypes);
    });

    const AxisSettings axis_CBEsum("#SigmaE_{kin}^{CB} / MeV", {300, 0, 1600});
    h_CBEsum_true = HistFac.makeTH1D("CB ESum (simple true)",axis_CBEsum, "h_CBEsum_true");
    h_CBEsum_rec = HistFac.makeTH1D("CB ESum (reconstructed)",axis_CBEsum, "h_CBEsum_rec");
}

void MCTrueOverview::perChannel_t::Fill(const TParticleTree_t& ptree, const TCandidateList& cands)
{
    // traverse through ptree in parallel to own tree histtree
    traverse_tree_and_fill(histtree, ptree);

    double CBEsum_true = 0;

    // really simplistic way of calculating the esum,
    // iterate over all leaf particles
    for(auto p : utils::ParticleTypeList::Make(ptree).GetAll()) {
        // super-simple CB geometry
        auto theta_deg = std_ext::radian_to_degree(p->Theta());
        if(!interval<double>(20,160).Contains(theta_deg))
            continue;
        // super-simple punch through for nucleons (realistic at least for protons?!)
        double Ek = p->Ek();
        if(p->Type() == ParticleTypeDatabase::Nucleon && Ek > 400)
            Ek = 400;
        CBEsum_true += Ek;
    }
    h_CBEsum_true->Fill(CBEsum_true);

    // more realistic way using reconstructed candidates (if any)
    double CBEsum_rec = 0;
    for(auto& cand : cands) {
        if(cand.Detector & Detector_t::Type_t::CB) {
            CBEsum_rec += cand.CaloEnergy;
        }
    }
    h_CBEsum_rec->Fill(CBEsum_rec);

    mult.Fill(ptree);
}

void MCTrueOverview::perChannel_t::Show(canvas& c) const
{
    c << h_CBEsum_true << h_CBEsum_rec;
    c << mult.h_CBTAPS;
    histtree->Map([&c] (const histnode_t& n) {
        n.Show(c);
    });
}

void MCTrueOverview::perChannel_t::traverse_tree_and_fill(const histtree_t& histtree,
                                                          const TParticleTree_t& ptree) const
{
    if(histtree->Daughters().size() != ptree->Daughters().size()) {
        throw runtime_error("Number of Daughters mismatch");
    }

    if(histtree->IsLeaf()) {
        return;
    }

    auto& d1 = histtree->Daughters();
    auto& d2 = ptree->Daughters();
    auto it_d1 = d1.begin();
    auto it_d2 = d2.begin();

    for(; it_d1 != d1.end() && it_d2 != d2.end(); ++it_d1, ++it_d2 ) {
        auto& node1 = *it_d1;
        auto& node2 = *it_d2;
        if(node1->IsLeaf()) {
            const auto& p = node2->Get();
            histtree->Get().Fill(*p);
        }
        else {
            // traverse in parallel
            traverse_tree_and_fill(node1, node2);
        }
    }
}

MCTrueOverview::perChannel_t::histnode_t::histnode_t(std::unique_ptr<const HistogramFactory> histFacPtr,
                                                     const vector<typeptr_t>& leafTypes)
{
    if(!histFacPtr)
        return;
    auto& HistFac = *histFacPtr;
    for(auto typeptr : leafTypes) {
        if(hists.find(typeptr) == hists.end())
            hists.emplace(make_pair(typeptr, HistogramFactory(typeptr->Name(), HistFac, typeptr->PrintName())));
    }
}

void MCTrueOverview::perChannel_t::histnode_t::Fill(const TParticle& p)
{
    auto& h  = hists.at(addressof(p.Type()));
    h.h_EkTheta->Fill(p.Ek(), std_ext::radian_to_degree(p.Theta()));
}

void MCTrueOverview::perChannel_t::histnode_t::Show(canvas& c) const
{
    for(auto& it_h : hists) {
        auto& h = it_h.second;
        c << drawoption("colz") << h.h_EkTheta;
    }
}



MCTrueOverview::perChannel_t::histnode_t::perType_t::perType_t(const HistogramFactory& HistFac)
{
    const AxisSettings axis_Theta("#theta / #circ", {50, 0, 180});
    const AxisSettings axis_Ek("E_{k} / MeV", {100, 0, 1000});
    h_EkTheta = HistFac.makeTH2D("E_{k} vs. #theta", axis_Ek, axis_Theta, "h_EkTheta");
}



MCTrueOverview::CBTAPS_Multiplicity::CBTAPS_Multiplicity(const HistogramFactory& histFac, const int nParticles):
    h_CBTAPS(histFac.makeTH1D("CB/TAPS Multiplicity","CB / TAPS","Particles/Event", BinSettings(unsigned(nParticles))))
{
    for(int i=1;i<=nParticles; ++i) {
        const string label = formatter() << nParticles-(i-1) << "/" << i-1;
        h_CBTAPS->GetXaxis()->SetBinLabel(i, label.c_str());
    }
}


void MCTrueOverview::CBTAPS_Multiplicity::Fill(const TParticleTree_t& particles)
{

    const auto taps = [particles] () {
        int n = 0;
        particles->Map_nodes([&n] (const TParticleTree_t& node) {
            if(node->IsLeaf() && node->Get()->Theta() < degree_to_radian(20.0))
                ++n;
        });
        return n;
    }();

    h_CBTAPS->Fill(taps);
}

AUTO_REGISTER_PHYSICS(MCTrueOverview)
