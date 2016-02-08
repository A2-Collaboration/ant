#include "FindProton.h"

#include "TTree.h"

using namespace ant;
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
}



FindProton::FindProton(const string& name, OptionsPtr opts):
    Physics(name,opts)
{
    auto tree = HistFac.makeTTree("tree");
    branches.SetBranchtes(tree);
}

FindProton::~FindProton()
{

}

void FindProton::ProcessEvent(const TEvent& event, Physics::manager_t& manager)
{



}


AUTO_REGISTER_PHYSICS(FindProton)
