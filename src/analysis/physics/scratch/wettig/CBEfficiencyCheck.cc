#include "CBEfficiencyCheck.h"
#include "root-addons/cbtaps_display/TH2CB.h"
#include "TCanvas.h"
#include "TColor.h"
#include "base/std_ext/memory.h"
#include "base/std_ext/string.h"
#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


CBEfficiencyCheck::CBEfficiencyCheck(const std::string& name, OptionsPtr opts):
    Physics(name,opts),
    w_px(opts->Get<int>("x_px", 1920))
{
    canv  = new TCanvas("CBeffCheck", "CBeffCheck");
    hist_hits = std_ext::make_unique<TH2CB>("", "", false);
    hist_hits->SetTitle("All Hits of all Clusters");

    hist_centClust = std_ext::make_unique<TH2CB>("", "", false);
    hist_centClust->SetTitle("Central Hits of all Clusters");

    grid = std_ext::make_unique<TH2CB>("", "", false);
    grid->SetTitle("");
}

void CBEfficiencyCheck::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& cands = event.Reconstructed().Candidates;

    for(const auto& c : cands) {

        if(c.Detector & Detector_t::Type_t::CB) {
            const auto& cluster  = c.FindCaloCluster();
            const auto centEl = cluster->CentralElement;

            hist_centClust->SetElement(centEl,hist_centClust->GetElement(centEl)+1);

            for(const auto& hit : cluster->Hits) {
//                hist_hits->SetElement(hit.Channel, hist_hits->GetElement(hit.Channel)+hit.Energy);
                hist_hits->SetElement(hit.Channel, hist_hits->GetElement(hit.Channel)+1);
            }
        }
    }
}

void CBEfficiencyCheck::Finish()
{
    canv->cd();
    canv->SetLogz();
    canv->SetMargin(0,0,0,5);

    const double ratio = 2;

    canv->SetCanvasSize(UInt_t(w_px), unsigned(w_px/ratio));

    // to get rid of the ugly black frame around the plot
    canv->cd(1);
    gPad->SetFrameLineColor(0);
    hist_hits->SetAxisColor(0);
    hist_hits->GetYaxis()->SetAxisColor(0);
    hist_centClust->SetAxisColor(0);
    hist_centClust->GetYaxis()->SetAxisColor(0);



    hist_hits->Draw("col");
    grid->Draw("text same");

    std_ext::formatter f, f2;
    f << "cbeff_allhits.pdf";
    canv->SaveAs(f.str().c_str());

    hist_centClust->Draw("col");
    grid->Draw("text same");

    f2 << "cbeff_centralclusters.pdf";
    canv->SaveAs(f2.str().c_str());
}

void CBEfficiencyCheck::ShowResult()
{

}


AUTO_REGISTER_PHYSICS(CBEfficiencyCheck)
