#include "XMasCB.h"
#include "base/cbtaps_display/TH2CB.h"
#include "base/std_ext/memory.h"
#include "base/std_ext/string.h"
#include "TCanvas.h"

using namespace ant;
using namespace ant::analysis::physics;
using namespace ant::analysis;
using namespace ant::analysis::data;
using namespace std;

XMasCB::XMasCB(const std::string& name, PhysOptPtr opts):
    Physics("XMasCB", opts)
{
    hist = std_ext::make_unique<TH2CB>("","");
    c    = new TCanvas("xmascb", "XMasCB");
    ext  = opts->Get<string>("FileType", "pdf");
}

XMasCB::~XMasCB()
{}

void XMasCB::ProcessEvent(const Event& event)
{
    hist->ResetElements(0.0);

    for(const auto& c : event.Reconstructed().AllClusters()) {
        if(c.Detector & Detector_t::Type_t::CB) {
            for(const auto& hit : c.Hits) {
                for(const auto& datum : hit.Data) {
                    if(datum.Type == Channel_t::Type_t::Integral) {
                        hist->SetElement(hit.Channel, hist->GetElement(hit.Channel)+datum.Value);
                    }
                }
            }
        }
    }

    c->cd();
    c->SetLogz();

    hist->Draw("col");

    std_ext::formatter f;
    f << "xmas_cb_" << n++ << "." << ext;

    c->SaveAs(f.str().c_str());

}

AUTO_REGISTER_PHYSICS(XMasCB)
