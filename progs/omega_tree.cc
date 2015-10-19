#include "base/WrapTFile.h"
#include "base/Logger.h"
#include "base/std_ext/math.h"

#include "TRint.h"
#include "analysis/plot/root_draw.h"
#include <string>
#include <list>
#include "TH1D.h"
#include "TTree.h"
#include "base/Logger.h"
#include "analysis/plot/HistogramFactories.h"
#include "TCutG.h"
#include "TEventList.h"

using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace std;

template <typename T>
bool GetObjectFromFile(const std::string& filename, const std::string& objname, T* obj) {
    return true;
}

void CopyBins(const TH1* src, TH1* dest) {
    for(int i=0; i<src->GetNbinsX(); ++i) {
        dest->SetBinContent(i,src->GetBinContent(i));
    }
}

void FillggIM(TH1* hist, const double v[3]) {
    hist->Fill(v[0]);
    hist->Fill(v[1]);
    hist->Fill(v[2]);
}

int main(int argc, char** argv) {

    if(argc != 2)
        return EXIT_FAILURE;

    WrapTFileInput infile(argv[1]);
    try {
        infile.OpenFile("eventlist.root");
    } catch (...) {}

    try {
    infile.OpenFile("cut1.root");
    } catch (...) {}

    try {
    infile.OpenFile("cut2.root");
    } catch (...) {}

    TTree* tree = nullptr;
    if (!infile.GetObject("omegaetag2", tree)) {
        LOG(ERROR) << "Tree not found" << endl;
    }

    TEventList* eventlist = nullptr;
    if( infile.GetObject("eventlist", eventlist)) {
        tree->SetEventList(eventlist);
        LOG(INFO) << "Event list set";
    }

    TCutG* cut1 = nullptr;
    if( !infile.GetObject("cut1", cut1)) {
        LOG(INFO) << "Cut1 not found";
    }

    TCutG* cut2 = nullptr;
    if( !infile.GetObject("cut2", cut2)) {
        LOG(INFO) << "Cut2 not found";
    }


    double calcpIM = 0.0;
    double gggIM = 0.0;
    double ggIM[3] = {};
    unsigned rf = 2;

    tree->SetBranchAddress("calcpIM", &calcpIM);
    tree->SetBranchAddress("gggIM",   &gggIM);
    tree->SetBranchAddress("ggIM",     ggIM);
    tree->SetBranchAddress("rf",      &rf);

    WrapTFileOutput("out.root", WrapTFileOutput::mode_t::recreate, true);

    SmartHistFactory hf("omega_tree");

    const BinSettings gggIMcalcPIMx(350,600,850);
    const BinSettings gggIMcalcPIMy(350,850,1200);

    auto gggIM_calcpIM = hf.makeTH2D(
                "3 #gamma IM vs calculated proton IM",
                "3#gamma IM [MeV]",
                "calculated proton IM [MeV]",
                gggIMcalcPIMx,
                gggIMcalcPIMy
                );
    const BinSettings bins_ggIM(900);

    auto h_ggIM = hf.makeTH1D(
                "2 #gamma IM (no Cuts)",
                "2#gamma IM [MeV]",
                "",
                bins_ggIM
                );

    auto h_ggIM_cut1 = hf.makeTH1D(
                "2 #gamma IM (cut1)",
                "2#gamma IM [MeV]",
                "",
                bins_ggIM
                );

    auto h_ggIM_cut2 = hf.makeTH1D(
                "2 #gamma IM (cut2)",
                "2#gamma IM [MeV]",
                "",
                bins_ggIM
                );

    auto h_ggIM_cut1_b = hf.makeTH1D(
                "2 #gamma IM (cut1, bg)",
                "2#gamma IM [MeV]",
                "",
                bins_ggIM
                );
    h_ggIM_cut1_b->SetLineColor(kRed);

    auto h_ggIM_cut2_b = hf.makeTH1D(
                "2 #gamma IM (cut2, bg)",
                "2#gamma IM [MeV]",
                "",
                bins_ggIM
                );
    auto h_ggIM_cut1_s = hf.makeTH1D(
                "2 #gamma IM (cut1, sig)",
                "2#gamma IM [MeV]",
                "",
                bins_ggIM
                );
    h_ggIM_cut1_s->SetLineColor(kGreen);

    auto h_ggIM_cut2_s = hf.makeTH1D(
                "2 #gamma IM (cut2, sig)",
                "2#gamma IM [MeV]",
                "",
                bins_ggIM
                );

    TH1D* h_ggIM_scaledbg = hf.makeTH1D(
                "2 #gamma IM (cut2, bg, scaled)",
                "2#gamma IM [MeV]",
                "",
                bins_ggIM
                );


    TH1D* h_ggIM_cleaned = hf.makeTH1D(
                "2 #gamma IM (cut1, cleaned)",
                "2#gamma IM [MeV]",
                "",
                bins_ggIM
                );
    const TVector2 center(735, 1025);

    const BinSettings gggIMcalcpIMrotx(gggIMcalcPIMx.Bins(), interval<double>(gggIMcalcPIMx) - center.X());
    const BinSettings gggIMcalcpIMroty(gggIMcalcPIMy.Bins(), interval<double>(gggIMcalcPIMy) - center.Y());

    auto gggIM_calcpIM_rotated = hf.makeTH2D(
                "3 #gamma IM vs calculated proton IM (45#circ)",
                "x",
                "y",
                gggIMcalcpIMrotx,
                gggIMcalcpIMroty
                );

    const auto entries = eventlist ? eventlist->GetN() : tree->GetEntries();

    LOG(INFO) << entries << " Entries";

    for(Long64_t entry=0; entry<entries; ++entry) {

        if(eventlist) {
            tree->GetEntry(eventlist->GetEntry(entry));
        } else {
            tree->GetEntry(entry);
        }

        gggIM_calcpIM->Fill(gggIM, calcpIM);

        TVector2 x(gggIM, calcpIM);
        x -= center;
        auto r = x.Rotate(-degree_to_radian(45.0));

        gggIM_calcpIM_rotated->Fill(r.X(),r.Y());

        FillggIM(h_ggIM, ggIM);

        if(cut1 && cut2) {
            const auto cut1_pass = cut1->IsInside(gggIM, calcpIM);
            const auto cut2_pass = cut2->IsInside(gggIM, calcpIM);
            const auto cut2not1  = cut2_pass && !cut1_pass;
            const auto is_sigerf = rf != 2;

            if(cut1_pass) {
                FillggIM(h_ggIM_cut1,ggIM);
                FillggIM( is_sigerf ? h_ggIM_cut1_s : h_ggIM_cut1_b ,ggIM);
            }

            if(cut2not1) {
                FillggIM(h_ggIM_cut2,ggIM);
                FillggIM( is_sigerf ? h_ggIM_cut2_s : h_ggIM_cut2_b ,ggIM);
            }
        }

    }

    const auto area1 = cut1->Area();
    const auto area2 = cut2->Area();
    const auto area2not1 = area2 - area1;
    const auto ratio = area1 / area2not1;
    LOG(INFO) << "Ratio " << ratio;

    CopyBins(h_ggIM_cut2_b, h_ggIM_scaledbg);
    h_ggIM_scaledbg->Scale(ratio);
    h_ggIM_scaledbg->SetLineColor(kBlue);

    CopyBins(h_ggIM_cut1, h_ggIM_cleaned);
    h_ggIM_cleaned->Add(h_ggIM_scaledbg, -1);

    int fake_argc=0;
    char** fake_argv=nullptr;
    TRint app("Ant",&fake_argc,fake_argv,nullptr,0,true);

    canvas c("Omega Tree");
    c << drawoption("colz")
      << gggIM_calcpIM
      << endc;

    if(cut1) cut1->Draw("same");
    if(cut2) cut2->Draw("same");

    canvas cr("Omega Tree: Rotated");
    cr << drawoption("colz")
      << gggIM_calcpIM_rotated
      << endc;


    canvas c_ggIM("Omega Tree: ggIM");

    c_ggIM << h_ggIM << h_ggIM_cut1 << h_ggIM_cut2 << endc;

    canvas c_ggIM2("Omega Tree: ggIM 2");

    hstack s1("ggIMs", "ggIMs");
    s1 << h_ggIM_scaledbg << h_ggIM_cut1_b << h_ggIM_cut1_s << h_ggIM_cut1;

    c_ggIM2 << h_ggIM_cut1_b << h_ggIM_cut2_b
            << h_ggIM_cut1_s << h_ggIM_cut2_s
            << drawoption("nostack") << s1
            << h_ggIM_cleaned << endc;




    app.Run(kTRUE);

    return EXIT_SUCCESS;
}
