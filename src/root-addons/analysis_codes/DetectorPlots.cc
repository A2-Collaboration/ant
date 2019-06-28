#include "DetectorPlots.h"
#include "analysis/plot/HistogramFactory.h"
#include "../cbtaps_display/TH2CB.h"
#include "../cbtaps_display/TH2TAPS.h"
#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/CB.h"
#include "base/std_ext/math.h"
#include <iostream>

#include "TCanvas.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"

using namespace std;
using namespace ant;


/**
 * @brief create an int representation from a Detector_t::ElementFlags_t
 * @param f
 * @return
 */
inline int getIntFromFlags(const Detector_t::ElementFlags_t& f) noexcept {
    int i = 0;
    if(f.test(Detector_t::ElementFlag_t::BadTDC)) i+=1;
    if(f.test(Detector_t::ElementFlag_t::Broken)) i+=2;
    if(f.test(Detector_t::ElementFlag_t::Missing)) i+=4;
    if(f.test(Detector_t::ElementFlag_t::NoCalibFill)) i+=8;
    if(f.test(Detector_t::ElementFlag_t::NoCalibUseDefault)) i+=16;
    return i;
};

TH2TAPS* DetectorPlots::MakeTAPSTheta(const string& setup_name)
{
    ExpConfig::Setup::SetByName(setup_name);
    const auto det = ExpConfig::Setup::GetDetector(Detector_t::Type_t::TAPS);

    auto taps = new TH2TAPS();

    for(unsigned ch = 0; ch < det->GetNChannels(); ++ch) {
        taps->SetElement(ch, std_ext::radian_to_degree(det->GetPosition(ch).Theta()));
    }

    taps->SetZTitle("Element #theta [#circ]");

    return taps;

}

TH2TAPS* DetectorPlots::MakeTAPSPhi(const string& setup_name)
{
    ExpConfig::Setup::SetByName(setup_name);
    const auto det = ExpConfig::Setup::GetDetector(Detector_t::Type_t::TAPS);

    auto taps = new TH2TAPS();

    for(unsigned ch = 0; ch < det->GetNChannels(); ++ch) {
        taps->SetElement(ch, std_ext::radian_to_degree(det->GetPosition(ch).Phi()));
    }

    taps->SetZTitle("Element #phi [#circ]");

    return taps;

}



TH2CB* DetectorPlots::MakeCBTheta(const string& setup_name)
{
    ExpConfig::Setup::SetByName(setup_name);
    const auto det = ExpConfig::Setup::GetDetector(Detector_t::Type_t::CB);

    auto cb = new TH2CB();

    for(unsigned ch = 0; ch < det->GetNChannels(); ++ch) {
        cb->SetElement(ch, std_ext::radian_to_degree(det->GetPosition(ch).Theta()));
    }

    cb->SetZTitle("Element #theta [#circ]");

    return cb;

}

TH2CB* DetectorPlots::MakeCBPhi(const string& setup_name)
{
    ExpConfig::Setup::SetByName(setup_name);
    const auto det = ExpConfig::Setup::GetDetector(Detector_t::Type_t::CB);

    auto cb = new TH2CB();

    for(unsigned ch = 0; ch < det->GetNChannels(); ++ch) {
        cb->SetElement(ch, std_ext::radian_to_degree(det->GetPosition(ch).Phi()));
    }

    cb->SetZTitle("Element #phi [#circ]");

    return cb;

}

void DetectorPlots::PlotCBTheta(const string& setup_name)
{
    new TCanvas();

    auto cb   = MakeCBTheta(setup_name);
    auto grid = new TH2CB();
    grid->FillElementNumbers();

    cb->Draw("colz");
    grid->Draw("same text");

}

void DetectorPlots::PlotCBPhi(const string& setup_name)
{
    new TCanvas();

    auto cb   = MakeCBPhi(setup_name);
    auto grid = new TH2CB();
    grid->FillElementNumbers();

    cb->Draw("colz");
    grid->Draw("same text");
}

void DetectorPlots::PlotCBIgnored(const string& setup_name, bool draw_element_numbers)
{
    ExpConfig::Setup::SetByName(setup_name);
    const auto det = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();

    auto cb   = new TH2CB();
    unsigned total = 0;

    for(unsigned ch=0;ch<det->GetNChannels();ch++) {
        cb->SetElement(ch, getIntFromFlags(det->GetElementFlags(ch)));
        if(det->IsIgnored(ch) && !det->IsHole(ch)) {
            total++;
        }
    }
    cout << "CB total ignored: " << total << endl;




    new TCanvas();
    cb->Draw("col");
    if(draw_element_numbers) {
        auto grid = new TH2CB();
        grid->FillElementNumbers();
        grid->Draw("same text");
    }
}

void DetectorPlots::PlotCBNearestAngles(const string& setup_name)
{
    ExpConfig::Setup::SetByName(setup_name);
    const auto det = ExpConfig::Setup::GetDetector(Detector_t::Type_t::CB);

    ant::analysis::HistogramFactory HistFac("DetectorPlots");
    auto h = HistFac.makeTH1D("Angle between pairs of channels", {"#Delta#alpha / #circ",{30,{0,30}}}, "h_CB_OpAngle");

    for(unsigned ch1=0;ch1<det->GetNChannels();ch1++) {
        for(unsigned ch2=ch1+1;ch2<det->GetNChannels();ch2++) {
            auto pos1 = det->GetPosition(ch1);
            auto pos2 = det->GetPosition(ch2);
            h->Fill(std_ext::radian_to_degree(pos1.Angle(pos2)));
        }
    }

    h->Draw();
}

void DetectorPlots::PlotTAPSTheta(const string& setup_name)
{
    new TCanvas();

    auto taps   = MakeTAPSTheta(setup_name);
    auto grid = new TH2TAPS();
    grid->FillElementNumbers();

    taps->Draw("colz");
    grid->Draw("same text");

}

void DetectorPlots::PlotTAPSPhi(const string& setup_name)
{
    new TCanvas();

    auto taps   = MakeTAPSPhi(setup_name);
    auto grid = new TH2TAPS();
    grid->FillElementNumbers();

    taps->Draw("colz");
    grid->Draw("same text");
}

void DetectorPlots::PlotTAPSIgnored(const string& setup_name, bool draw_element_numbers)
{
    ExpConfig::Setup::SetByName(setup_name);
    const auto det = ExpConfig::Setup::GetDetector(Detector_t::Type_t::TAPS);

    auto taps   = new TH2TAPS();
    unsigned total = 0;

    for(unsigned ch=0;ch<det->GetNChannels();ch++) {
        taps->SetElement(int(ch), getIntFromFlags(det->GetElementFlags(ch)));
        if(det->IsIgnored(ch)) {
            total++;
        }
    }
    cout << "TAPS total ignored: " << total << endl;



    new TCanvas();
    taps->Draw("col");
    if(draw_element_numbers) {
        auto grid = new TH2TAPS();
        grid->FillElementNumbers();
        grid->Draw("same text");
    }
}

void DetectorPlots::PlotTAPSNearestAngles(const string& setup_name)
{
    ExpConfig::Setup::SetByName(setup_name);
    const auto det = ExpConfig::Setup::GetDetector(Detector_t::Type_t::TAPS);

    ant::analysis::HistogramFactory HistFac("DetectorPlots");
    auto h = HistFac.makeTH1D("Angle between pairs of channels", {"#Delta#alpha / #circ",{100,{0,50}}}, "h_TAPS_OpAngle");

    for(unsigned ch1=0;ch1<det->GetNChannels();ch1++) {
        for(unsigned ch2=ch1+1;ch2<det->GetNChannels();ch2++) {
            auto pos1 = det->GetPosition(ch1);
            auto pos2 = det->GetPosition(ch2);
            h->Fill(std_ext::radian_to_degree(pos1.Angle(pos2)));
        }
    }

    h->Draw();
}

void DetectorPlots::PlotCBTAPSDetectorPositions(const string& setup_name, double CB_gap)
{
    ExpConfig::Setup::SetByName(setup_name);
    const auto cb = ExpConfig::Setup::GetDetector(Detector_t::Type_t::CB);
    const auto taps = ExpConfig::Setup::GetDetector(Detector_t::Type_t::TAPS);

    auto graph = new TGraph2D(cb->GetNChannels()+taps->GetNChannels());

    for(unsigned ch=0;ch<cb->GetNChannels();ch++) {
        auto pos = cb->GetPosition(ch);
        if(pos.y<0)
            pos.y -= CB_gap/2;
        else
            pos.y += CB_gap/2;
        graph->SetPoint(ch,pos.x,pos.y,pos.z);
    }

    for(unsigned ch=0;ch<taps->GetNChannels();ch++) {
        const auto& pos = taps->GetPosition(ch);
        graph->SetPoint(cb->GetNChannels()+ch,pos.x,pos.y,pos.z);
    }


    graph->GetXaxis()->SetTitle("X");
    graph->GetYaxis()->SetTitle("Y");
    graph->GetZaxis()->SetTitle("Z");
    graph->Draw("P0");
}

void DetectorPlots::PlotTaggerChannelEnergy(const string& setup_name)
{
    ExpConfig::Setup::SetByName(setup_name);
    auto tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();

    const auto n = tagger->GetNChannels();
    auto g = new TGraphErrors(int(n));

    for(unsigned i=0; i<n; ++i) {
        g->SetPoint(     int(i), i, tagger->GetPhotonEnergy(i));
        g->SetPointError(int(i), 0, tagger->GetPhotonEnergyWidth(i));
    }

    new TCanvas();
    g->Draw("AP");
}

void DetectorPlots::PlotTaggerChannelEnergyWidth(const string& setup_name)
{
    ExpConfig::Setup::SetByName(setup_name);
    auto tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();

    const auto n = tagger->GetNChannels();
    auto g = new TGraph(int(n));

    for(unsigned i=0; i<n; ++i) {
        g->SetPoint(int(i), i, tagger->GetPhotonEnergyWidth(i));
    }

    new TCanvas();
    g->Draw("AP");
}
