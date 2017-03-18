#include "DetectorPlots.h"
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

void DetectorPlots::PlotCBIgnored(const string& setup_name)
{
    ExpConfig::Setup::SetByName(setup_name);
    const auto det = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();

    auto cb   = new TH2CB();
    unsigned total = 0;
    for(unsigned ch=0;ch<det->GetNChannels();ch++) {
        if(det->IsIgnored(ch) && !det->IsHole(ch)) {
            cb->SetElement(ch, 1.0);
            total++;
        }
    }
    cout << "CB total ignored: " << total << endl;


    auto grid = new TH2CB();
    grid->FillElementNumbers();

    new TCanvas();
    cb->Draw("col");
    grid->Draw("same text");
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

void DetectorPlots::PlotTAPSIgnored(const string& setup_name)
{
    ExpConfig::Setup::SetByName(setup_name);
    const auto det = ExpConfig::Setup::GetDetector(Detector_t::Type_t::TAPS);

    auto taps   = new TH2TAPS();
    unsigned total = 0;
    for(unsigned ch=0;ch<det->GetNChannels();ch++) {
        if(det->IsIgnored(ch)) {
            taps->SetElement(ch, 1.0);
            total++;
        }
    }
    cout << "TAPS total ignored: " << total << endl;

    auto grid = new TH2TAPS();
    grid->FillElementNumbers();

    new TCanvas();
    taps->Draw("col");
    grid->Draw("same text");
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
