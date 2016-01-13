#include "DetectorPlots.h"
#include "../cbtaps_display/TH2CB.h"
#include "../cbtaps_display/TH2TAPS.h"
#include "expconfig/ExpConfig.h"
#include "base/std_ext/math.h"
#include <iostream>

#include "TCanvas.h"

using namespace std;
using namespace ant;


TH2TAPS* DetectorPlots::MakeTAPSTheta(const string& setup_name)
{
    auto taps = new TH2TAPS();

    ExpConfig::Setup::ManualName = setup_name;

    const auto setup = ExpConfig::Setup::GetLastFound();

    if(!setup) {
        cerr << "Setup \"" << setup_name << "\" not found!" << endl;
        return nullptr;
    }

    const auto det = setup->GetDetector(Detector_t::Type_t::TAPS);

    if(!det) {
        cerr << "No TAPS detector definied in Setup!" << endl;
        return nullptr;
    }

    for(unsigned ch = 0; ch < det->GetNChannels(); ++ch) {
        taps->SetElement(ch, std_ext::radian_to_degree(det->GetPosition(ch).Theta()));
    }

    taps->SetZTitle("Element #theta [#circ]");

    return taps;

}

TH2TAPS* DetectorPlots::MakeTAPSPhi(const string& setup_name)
{
    auto taps = new TH2TAPS();

    ExpConfig::Setup::ManualName = setup_name;

    const auto setup = ExpConfig::Setup::GetLastFound();

    if(!setup) {
        cerr << "Setup \"" << setup_name << "\" not found!" << endl;
        return nullptr;
    }

    const auto det = setup->GetDetector(Detector_t::Type_t::TAPS);

    if(!det) {
        cerr << "No TAPS detector definied in Setup!" << endl;
        return nullptr;
    }

    for(unsigned ch = 0; ch < det->GetNChannels(); ++ch) {
        taps->SetElement(ch, std_ext::radian_to_degree(det->GetPosition(ch).Phi()));
    }

    taps->SetZTitle("Element #phi [#circ]");

    return taps;

}



TH2CB* DetectorPlots::MakeCBTheta(const string& setup_name)
{
    auto cb = new TH2CB();

    ExpConfig::Setup::ManualName = setup_name;

    const auto setup = ExpConfig::Setup::GetLastFound();

    if(!setup) {
        cerr << "Setup \"" << setup_name << "\" not found!" << endl;
        return nullptr;
    }

    const auto det = setup->GetDetector(Detector_t::Type_t::CB);

    if(!det) {
        cerr << "No CB detector definied in Setup!" << endl;
        return nullptr;
    }

    for(unsigned ch = 0; ch < det->GetNChannels(); ++ch) {
        cb->SetElement(ch, std_ext::radian_to_degree(det->GetPosition(ch).Theta()));
    }

    cb->SetZTitle("Element #theta [#circ]");

    return cb;

}

TH2CB* DetectorPlots::MakeCBPhi(const string& setup_name)
{
    auto cb = new TH2CB();

    ExpConfig::Setup::ManualName = setup_name;

    const auto setup = ExpConfig::Setup::GetLastFound();

    if(!setup) {
        cerr << "Setup \"" << setup_name << "\" not found!" << endl;
        return nullptr;
    }

    const auto det = setup->GetDetector(Detector_t::Type_t::CB);

    if(!det) {
        cerr << "No CB detector definied in Setup!" << endl;
        return nullptr;
    }

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
