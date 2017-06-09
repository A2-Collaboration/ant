#include "catch.hpp"

#include "root-addons/analysis_codes/hadd.h"
#include "base/WrapTFile.h"
#include "base/tmpfile_t.h"
#include "base/std_ext/memory.h"

#include "TH1D.h"

using namespace std;
using namespace ant;

TEST_CASE("Hadd: Two simple hists", "[root-addons]") {

    // create two input files with histograms
    tmpfile_t in_file1;
    {
        WrapTFileOutput out(in_file1.filename);
        {
            auto h = out.CreateInside<TH1D>("h","",100,0,1);
            h->Sumw2();
            h->Fill(0.0, 3.0);
        }
        {
            auto h = out.CreateInside<TH1D>("h_avg","",100,0,1);
            h->SetBit(TH1::kIsAverage);
            h->SetBinContent(10, 3.0);
            h->SetBinError(10, 1.0);
            h->SetBinContent(20, 5.0);
            h->SetBinError(20, 1.0);
        }
        {
            auto h = out.CreateInside<TH1D>("h_lbl","",1,0,1);
            h->Fill("a", 3.0);
            h->Fill("b", 7.0);
        }
    }
    tmpfile_t in_file2;
    {
        WrapTFileOutput out(in_file2.filename);
        {
            auto h = out.CreateInside<TH1D>("h","",100,0,1);
            h->Sumw2();
            h->Fill(0.0, 6.0);
        }
        {
            auto h = out.CreateInside<TH1D>("h_avg","",100,0,1);
            h->SetBit(TH1::kIsAverage);
            h->SetBinContent(10, 6.0);
            h->SetBinError(10, 1.0);
            h->SetBinContent(20, 1.0);
            h->SetBinError(20, 1.0);
        }
        {
            auto h = out.CreateInside<TH1D>("h_lbl","",1,0,1);
            h->Fill("b", 10.0);
            h->Fill("a", 17.0);
        }
        {
            auto h = out.CreateInside<TH1D>("h_single","",100,0,1);
            h->Fill(1.0, 6.0); // fill overflow
        }

    }

    // merge files
    tmpfile_t tmp_outfile;
    {
        auto outputfile = std_ext::make_unique<TFile>(tmp_outfile.filename.c_str(), "RECREATE");
        hadd::sources_t sources;
        sources.emplace_back(std_ext::make_unique<TFile>(in_file1.filename.c_str(), "READ"));
        sources.emplace_back(std_ext::make_unique<TFile>(in_file2.filename.c_str(), "READ"));
        unsigned nPaths = 0;
        hadd::MergeRecursive(*outputfile, sources, nPaths);
        outputfile->Write();
        CHECK(nPaths == 1);
    }


    // test correct merging
    {
        WrapTFileInput input(tmp_outfile.filename);
        {
            auto h = input.GetSharedHist<TH1D>("h");
            CHECK(h->GetNbinsX() == 100);
            CHECK(h->GetBinContent(1) == Approx(9.0));
        }
        {
            auto h = input.GetSharedHist<TH1D>("h_avg");
            CHECK(h->GetNbinsX() == 100);
            CHECK(h->GetBinContent(10) == Approx(4.5));
            CHECK(h->GetBinContent(20) == Approx(3.0));
        }
        {
            auto h = input.GetSharedHist<TH1D>("h_lbl");
            CHECK(h->GetNbinsX() == 2);
            CHECK(h->GetBinContent(h->GetXaxis()->FindBin("a")) == Approx(20.0));
            CHECK(h->GetBinContent(h->GetXaxis()->FindBin("b")) == Approx(17.0));
        }
        {
            auto h = input.GetSharedHist<TH1D>("h_single");
            CHECK(h->GetNbinsX() == 100);
            CHECK(h->GetBinContent(101) == Approx(6.0));
        }
    }

}