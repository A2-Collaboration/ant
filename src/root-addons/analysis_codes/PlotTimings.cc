
#include "PlotTimings.h"

#include "analysis/plot/root_draw.h"

#include "TTree.h"
#include "TH1D.h"

#include <vector>
#include <iostream>
#include <limits>

using namespace std;
using namespace ant;

void PlotTimings::Plot(TTree* tree)
{
    if(!tree)
        return;

    vector<unsigned>* ADC_index = nullptr;
    vector<unsigned>* ADC_value = nullptr;

    tree->SetBranchAddress("ADC_index",addressof(ADC_index));
    tree->SetBranchAddress("ADC_value",addressof(ADC_value));

    cout << "Tree Entries: " << tree->GetEntries() << endl;

    TH1D* h =     new TH1D("h","",24000,-12000,12000);
    TH1D* h_raw = new TH1D("h_raw","",24000,-12000,12000);

    unsigned multiple_refhits = 0;

    for(long long entry=0;entry<tree->GetEntries();entry++) {
        tree->GetEntry(entry);
        // search for ref hit
        const unsigned refhit = 29192;
        double refhit_val = numeric_limits<double>::quiet_NaN();
        for(size_t i=0;i<ADC_index->size();i++) {
            if(ADC_index->at(i) == refhit) {
                if(isfinite(refhit_val))
                    multiple_refhits++;
                refhit_val = ADC_value->at(i);
            }
        }
        if(!isfinite(refhit_val)) {
            continue;
        }
        for(size_t i=0;i<ADC_index->size();i++) {
            if(ADC_index->at(i) == refhit)
                continue;

            double val = ADC_value->at(i);
            h->Fill(val-refhit_val);
            h_raw->Fill(val);

        }
    }
    cout << "MultiRefHits: " << multiple_refhits << endl;

    canvas() << h << h_raw << endc;

}
