#ifndef DELTAPLUSPHYSICS_H
#define DELTAPLUSPHYSICS_H

#include "AntPhysics.h"
#include "plot/Histogram.h"
#include "base/interval.h"
#include <map>
#include <string>
#include "TLorentzVector.h"
class TH1;

namespace ant {
namespace analysis {

class DeltaPlusPhysics: public ant::Physics {

protected:

    // Class that groups histograms together
    class Histogm {

    public:
        std::string  pref;  // prefix to label whole group of histograms
        mutable std::map<std::string, TH1* > h; // container for histograms by name (without prefix)
        std::map<std::string, std::string> h_title; // container for histogram titles by name (without prefix)

        // Add 1D histogram
        void AddHistogram(const std::string& name,       // short specifier for histogram
                          const std::string& title,      // descriptive title for histogram
                          const std::string& x_label,    // x axis label
                          const std::string& y_label,    // y axis label
                          const int x_bins_n,       // number of bins in x
                          const double x_bins_low,  // lower bound of x axis
                          const double x_bins_up    // upper bound of x axis
                          );

        // Add 2D histogram
        void AddHistogram(const std::string& name,       // short specifier for histogram
                          const std::string& title,      // descriptive title for histogram
                          const std::string& x_label,    // x axis label
                          const std::string& y_label,    // y axis label
                          const int x_bins_n,       // number of bins in x
                          const double x_bins_low,  // lower bound of x axis
                          const double x_bins_up,   // upper bound of y axis
                          const int y_bins_n,       // number of bins in y
                          const double y_bins_low,  // lower bound of y axis
                          const double y_bins_up    // upper bound of y axis
                          );

        Histogm( const std::string& prefix );

        void Draw();

        Histogm& operator*= ( const Double_t factor );

        Histogm operator= (const Histogm& other);

        void AddScaled( const Histogm& h2, const Double_t f=1.0 );

        TH1* operator[] (const std::string& key) {
            return h[key];
        }

        const TH1* operator[] (const std::string& key) const {
            return h[key];
        }
    };

    Histogm prompt;
    Histogm random;
    Histogm diff;

    interval<double> pi0_cut;
    interval<double> prompt_window;
    interval<double> random_window;

    const TLorentzVector target;

public:
    DeltaPlusPhysics(const std::string& name="DeltaPlusPhysics");
    virtual ~DeltaPlusPhysics() {}


    // Physics interface
public:
    void ProcessEvent(const ant::Event &event);
    void Finish();
    void ShowResult();
};

}
}
#endif
