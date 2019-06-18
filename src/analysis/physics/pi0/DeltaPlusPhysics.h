#pragma once

#include "physics/Physics.h"
#include "base/interval.h"
#include <map>
#include <string>

class TH1;

namespace ant {
namespace analysis {
namespace physics {

class DeltaPlusPhysics: public Physics {

protected:

    // Class that groups histograms together
    class Histogm {

    public:
        std::string  pref;  // prefix to label whole group of histograms
        mutable std::map<std::string, TH1* > h; // container for histograms by name (without prefix)

        Histogm(const Histogm&) = default;

        Histogm(HistogramFactory HistFac);

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

    const LorentzVec target;

public:
    DeltaPlusPhysics(const std::string& name, OptionsPtr opts);
    virtual ~DeltaPlusPhysics() {}


    // Physics interface
public:
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};

}
}
}
