#include "TaggerScalers.h"

#include <numeric>

#include "SlowControlProcessors.h"

#include "expconfig/ExpConfig.h"


using namespace std;
using namespace ant::analysis::slowcontrol;
using namespace ant::analysis::slowcontrol::variable;

void TaggerScalers::Init()
{
    auto taggerdetector = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    if(!taggerdetector)
        return;
    nChannels = taggerdetector->GetNChannels();

    if(taggerdetector->Type == Detector_t::Type_t::EPT) {
        mode = mode_t::EPT_2014;
    }
    else if(taggerdetector->Type == Detector_t::Type_t::Tagger) {
        mode = mode_t::Tagger;
    }
}

list<Variable::ProcessorPtr> TaggerScalers::GetNeededProcessors() const
{
    if(mode == mode_t::EPT_2014) {
        /// \todo check how EPT_2012 scalers were recorded
        /// for 2014, we know that EPT_Scalers are in Beampolmon VUPROMs
        return {Processors::EPT_Scalers, Processors::Beampolmon, Processors::EPT_Or};
    }

    if(mode == mode_t::Tagger) {
        return {Processors::Tagger_Scalers, Processors::Beampolmon, Processors::Tagger_Or};
    }
    return {};
}

std::vector<double> TaggerScalers::Get() const
{
    vector<double> scalers(nChannels, std::numeric_limits<double>::quiet_NaN());
    if(mode == mode_t::EPT_2014) {
        const double reference = Processors::Beampolmon->Reference_1MHz.Get();
        for(const auto& kv : Processors::EPT_Scalers->Get()) {
            if(kv.Key<scalers.size())
                scalers[kv.Key] = 1.0e6*kv.Value/reference;
        }
    }

    // Here I would add in the access to Tagger_scalers...

    return scalers;
}

std::vector<int64_t> TaggerScalers::GetCounts() const
{
    vector<int64_t> counts(nChannels, std::numeric_limits<int64_t>::quiet_NaN());
    if (mode == mode_t::EPT_2014)
    {
        for (const auto& kv: Processors::EPT_Scalers->Get())
        {
            if (kv.Key < counts.size())
                counts[kv.Key] = kv.Value;
        }
    }

    // Here I would add in the access to Tagger_scalers...

    return counts;
}

double TaggerScalers::GetTaggerOr() const
{
    return Processors::EPT_Or->Get() * 1.0e6 / Processors::Beampolmon->Reference_1MHz.Get();

}

double TaggerScalers::GetTaggerOrAsSum() const
{
    auto scalers = Get();
    return accumulate(scalers.begin(),scalers.end(),0);
}



