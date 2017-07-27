#include "TaggerScalers.h"

#include <numeric>

#include "SlowControlProcessors.h"

#include "expconfig/ExpConfig.h"


using namespace std;
using namespace ant::analysis::slowcontrol;
using namespace ant::analysis::slowcontrol::variable;

void TaggerScalers::Init(const input::reader_flags_t& reader_flags)
{
    Variable::Init(reader_flags);

    auto taggerdetector = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
    nChannels = taggerdetector->GetNChannels();

    if(taggerdetector->Type == Detector_t::Type_t::EPT) {
        mode = mode_t::EPT_2014;
    }
    else if(taggerdetector->Type == Detector_t::Type_t::Tagger) {
        mode = mode_t::Tagger;
    }
    else {
        throw runtime_error("Tagger type not implemented for TaggerScalers Slowcontrol variable");
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

std::vector<double> TaggerScalers::GetRates() const
{
    if(!slowcontrol_provided)
        return vector<double>(nChannels, 1.0);

    const auto counts = GetCounts();
    vector<double> rates(counts.size());
    std::transform(counts.begin(), counts.end(), rates.begin(), [] (double v) {
        const double reference = Processors::Beampolmon->Reference_1MHz.Get();
        return 1.0e6*v/reference;
    });
    return rates;
}

std::vector<int64_t> TaggerScalers::GetCounts() const
{
    if(!slowcontrol_provided)
        return vector<int64_t>(nChannels, 1.0);

    vector<int64_t> counts(nChannels, std::numeric_limits<int64_t>::quiet_NaN());
    if (mode == mode_t::EPT_2014)
    {
        for (const auto& kv: Processors::EPT_Scalers->Get())
        {
            if (kv.Key < counts.size())
                counts[kv.Key] = kv.Value;
        }
    }

    if (mode == mode_t::Tagger)
    {
        for (const auto& kv: Processors::Tagger_Scalers->Get())
        {
            if (kv.Key < counts.size())
                counts[kv.Key] = kv.Value;
        }
    }

    return counts;
}

double TaggerScalers::GetTaggerOr() const
{
    if(!slowcontrol_provided)
        return 1.0;

    if (mode == mode_t::EPT_2014)
        return Processors::EPT_Or->Get() * 1.0e6 / Processors::Beampolmon->Reference_1MHz.Get();

    if (mode == mode_t::Tagger)
        return Processors::Tagger_Or->Get() * 1.0e6 / Processors::Beampolmon->Reference_1MHz.Get();

    return 0;

}

double TaggerScalers::GetTaggerOrAsSum() const
{
    const auto scalers = GetRates();
    return accumulate(scalers.begin(), scalers.end(), 0);
}



