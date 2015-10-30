#include "Energy.h"

#include "calibration/DataManager.h"
#include "calibration/fitfunctions/FitGausPol0.h"
#include "calibration/gui/CalCanvas.h"

#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"

#include "tree/TCalibrationData.h"
#include "tree/TDetectorRead.h"

#include "base/Logger.h"

#include <cstdint>
#include <ctime>
#include <algorithm>
#include <sstream>
#include <vector>
#include <list>

using namespace std;
using namespace ant;
using namespace ant::calibration;

Energy::Energy(Detector_t::Type_t detectorType,
               std::shared_ptr<DataManager> calmgr,
               Calibration::Converter::ptr_t converter,
               double defaultPedestal,
               double defaultGain,
               double defaultThreshold,
               double defaultRelativeGain,
               Channel_t::Type_t channelType) :
    Calibration::Module(
        std_ext::formatter()
        << Detector_t::ToString(detectorType)
        << "_"
        << ( channelType == Channel_t::Type_t::IntegralShort ? "Short" : "" )
        << "Energy"
           ),
    DetectorType(detectorType),
    ChannelType(channelType),
    calibrationManager(calmgr),
    Converter(move(converter)),
    Pedestals(defaultPedestal,"Pedestals", true), // pedestals are always extendable
    Gains(defaultGain,"Gains"),
    Thresholds(defaultThreshold,"Thresholds"),
    RelativeGains(defaultRelativeGain,"RelativeGains")
{
    if(Converter==nullptr)
        throw std::runtime_error("Given converter should not be nullptr");
}

Energy::~Energy()
{
}

void Energy::ApplyTo(const readhits_t& hits, extrahits_t& extrahits)
{
    const auto& dethits = hits.get_item(DetectorType);

    // now calibrate the Energies (ignore any other kind of hits)
    for(TDetectorReadHit* dethit : dethits) {
        if(dethit->GetChannelType() != ChannelType)
            continue;



        // Values might already be filled
        // (for example by previous calibration run, or A2Geant unpacker),
        // then we apply the threshold and the relative gain only
        std::vector<double> values;

        // prefer RawData if available
        if(!dethit->RawData.empty()) {
            // convert to not-so-raw values (still not MeV scale)
            values = Converter->Convert(dethit->RawData);

            // for pedestal calibration, we insert extra hits here
            // containing the raw values
            if(NeedsPedestals()) {
                extrahits.emplace_back(
                            LogicalChannel_t{
                                dethit->GetDetectorType(),
                                ChannelType == Channel_t::Type_t::IntegralShort ?
                                Channel_t::Type_t::PedestalShort : Channel_t::Type_t::Pedestal,
                                dethit->Channel
                            },
                            values
                            );
            }

            // apply pedestal/gain to each of the values (might be multihit)
            for(double& value : values) {
                if(NeedsPedestals()) {
                    if(Pedestals.Values.empty())
                        value -= Pedestals.DefaultValue;
                    else
                        value -= Pedestals.Values[dethit->Channel];
                }

                if(Gains.Values.empty())
                    value *= Gains.DefaultValue;
                else
                    value *= Gains.Values[dethit->Channel];
            }

        }
        else {
            // maybe the values are already filled
            values = dethit->Values;
            dethit->Values.resize(0);
        }

        // always apply the threshold cut and the relative gains
        dethit->Values.reserve(values.size());

        for(double value : values) {
            if(RelativeGains.Values.empty())
                value *= RelativeGains.DefaultValue;
            else
                value *= RelativeGains.Values[dethit->Channel];

            const double threshold = Thresholds.Values.empty()
                                     ? Thresholds.DefaultValue
                                     : Thresholds.Values[dethit->Channel];
            if(value<threshold)
                continue;

            // only add if it passes the threshold
            dethit->Values.push_back(value);
        }

    }
}

std::vector<std::list<TID> > Energy::GetChangePoints() const
{
    vector<list<TID>> changePointLists;

    for (auto calibration: AllCalibrations) {
        changePointLists.push_back(
                    calibrationManager->GetChangePoints(
                        GUI_CalibType::ConstructName(GetName(), calibration->Name)
                        )
                    );
    }
    return changePointLists;
}
void Energy::Update(size_t index, const TID& tid)
{
    auto calibration = AllCalibrations[index];

    TCalibrationData cdata;
    if(calibrationManager->GetData(
           GUI_CalibType::ConstructName(GetName(), calibration->Name),
           tid, cdata))
    {
        auto& values = calibration->Values;
        for (const auto& val: cdata.Data) {
            if(values.size()<val.Key+1)
                values.resize(val.Key+1);
            values[val.Key] = val.Value;
        }
    }
    else {
        LOG_IF(!calibration->Values.empty(), WARNING)
                << "No calibration data found for " << calibration->Name
                << " at changepoint TID=" << tid << ", using default values";
        calibration->Values.resize(0);
    }
}

Energy::GUI_CalibType::GUI_CalibType(const string& basename, CalibType& type,
                                     const shared_ptr<DataManager>& calmgr,
                                     const shared_ptr<Detector_t>& detector_) :
    gui::Manager_traits(basename),
    calibType(type),
    calibrationManager(calmgr),
    detector(detector_)
{}

string Energy::GUI_CalibType::GetName() const
{
    // serves as the CalibrationID for the manager,
    // and as the histogram name
    return ConstructName(Manager_traits::GetName(), calibType.Name);
}

string Energy::GUI_CalibType::GetHistogramName() const
{
    return GetName();
}

unsigned Energy::GUI_CalibType::GetNumberOfChannels() const
{
    return detector->GetNChannels();
}

void Energy::GUI_CalibType::StartRange(const interval<TID>& range)
{
    // always make sure the values are large enough
    std::vector<double>& values = calibType.Values;
    values.resize(GetNumberOfChannels(), calibType.DefaultValue);

    TCalibrationData cdata;
    if(calibrationManager->GetData(GetName(), range.Start(), cdata)) {
        for(const TKeyValue<double>& kv : cdata.Data) {
            values[kv.Key] = kv.Value;
        }
        for(const TKeyValue<vector<double>>& kv : cdata.FitParameters) {
            fitParameters.insert(make_pair(kv.Key, kv.Value));
        }
        LOG(INFO) << GetName() << ": Loaded previous values from database";
    }
    else {
        LOG(INFO) << GetName() << ": No previous values found, built from default value";
    }

    // save a copy for comparison at finish stage
    previousValues = calibType.Values;

}

void Energy::GUI_CalibType::StoreFinishRange(const interval<TID>& range)
{
    TCalibrationData cdata(
                GetName(),
                range.Start(),
                range.Stop(),
                calibType.Extendable
                );

    std::vector<double>& values = calibType.Values;

    // fill data
    for(unsigned ch=0;ch<values.size();ch++) {
        cdata.Data.emplace_back(ch, values[ch]);
    }

    // fill fit parameters (if any)
    for(const auto& it_map : fitParameters) {
        const unsigned ch = it_map.first;
        const vector<double>& params = it_map.second;
        cdata.FitParameters.emplace_back(ch, params);
    }

    calibrationManager->Add(cdata);
}

Energy::GUI_Pedestals::GUI_Pedestals(const string& basename,
                          CalibType& type,
                          const std::shared_ptr<DataManager>& calmgr,
                          const std::shared_ptr<Detector_t>& detector) :
    GUI_CalibType(basename, type, calmgr, detector),
    func(make_shared<gui::FitGausPol0>())
{

}

void Energy::GUI_Pedestals::InitGUI(gui::ManagerWindow_traits* window)
{
    canvas = window->AddCalCanvas();
}

gui::Manager_traits::DoFitReturn_t Energy::GUI_Pedestals::DoFit(TH1* hist, unsigned channel,
                                                                    const Manager_traits::DoFitOptions_t& options)
{
    if(detector->IsIgnored(channel))
        return DoFitReturn_t::Skip;

    TH2* hist2 = dynamic_cast<TH2*>(hist);

    h_projection = hist2->ProjectionX("",channel+1,channel+1);

    func->SetDefaults(h_projection);
    const auto it_fit_param = fitParameters.find(channel);
    if(it_fit_param != fitParameters.end()
       && !options.IgnorePreviousFitParameters) {
        VLOG(5) << "Loading previous fit parameters for channel " << channel;
        func->Load(it_fit_param->second);
    }

    func->Fit(h_projection);

    /// \todo implement automatic stop if fit failed?

    // goto next channel
    return DoFitReturn_t::Next;
}

void Energy::GUI_Pedestals::DisplayFit()
{
    canvas->Show(h_projection, func.get());
}

void Energy::GUI_Pedestals::StoreFit(unsigned channel)
{

    const double oldValue = previousValues[channel];
    const double newValue = func->GetPeakPosition();

    calibType.Values[channel] = newValue;

    const double relative_change = 100*(newValue/oldValue-1);

    LOG(INFO) << "Stored Ch=" << channel << ":  "
              <<" Pedestal changed " << oldValue << " -> " << newValue
              << " (" << relative_change << " %)";


    // don't forget the fit parameters
    fitParameters[channel] = func->Save();
}

bool Energy::GUI_Pedestals::FinishRange()
{
    return false;
}
