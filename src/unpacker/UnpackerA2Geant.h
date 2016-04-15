#pragma once

#include "Unpacker.h"

#include "expconfig/ExpConfig.h"
#include "base/Detector_t.h"

#include "Rtypes.h"

#include <memory>
#include <list>
#include <map>
#include <vector>
#include <cstdint>

class TTree;

namespace ant {

struct TID;
class WrapTFileInput;

namespace unpacker {
namespace geant {
struct promptrandom_t;
}}

/**
 * @brief The UnpackerA2Geant class
 */
class UnpackerA2Geant : public Unpacker::Module
{
public:
    UnpackerA2Geant();
    virtual ~UnpackerA2Geant();
    virtual bool OpenFile(const std::string& filename) override;
    virtual TEvent NextEvent() noexcept override;

    class Exception : public Unpacker::Exception {
        using Unpacker::Exception::Exception; // use base class constructor
    };

    virtual double PercentDone() const override;

private:
    TID* id = nullptr;

    std::unique_ptr<WrapTFileInput> inputfile;
    TTree* geant;

    std::shared_ptr<TaggerDetector_t> taggerdetector;
    std::shared_ptr<Detector_t> cb_detector;
    std::shared_ptr<Detector_t> pid_detector;
    std::shared_ptr<Detector_t> taps_detector;
    std::shared_ptr<Detector_t> tapsveto_detector;

    std::unique_ptr<unpacker::geant::promptrandom_t> promptrandom;

    // keep in syn with A2CBoutput.h in a2geant
    static constexpr int GEANT_MAX_TAPSHITS = 438;
    static constexpr int GEANT_MAX_CBHITS   = 720;
    static constexpr int GEANT_MAX_PIDHITS  =  24;
    static constexpr int GEANT_MAX_MWPCHITS = 400;
    static constexpr int GEANT_MAX_PART     = 100;

    Long64_t current_entry = -1;

    // Brach memories
    Int_t           fnhits = 0;
    Int_t           fnpart = 0;
    Int_t           fntaps = 0;
    Int_t           fnvtaps = 0;
    Int_t           fvhits = 0;
    Float_t         plab[GEANT_MAX_PART] = {};
    Float_t         tctaps[GEANT_MAX_TAPSHITS] = {};
    Float_t         fvertex[3] = {};
    Float_t         fbeam[5] = {};
    Float_t         dircos[GEANT_MAX_PART][3] = {};
    Float_t         ecryst[GEANT_MAX_CBHITS] = {};
    Float_t         tcryst[GEANT_MAX_CBHITS] = {};
    Float_t         ectapfs[GEANT_MAX_TAPSHITS] = {};
    Float_t         ectapsl[GEANT_MAX_TAPSHITS] = {};
    Float_t         elab[GEANT_MAX_PART] = {};
    Float_t         feleak = 0;
    Float_t         fenai = 0;
    Float_t         fetot = 0;
    Float_t         eveto[GEANT_MAX_PIDHITS] = {};
    Float_t         tveto[GEANT_MAX_PIDHITS] = {};
    Float_t         evtaps[GEANT_MAX_TAPSHITS] = {};
    Int_t           icryst[GEANT_MAX_CBHITS] = {};
    Int_t           ictaps[GEANT_MAX_TAPSHITS] = {};
    Int_t           ivtaps[GEANT_MAX_TAPSHITS] = {};
    Int_t           idpart[GEANT_MAX_PART] = {};
    Int_t           iveto[GEANT_MAX_PIDHITS] = {};
    Int_t           fnmwpc = 0;
    Int_t           imwpc[GEANT_MAX_MWPCHITS] = {};
    Float_t         mposx[GEANT_MAX_MWPCHITS] = {};
    Float_t         mposy[GEANT_MAX_MWPCHITS] = {};
    Float_t         mposz[GEANT_MAX_MWPCHITS] = {};
    Float_t         emwpc[GEANT_MAX_MWPCHITS] = {};

    bool tid_from_file = false;

};

/**
 * @brief The UnpackerA2GeantConfig class provides stuff the A2Geant unpacker needs
 */
class UnpackerA2GeantConfig {
public:
    virtual std::list< std::shared_ptr< Detector_t > > GetDetectors() const = 0;

    /**
     * @brief The promptrandom_config_t struct
     *
     * Contains information about the timing spectrum of the TaggerDetector_t
     *
     * Fit the timing spectrum of the tagger with gaus(0)+pol0(3) with ROOT's FitPanel,
     * then the following parameters are used
     * p0 = height (not integral)
     * p1 = position
     * p2 = sigma
     * p3 = offset
     */
    struct promptrandom_config_t {
        /**
         * @brief PromptSigma the sigma of the fit
         */
        double PromptSigma = 0;
        /**
         * @brief RandomPromptRatio per unit time interval (in ns)
         *
         * can be calculated by
         * offset/(sqrt(2pi)*sigma*height)
         */
        double RandomPromptRatio = 0;

        /**
         * @brief TimeWindowLength is by default the maximum window possible, in ns
         */
        interval<double> TimeWindow = interval<double>{-1000, 1000};
    };

    virtual promptrandom_config_t GetPromptRandomConfig() const {
        return {}; // use default
    }
protected:
    ~UnpackerA2GeantConfig() = default;
};

} // namespace ant
