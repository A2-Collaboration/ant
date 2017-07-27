#pragma once

#include "Unpacker.h"

#include "expconfig/ExpConfig.h"
#include "base/Detector_t.h"
#include "base/WrapTTree.h"

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
    virtual TEvent NextEvent() override;
    virtual bool ProvidesSlowControl() const override { return false; }

    class Exception : public Unpacker::Exception {
        using Unpacker::Exception::Exception; // use base class constructor
    };

    virtual double PercentDone() const override;

private:
    // important to declare inputfile before WrapTTree
    std::unique_ptr<WrapTFileInput> inputfile;

    std::shared_ptr<TaggerDetector_t> taggerdetector;
    std::shared_ptr<Detector_t> cb_detector;
    std::shared_ptr<Detector_t> pid_detector;
    std::shared_ptr<Detector_t> taps_detector;
    std::shared_ptr<Detector_t> tapsveto_detector;

    std::unique_ptr<unpacker::geant::promptrandom_t> promptrandom;

    long long current_entry = -1;

    struct TIDTree_t : WrapTTree {
        ADD_BRANCH_T(TID, tid)
    };

    TIDTree_t tidTree;

    struct GeantTree_t : WrapTTree {
        ADD_BRANCH_T(ROOTArray<Float_t>, plab)
        ADD_BRANCH_T(ROOTArray<Float_t>, tctaps)
        ADD_BRANCH_T(ROOTArray<Float_t>, vertex)
        ADD_BRANCH_T(ROOTArray<Float_t>, beam)
        ADD_BRANCH_T(ROOTArray_Float<3>, dircos)
        ADD_BRANCH_T(ROOTArray<Float_t>, ecryst)
        ADD_BRANCH_OPT_T(ROOTArray<Float_t>, tcryst)
        ADD_BRANCH_T(ROOTArray<Float_t>, ectapfs)
        ADD_BRANCH_T(ROOTArray<Float_t>, ectapsl)
        ADD_BRANCH_T(ROOTArray<Float_t>, elab)
        ADD_BRANCH_T(ROOTArray<Float_t>, eveto)
        ADD_BRANCH_OPT_T(ROOTArray<Float_t>, tveto)
        ADD_BRANCH_T(ROOTArray<Float_t>, evtaps)
        ADD_BRANCH_T(ROOTArray<Int_t>,   icryst)
        ADD_BRANCH_T(ROOTArray<Int_t>,   ictaps)
        ADD_BRANCH_OPT_T(ROOTArray<Int_t>,   ivtaps)
        ADD_BRANCH_T(ROOTArray<Int_t>,   idpart)
        ADD_BRANCH_T(ROOTArray<Int_t>,   iveto)
        ADD_BRANCH_OPT_T(ROOTArray<Int_t>,   imwpc)
        ADD_BRANCH_OPT_T(ROOTArray<Float_t>, mposx)
        ADD_BRANCH_OPT_T(ROOTArray<Float_t>, mposy)
        ADD_BRANCH_OPT_T(ROOTArray<Float_t>, mposz)
        ADD_BRANCH_OPT_T(ROOTArray<Float_t>, emwpc)
    };

    GeantTree_t geantTree;

    bool oldTreeFormat = false;

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

        /**
         * @brief TimeOffset if the tagger time is determined with some offsetting reference
         */
        double PromptOffset = 0;
    };

    virtual promptrandom_config_t GetPromptRandomConfig() const {
        return {}; // use default
    }
protected:
    ~UnpackerA2GeantConfig() = default;
};

} // namespace ant
