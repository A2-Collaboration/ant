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
#include <limits>

class TTree;

namespace ant {

class TID;
class WrapTFileInput;
class UnpackerA2GeantConfig; // see below

class UnpackerA2Geant : public Unpacker::Module
{
public:
    UnpackerA2Geant();
    virtual ~UnpackerA2Geant();
    virtual bool OpenFile(const std::string& filename) override;
    virtual std::unique_ptr<TDataRecord> NextItem() noexcept override;

    class Exception : public Unpacker::Exception {
        using Unpacker::Exception::Exception; // use base class constructor
    };

    double PercentDone() const;

private:
    TID* id = nullptr;

    std::unique_ptr<WrapTFileInput> filemanager;
    TTree* geant;

    std::unique_ptr<THeaderInfo> headerInfo;
    std::list< std::shared_ptr<TaggerDetector_t> >  taggerdetectors;

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

// we define some methods here which
// the configs are required to implement
class UnpackerA2GeantConfig : public ExpConfig::Unpacker<UnpackerA2GeantConfig> {
public:
    virtual std::list< std::shared_ptr< Detector_t > > GetDetectors() const = 0;
};

} // namespace ant
