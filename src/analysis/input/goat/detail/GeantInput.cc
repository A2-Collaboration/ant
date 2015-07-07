#include "GeantInput.h"

#include "TTree.h"
#include <iostream>

using namespace ant;
using namespace input;


GeantInput::GeantInput()
{
}

GeantInput::~GeantInput()
{
}

void GeantInput::GetEntry()
{
}

bool GeantInput::SetupBranches(TreeRequestManager &&input_files) {

    TTree* geant = input_files.GetTree("h12");

    if(geant==nullptr)
        return false;

    geant->SetBranchAddress("nhits",&fnhits);
    geant->SetBranchAddress("npart",&fnpart);
    geant->SetBranchAddress("ntaps",&fntaps);
    geant->SetBranchAddress("nvtaps",&fnvtaps);
    geant->SetBranchAddress("vhits",&fvhits);
    geant->SetBranchAddress("plab",plab);
    geant->SetBranchAddress("tctaps",tctaps);
    geant->SetBranchAddress("vertex",fvertex);
    geant->SetBranchAddress("beam",fbeam);
    geant->SetBranchAddress("dircos",dircos);
    geant->SetBranchAddress("ecryst",ecryst);
    geant->SetBranchAddress("tcryst",tcryst);
    geant->SetBranchAddress("ectapfs",ectapfs);
    geant->SetBranchAddress("ectapsl",ectapsl);
    geant->SetBranchAddress("elab",elab);
    geant->SetBranchAddress("eleak",&feleak);
    geant->SetBranchAddress("enai",&fenai);
    geant->SetBranchAddress("etot",&fetot);
    geant->SetBranchAddress("eveto",eveto);
    geant->SetBranchAddress("tveto",tveto);
    geant->SetBranchAddress("evtaps",evtaps);
    geant->SetBranchAddress("icryst",icryst);
    geant->SetBranchAddress("ictaps",ictaps);
    geant->SetBranchAddress("ivtaps",ivtaps);
    geant->SetBranchAddress("idpart",idpart);
    geant->SetBranchAddress("iveto",iveto);
    geant->SetBranchAddress("nmwpc",&fnmwpc);
    geant->SetBranchAddress("imwpc",imwpc);
    geant->SetBranchAddress("mposx",mposx);
    geant->SetBranchAddress("mposy",mposy);
    geant->SetBranchAddress("mposz",mposz);
    geant->SetBranchAddress("emwpc",emwpc);
    geant->SetBranchAddress("mc_evt_id",&mc_evt_id);
    geant->SetBranchAddress("mc_rnd_id",&mc_evt_id);

    return true;
}

TLorentzVector GeantInput::GetBeam() const
{
    double x = fbeam[0], y = fbeam[1], z = fbeam[2], t = fbeam[3];

    return TLorentzVector(t*x, t*y, t*z, t);
}

TVector3 GeantInput::GetVertex() const
{
    return TVector3(fvertex);
}

Double_t GeantInput::GetCBESum() const
{
    return fetot;
}

Double_t GeantInput::GetPIDESum() const
{
    return sumArray(GetPIDHitEnergies(), GetNPIDHits());
}

Double_t GeantInput::GetTAPSESum() const
{
    return sumArray(ectapsl, fntaps);
}

Double_t GeantInput::GetTAPSVetoESum() const
{
    return sumArray(evtaps, fnvtaps);
}


Int_t GeantInput::GetNCBHits() const
{
    return fnhits;
}

Int_t GeantInput::GetCBHitIndex(const UInt_t n) const throw (std::out_of_range)
{
    if( n >= (UInt_t)GetNCBHits() )
        throw std::out_of_range(Form("%s: CB Hit index out of bounds. (%d/%d)", __func__, n, GetNCBHits()));

    return icryst[n];
}

Float_t GeantInput::GetCBHitEnergy(const UInt_t n) const throw (std::out_of_range)
{
    if( n >= (UInt_t)GetNCBHits() )
        throw std::out_of_range(Form("%s: CB Hit index out of bounds. (%d/%d)", __func__, n, GetNCBHits()));

    return ecryst[n];
}

const Int_t * const GeantInput::GetCBHitIndices() const
{
    return icryst;
}

const Float_t * const GeantInput::GetCBHitEnergies() const
{
    return ecryst;
}


Int_t GeantInput::GetNTAPSHits() const
{
    return fntaps;
}

const Int_t * const GeantInput::GetTAPSHitIndices() const
{
    return ictaps;
}

Int_t GeantInput::GetTAPSHitIndex(const UInt_t n) const throw (std::out_of_range)
{
    if( n >= (UInt_t)GetNTAPSHits() )
        throw std::out_of_range(Form("%s: TAPS Hit index out of bounds. (%d/%d)", __func__, n, GetNTAPSHits()));

    return ictaps[n];
}

Float_t GeantInput::GetTAPSHitEnergyLong(const UInt_t n) const throw (std::out_of_range)
{
    if( n >= (UInt_t)GetNTAPSHits() )
        throw std::out_of_range(Form("%s: TAPS Hit index out of bounds. (%d/%d)", __func__, n, GetNTAPSHits()));

    return ectapsl[n];
}

Float_t GeantInput::GetTAPSHitEnergyShort(const UInt_t n) const throw (std::out_of_range)
{
    if( n >= (UInt_t)GetNTAPSHits() )
        throw std::out_of_range(Form("%s: TAPS Hit index out of bounds. (%d/%d)", __func__, n, GetNTAPSHits()));

    return ectapfs[n];
}

const Float_t * const GeantInput::GetTAPSHitEnergiesLong() const
{
    return ectapsl;
}

const Float_t * const GeantInput::GetTAPSHitEnergiesShort() const
{
    return ectapfs;
}

Float_t GeantInput::GetTAPSHitTime(const UInt_t n) const throw (std::out_of_range)
{
    if( n >= (UInt_t)GetNTAPSHits() )
        throw std::out_of_range(Form("%s: TAPS Hit index out of bounds. (%d/%d)", __func__, n, GetNTAPSHits()));

    return tctaps[n];
}

const Float_t * const GeantInput::GetTAPSHitTimes() const
{
    return tctaps;
}

Int_t GeantInput::GetNTAPSVetoHits() const
{
    return fnvtaps;
}

Int_t GeantInput::GetTAPSVetoHitIndex(const UInt_t n) const throw (std::out_of_range)
{
    if( n >= (UInt_t)GetNTAPSVetoHits() )
        throw std::out_of_range(Form("%s: TAPS Hit index out of bounds. (%d/%d)", __func__, n, GetNTAPSVetoHits()));

    return ivtaps[n];
}

const Int_t * const GeantInput::GetTAPSVetoHitIndices() const
{
    return ivtaps;
}

Float_t GeantInput::GetTAPSVetoHitEnergy(const UInt_t n) const throw (std::out_of_range)
{
    if( n >= (UInt_t)GetNTAPSVetoHits() )
        throw std::out_of_range(Form("%s: TAPS Hit index out of bounds. (%d/%d)", __func__, n, GetNTAPSVetoHits()));

    return evtaps[n];
}

const Float_t * const GeantInput::GetTAPSVetoHitEnergies() const
{
    return evtaps;
}

Int_t GeantInput::GetNPIDHits() const
{
    return fvhits;
}

Int_t GeantInput::GetPIDHitIndex(const UInt_t n) const throw (std::out_of_range)
{
    if( n >= (UInt_t)GetNPIDHits() )
        throw std::out_of_range(Form("%s: PID Hit index out of bounds. (%d/%d)", __func__, n, GetNPIDHits()));

    return iveto[n];
}

Float_t GeantInput::GetPIDHitEnergy(const UInt_t n) const throw (std::out_of_range)
{
    if( n >= (UInt_t)GetNPIDHits() )
        throw std::out_of_range(Form("%s: PID Hit index out of bounds. (%d/%d)", __func__, n, GetNPIDHits()));

    return eveto[n];
}

const Int_t * const GeantInput::GetPIDHitIndices() const
{
    return iveto;
}

const Float_t * const GeantInput::GetPIDHitEnergies() const
{
    return eveto;
}

Float_t GeantInput::GetPIDHitTime(const UInt_t n) const throw (std::out_of_range)
{
    if( n >= (UInt_t)GetNPIDHits() )
        throw std::out_of_range(Form("%s: PID Hit index out of bounds. (%d/%d)", __func__, n, GetNPIDHits()));

    return tveto[n];
}

const Float_t * const GeantInput::GetPIDHitTimes() const
{
    return tveto;
}

Int_t GeantInput::GetNMWPCHits() const
{
    return fnmwpc;
}

const Int_t * const GeantInput::GetMWPCHitIndices() const
{
    return imwpc;
}

const Float_t * const GeantInput::GetMWPCHitEnergies() const
{
    return emwpc;
}

const Float_t * const GeantInput::GetMWPCHitPosX() const
{
    return mposx;
}

const Float_t * const GeantInput::GetMWPCHitPosY() const
{
    return mposy;
}

const Float_t * const GeantInput::GetMWPCHitPosZ() const
{
    return mposz;
}

TVector3 GeantInput::GetMWPCVector(const UInt_t n) const throw(std::out_of_range)
{
    if( n >= (UInt_t)GetNMWPCHits() )
        throw std::out_of_range(Form("%s: MWPC Particle index out of bounds. (%d/%d)", __func__, n, GetNMWPCHits()));

    return TVector3(mposx[n], mposy[n], mposz[n]);
}

Float_t GeantInput::GetMWPCEnergy(const UInt_t n) const throw(std::out_of_range)
{
    if( n >= (UInt_t)GetNMWPCHits() )
        throw std::out_of_range(Form("%s: MWPC Particle index out of bounds. (%d/%d)", __func__, n, GetNMWPCHits()));

    return emwpc[n];
}

Int_t GeantInput::GetMWPCIndex(const UInt_t n) const throw(std::out_of_range)
{
    if( n >= (UInt_t)GetNMWPCHits() )
        throw std::out_of_range(Form("%s: MWPC Particle index out of bounds. (%d/%d)", __func__, n, GetNMWPCHits()));

    return imwpc[n];
}

void GeantInput::BuildCBHitPattern(GeantInput::hitvector &pattern) const
{
    try {
        buildPattern(pattern, icryst, ecryst, GetNCBHits(), 720);
    } catch (const std::out_of_range& oor) {
        std::cerr << "CB crystal index out of range " << oor.what() << std::endl;
    }
}

void GeantInput::BuildTAPSHitPattern(GeantInput::hitvector &pattern) const
{
    try {
        buildPattern(pattern, ictaps, ectapsl, GetNTAPSHits(), 438);
    } catch (const std::out_of_range& oor) {
        std::cerr << "TAPS crystal index out of range " << oor.what() << std::endl;
    }
}

UInt_t GeantInput::GetNTrueParticles() const
{
    return fnpart;
}

UInt_t GeantInput::GetTrueID(const UInt_t n) const throw(std::out_of_range)
{
    if( n >= GetNTrueParticles() )
        throw std::out_of_range(Form("%s: MC-True Particle index out of bounds. (%d/%d)", __func__, n, GetNTrueParticles()));

    return idpart[n];

}

TLorentzVector GeantInput::GetTrueVector(const UInt_t n) const throw(std::out_of_range)
{
    if( n >= GetNTrueParticles() )
        throw std::out_of_range(Form("%s: MC-True Particle index out of bounds. (%d/%d)", __func__, n, GetNTrueParticles()));

    TVector3 p(dircos[n][0], dircos[n][1], dircos[n][2]);
    p *= plab[n];

    TLorentzVector lv( p, elab[n]);

    return lv;
}


Double_t GeantInput::sumArray(const Float_t * const data, const Int_t size)
{
    Double_t sum(0.0);

    for( Int_t i=0; i<size; ++i )
        sum += data[i];

    return sum;
}

void GeantInput::buildPattern(hitvector &pattern, const Int_t * const indices, const Float_t * const data, const Int_t nhits, const Int_t patternsize)
{

    pattern.assign(patternsize, 0.0);

    for(Int_t i=0; i < nhits; ++i ) {
        pattern.at(i) = data[ indices[i] ];
    }
}
