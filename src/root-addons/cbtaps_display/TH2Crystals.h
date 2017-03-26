#pragma once

#include "Rtypes.h"
#include "TH2Poly.h"
#include "TVector2.h"

#include <string>

// disable the override keyword
#ifdef __CINT__
#define override
#endif

#define ANT_CBTAPS_DISPLAY_VERSION 2

class TMarker;
class TGraph;
class TVector2;

namespace ant {

/**
 * @brief Base class for A2 Calorimeter display classes.
 */

class TH2Crystals: public TH2Poly {

protected:
    using TH2Poly::SetBinContent;
    using TH2Poly::GetBinContent;

    bool draw_glue_pads = false;

    /**
     * @brief Add the content of another TH2Crystal object, if number of bins match
     * @param h The other TH2Crystals object
     *
     * This method is kept protected in the base class and gets used by defived classes.
     * This avoids the possibility to add a TAPS to a CB or something like that.
     *
     * This method is efficient (use of iterators, interating once).
     */
    virtual void FillElements( const TH2Crystals& h);

    /**
     * @brief Replace all values with the content of another TH2Crystal object, if number of bins match
     * @param h The other TH2Crystals object
     *
     * This method is kept protected in the base class and gets used by defived classes.
     * This avoids the possibility to add a TAPS to a CB or something like that.
     *
     * This method is efficient (use of iterators, interating once).
     */
    virtual void SetElements( const TH2Crystals& h);

    void SetMarkerOnBin(Int_t bin);

    /**
     * @brief calcCOG calculates center of gravity
     * @param g
     * @param x
     * @param y
     */
    void calcCOG(TGraph* g, double& x, double& y) const;

public:

    /**
     * @brief Base Constructior. Sets name and title, and links the object into the
     *        current ROOT directory (fix for a Bug in TH2Poly)
     * @param name Name for the new object
     * @param title Title for the new object
     */
    TH2Crystals(const std::string &name="", const std::string &title="");

    virtual ~TH2Crystals() {}

    virtual Double_t GetElement( const UInt_t element ) const;
    virtual void SetElement( const UInt_t element, Double_t value );

    /**
     * @brief Fill the number of the TH2Poly bin. Useful for debugging only, I guess.
     */
    virtual void FillBinNumbers();

    /**
     * @brief Fill the element numbers.
     */
    virtual void FillElementNumbers();  //*MENU*

    /**
     * @brief Set values of all elements to the ones stored in a std::vector (values ordered by element number)
     * @param pattern Values to set. size() has to be equal to number of elements of the detector.
     */
    virtual void SetElements( const std::vector<Double_t>& pattern );

    /**
     * @brief Set values of all elements to the ones stored in a TH1 (values ordered by element number)
     * @param h The histogram with the values to set. Number of x-bins (without under-/overflow) has to be equal to number of elements of the detector.
     */
    virtual void SetElements( const TH1& h );

    virtual void FillElement( const UInt_t element, const Double_t w );
    virtual void FillElements( const std::vector<Double_t>& pattern );
    virtual void FillElements( const TH1& h );

    /**
     * @brief Get the number of elements
     *        Usually this is the number of bins in the TH2Poly, but not always.
     * @return number of elements
     */
    virtual Int_t GetNumberOfElements() const;

    /**
     * @brief Reset Elements (clear everything)
     * @param value Value to set the elements to (default 0.0)
     */
    virtual void ResetElements( const Double_t value=0.0 ); //*MENU*

    /**
     * @brief CreateMarker creates marker at given element
     * @param element
     * @return
     */
    virtual void CreateMarker(UInt_t element); // *MENU*

    /**
     * @brief CreateMarker black and white markers at position p with two (different) styles
     * @param p position
     * @param marker_style1 style of black marker
     * @param marker_style2 style of white marker
     */
    virtual void CreateMarker(const TVector2& p, const int marker_style_black, const int marker_style_white);

    /**
     * @brief ClearMarkers removes all markers created by CreateMarker
     */
    virtual void ClearMarkers(); // *MENU*



    ClassDef(TH2Crystals,ANT_CBTAPS_DISPLAY_VERSION)
};

}
