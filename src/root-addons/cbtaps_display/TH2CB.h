#pragma once

#include <set>
#include <vector>
#include <string>

#include "TH2Crystals.h"

#include "TH2Poly.h"
#include "Rtypes.h"


/**
 * @brief The TH2CB class
 *
 * A TH2Poly of the Crystal Ball mesh.
 *
 *  \image html cb.png
 *
 * bin number (bin): Index of a bin in the TH2Poly. [1..n], used to access the individual bins with SetBinContent(bin, value), etc.
 *  A bin exists for every physical crystal. There are 672 crystals.
 *
 * virtual bin number (vbin): A virtual bin number including non-existant crystals in the hole regions, [1..720].
 *
 * Major-Minior-Crystal Number (MMC): The Crystal Ball is divided into major triganles (20) which again are subdivided into minor triangles (4).
 *  Each mintor triangle consists of 9 crystals. The MMC number addesses crystals by their Major-Minor-Crystal indices.
 *
 * Element Number (element): The crystals are not cabled in the ordering of the crystals. Acqu uses elements.
 */

namespace ant {

class TH2DrawTool;

class TH2CB: public TH2Crystals {

protected:

    virtual void Build();
    void MakeLevel(TH2DrawTool& c, const UInt_t n, std::set<Int_t>::const_iterator& nexthole, Int_t& vbins);

    bool draw_glue_pads = false;

public:
    TH2CB( const std::string& name="", const std::string& title="", bool glue_pads=false);
    virtual ~TH2CB() {}

    static Int_t GetBinOfMMC(const UChar_t major, const UChar_t minor, const UChar_t crystal);
    static Int_t GetVBinOfMMC(const UChar_t major, const UChar_t minor, const UChar_t crystal);
    static Int_t GetBinOfVBin( const Int_t vbin );

    static bool IsInHole( const UChar_t a, const UChar_t b, const UChar_t c);
    static bool IsInHole( const Int_t vbin );

    // pull in methods from base class
    using TH2Crystals::SetElements;
    using TH2Crystals::FillElements;

    virtual void SetElements( const TH2CB& h);
    virtual void FillElements( const TH2CB& h);

    /**
     * @brief Fill in the crystal numbers.
     *
     * Crystal numers count in the same way as the element numbers (ex.: 1/1/1), stating with 0.
     */
    virtual void FillCrystalNumbers();     //*MENU*

    /**
     * @brief Fill the Major-Minor-Crystal numbers in. ex: 1/1/1 -> 111 and 14/2/8 -> 1428.
     */
    void FillMMCNumbers();     //*MENU*

    /**
     * @brief Fill in the element numbers. Crystal numers are mapped to element numbers.
     */
    void FillElementNumbers() override;     //*MENU*

    /**
     * @brief Fill a hit pattern (unmapped), only containing existing crystals
     * @param pattern Vector of the lenfth 672 (is checked!)
     * @see FillCrystals720()
     */
    void FillCrystals672( const std::vector<Double_t>& pattern );

    /**
     * @brief Fil a hit pattern (unmapped), containing crystals in holes
     * @param pattern Vector of the lenfth 720 (is checked!)
     * @see FillCrystals672()
     */
    void FillCrystals720( const std::vector<Double_t>& pattern );

    /**
     * @brief Get value of a crystal (unmapped), only counting exising crystals (no holes)
     * @param i Crystal number [0..671]
     * @return Content of the crystal
     */
    Double_t GetCrystal672(const UInt_t i) const;

    /**
     * @brief Get value of a crystal (unmapped), counting also crystals in holes
     * @param i Crystal number, might be in a hole, [0..719]
     * @return Content of the crystal. If the crystal is in a hole, 0.0 ist returend
     */
    Double_t GetCrystal720(const UInt_t i) const;

    /**
     * @brief Set the value of a crystal (unmapped), counting only existing crystals (no holes)
     * @param i Crystal number [0..671]
     * @param value Value to set it it
     * @see SetCrystal720()
     */
    void SetCrystal672(const UInt_t i, Double_t value);

    /**
     * @brief Get value of a crystal (unmapped), counting also crystals in holes
     * @param i Crystal number, might be in a hole, [0..719]
     * @param value Value to set it it
     * @see SetCrystal672()
     */
    void SetCrystal720(const UInt_t i, Double_t value);

    /**
     * @brief Get value of an element (mapped)
     * @param element Element number [0..720]
     * @return Content of the crystal, 0 if is inside a hole
     */
    virtual Double_t GetElement(const UInt_t element) const override;

    /**
     * @brief Set the value of an element (mapped)
     * @param element Element number [0..720]
     * @param value Value to set it it
     * @see SetCrystal720()
     */
    virtual void SetElement(const UInt_t element, Double_t value) override; 

    /**
     * @brief Fill a hit pattern (mapped), ordered by element numers
     * @param pattern Vector of the lenfth 720 (is checked!)
     * @see FillCrystals720()
     * @see FillCrystals672()
     */
    virtual void SetElements(const std::vector<Double_t>& pattern) override;
    virtual void SetElements(const TH1& h) override;

    virtual void FillElements(const TH1& h) override;

    /**
     * @brief Get the crystal number for an element number
     * @param element The element number
     * @return crystal number
     */
    static UInt_t GetCrystalOfElement( const UInt_t element );

    /**
     * @brief Get the element number for a crystal number
     * @param crystal The crystal number
     * @return element number
     */
    static UInt_t GetElementOfCrystal(const UInt_t crystal );

    /**
     * @brief Get the number of elements
     *        Crystal Ball also counts crystal positions that are in the hole regions, so there are element indices that do not have a crystal.
     * @return Number of elements (720)
     */
    Int_t GetNumberOfElements() const override { return 720; }

    virtual void CreateMarker(UInt_t element) override; // *MENU*
    
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Winconsistent-missing-override"
#endif
    ClassDef(TH2CB,ANT_CBTAPS_DISPLAY_VERSION)
#ifdef __clang__
#pragma clang diagnostic pop
#endif
};


}
