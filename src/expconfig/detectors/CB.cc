#include "CB.h"
#include <cassert>
#include "base/std_ext.h"

using namespace std;
using namespace ant;
using namespace ant::expconfig::detector;

void CB::BuildMappings(vector<UnpackerAcquConfig::hit_mapping_t> &hit_mappings,
                       vector<UnpackerAcquConfig::scaler_mapping_t>&) const {
  // CB has only hit_mappings to add,
  // no scalers
  for(const CBElement_t& element : elements)  {
    hit_mapping_t m_adc;
    m_adc.LogicalChannel = {Type, Channel_t::Type_t::Integral, element.Channel};
    m_adc.RawChannels.push_back(element.ADC);
    hit_mappings.emplace_back(move(m_adc));

    hit_mapping_t m_tdc;
    m_tdc.LogicalChannel = {Type, Channel_t::Type_t::Timing, element.Channel};
    m_tdc.RawChannels.push_back(element.TDC);
    hit_mappings.emplace_back(move(m_tdc));
  }
}

vector<CB::CBElement_t> CB::initElements() {
  // hard-code the holes
  // we even use insert range for single elements
  // for the sake of code cleanness
  vector<unsigned> holes;
  std_ext::insertRange(holes,  26,  26);
  std_ext::insertRange(holes,  29,  38);
  std_ext::insertRange(holes,  40,  40);
  std_ext::insertRange(holes, 311, 311);
  std_ext::insertRange(holes, 315, 316);
  std_ext::insertRange(holes, 318, 319);
  std_ext::insertRange(holes, 353, 366);
  std_ext::insertRange(holes, 400, 402);
  std_ext::insertRange(holes, 405, 405);
  std_ext::insertRange(holes, 408, 408);
  std_ext::insertRange(holes, 679, 679);
  std_ext::insertRange(holes, 681, 689);
  std_ext::insertRange(holes, 691, 692);

#include "detail/CB_elements.h" // defines all_elements

  // build the really existing elements
  vector<CB::CBElement_t> elements;
  for(const CB::CBElement_t& elem : all_elements) {
    // if the element is a hole, don't add it
    if(find(holes.cbegin(), holes.cend(), elem.Channel) != holes.cend())
      continue;
    elements.push_back(elem);
  }
  assert(elements.size() == 672);
  return elements;
}


const vector<CB::CBElement_t> CB::elements = CB::initElements();



