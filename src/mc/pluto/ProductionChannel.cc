#include "ProductionChannel.h"

#include "database/gp_ppi0.h"


using namespace std;
using namespace ant::mc::pluto;


const ChannelDataBase::XSections_t ChannelDataBase::XSections = MakeXSections();



ChannelDataBase::XSections_t ChannelDataBase::MakeXSections()
{
    XSections_t x;

    x.insert(gp_ppi0);

    return x;
}
