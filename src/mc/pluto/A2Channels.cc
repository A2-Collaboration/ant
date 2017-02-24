#include "A2Channels.h"
#include "mc/pluto/chgen.h"

#include <list>
#include <iostream>
#include <sstream>
#include <fstream>

#include "TCanvas.h"
#include "TH1D.h"

#include "base/Logger.h"
#include "analysis/utils/particle_tools.h"

// Pluto
#include "PDecayChannel.h"

using namespace std;
using namespace ant::mc::pluto;


double A2ChannelManager::Xsection(const string &name, const double Egamma) const
{
    auto ch = _XList.find(name);
    if ( ch == _XList.end())
    {
        cerr << "  Warning:  '" << name << "' does not exist, returning 0 for Xsection." << endl;
        return 0;
    }

    // don't crash for unknown energies, just use last known value
    if (Egamma < ch->second.Energies.front()){
        return  ch->second.Xsections.front();
    }
    if (Egamma > ch->second.Energies.back()){
        return  ch->second.Xsections.back();
    }

    // interpolate:
    double xsec = ROOT::Math::Interpolator(ch->second.Energies,ch->second.Xsections).Eval(Egamma);

    //quickfix for Interpolator smoothing into negative numbers
    if (xsec < 0 ) xsec = 0;

    return xsec;
}

double A2ChannelManager::TotalXsection(const double& Egamma) const
{
    double sigmatot = 0;

    for ( auto ch: _XList)
        sigmatot += Xsection(ch.first, Egamma);

    return sigmatot;
}

void A2ChannelManager::unifyDecayName(string &decay) const
{
    const string wspace = " ";

    //whitespace at the beginning
    auto start = decay.find_first_not_of(wspace);
    if ( start == string::npos)
    {
        decay="";
        return;
    }


    //... and end
    auto stop = decay.find_last_not_of(wspace);

    decay = decay.substr(start,stop - start + 1);

    //... and  reduce in middle to one
    start = decay.find_first_of(wspace);
    while (start != string::npos)
    {
        stop = decay.find_first_not_of(wspace, start);
        decay.replace(start, stop - start, wspace);
        start = decay.find_first_of(wspace, start + wspace.length());
    }

}



A2ChannelManager::A2ChannelManager()
{

    _XList = LoadStdDecays();

}

vector<string> A2ChannelManager::GetChannels() const
{
    vector<string> names;
    for (const auto& ch: _XList) names.push_back(ch.first);
    return names;
}

