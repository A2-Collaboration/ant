#include "A2Channels.h"
#include "chgen.h"

#include <list>
#include <iostream>
#include <sstream>
#include <fstream>

#include "TCanvas.h"
#include "TH1D.h"

using namespace std;


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

// ugly code, sorry...
const bool A2ChannelManager::ParseFile(const string &filename)
{
    const string decayDescriptor("Particle");
    ifstream dataFile;
    dataFile.open(filename);

    string decay;
    ParticleData pdata;

    // data file readable?
    if (!dataFile.is_open())
    {
        cerr << "  Warning:  " << filename << " not readable, skipping this file!" << endl;
        return false;
    }

    string line;
    unsigned int n_lines = 0;

    while(getline(dataFile,line))
    {
        // cut away comments
        line = line.substr(0,line.find("#"));
        // skip empty lines
        if (line.empty())
            continue;

        n_lines++;

        //first line should be a decay-particle description:
        if ((n_lines ==1 && !(line.find(decayDescriptor) == 0)))
        {
            cerr << "  Warning:  no parrticle description found skipping " << filename << endl;
            return false;
        }


        stringstream s_line_data(line);

        double parameter;
        vector<double> datapair;
        string descriptor;

        while( s_line_data >> parameter )
            datapair.push_back(parameter);

        //particle data should be a pair of two doubles (energy & xsection)
        if(datapair.size() == 2){
            pdata.Energies.push_back(datapair[0]);
            pdata.Xsections.push_back(datapair[1]);
        } else if (datapair.size() == 0)  // if not, look for a decay-name
        {
            stringstream s_line_text(line);
            s_line_text >> descriptor;
            if ( descriptor.find(decayDescriptor) == 0 ){
                if (n_lines!=1)
                {
                    _XList[decay] = pdata;
                    pdata = ParticleData();
                }
                decay = (s_line_text.str()).substr(s_line_text.tellg());
                unifyDecayName(decay);
            }

        } else
        {
            cerr << "  Warning:  skipping inputline with wrong number of parameters in" << filename << endl;
            continue;
        }
    }
    //add last dataset
    _XList[decay] = pdata;

    return true;
}

A2ChannelManager::A2ChannelManager(vector<string> dataFiles)
{
    unsigned int foundData = 0;

    _XList = XsecList();

    // if no external decays Provided, use precompiled ones:
    if ( dataFiles.empty() )
    {
        _XList = LoadStdDecays();
        return;
    }

    for ( auto& filename: dataFiles)
        foundData += ParseFile(filename);

    if ( foundData == 0 )
    {
        cerr << "  Warning:  No readable files found, falling back to precompiled data" << endl;
        _XList = LoadStdDecays();
    }
}

vector<string> A2ChannelManager::GetChannels() const
{
    vector<string> names;
    for ( auto& ch: _XList) names.push_back(ch.first);
    return names;
}



PDecayChannel* A2ChannelManager::GenerateDecays(const double &Energy)
{
    PDecayChannel* primaries = new PDecayChannel();
    for ( auto& ch: _XList){
        istringstream ss(ch.first);
        vector<string> dparticles{ istream_iterator<string>{ss},
                                   istream_iterator<string>{}};
        double xsec = Xsection(ch.first,     Energy);
        if ( xsec > 0 ){
            if ( dparticles.size() == 1 ){
                primaries->AddChannel(xsec,// * TotalXsection(Energy),
                                      (char *)"p",
                                      const_cast<char*>(dparticles[0].c_str()) );
            }
            if ( dparticles.size() == 2 ){
                primaries->AddChannel(xsec ,//* TotalXsection(Energy),
                                      (char *)"p",
                                      const_cast<char*>(dparticles[0].c_str()),
                        const_cast<char*>(dparticles[1].c_str()) );
            }
            if ( dparticles.size() == 3 ){
                primaries->AddChannel(xsec ,//* TotalXsection(Energy),
                                      (char *)"p",
                                      const_cast<char*>(dparticles[0].c_str()),
                        const_cast<char*>(dparticles[1].c_str()),
                        const_cast<char*>(dparticles[2].c_str()) );
            }
        }
    }

    return primaries;
}
