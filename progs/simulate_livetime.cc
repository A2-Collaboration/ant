#include "base/Logger.h"
#include "base/CmdLine.h"
#include "base/std_ext/math.h"

#include <random>
#include <cmath>
#include <ostream>

using namespace std;
using namespace ant;

default_random_engine r_gen;
uniform_real_distribution<double> r_zero_one;

double get_poisson_time(const double lambda) {
    return -std::log(r_zero_one(r_gen))/lambda;
}

ostream& operator<<(ostream& s, const std_ext::RMS& rms) {
    return s << std::fixed << rms.GetMean() << " +/- " << rms.GetRMS();
}

int main( int argc, char** argv )
{
    SetupLogger();

    TCLAP::CmdLine cmd("simulate_livetime", ' ', "0.1");

    cmd.parse(argc, argv);

    const double PbGlassRate = 10000;

    const double LiveTime = 0.8;
    const double ExpTriggerRate = 4000;
    const double DeadTime = (1-LiveTime)/ExpTriggerRate;

    const double MaxTotalTime = 200;

    LOG(INFO) << "Dead time: " << DeadTime;

    std_ext::RMS RMS_poisson_time;
    unsigned nHits = 0;
    unsigned nEvents = 0;
    {
        double total_time = 0.0;
        double system_dead_until = 0.0;

        while(total_time<MaxTotalTime) {
            auto poisson_time = get_poisson_time(PbGlassRate);
            RMS_poisson_time.Add(poisson_time);

            if(total_time>=system_dead_until) {
                nEvents++;

                system_dead_until = total_time + DeadTime;
            }

            total_time += poisson_time;
            nHits++;
        }
    }

    LOG(INFO) << RMS_poisson_time;
    LOG(INFO) << "Hits=" << nHits << " Events=" << nEvents;
    const double EffectiveLiveTime = (double)nEvents/nHits;
    LOG(INFO) << "Effective Livetime " << EffectiveLiveTime;

    return EXIT_SUCCESS;
}