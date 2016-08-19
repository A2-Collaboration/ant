#include "base/Logger.h"
#include "base/CmdLine.h"
#include "base/std_ext/math.h"

#include <random>
#include <cmath>
#include <ostream>


#include "TH1D.h"
#include "TRint.h"
#include "analysis/plot/root_draw.h"



using namespace std;
using namespace ant;

static default_random_engine r_gen;
static uniform_real_distribution<double> r_zero_one;

double get_next_time(const double rate) {
//    uniform_real_distribution<double> r(0,2.0/rate);
//    normal_distribution<double> r(1.0/rate,2.0 / rate);
//    return r(r_gen);

    return -std::log(1.0-r_zero_one(r_gen))/rate;
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

    std_ext::RMS RMS_next_time;
    unsigned nHits = 0;
    unsigned nEvents = 0;

    double total_time = 0.0;
    double system_dead_until = 0.0;

    int fake_argc=1;
    char* fake_argv[2];
    fake_argv[0] = argv[0];
    auto app = new TRint("Ant-simulate_livetime",&fake_argc,fake_argv);
    TH1D* h = new TH1D("I1","I1(t)",50,0,0);
    TH1D* g = new TH1D("t","t",1000,0,50.0/PbGlassRate);
    TH1D* i = new TH1D("tm","t_measure",1000,0,50.0/PbGlassRate);

    while(total_time<MaxTotalTime) {
        auto next_time = get_next_time(PbGlassRate);
        h->Fill(next_time);
        RMS_next_time.Add(next_time);

        if(total_time>=system_dead_until) {
            nEvents++;
            i->Fill(total_time);
            system_dead_until = total_time + DeadTime;
        }

        total_time += next_time;
        g->Fill(total_time);
        nHits++;
    }

    LOG(INFO) << "Average Pb-glass rate: " << 1.0/RMS_next_time.GetMean();
//    LOG(INFO) << "Hits=" << nHits << " Events=" << nEvents;
    const double EffectiveRate = nEvents/total_time;
    LOG(INFO) << "Trigger rate " << EffectiveRate;

    i->SetLineColor(kRed);
    canvas("check") << h << g << samepad << i << endc;
    app->Run(kTRUE);

    return EXIT_SUCCESS;
}
