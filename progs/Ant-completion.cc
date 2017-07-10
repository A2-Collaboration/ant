#include "expconfig/ExpConfig.h"
#include "analysis/physics/Physics.h"

#include <iostream>

using namespace std;
using namespace ant;

int doAnt(const vector<string>& args) {
    if(args.front() == "--physics") {
        for(const auto& name : analysis::PhysicsRegistry::GetList()) {
            cout << name << endl;
        }
        return EXIT_SUCCESS;
    }
    else if(args.front() == "--setup") {
        for(const auto& name : ExpConfig::Setup::GetNames()) {
            cout << name << endl;
        }
        return EXIT_SUCCESS;
    }
    else if(args.front() == "--calibration") {
        if(args.size()==2) {
            ExpConfig::Setup::SetByName(args.at(1));
            auto& setup = ExpConfig::Setup::Get();
            for(const auto& calibration : setup.GetCalibrations()) {
                cout << calibration->GetName() << endl;
            }
            return EXIT_SUCCESS;
        }
    }
    return EXIT_FAILURE;
}

int main(int argc, char** argv) {

    if(argc < 3)
        return EXIT_FAILURE;

    const string progname(argv[1]);
    const vector<string> args(argv+2, argv+argc);

    if(progname == "Ant")
        return doAnt(args);
    // may add more completion programs here

    return EXIT_FAILURE;

}
