#include "Setup.h"

#include "detectors/CB.h"


using namespace std;

namespace ant {
namespace expconfig {
namespace setup {


class Setup_CBSourceCalib : public Setup
{
public:

    Setup_CBSourceCalib(const std::string& name, OptionsPtr opt) : Setup(name, opt)
    {
        auto cb = make_shared<detector::CB>();
        AddDetector(cb);
    }
    bool Matches(const TID&) const override {
        return false;
    }
};

AUTO_REGISTER_SETUP(Setup_CBSourceCalib)
}
}
}
