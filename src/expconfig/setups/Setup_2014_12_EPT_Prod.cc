#include "Setup_2014_EPT.h"

using namespace std;

namespace ant {
namespace expconfig {
namespace setup {

class Setup_2014_12_EPT_Prod : public Setup_2014_EPT
{
public:

    Setup_2014_12_EPT_Prod(const std::string& name, SetupOptPtr opt)
        : Setup_2014_EPT(name, opt)
    {
        /// \todo add ignored elements
    }


    bool Matches(const THeaderInfo& header) const override {
        if(!Setup_2014_EPT::Matches(header))
            return false;
        if(!std_ext::time_between(header.Timestamp, "2014-12-02", "2014-12-22"))
            return false;
        return true;
    }
};

// don't forget registration
AUTO_REGISTER_SETUP(Setup_2014_12_EPT_Prod)

}}} // namespace ant::expconfig::setup
