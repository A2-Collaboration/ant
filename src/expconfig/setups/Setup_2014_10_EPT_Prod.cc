#include "Setup_2014_EPT.h"

using namespace std;

namespace ant {
namespace expconfig {
namespace setup {

class Setup_2014_10_EPT_Prod : public Setup_2014_EPT
{
public:
    virtual std::string GetName() const override {
        return "Setup_2014-10_EPT_Prod";
    }

    Setup_2014_10_EPT_Prod()
        : Setup_2014_EPT(GetName())
    {
        /// \todo add ignored elements
    }


    bool Matches(const THeaderInfo& header) const override {
        if(!Setup_2014_EPT::Matches(header))
            return false;
        if(!std_ext::time_between(header.Timestamp, "2014-10-14", "2014-11-03"))
            return false;
        return true;
    }
};

// don't forget registration
AUTO_REGISTER_SETUP(Setup_2014_10_EPT_Prod)

}}} // namespace ant::expconfig::setup
