#include "Setup_2014_EPT.h"

using namespace std;

namespace ant {
namespace expconfig {
namespace setup {

class Setup_2014_07_EPT_Prod : public Setup_2014_EPT
{
public:
    virtual std::string GetName() const override {
        return "Setup_2014-07_EPT_Prod";
    }

    Setup_2014_07_EPT_Prod()
        : Setup_2014_EPT(GetName())
    {
        /// \todo add ignored elements
    }


    bool Matches(const THeaderInfo& header) const override {
        if(!Setup_2014_EPT::Matches(header))
            return false;
        if(!std_ext::time_between(header.Timestamp, "2014-07-29", "2014-08-25"))
            return false;
        return true;
    }
};

// don't forget registration
AUTO_REGISTER_SETUP(Setup_2014_07_EPT_Prod)

}}} // namespace ant::expconfig::setup
