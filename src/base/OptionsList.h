#pragma once

#include <string>
#include <memory>
#include <map>

namespace ant {

class OptionsList {
protected:

    std::shared_ptr<const OptionsList> parent;

    std::map<std::string, std::string> options;

public:
    OptionsList(std::shared_ptr<const OptionsList> Parent=nullptr);

    void SetOption(const std::string& str, const std::string delim="=");
    void SetOptions(const std::string& str, const std::string optdelim=",", const std::string valdelim="=");

    std::string GetOption(const std::string& key) const;
    bool IsFlagSet(const std::string& key) const;


};

}