#pragma once

#include <string>
#include <memory>
#include <map>
#include <sstream>

namespace ant {

class OptionsList {
protected:

    std::shared_ptr<const OptionsList> parent;

    std::map<std::string, std::string> options;
    std::string GetOption(const std::string& key) const;

public:
    OptionsList(std::shared_ptr<const OptionsList> Parent=nullptr);

    void SetOption(const std::string& str, const std::string delim="=");
    void SetOptions(const std::string& str, const std::string optdelim=",", const std::string valdelim="=");

    std::string Flatten() const;

    template <typename T>
    T Get(const std::string& key, const T& def_value = T()) const {
        std::stringstream ss(GetOption(key));
        T ret;
        if(ss >> ret) {
            return ret;
        }
        return def_value;
    }
};

template<>
bool OptionsList::Get<bool>(const std::string& key, const bool& def_value) const;

}
