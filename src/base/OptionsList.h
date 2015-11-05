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

    template <typename T>
    T Get(const std::string& key, const T& def_value) const {
        const auto& v = GetOption(key);
        if(v.empty()) {
            return def_value;
        }
        return T(v);
    }


};

template<>
double OptionsList::Get<double>(const std::string& key, const double& def_value) const;

template<>
int OptionsList::Get<int>(const std::string& key, const int& def_value) const;

}
