#pragma once

#include <string>
#include <memory>
#include <map>
#include <set>
#include <sstream>

namespace ant {

class OptionsList;

using OptionsPtr = std::shared_ptr<const OptionsList>;

class OptionsList {
protected:

    OptionsPtr parent;

    // options are stored as key=value, with additional bool
    // to indicate that this Option was used by some client
    // (physics classes or setups)
    // see GetUnused() and GetNotFound()
    // this makes debugging missing or wrongly specified options
    // much easier
    struct option_t {
        std::string Value;
        bool Used = false;
    };

    mutable std::map<std::string, option_t > options;
    mutable std::set<std::string> notfound;

    std::string GetOption(const std::string& key) const;

public:
    OptionsList(std::shared_ptr<const OptionsList> Parent=nullptr);

    void SetOption(const std::string& str, const std::string delim="=");
    void SetOptions(const std::string& str, const std::string optdelim=",", const std::string valdelim="=");

    /**
     * @brief Flatten
     * @return the instance as string representation
     */
    std::string Flatten() const;
    /**
     * @brief GetUnconsumed
     * @return string representations of options which were not used by GetOption
     */
    std::set<std::string> GetUnused() const;
    std::set<std::string> GetNotFound() const;

    /**
     * @brief main method to retrieve options, already parsed to requested type
     * @return the parsed value on success, or the default value
     */
    template <typename T>
    T Get(const std::string& key, const T& def_value = T()) const {
        std::stringstream ss(GetOption(key));
        T ret(def_value);
        if(ss >> ret) {
            return ret;
        }
        return def_value;
    }
};

// bool has some specialization
template<>
bool OptionsList::Get<bool>(const std::string& key, const bool& def_value) const;

}
