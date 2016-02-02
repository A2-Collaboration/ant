#pragma once

#include <string>
#include <memory>
#include <map>
#include <vector>
#include <sstream>

namespace ant {

class OptionsList;

using OptionsPtr = std::shared_ptr<const OptionsList>;

class OptionsList {
protected:

    OptionsPtr parent;

    // options are stored as key=value, with additional bool
    // to indicate that this Option was used by someones
    // see GetUnused()
    mutable std::map<std::string, std::pair<bool, std::string> > options;
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
    std::vector<std::string> GetUnused() const;

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
