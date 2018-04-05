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

    using options_t = std::map<std::string, option_t>;
    std::unique_ptr< options_t > options;
    using notfound_t = std::set<std::string>;
    std::unique_ptr< notfound_t > notfound;

    std::string GetOption(const std::string& key) const;

public:
    OptionsList(std::shared_ptr<const OptionsList> Parent=nullptr);

    void SetOption(const std::string& str, const std::string delim="=");
    void SetOptions(const std::string& str, const std::string optdelim=",", const std::string valdelim="=");

    /**
     * @brief HasOption
     * @param key
     * @return bool if the option is set
     */
    bool HasOption(const std::string& key) const;
    /**
     * @brief HasOptionStartsWith
     * @param key
     * @return bool if the option starts with the given string parameter key
     */
    bool HasOptionStartsWith(const std::string& key) const;
    /**
     * @brief OptionStartsWith
     * @param key
     * @return string with option key if it starts with the given string parameter key
     */
    std::string OptionStartsWith(const std::string& key) const;

    /**
     * @brief HasUnusedOptionStartsWith
     * @param key
     * @return bool if there's an unused option starting with the given string parameter key
     */
    bool HasUnusedOptionStartsWith(const std::string& key) const;
    /**
     * @brief UnusedOptionStartsWith
     * @param key
     * @return string unused option key if it starts with the given string parameter key
     */
    std::string UnusedOptionStartsWith(const std::string& key) const;

    /**
     * @brief Flatten
     * @return the instance as string representation
     */
    std::string Flatten() const;
    /**
     * @brief GetUnused
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
        typename std::remove_const<T>::type ret(def_value);
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
