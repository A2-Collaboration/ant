#include "OptionsList.h"

#include "Logger.h"
#include <cstdlib>

using namespace std;
using namespace ant;

OptionsList::OptionsList(std::shared_ptr<const OptionsList> Parent):
    parent(Parent)
{
}

void OptionsList::SetOption(const string& str, const string delim)
{
    const auto delimiter_pos = str.find(delim);
    if( delimiter_pos != str.npos) {
        const std::string key = str.substr(0, delimiter_pos);
        const std::string val = str.substr(delimiter_pos + delim.length(), str.npos);
        options.insert({key,val});
    } else {
        LOG(WARNING) << "Can't parse option string \"" << str << "\"";
    }
}

void OptionsList::SetOptions(const string& str,const string optdelim, const string valdelim)
{
    string::size_type p = 0;
    string::size_type np = 0;

    do {

        np = str.find(optdelim, p);

        SetOption(str.substr(p,np), valdelim);

        p = np+optdelim.size();

    } while(np != str.npos);
}

string OptionsList::GetOption(const string& key) const
{
    const auto entry = options.find(key);

    if(entry == options.end()) {
        if(parent) {
            return parent->GetOption(key);
        }
        return "";
    }

    return entry->second;

}

string OptionsList::Flatten() const
{
    auto strptrcmp = [] (const string* s1, const string* s2) { return *s1 == *s2; };
    std::map<const string*, const string*, decltype (strptrcmp)> m(strptrcmp);

    for(const auto& entry : options) {
        m[addressof(entry.first)] = addressof(entry.second);
    }

    stringstream s;
    for(const auto& e : m) {
        s << *e.first << "=" << *e.second <<":";
    }
    return s.str();
}

template<>
bool ant::OptionsList::Get<bool>(const string& key, const bool& def_value) const
{
    auto v = GetOption(key);
    std::transform(v.begin(), v.end(), v.begin(), ::tolower);

    if(v=="on" || v=="yes" || v=="true" || v=="1") {
        return true;
    }

    if(v=="off" || v=="no" || v=="false" || v=="0") {
        return false;
    }

    return def_value;
}
