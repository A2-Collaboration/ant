#include "GitInfo.h"

#include "std_ext/system.h"
#include "std_ext/string.h"

#include "base/Paths.h"

#include <cstdio>
#include <sstream>

using namespace std;
using namespace ant;

GitInfo::GitInfo(const string& repofolder_) :
    repofolder(repofolder_.empty() ? ANT_PATH_GITREPO : repofolder_)
{

}

string GitInfo::GetUser()
{
    const string username = exec_git("config user.name");
    if(username.empty())
        return "";
    const string email = exec_git("config user.email");
    if(email.empty())
        return username;
    return std_ext::formatter() << username << " <" << email << ">";
}

string GitInfo::GetDescription()
{
    return exec_git("describe --always --dirty");
}



string GitInfo::exec_git(const string& args)
{
    if(!std_ext::system::testopen(repofolder))
        return "";
    stringstream ss_cmd;
    ss_cmd << "cd " << repofolder << " && git " << args;
    return std_ext::system::exec(ss_cmd.str());
}
