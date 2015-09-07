#include "system.h"

#include "string.h"

#include <unistd.h>
#include <dirent.h>
#include <cstdio>
#include <cstring>



using namespace ant::std_ext;
using namespace std;

bool system::isInteractive()
{
    return isatty(fileno(stdin));
}

list<string> system::lsFiles(const string& folder, const string& extension)
{
    list<string> filenames;

    DIR* dir = opendir (folder.c_str());
    if(dir) {
        dirent* ent;
        while((ent = readdir(dir))) {
            string filename = ent->d_name;
            if(!std_ext::string_ends_with(filename,extension))
                continue;
            filenames.emplace_back(folder+"/"+filename);
        }
        closedir(dir);
    }
    return filenames;
}

string system::exec(const string& cmd)
{
    string cmd_silent = cmd + " 2>/dev/null";
    FILE* pipe = popen(cmd_silent.c_str(), "r");
    if (!pipe) return "";
    char buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
        if(fgets(buffer, 128, pipe) != NULL)
            result += buffer;
    }
    pclose(pipe);
    result = std_ext::string_sanitize(result.c_str());
    return result;
}


bool ant::std_ext::system::testopen(const std::string& filename, std::string errmsg)
{
    FILE* fp = fopen(filename.c_str(),"r");
    if(fp==NULL) {
        errmsg = strerror(errno);
        return false;
    }
    fclose(fp);
    errmsg = "";
    return true;
}
