#include "system.h"

#include "string.h"

#include <unistd.h>
#include <dirent.h>
#include <cstdio>
#include <cstring>
#include <stdexcept>



using namespace ant::std_ext;
using namespace std;

bool system::isInteractive()
{
    return isatty(fileno(stdin));
}

list<string> system::lsFiles(const string& folder, const string& extension,
                             bool ignoreDotDirs,
                             bool doNotPrependFolder)
{
    list<string> filenames;

    DIR* dir = opendir (folder.c_str());
    if(dir) {
        dirent* ent;
        while((ent = readdir(dir))) {
            string filename = ent->d_name;
            if(!std_ext::string_ends_with(filename,extension))
                continue;
            if(ignoreDotDirs && (filename == "." || filename == ".."))
                continue;
            filenames.emplace_back(doNotPrependFolder ? filename : folder+"/"+filename);
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

string system::absolutePath(const string& path, string cwd)
{
    if(path.empty())
        return "";

    if(path[0] == '/')
        return path;

    if(cwd.empty()) {
        char buf[PATH_MAX];
        if(getcwd(buf,sizeof(buf)) == NULL)
            throw runtime_error("Cannot get current working dir");
        cwd = buf;
    }

    return cwd+"/"+path;
}


bool ant::std_ext::system::testopen(const std::string& filename, std::string& errmsg)
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
