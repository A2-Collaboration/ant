#include "system.h"

#include "string.h"

#include <unistd.h>
#include <dirent.h>
#include <cstdio>



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
