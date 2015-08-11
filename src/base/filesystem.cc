#include "filesystem.h"

#include "TSystemDirectory.h"
#include "TList.h"
#include "TCollection.h"

using namespace ant;
using namespace std;

list<std::string> filesystem::lsFiles(const std::string& path, const std::string& extension)
{
    list<string> filenames;

    TSystemDirectory dir(path.c_str(), path.c_str());
    TList *files = dir.GetListOfFiles();
    if (files) {
        TSystemFile *file;
        TString fname;
        TIter next(files);
        while ((file=(TSystemFile*)next())) {
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(extension.c_str())) {
                filenames.emplace_back(path+"/"+fname);
            }
        }
    }
    return filenames;
}
