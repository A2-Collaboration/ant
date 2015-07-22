#ifndef ANT_TFILEWRAPPER_H
#define ANT_TFILEWRAPPER_H

#include <memory>
#include <string>

class TFile;

namespace ant {

class TFileWrapper {
protected:
    std::unique_ptr<TFile> file;

public:
    TFileWrapper(const std::string& filename);

    ///@fixme is this required?? remove if possible
    TFile* operator* () { return file.get(); }

    bool isOpen() const;
    void cd();

    ~TFileWrapper();
    TFileWrapper(const TFileWrapper&) = delete;
    TFileWrapper& operator= (const TFileWrapper&) = delete;
};
}

#endif
