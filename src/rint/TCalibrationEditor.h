#ifndef ANT_TCALIBRATIONEDITOR_H
#define ANT_TCALIBRATIONEDITOR_H
#include <string>
#include <map>
#include <functional>
#include "Rtypes.h"

namespace ant
{
namespace calibration
{
class Editor;
}
}

#define ANT_RINT_VERSION 1

namespace ant {

class TCalibrationEditor
{
#ifndef __CINT__

    // This Map defines the functions the Editor should have
    using CommandMap = std::map<std::string,std::function<int()>>;

private:
    calibration::Editor* ed;
    std::string currentCalibration;

    CommandMap cmds;

    int chcal();
    int listcal();
    int show();
    int remove();
    int removeRange();
    int exit();
    int addd();


#endif


public:

    TCalibrationEditor();

    /// ONLY FOR DEBUG!!!!!
    void AddSomeRandomData();

    void ListCalibrations() const;
    void ListCommands();

    void AddFromFile(const std::string& fileName);
    void SaveToFile(const std::string& fileName);

    void Remove(const std::string& calibrationID, const std::uint32_t& index);
    void Remove(const std::string& calibrationID, const std::uint32_t& index1, const std::uint32_t& index2);
    void ReduceToValid(const std::string& calibrationID);

    void ShowHistory(const std::string& calibrationID) const;
    /**
     * @brief ShowValid shows only calibration steps which are used
     * @param calibrationID
     */
    void ShowValid(const std::string&  calibrationID) const;

    int Loop();

    int Execute(const std::string& command);

    virtual ~TCalibrationEditor();

    ClassDef(TCalibrationEditor, ANT_RINT_VERSION)
};

}

#endif
