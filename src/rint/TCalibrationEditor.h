#ifndef ANT_TCALIBRATIONEDITOR_H
#define ANT_TCALIBRATIONEDITOR_H
#include <string>
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
private:
    calibration::Editor* ed;
public:

    TCalibrationEditor();

    /// ONLY FOR DEBUG!!!!!
    void AddSomeRandomData();

    void ListCalibrations() const;

    void AddFromFile(const std::string& fileName);
    void SaveToFile(const std::string& fileName);

    void Remove(const std::string& calibrationID, const uint32_t& index);

    void ShowHistory(const std::string& calibrationID) const;
    /**
     * @brief ShowValid shows only calibration steps which are used
     * @param calibrationID
     */
    void ShowValid(const std::string&  calibrationID) const;

    virtual ~TCalibrationEditor();

    ClassDef(TCalibrationEditor, ANT_RINT_VERSION)
};

}

#endif
