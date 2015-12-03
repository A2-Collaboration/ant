#pragma once

#include "tree/TCalibrationData.h"
#include <string>

class TH2D;

namespace ant
{
namespace calibration
{

struct Editor
{
    const std::string filename;
    ant::TCalibrationData cdata;

    Editor( const std::string& fileName);

    void ResetData();

    void Save() const;
    void SaveAs(const std::string& currentFileName) const;
};

}
}
