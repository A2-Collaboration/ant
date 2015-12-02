#pragma once
#include "analysis/plot/root_draw.h"
#include "base/interval.h"

#include <memory>
#include <string>
#include <vector>

class TH2D;

namespace ant
{

class TCalibrationData;

namespace calibration
{

class Editor
{

private:
    const std::string filename;

public:
    std::shared_ptr<ant::TCalibrationData> cdata;


    Editor( const std::string& fileName);

    //File Operations
    void Save() const;

};

}
}
