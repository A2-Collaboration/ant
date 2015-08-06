#pragma once

#include "CalibrationDataManager.h"

#include <memory>
#include <string>

namespace ant
{
namespace calibration
{

class Editor
{
private:
    DataManager calman;
public:
    Editor(const std::string& filename): calman(filename){}
};

}
}
