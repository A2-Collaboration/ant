#pragma once

#include "analysis/physics/Physics.h"
#include "reconstruct/Reconstruct_traits.h"

namespace ant {

class Calibration {
public:

    struct Converter {
        using ptr_t = std::shared_ptr<Converter>;

        virtual vector<double> Convert(const vector<uint8_t>& rawData) const = 0;
        virtual ~Converter() = default;
    };

    class SimpleModule :
            public CalibrationApply_traits
    {
    public:
        std::string GetName() const { return name; }
        virtual ~SimpleModule() = default;
    protected:
        SimpleModule(const std::string& name_) :
            name(std::string("Calibration_")+name_)
        {}
    private:
        const std::string name;
    };

    /**
   * @brief The Calibration::Module class
   * A calibration module is two things:
   * * A physics class to make histograms and determine calibration factors
   * * and a CalibrationApply class that can apply those factors to data
   * * furthermore, it may update its constants during the analysis
   */
    class Module :
            public SimpleModule,
            public Physics,
            public Updateable_traits
    {
    protected:
        Module(const std::string& name_) :
            SimpleModule(name_),
            Physics(GetName()) /// \todo put calibration histograms in subdir?
        {}
    };

};

} // namespace ant

