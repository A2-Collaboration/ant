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

    /**
     * @brief The BaseModule class has just a name
     */
    class BaseModule
    {
    public:
        std::string GetName() const { return name; }
        virtual ~BaseModule() = default;
    protected:
        BaseModule(const std::string& name_) :
            name(std::string("Calibration_")+name_)
        {}
    private:
        const std::string name;
    };

    /**
     * @brief The Module class is typically used for calibrations
     *
     * Note that each module has the freedom to implement ReconstructHook::DetectorReadHits
     * or ReconstructHook::Clusters or both.
     */
    class Module :
            public BaseModule,
            public Physics,
            public Updateable_traits
    {
    protected:
        Module(const std::string& name_) :
            BaseModule(name_),
            Physics(GetName()) /// \todo put calibration histograms in subdir?
        {}
    };

};

} // namespace ant

