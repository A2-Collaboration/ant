#pragma once

#include "analysis/physics/Physics.h"
#include "reconstruct/Reconstruct_traits.h"

#include <vector>

namespace ant {

class Physics;

class Calibration {
public:

    struct Converter {
        using ptr_t = std::shared_ptr<Converter>;

        virtual std::vector<double> Convert(const std::vector<uint8_t>& rawData) const = 0;
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
            public Updateable_traits
    {
    protected:
        Module(const std::string& name_) :
            BaseModule(name_)
        {}
    };

    class PhysicsModule : public Module
    {
    protected:
        using Module::Module;
    public:
        virtual std::unique_ptr<Physics> GetPhysicsModule() = 0;
    };

};

} // namespace ant

