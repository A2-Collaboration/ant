#pragma once


#include "base/Format.h"
#include "base/printable.h"

namespace ant {

struct TUnpackerMessage : printable_traits
{
    enum class Level_t : std::uint8_t {
        Info, Warn, DataError, DataDiscard, HardwareError
    };

    Level_t Level;
    std::string Message;
    std::vector<double> Payload;

    TUnpackerMessage(Level_t level,
                     const std::string& message) :
        Level(level),
        Message(message)
    {}
    TUnpackerMessage() {}
    virtual ~TUnpackerMessage() {}

    template<class Archive>
    void serialize(Archive& archive) {
        archive(Level, Message, Payload);
    }

    const char* LevelToString() const {
        return LevelToString(Level);
    }

    static const char* LevelToString(const Level_t& level) {
        switch(level) {
        case Level_t::Info:
            return "Info";
        case Level_t::Warn:
            return "Warn";
        case Level_t::DataError:
            return "DataError";
        case Level_t::DataDiscard:
            return "DataDiscard";
        case Level_t::HardwareError:
            return "HardwareError";
        }
        throw std::runtime_error("Not implemented");
    }


    std::string FormattedMessage() const {
        return fmt::format_vector(Message, Payload);
    }

    virtual std::ostream& Print( std::ostream& s) const override {
        return s << "TUnpackerMessage"
                 << " Level=" << LevelToString()
                 << " Msg='" << FormattedMessage() << "'";
    }


};

} // namespace ant
