#pragma once

#include <type_traits>

namespace ant {

template <typename EnumClass, typename FieldType=unsigned>
struct enumfield {
private:
    FieldType field = 0;

    constexpr enumfield(const FieldType f) noexcept: field(f) {}
public:

    using EnumType = typename std::underlying_type<EnumClass>::type;

    constexpr enumfield(EnumClass e) noexcept: field(bit(e)) {}
    constexpr enumfield() noexcept {}

    enumfield(const enumfield&) = default;
    enumfield& operator=(const enumfield&) = default;

    enumfield(enumfield&&) = default;
    enumfield& operator=(enumfield&&) = default;

    ~enumfield() = default;

    constexpr bool isChecked(const EnumClass e) const noexcept {
        return field & bit(e);
    }

    static constexpr FieldType bit(const EnumClass e) noexcept {
        return 1 << static_cast<EnumType>(e);
    }

    void Set(const EnumClass e) noexcept {
        field |= bit(e);
    }

    void Reset(const EnumClass e) noexcept {
        field &= ~bit(e);
    }

    enumfield& operator&= (const enumfield e) {
        field &= e.field;
        return *this;
    }

    enumfield& operator|= (const enumfield e) {
        field |= e.field;
        return *this;
    }

    enumfield& operator^= (const enumfield e) {
        field ^= e.field;
        return *this;
    }

    constexpr bool operator== (const EnumClass e) const noexcept {
        return field == bit(e);
    }

    constexpr bool operator== (const enumfield& rhs) const noexcept {
        return field == rhs.field;
    }

    constexpr enumfield operator&(const enumfield& r) const noexcept {
        return enumfield(field & r.field);
    }

    constexpr enumfield operator|(const enumfield& r) const noexcept {
        return enumfield(field | r.field);
    }

    constexpr enumfield operator^(const enumfield& r) const noexcept {
        return enumfield(field ^ r.field);
    }

};

}
