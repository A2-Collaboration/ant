#pragma once

#include <bitset>

namespace ant {

// heavily based on http://softwareengineering.stackexchange.com/a/338472

template<typename Enum, std::size_t number_of_bits = 32>
class bitflag
{
public:
    constexpr bitflag() = default;
    constexpr bitflag(Enum value) : bits(1 << static_cast<std::size_t>(value)) {}
    constexpr bitflag(const bitflag& other) : bits(other.bits) {}

    bool operator==(const bitflag& o) const { return bits == o.bits; }
    bool operator!=(const bitflag& o) const { return bits != o.bits; }

    bitflag operator|(Enum value) const { bitflag result = *this; result.bits |= 1 << static_cast<std::size_t>(value); return result; }
    bitflag operator&(Enum value) const { bitflag result = *this; result.bits &= 1 << static_cast<std::size_t>(value); return result; }
    bitflag operator^(Enum value) const { bitflag result = *this; result.bits ^= 1 << static_cast<std::size_t>(value); return result; }
    bitflag operator~() const { bitflag result = *this; result.bits.flip(); return result; }

    bitflag& operator|=(Enum value) { bits |= 1 << static_cast<std::size_t>(value); return *this; }
    bitflag& operator&=(Enum value) { bits &= 1 << static_cast<std::size_t>(value); return *this; }
    bitflag& operator^=(Enum value) { bits ^= 1 << static_cast<std::size_t>(value); return *this; }

    bitflag operator|(const bitflag& o) const { bitflag result = *this; result.bits |= o.bits; return result; }
    bitflag operator&(const bitflag& o) const { bitflag result = *this; result.bits &= o.bits; return result; }
    bitflag operator^(const bitflag& o) const { bitflag result = *this; result.bits ^= o.bits; return result; }

    bitflag& operator|=(const bitflag& o) { bits |= o.bits; return *this; }
    bitflag& operator&=(const bitflag& o) { bits &= o.bits; return *this; }
    bitflag& operator^=(const bitflag& o) { bits ^= o.bits; return *this; }

    constexpr bool any() const { return bits.any(); }
    constexpr bool all() const { return bits.all(); }
    constexpr bool none() const { return bits.none(); }
    constexpr explicit operator bool() const { return any(); }

    constexpr bool test(Enum value) const { return bits.test(static_cast<std::size_t>(value)); }
    void set(Enum value) { bits.set(static_cast<std::size_t>(value)); }
    void unset(Enum value) { bits.reset(static_cast<std::size_t>(value)); }

protected:
    std::bitset<number_of_bits> bits;
};

template<typename Enum, std::size_t n>
constexpr typename std::enable_if<std::is_enum<Enum>::value, bitflag<Enum, n> >::type
operator|(Enum left, const bitflag<Enum, n>& right)
{
    return bitflag<Enum, n>(left) | right;
}
template<typename Enum, std::size_t n>
constexpr typename std::enable_if<std::is_enum<Enum>::value, bitflag<Enum, n>>::type
operator&(Enum left, const bitflag<Enum, n>& right)
{
    return bitflag<Enum, n>(left) & right;
}
template<typename Enum, std::size_t n>
constexpr typename std::enable_if<std::is_enum<Enum>::value, bitflag<Enum, n>>::type
operator^(Enum left, const bitflag<Enum, n>& right)
{
    return bitflag<Enum, n>(left) ^ right;
}

}