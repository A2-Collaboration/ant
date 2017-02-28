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

    constexpr bitflag operator|(Enum value) const { bitflag result = *this; result.bits |= 1 << static_cast<std::size_t>(value); return result; }
    constexpr bitflag operator&(Enum value) const { bitflag result = *this; result.bits &= 1 << static_cast<std::size_t>(value); return result; }
    constexpr bitflag operator^(Enum value) const { bitflag result = *this; result.bits ^= 1 << static_cast<std::size_t>(value); return result; }
    constexpr bitflag operator~() const { bitflag result = *this; result.bits.flip(); return result; }

    constexpr bitflag& operator|=(Enum value) { bits |= 1 << static_cast<std::size_t>(value); return *this; }
    constexpr bitflag& operator&=(Enum value) { bits &= 1 << static_cast<std::size_t>(value); return *this; }
    constexpr bitflag& operator^=(Enum value) { bits ^= 1 << static_cast<std::size_t>(value); return *this; }

    constexpr bool any() const { return bits.any(); }
    constexpr bool all() const { return bits.all(); }
    constexpr bool none() const { return bits.none(); }
    constexpr operator bool() { return any(); }

    constexpr bool test(Enum value) const { return bits.test(1 << static_cast<std::size_t>(value)); }
    constexpr void set(Enum value) { bits.set(1 << static_cast<std::size_t>(value)); }
    constexpr void unset(Enum value) { bits.reset(1 << static_cast<std::size_t>(value)); }

private:
    std::bitset<number_of_bits> bits;
};

template<typename Enum>
constexpr typename std::enable_if<std::is_enum<Enum>::value, bitflag<Enum>>::type
operator|(Enum left, Enum right)
{
    return bitflag<Enum>(left) | right;
}
template<typename Enum>
constexpr typename std::enable_if<std::is_enum<Enum>::value, bitflag<Enum>>::type
operator&(Enum left, Enum right)
{
    return bitflag<Enum>(left) & right;
}
template<typename Enum>
constexpr typename std::enable_if_t<std::is_enum<Enum>::value, bitflag<Enum>>::type
operator^(Enum left, Enum right)
{
    return bitflag<Enum>(left) ^ right;
}


}