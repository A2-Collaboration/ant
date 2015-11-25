#pragma once

namespace ant {
namespace std_ext {

template <typename To, typename From>
To clipNegative(const From& x) {
    return To( x >= 0 ? x : 0);
}

}
}
