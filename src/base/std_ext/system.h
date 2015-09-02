#pragma once

namespace ant {
namespace std_ext {

struct system {

    /**
     * @brief Check if running an interactive shell / if a TTY is attached to stdin
     * @return
     */
    static bool isInteractive();
};
}
}
