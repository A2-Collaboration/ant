#include "system.h"

#include <unistd.h>
#include <stdio.h>

bool ant::std_ext::system::isInteractive()
{
        return isatty(fileno(stdin));
}
