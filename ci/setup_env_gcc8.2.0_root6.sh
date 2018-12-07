#!/bin/bash

export COVERAGE=false
export CURLCMD="curl -SLs --retry 10"
source <($CURLCMD https://raw.githubusercontent.com/A2-Collaboration/travis-container-packets/cmake/setup_3.13.sh)
source <($CURLCMD https://raw.githubusercontent.com/A2-Collaboration/travis-container-packets/gcc/setup_8.2.sh)
source <($CURLCMD https://raw.githubusercontent.com/A2-Collaboration/travis-container-packets/cern-root6/setup.sh)
source <($CURLCMD https://raw.githubusercontent.com/A2-Collaboration/travis-container-packets/APLCONpp/setup_gcc8.2.sh)
source <($CURLCMD https://raw.githubusercontent.com/A2-Collaboration/travis-container-packets/pluto6/setup.sh)
