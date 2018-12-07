#!/bin/bash

export COVERAGE=true
export CURLCMD="curl -SLs --retry 10"
source <($CURLCMD https://raw.githubusercontent.com/A2-Collaboration/travis-container-packets/cmake/setup.sh)
source <($CURLCMD https://raw.githubusercontent.com/A2-Collaboration/travis-container-packets/gcc/setup.sh)
source <($CURLCMD https://raw.githubusercontent.com/A2-Collaboration/travis-container-packets/cern-root/setup.sh)
source <($CURLCMD https://raw.githubusercontent.com/A2-Collaboration/travis-container-packets/APLCONpp/setup.sh)
source <($CURLCMD https://raw.githubusercontent.com/A2-Collaboration/travis-container-packets/gsi-pluto/setup.sh)
