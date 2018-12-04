# Get Doxygen
wget https://github.com/A2-Collaboration/travis-container-packets/blob/doxygen/doxygen-1.8.13.linux.bin.tar.gz?raw=true -O doxygen.tar.gz
mkdir doxygen
tar -xzf doxygen.tar.gz -C doxygen --strip-components=1
rm doxygen.tar.gz
export PATH=$(pwd)/doxygen/bin:$PATH
