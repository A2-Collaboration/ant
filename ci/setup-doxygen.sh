# Get Doxygen
wget http://ftp.stack.nl/pub/users/dimitri/doxygen-1.8.13.linux.bin.tar.gz -O doxygen.tar.gz
mkdir doxygen
tar -xzf doxygen.tar.gz -C doxygen --strip-components=1
#rm doxygen.tar.gz
export PATH=$(pwd)/doxygen/bin:$PATH
