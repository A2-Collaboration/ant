#include "DataReader.h"
#include "base/ReadTFiles.h"



ant::input::DataReader::DataReader(const std::shared_ptr<ReadTFiles>& rootfiles) :
    files(rootfiles)
{

}

ant::input::DataReader::~DataReader()
{

}

