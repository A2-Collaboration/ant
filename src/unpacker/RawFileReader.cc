#include "RawFileReader.h"

#include "base/Logger.h"
#include "base/std_ext/memory.h"

#include <cstdio> // for BUFSIZ
#include <cstring> // for strerror
#include <limits>
#include <iomanip>

extern "C" {
#include <lzma.h>
}

using namespace std;
using namespace ant;

ant::RawFileReader::~RawFileReader() {}

double RawFileReader::PercentDone() const
{
    return double(p->pos()) / double(p->filesize_total());
}

void RawFileReader::open(const string& filename, const size_t inbufsize) {
    // open it as plain raw file
    ifstream file(filename.c_str());

    // already check here if the open was ok
    if(!file)
        throw Exception(string("Error when opening file ")
                        +filename
                        +": "
                        +string(strerror(errno)));

    if(XZ::test(file)) {
        p = std_ext::make_unique<XZ>(filename, inbufsize);
    }
    else {
        p = std_ext::make_unique<PlainBase>(filename);
    }

    progress = MakeProgressCounter();
}

RawFileReader::progress_t RawFileReader::MakeProgressCounter()
{
    // in future, there might be more than one compressed reader
    // so create the progress counter only according to p->gcount_compressed()
    ProgressCounter::Updater_t updater;
    if(p->gcount_compressed()<0) {
        // uncompressed reader
        updater = [this] (std::chrono::duration<double> elapsed_seconds) {
            const double bytes_per_s = (totalBytesRead - last_totalBytesRead)/elapsed_seconds.count();
            LOG(INFO) << "Reading file with "
                      << std::fixed << setprecision(3) << bytes_per_s/(1<<20)
                      << " MB/s";
            last_totalBytesRead = totalBytesRead;
        };
    }
    else {
        // compressed reader
        updater = [this] (std::chrono::duration<double> elapsed_seconds) {
            const double bytes_per_s_compressed = (totalBytesRead_compressed - last_totalBytesRead_compressed)
                                                  / elapsed_seconds.count();
            const double bytes_per_s = (totalBytesRead - last_totalBytesRead)/elapsed_seconds.count();

            LOG(INFO) << "Reading compressed/uncompressed file with "
                      << std::fixed << setprecision(3) << bytes_per_s_compressed/(1<<20)
                      << "/"
                      << std::fixed << setprecision(3) << bytes_per_s/(1<<20)
                      << " MB/s, ratio="
                      << std::fixed << setprecision(2) << 100*bytes_per_s_compressed/bytes_per_s
                      << " %";
            last_totalBytesRead = totalBytesRead;
            last_totalBytesRead_compressed = totalBytesRead_compressed;
        };
    }

    return std_ext::make_unique<ProgressCounter>(updater);
}

struct RawFileReader::XZ::lzma_stream : ::lzma_stream {};

RawFileReader::XZ::XZ(const std::string &filename, const size_t inbufsize) :
    PlainBase(filename),
    inbuf(inbufsize),
    decompressFailed(false),
    gcount_(0),
    gcount_compressed_(0),
    eof_(false),
    strm(new lzma_stream(),
         [] (lzma_stream* strm) { lzma_end(strm); delete strm; })
{
    init_decoder();
}

RawFileReader::XZ::~XZ() {}

bool RawFileReader::XZ::test(ifstream& file) {
    // check some magic bytes in the beginning
    // to determine possible compression
    std::vector<char> file_bytes(6);
    file.seekg(0, file.beg); // ensure start of file
    file.read(file_bytes.data(), file_bytes.size());


    vector<char> magic_bytes_xz{ static_cast<char>(0xFD), '7', 'z', 'X', 'Z', 0x00 };

    return file_bytes == magic_bytes_xz;
}

void RawFileReader::XZ::init_decoder()
{
    // using C-style init is a bit messy in C++
    using lzma_stream_pod = ::lzma_stream;
    auto ptr = reinterpret_cast<lzma_stream_pod*>(strm.get());
    *ptr = LZMA_STREAM_INIT;

    lzma_ret ret = lzma_stream_decoder(strm.get(), UINT64_MAX, LZMA_CONCATENATED);

    // Return successfully if the initialization went fine.
    if (ret == LZMA_OK) {
        return;
    }

    // otherwise throw some exceptions
    switch (ret) {
    case LZMA_MEM_ERROR:
        throw Exception("Memory allocation failed");
    case LZMA_OPTIONS_ERROR:
        throw Exception("Unsupported decompressor flags");
    default:
        throw Exception("Unknown error, possibly a bug");
    }
}

void RawFileReader::XZ::read(char* s, streamsize n) {

    lzma_action action = PlainBase::eof() ? LZMA_FINISH : LZMA_RUN;

    strm->next_out = reinterpret_cast<uint8_t*>(s);
    strm->avail_out = n;

    gcount_compressed_ = 0;

    while (true) {

        if (strm->avail_in == 0 && !PlainBase::eof()) {
            // read from the underlying compressed file into buffer
            PlainBase::read(reinterpret_cast<char*>(inbuf.data()), inbuf.size());
            strm->next_in = inbuf.data();
            strm->avail_in = PlainBase::gcount();
            gcount_compressed_ += strm->avail_in;

            if(PlainBase::eof()) {
                action = LZMA_FINISH;
            }
            else if (!PlainBase::operator bool()) {
                throw Exception(string("Error while reading from compressed input file: ")
                                +string(strerror(errno)));
            }
        }

        lzma_ret ret = lzma_code(strm.get(), action);

        if(ret == LZMA_STREAM_END) {
            gcount_ = n - strm->avail_out; // number of decompressed bytes
            if(strm->avail_out > 0)
                eof_ = true;
            return;
        }

        if(strm->avail_out == 0) {
            gcount_ = n;
            return;
        }

        if(ret == LZMA_OK)
            continue;

        decompressFailed = true;

        // ret indicates some error condition
        switch (ret) {
        case LZMA_MEM_ERROR:
            throw Exception("Memory allocation failed");
            break;

        case LZMA_FORMAT_ERROR:
            throw Exception("The input is not in the .xz format");
            break;

        case LZMA_OPTIONS_ERROR:
            throw Exception("Unsupported compression options");
            break;

        case LZMA_DATA_ERROR:
            throw Exception("Compressed file is corrupt");
            break;

        case LZMA_BUF_ERROR:
            throw Exception("Compressed file is truncated or "
                            "otherwise corrupt");

        default:
            throw Exception("Unknown error, possibly a bug");

        }
    }
}


