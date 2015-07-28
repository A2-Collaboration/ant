#pragma once

#include <fstream>
#include <string>
#include <memory>
#include <cstdint>
#include <vector>
#include <chrono>

namespace ant {

/**
 * @brief The RawFileReader class
 * Super-simple wrapper for reading binary files,
 * even if they are compressed :)
 * Possible IO errors are propagated as exceptions
 */
class RawFileReader {
public:

    /**
     * @brief OutputPerformanceStats
     *
     * Set to finite number to output performance stats every
     * given seconds for all RawFileReader's in use
     */
    static double OutputPerformanceStats;

    RawFileReader() :
        performanceBytesRead(-1), // nothing read at all
        performanceBytesRead_compressed(-1) // needed if there's underlying compression
    {}

    virtual ~RawFileReader();

    /**
   * @brief open
   * @param filename
   * @param inbufsize
   *
   * Parameter inbufsize is ignored if non-compressed data is read
   */
    void open(const std::string& filename, const size_t inbufsize = BUFSIZ);

    /**
   * @brief operator bool
   *
   * is false if the reader has a problem
   */
    explicit operator bool() const {
        return p->operator bool();
    }

    /**
   * @brief read n bytes to buffer s
   * @param s buffer (must be allocated by caller)
   * @param n
   */
    void read(char* s, std::streamsize n) {
        p->read(s,n);
        // track how much has been read if required
        HandlePerformanceStats();
    }

    void read(std::uint32_t* s, std::streamsize n) {
        read(reinterpret_cast<char*>(s), n*uint32_t_factor);
    }

    /**
   * @brief gcount
   * @return number of bytes read
   */
    std::streamsize gcount() const {
        return p->gcount();
    }

    /**
   * @brief eof
   * @return true if last read went beyond end of file.
   *
   * Note that reading the exact number of bytes makes eof() stay false
   */
    bool eof() const {
        return p->eof();
    }

    void expand_buffer(std::vector<std::uint32_t>& buffer, size_t totalSize) {
        if(buffer.size()>=totalSize)
            return;
        const std::streamsize toBeRead = totalSize - buffer.size();
        buffer.resize(totalSize); // make space in buffer
        read(std::addressof(buffer[totalSize-toBeRead]), toBeRead);
        if(uint32_t_factor*toBeRead != gcount()) {
            throw Exception("Not enough bytes available from file to expand buffer");
        }
    }

    class Exception : public std::runtime_error {
        using std::runtime_error::runtime_error; // use base class constructor
    };

private:
    static constexpr std::streamsize uint32_t_factor = sizeof(std::uint32_t)/sizeof(char);

    /**
   * @brief The PlainBase class
   *
   * Plain and simple filereader, encapsulating an ifstream
   * and providing a base class for the more complicated decompression classes
   *
   * Only the really needed methods are exported
   */
    class PlainBase {
    public:
        explicit PlainBase(const std::string& filename)
            : file(filename.c_str(), std::ios::binary)
        {
            const std::streampos begin = file.tellg();
            file.seekg (0, std::ios::end);
            filesize_ = file.tellg() - begin;
            file.seekg(0, std::ios::beg);
        }

        virtual ~PlainBase() = default;

        virtual explicit operator bool() const {
            return !file.operator!(); // some older ifstream version don't implement "operator bool"
        }

        virtual void read(char* s, std::streamsize n) {
            file.read(s, n);
            gcount_total_ += file.gcount();
        }

        virtual bool eof() const {
            return file.eof();
        }

        virtual std::streamsize gcount() const {
            return file.gcount();
        }

        // PlainBase has no compression layer,
        // indicate by returning negative value
        virtual std::streamsize gcount_compressed() const {
            return -1;
        }

        // used for estimating ETA in HandlePerformanceStats()
        virtual std::streamsize filesize_remaining() const {
            return filesize_ - gcount_total_;
        }

    private:
        std::ifstream file;
        std::streamsize filesize_;
        std::streamsize gcount_total_;
    }; // class RawFileReader::Plain

    /**
   * @brief The XZ class
   * Decompresses the bytes before giving them to the user
   *
   * Adapted from
   * http://git.tukaani.org/?p=xz.git;a=blob_plain;f=doc/examples/02_decompress.c;hb=HEAD
   *
   */
    class XZ : public PlainBase {
    public:

        XZ(const std::string& filename, const size_t inbufsize);

        virtual ~XZ();

        static bool test(std::ifstream& file);

        virtual explicit operator bool() const override {
            return PlainBase::operator bool() && !decompressFailed;
        }

        virtual void read(char *s, std::streamsize n) override;

        virtual std::streamsize gcount() const override {
            return gcount_;
        }

        virtual std::streamsize gcount_compressed() const {
            return gcount_compressed_;
        }

        virtual bool eof() const override {
            return eof_;
        }

    private:
        std::vector<uint8_t> inbuf;
        bool decompressFailed;
        std::streamsize gcount_;
        std::streamsize gcount_compressed_;
        bool eof_;

        template<typename T>
        using deleted_unique_ptr = std::unique_ptr<T, std::function<void(T*)>>;

        struct lzma_stream;
        deleted_unique_ptr<lzma_stream> strm;
        void init_decoder();

    }; // class RawFileReader::XZ


    // private stuff for RawFileReader
    std::unique_ptr<PlainBase> p;


    std::chrono::time_point<std::chrono::system_clock> lastPerformanceOutput;
    std::streamsize performanceBytesRead;
    std::streamsize performanceBytesRead_compressed;
    void HandlePerformanceStats();



}; // class RawFileReader



} // namespace ant
