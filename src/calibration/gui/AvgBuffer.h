#pragma once

#include <memory>
#include <list>
#include <queue>
#include <cassert>

#include "AvgBuffer_traits.h"

#include "base/interval.h"
#include "base/SavitzkyGolay.h"
#include "tree/TID.h"

#include "TH1.h"
#include "TArray.h"

namespace ant {
namespace calibration {
namespace gui {

// using AvgBufferItem_traits allows using other types just
// behaving similar to a histogram in terms of an AvgBuffer
template<typename AvgBufferItem>
struct AvgBufferItem_traits {};

// provide commonly used specialization for TH1
template<>
struct AvgBufferItem_traits<TH1> {
    static std::unique_ptr<TH1> Clone(const TH1& h) { return std::unique_ptr<TH1>(dynamic_cast<TH1*>(h.Clone())); }
    static void   Add(TH1& dest, const TH1& src) { dest.Add(std::addressof(src)); };
    static int    GetNBins(const TH1& h) { return dynamic_cast<const TArray&>(h).GetSize(); };
    static double GetBin(const TH1& h, int bin) { return h.GetBinContent(bin); }
    static void   SetBin(TH1& h, int bin, double v) { h.SetBinContent(bin, v); }
};


template<typename AvgBufferItem>
class AvgBuffer_Sum : public AvgBuffer_traits<AvgBufferItem> {
protected:
    using Traits = AvgBufferItem_traits<AvgBufferItem>;

    std::unique_ptr<AvgBufferItem> m_totalsum; // the total sum of all histograms
    interval<TID> m_range{TID(), TID()};
    bool flushed = false;
public:

    AvgBuffer_Sum() = default;
    virtual ~AvgBuffer_Sum() = default;

    void Push(std::shared_ptr<AvgBufferItem> h, const interval<TID>& id) override
    {
        // use m_totalsum to detect if this is the first push
        if(m_totalsum == nullptr) {
            m_range = id;
            m_totalsum = std::unique_ptr<AvgBufferItem>(Traits::Clone(*h));
            return;
        }

        Traits::Add(*m_totalsum, *h);
        m_range.Extend(id);
    }

    void Flush() override {
        flushed = true;
    }

    // indicate empty until not flushed (slurp all input files via Push())
    // or indicate empty after first Next() call (makes it one slice buffer)
    bool Empty() const override {
        return !flushed || !m_totalsum;
    }

    const AvgBufferItem& CurrentItem() const override {
        return *m_totalsum;
    }

    const interval<TID>& CurrentRange() const override {
        return m_range;
    }

    // first Next call invalidates this buffer
    void Next() override {
        m_totalsum = nullptr;
    }
};

template<typename AvgBufferItem>
class AvgBuffer_SavitzkyGolay : public AvgBuffer_traits<AvgBufferItem> {
protected:
    using Traits = AvgBufferItem_traits<AvgBufferItem>;

    struct buffer_entry {
        buffer_entry(const std::shared_ptr<AvgBufferItem>& h, const interval<TID>& ID) : hist(h), id(ID) {}
        std::shared_ptr<AvgBufferItem> hist;
        interval<TID> id;
    };


    struct buffer_t : std::list<buffer_entry> {
        using std::list<buffer_entry>::list;
        using const_iterator = typename std::list<buffer_entry>::const_iterator;
        const_iterator middle() const {
            return std::next(this->begin(), this->size()/2);
        }
    };

    buffer_t m_buffer; // buffered histograms for smoothing

    std::queue<buffer_entry> worklist;

    bool startup_done = false;
    const std::size_t m_sum_length;

    std::shared_ptr<AvgBufferItem> GetSmoothedClone(typename buffer_t::const_iterator i) const {
        // normalize the bin contents to length of run

        double normalization = i->id.Stop().Lower - i->id.Start().Lower;
        normalization /= this->total_length/this->total_n;
        // expect at least one event in range and identical timestamps
        // (otherwise length is hard to estimate here)
        if(i->id.Start().Timestamp != i->id.Stop().Timestamp || !(normalization > 0)) {
            normalization = 1.0;
        }

        const auto h = std::shared_ptr<AvgBufferItem>(Traits::Clone(*i->hist));

        // to get the number of cells (or total number of all bins)
        // this cast is necessary, as GetNcells is not there in current ROOT5 branch?!
        const auto nBins = Traits::GetNBins(*h);;

        // h is the destination of the smoothing

        // range is relative to i and inclusive, so take distance-to-end-1
        const interval<int> range(-std::distance(m_buffer.begin(), i),
                                  std::distance(i, m_buffer.end())-1);

        for(auto bin=0;bin<nBins;bin++) {
            auto getY = [i,bin,normalization] (const int i_) {
                return Traits::GetBin(*std::next(i, i_)->hist, bin)/normalization;
            };
            auto setY = [h,bin] (const double v) {
                Traits::SetBin(*h, bin, v);
            };
            sg.Convolute(getY, setY, range);
        }

        return h;
    }

    const SavitzkyGolay sg;

public:

    AvgBuffer_SavitzkyGolay(std::size_t length, std::size_t polorder) :
        m_sum_length(length),
        // always setup SG filter symmetric
        sg(m_sum_length, polorder)
    {
        if(m_sum_length<1)
            throw std::runtime_error("Average window size must be at least 1");
    }
    virtual ~AvgBuffer_SavitzkyGolay() = default;

    void Push(std::shared_ptr<AvgBufferItem> h, const interval<TID>& id) override
    {
        // add the item to the buffer
        m_buffer.emplace_back(buffer_entry(h, id));


        // pop elements from buffer
        if(m_buffer.size() > m_sum_length)
        {
            m_buffer.pop_front();
        }

        // check if sufficient size of buffer is reached
        if(m_buffer.size() >= m_sum_length) {
            if(!startup_done) {
                for(auto i=m_buffer.begin(); i!=m_buffer.middle(); ++i) {
                    worklist.emplace(GetSmoothedClone(i), i->id);
                }
                startup_done = true;
            }
            worklist.emplace(GetSmoothedClone(m_buffer.middle()),
                             m_buffer.middle()->id);
        }

        assert(m_buffer.size() <= m_sum_length);
    }

    void Flush() override {
        if(m_buffer.empty())
            return;

        for(auto it = std::next(m_buffer.middle()); it != m_buffer.end(); ++it) {
            worklist.emplace(GetSmoothedClone(it), it->id);
        }

        m_buffer.clear();
    }

    const AvgBufferItem& CurrentItem() const override {
        return *worklist.front().hist;
    }

    const interval<TID>& CurrentRange() const override {
        return worklist.front().id;
    }

    void Next() override {
        worklist.pop();
    }
    bool Empty() const override {
        return worklist.empty();
    }
};

}
}
}
