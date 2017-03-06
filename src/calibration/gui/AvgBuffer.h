#pragma once

#include <memory>
#include <list>
#include <queue>
#include <cassert>

#include "base/interval.h"
#include "tree/TID.h"

#include "TH1.h"

namespace ant {
namespace calibration {
namespace gui {

class AvgBuffer_traits {
public:
    virtual void Push(std::shared_ptr<TH1> hist, const interval<TID>& range) =0;
    virtual bool Empty() const =0;
    virtual void Flush() =0;
    virtual void Next() =0;

    virtual const TH1& CurrentHist() const =0;
    virtual const interval<TID>& CurrentRange() const =0;

    virtual ~AvgBuffer_traits() = default;
};

class AvgBuffer_JustSum : public AvgBuffer_traits {
protected:
    std::unique_ptr<TH1> m_totalsum; // the total sum of all histograms
    interval<TID> m_range{TID(), TID()};
    bool flushed = false;
public:

    AvgBuffer_JustSum() = default;
    virtual ~AvgBuffer_JustSum() = default;

    void Push(std::shared_ptr<TH1> h, const interval<TID>& id) override
    {
        // use m_totalsum to detect if this is the first push
        if(m_totalsum == nullptr) {
            m_range = id;
            m_totalsum = std::unique_ptr<TH1>(dynamic_cast<TH1*>(h->Clone()));
            return;
        }

        m_totalsum->Add(h.get(), 1.0);
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

    const TH1& CurrentHist() const override {
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

class AvgBuffer_MovingWindow : public AvgBuffer_traits {
protected:

    struct buffer_entry {
        buffer_entry(const std::shared_ptr<TH1>& h, const interval<TID>& ID) : hist(h), id(ID) {}
        std::shared_ptr<TH1> hist;
        interval<TID> id;
    };


    struct buffer_t : std::list<buffer_entry> {
        using std::list<buffer_entry>::list;
        typename std::list<buffer_entry>::const_iterator middle() const {
            return std::next(this->begin(), this->size()/2);
        }
    };

    buffer_t m_buffer; // contains the histograms summed up in m_movingsum

    std::queue<buffer_entry> worklist;

    std::unique_ptr<TH1> m_movingsum; // the current moving sum of pushed histograms

    bool startup_done = false;
    const std::size_t m_sum_length;

    std::shared_ptr<TH1> GetMovingSumClone() const {
        return std::shared_ptr<TH1>(dynamic_cast<TH1*>(m_movingsum->Clone()));
    }

public:

    AvgBuffer_MovingWindow(std::size_t sum_length) :
        m_sum_length(sum_length) {
        if(sum_length<1)
            throw std::runtime_error("Average window size must be at least 1");
    }
    virtual ~AvgBuffer_MovingWindow() = default;

    void Push(std::shared_ptr<TH1> h, const interval<TID>& id) override
    {
        if(m_movingsum == nullptr) {
            m_movingsum = std::unique_ptr<TH1>(dynamic_cast<TH1*>(h->Clone()));
        } else {
            m_movingsum->Add(h.get(), 1.0);
        }

        // add the item to the buffer
        m_buffer.emplace_back(buffer_entry(h, id));


        // pop elements from buffer
        if(m_buffer.size() > m_sum_length)
        {
            const auto& last = m_buffer.front().hist;
            m_movingsum->Add(last.get(), -1.0);
            m_buffer.pop_front();
        }

        // check if max_size is reached
        if(m_buffer.size() >= m_sum_length) {
            auto movingsum = GetMovingSumClone();
            if(!startup_done) {
                for(auto i=m_buffer.begin(); i!=m_buffer.middle(); ++i) {
                    worklist.emplace(movingsum, i->id);
                }
                startup_done = true;
            }
            worklist.emplace(movingsum, m_buffer.middle()->id);
        }

        assert(m_buffer.size() <= m_sum_length);
    }

    void Flush() override {
        if(m_buffer.empty())
            return;

        auto movingsum = GetMovingSumClone();

        for(auto it = std::next(m_buffer.middle()); it != m_buffer.end(); ++it) {
            worklist.emplace(movingsum, it->id);
        }

        m_buffer.clear();
    }

    const TH1& CurrentHist() const override {
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
