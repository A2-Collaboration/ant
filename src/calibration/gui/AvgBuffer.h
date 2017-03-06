#pragma once

#include <memory>
#include <list>
#include <queue>
#include <cassert>

namespace ant {
namespace calibration {
namespace gui {

template <typename HistType, typename IDType>
class AvgBuffer {
protected:

    struct buffer_entry {
        buffer_entry(const std::shared_ptr<HistType>& h, const IDType& ID) : hist(h), id(ID) {}
        std::shared_ptr<HistType> hist;
        IDType id;
    };


    struct buffer_t : std::list<buffer_entry> {
        using std::list<buffer_entry>::list;
        typename std::list<buffer_entry>::const_iterator middle() const {
            return std::next(this->begin(), this->size()/2);
        }
    };

    buffer_t m_buffer; // contains the histograms summed up in m_movingsum

    std::queue<buffer_entry> worklist;

    std::unique_ptr<HistType> m_movingsum; // the current moving sum of pushed histograms

    bool startup_done = false;
    const std::size_t m_sum_length;

    std::shared_ptr<HistType> GetMovingSumClone() const {
        return std::shared_ptr<HistType>(dynamic_cast<HistType*>(m_movingsum->Clone()));
    }

public:

    AvgBuffer(std::size_t sum_length) :  m_sum_length(sum_length) {}

    void Push(std::shared_ptr<HistType> h, const IDType& id)
    {
        if(m_movingsum == nullptr) {
            m_movingsum = std::unique_ptr<HistType>(dynamic_cast<HistType*>(h->Clone()));
        } else {
            m_movingsum->Add(h.get(), 1.0);
        }

        // special mode when m_sum_length==0 (no truely moving sum required)
        // just sum up the histograms, and only remember the ID span
        if(m_sum_length == 0) {
            if(m_buffer.empty()) {
                // we're just interested in the id
                m_buffer.emplace_back(buffer_entry(nullptr, id));
                startup_done = true; // for consistency
            }
            else {
                // extend the existing id
                m_buffer.front().id.Extend(id);
            }
            return;
        }

        // in normal mode, always add the item to the buffer
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

    void Flush() {
        if(m_buffer.empty())
            return;

        auto movingsum = GetMovingSumClone();

        if(m_sum_length == 0) {
            worklist.emplace(movingsum, m_buffer.front().id);
        }
        else {
            for(auto it = std::next(m_buffer.middle()); it != m_buffer.end(); ++it) {
                worklist.emplace(movingsum, it->id);
            }
        }
        m_buffer.clear();
    }

    const HistType& CurrentSum() { return *worklist.front().hist; }
    const IDType& CurrentID() const { return worklist.front().id; }

    void GotoNextID() { worklist.pop(); }
    bool Empty() const {return worklist.empty(); }

    std::size_t GetSumLength() const { return m_sum_length; }


};

}
}
}
