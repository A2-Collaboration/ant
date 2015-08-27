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
public:
    struct buffer_entry {
        buffer_entry(const std::shared_ptr<HistType>& h, const IDType& ID): hist(h), id(ID) {}
        std::shared_ptr<HistType> hist;
        IDType id;
    };

protected:

    using shpBuffer = std::list<buffer_entry>;

    shpBuffer           m_buffer;
    typename shpBuffer::const_iterator midpos;

    std::size_t         m_sum_length;

    std::unique_ptr<HistType> m_movingsum;

    bool startup_done = false;

    std::queue<IDType> worklist;

public:

    using const_iterator =  typename shpBuffer::const_iterator;

    AvgBuffer(const std::size_t sum_length):  midpos(m_buffer.begin()), m_sum_length(sum_length) {}

    void Push(std::shared_ptr<HistType> h, const IDType& id)
    {
        if(m_movingsum == nullptr) {
            m_movingsum = std::unique_ptr<HistType>(dynamic_cast<HistType*>(h->Clone()));
        } else {
            m_movingsum->Add(h.get(), 1.0);
        }

        // special mode when m_max_size==0 (no moving sum required)
        // just sum up the histograms, and only remember the ID span
        if(m_sum_length == 0) {
            if(m_buffer.empty()) {
                // we're just interested in the id
                m_buffer.emplace_back(buffer_entry(nullptr, id));
                midpos = m_buffer.begin();
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

        // check if max_size is reached
        if(m_sum_length <= m_buffer.size()) {
            if(!startup_done) {
                auto pos = m_buffer.begin();
                for(unsigned i=0;i<m_sum_length/2;++i) {
                    worklist.emplace(pos->id);
                    ++pos;
                }
                midpos = pos;
                startup_done = true;
            }
            else {
                worklist.emplace(midpos->id);
                ++midpos;
            }
        }
        else if(!startup_done) {
            // if startup not done yet, keep midpos at start of buffer
            midpos = m_buffer.begin();
        }


        // pop elements from buffer
        if(m_buffer.size() > m_sum_length)
        {
            const std::shared_ptr<HistType>& last = m_buffer.front().hist;
            m_movingsum->Add(last.get(), -1.0);
            m_buffer.pop_front();
        }

        assert(m_buffer.size() <= m_sum_length);
    }

    void Finish() {
        while(midpos != m_buffer.end()) {
            worklist.emplace(midpos->id);
            ++midpos;
        }
    }

    HistType* CurrentSum() { return m_movingsum.get(); }

    const IDType& CurrentID() const { return worklist.front(); }
    void GotoNextID() { worklist.pop(); }
    bool Empty() const {return worklist.empty(); }
    std::size_t GetSumLength() const { return m_sum_length; }


};

}
}
}
