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
        buffer_entry(std::shared_ptr<HistType>& h, const IDType& ID): hist(h), id(ID) {}
        std::shared_ptr<HistType> hist;
        IDType id;
    };

protected:

    using shpBuffer = std::list<buffer_entry>;

    shpBuffer           m_buffer;
    typename shpBuffer::const_iterator midpos;

    std::size_t         m_max_size;

    std::unique_ptr<HistType> m_movingsum;

    bool startup_done = false;

    std::queue<IDType> worklist;

public:

    using const_iterator =  typename shpBuffer::const_iterator;

    AvgBuffer(const std::size_t size):  midpos(m_buffer.begin()), m_max_size(size) {}

    void Push(std::shared_ptr<HistType> h, const IDType& id)
    {
        if(m_movingsum == nullptr) {
            m_movingsum = std::unique_ptr<HistType>(dynamic_cast<HistType*>(h->Clone()));
        } else {
            m_movingsum->Add(h.get(), 1.0);
        }

        m_buffer.emplace_back(buffer_entry(h, id));

        // check if max_size is reached
        if(m_max_size <= m_buffer.size()) {
            if(!startup_done) {
                auto pos = m_buffer.begin();
                for(unsigned i=0;i<m_max_size/2;++i) {
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
        if(m_buffer.size() > m_max_size)
        {
            const std::shared_ptr<HistType>& last = m_buffer.front().hist;
            m_movingsum->Add(last.get(), -1.0);
            m_buffer.pop_front();
        }

        assert(m_buffer.size() <= m_max_size);
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



};

}
}
}
