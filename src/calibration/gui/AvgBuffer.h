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

    std::unique_ptr<HistType> m_average;

    bool startup_done = false;

    std::queue<IDType> worklist;

public:

    using const_iterator =  typename shpBuffer::const_iterator;

    AvgBuffer(const std::size_t size):  midpos(m_buffer.begin()), m_max_size(size) {}

    void Push(std::shared_ptr<HistType> h, const IDType& id)
    {
        if(m_average == nullptr) {
            m_average = std::unique_ptr<HistType>(dynamic_cast<HistType*>(h->Clone()));
        } else {
            m_average->Add(h.get(), 1.0);
        }

        m_buffer.emplace_back(buffer_entry(h, id));

        if(isFull()) {
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
            midpos = m_buffer.begin();
        }



        if(m_buffer.size() > m_max_size)
        {
            Pop();
        }
        assert(m_buffer.size() <= m_max_size);
    }

    void Pop() {
        if(!m_buffer.empty()) {
            const std::shared_ptr<HistType>& last = m_buffer.front().hist;
            m_average->Add(last.get(), -1.0);
            m_buffer.pop_front();
        }
    }

    std::size_t size() const { return m_buffer.size(); }
    std::size_t max_size() const { return m_max_size; }

    bool isFull() const { return m_max_size <= size(); }

    HistType* Average() { return m_average.get(); }

    const_iterator begin() const { return m_buffer.begin(); }
    const_iterator end() const { return m_buffer.end(); }
    const_iterator mid() const { return midpos; }

    std::queue<IDType>& Worklist() { return worklist; }

    void PushRestToWorklist() {
        while(midpos != m_buffer.end()) {
            worklist.emplace(midpos->id);
            ++midpos;
        }
    }

};

}
}
}
