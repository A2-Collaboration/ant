#pragma once

#include <memory>
#include <list>
#include <stack>

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

    std::stack<IDType> worklist;

public:

    using const_iterator =  typename shpBuffer::const_iterator;


    AvgBuffer(const std::size_t size):  midpos(m_buffer.end()), m_max_size(size) {}

    void Push(std::shared_ptr<HistType> h, const IDType& id)
    {
        if(m_average == nullptr) {
            m_average = std::unique_ptr<HistType>(dynamic_cast<HistType*>(h->Clone()));
        } else {
            auto x = h.get();
            m_average->Add(x, 1.0);
        }

        m_buffer.emplace_back(buffer_entry(h, id));

        if(isFull()) {
            if(!startup_done) {

                auto pos = m_buffer.begin();
                for(unsigned i=0;i<m_max_size/2;++i) {
                    worklist.emplace(pos->id);
                }

                midpos = pos;

                startup_done = true;

            } else if(startup_done) {
                worklist.emplace(midpos->id);
                ++midpos;
            }
        }

        if(m_buffer.size() > m_max_size)
        {
            Pop();
        }

    }

    void Pop() {
        if(!m_buffer.empty()) {

            std::shared_ptr<HistType>& last = m_buffer.front().hist;
            m_average->Add(last.get(), -1.0);
            m_buffer.pop_front();
        }
    }

    std::size_t size() const { return m_buffer.size(); }
    std::size_t max_size() const { return m_max_size; }

    bool isFull() const { return m_max_size >= size(); }

    HistType* Average() { return m_average.get(); }

    const_iterator begin() const { return m_buffer.begin(); }
    const_iterator end() const { return m_buffer.end(); }
    const_iterator mid() const { return midpos; }

    std::stack<IDType>& Worklist() { return worklist; }

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
