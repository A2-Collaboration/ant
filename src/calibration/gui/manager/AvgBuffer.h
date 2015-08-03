#pragma once

#include <memory>
#include <list>

namespace ant {
namespace calibration {
namespace gui {

template <typename HistType, typename IDType>
class AvgBuffer {
protected:

    struct buffer_entry {
        buffer_entry(std::shared_ptr<HistType>& h, const IDType& ID): hist(h), id(ID) {}
        std::shared_ptr<HistType> hist;
        IDType id;
    };

    using shpBuffer = std::list<buffer_entry>;

    shpBuffer           m_buffer;

    std::size_t         m_max_size;

    std::unique_ptr<HistType> m_average;

public:

    using const_iterator =  typename shpBuffer::const_iterator;


    AvgBuffer(const std::size_t size): m_max_size(size) {}

    void Push(std::shared_ptr<HistType> h, const IDType& id)
    {
        if(m_average.get() == nullptr) {
            m_average = std::unique_ptr<HistType>( static_cast<HistType*>(h->Clone()) );
        } else {
            m_average->Add(h.get(), 1.0);
        }

        m_buffer.emplace_back(buffer_entry(h, id));

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

    bool isFull() const { return m_max_size == size(); }

    HistType* Average() { return m_average.get(); }

    const_iterator begin() const { return m_buffer.begin(); }
    const_iterator end() const { return m_buffer.end(); }

};

}
}
}
