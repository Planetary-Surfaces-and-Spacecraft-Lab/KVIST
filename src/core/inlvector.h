#include <vector>
#include <exception>

// Ian Numerics Library
// Based on https://codereview.stackexchange.com/questions/211241/simple-re-implementation-of-stdvector
// The intention behind this vector is that it allows a scalar function to
// initialize a vector better than std::vector
// Note: Because c++ templates are dumb - this has to be in it's own header

template <class Ty>
class memory_block
{
public:
    using pointer = Ty *;
    using const_pointer = Ty const *;

    memory_block() : m_size(std::size_t{0}),
                     m_buffer(nullptr) {}

    explicit memory_block(std::size_t count) : m_size(count),
                                               m_buffer(allocate(m_size)) {}

    memory_block(memory_block const &other) : m_size(other.m_size),
                                              m_buffer(allcoate(m_size)) {}

    memory_block(memory_block &&other) : m_size(std::move(other.m_size)),
                                         m_buffer(std::move(other.m_buffer))
    {
        other.m_size = std::size_t{0};
        other.m_buffer = nullptr;
    }

    ~memory_block()
    {
        deallocate(m_buffer);
    }

    void swap(memory_block &other);

    memory_block &operator=(memory_block const &other);
    memory_block &operator=(memory_block &&other);

    std::size_t size() const;

    pointer data()
    {
        return m_buffer;
    }

    const_pointer data() const
    {
        return m_buffer;
    }

    pointer at(std::size_t index)
    {
        assert(index < m_size); // maybe throw instead
        return m_buffer + index;
    }

    const_pointer at(std::size_t index) const
    {
        assert(index < m_size); // maybe throw instead
        return m_buffer + index;
    }

private:
    static pointer allocate(std::size_t size)
    {
        return static_cast<pointer>(::operator new(sizeof(Ty) * size));
    }

    static void deallocate(pointer buffer)
    {
        delete buffer; // static_cast<void*>(buffer);
    }

    std::size_t m_size;
    Ty *m_buffer;
};

template <class Ty>
class inlvector
{
public:
    inlvector<Ty>();
    explicit inlvector<Ty>(std::size_t count);
    inlvector<Ty>(std::size_t count, const Ty &val);

    inlvector<Ty>(const inlvector &other);
    inlvector<Ty>(inlvector &&other);

    ~inlvector<Ty>();

    void swap(inlvector &other);

    inlvector<Ty> &operator=(const inlvector<Ty> &other);
    inlvector<Ty> &operator=(inlvector<Ty> &&other);

    Ty &operator[](const std::size_t idx);
    const Ty &operator[](const std::size_t idx) const;

    std::vector<Ty> stdvec();

    std::size_t size() const;
    std::size_t capacity() const;

    void push_back(const Ty &val);
    void pop_back();

private:
    static const std::size_t M_INITIAL_SIZE = std::size_t{10};

    std::size_t m_size;
    memory_block<Ty> m_buffer;

    void grow(std::size_t amount);
    void reallocate(std::size_t min_size);

    void construct(std::size_t index, const Ty &value);
    void destruct(std::size_t index);
    void destruct_all();
    Ty &get(std::size_t idx);
};



template <class Ty>
inline void swap(memory_block<Ty> &a, memory_block<Ty> &b)
{
    a.swap(b);
}

template <class Ty>
inline void swap(inlvector<Ty> &a, inlvector<Ty> &b)
{
    a.swap(b);
}

template <class Ty>
void memory_block<Ty>::swap(memory_block<Ty> &other)
{
    std::swap(m_size, other.m_size);
    std::swap(m_buffer, other.m_buffer);
}

template <class Ty>
memory_block<Ty>& memory_block<Ty>::operator=(memory_block<Ty> const &other)
{
    auto temp = memory_block(other);
    swap(temp);
    return *this;
}

template <class Ty>
memory_block<Ty>& memory_block<Ty>::operator=(memory_block<Ty> &&other)
{
    auto temp = memory_block(std::move(other));
    swap(temp);
    return *this;
}

template <class Ty>
std::size_t memory_block<Ty>::size() const
{
    return m_size;
}


// Ian Numerics Library
// Based on https://codereview.stackexchange.com/questions/211241/simple-re-implementation-of-stdvector
template <class Ty>
inlvector<Ty>::inlvector() : m_size(0u),
                            m_buffer(M_INITIAL_SIZE)
{
}

template <class Ty>
inlvector<Ty>::inlvector(std::size_t count) : m_size(count),
                                             m_buffer(m_size)
{
    std::uninitialized_value_construct_n(m_buffer.data(), m_size); // value construct each element w/ placement new (C++17)
}

template <class Ty>
inlvector<Ty>::inlvector(std::size_t count, const Ty &value) : m_size(count),
                                                              m_buffer(m_size)
{
    std::uninitialized_fill_n(m_buffer.data(), m_size, value); // copy construct each element w/ placement new
}

template <class Ty>
inlvector<Ty>::inlvector(const inlvector &other) : m_size(other.m_size),
                                               m_buffer(m_size) // note: allocates only what we need to contain the existing elements, not the same as the capacity of the other buffer
{
    std::uninitialized_copy_n(other.m_buffer.data(), other.m_size, m_buffer.data()); // copy construct each element from old buffer to new buffer w/ placement new
}

template <class Ty>
inlvector<Ty>::inlvector(inlvector &&other) : m_size(std::move(other.m_size)),
                                          m_buffer(std::move(other.m_buffer)) // take ownership of the buffer
{
    other.m_size = std::size_t{0}; // other vector is now empty (nothing needs to be constructed / destructed)
}

template <class Ty>
inlvector<Ty>::~inlvector()
{
    destruct_all();
}

template <class Ty>
void inlvector<Ty>::swap(inlvector &other)
{
    using std::swap;
    swap(m_size, other.m_size);
    swap(m_buffer, other.m_buffer);
}

template <class Ty>
Ty& inlvector<Ty>::get(std::size_t idx)
{
    auto b = m_buffer.data();
    return b[idx];
}

template<class Ty>
std::vector<Ty> inlvector<Ty>::stdvec()
{
    std::vector<Ty> ret(size());
    for(int i=0; i<size(); ++i) {
        Ty val = this->get(i);
        ret.push_back(val);
    }
    return ret;
}

template <class Ty>
inlvector<Ty> &inlvector<Ty>::operator=(const inlvector &other)
{
    auto temp = inlvector(other);
    swap(temp);
    return *this;
}

template <class Ty>
inlvector<Ty> &inlvector<Ty>::operator=(inlvector<Ty> &&other)
{
    auto temp = inlvector(std::move(other));
    swap(temp);
    return *this;
}

template <class Ty>
std::size_t inlvector<Ty>::size() const
{
    return m_size;
}

template <class Ty>
Ty &inlvector<Ty>::operator[](const std::size_t idx)
{
    if (idx >= size())
    {
        throw std::range_error("index out of range");
    }
    return this->get(idx);
}

template <class Ty>
const Ty &inlvector<Ty>::operator[](const std::size_t idx) const
{
    if (idx >= size())
    {
        throw std::range_error("index out of range");
    }
    return this->get(idx);
}

template <class Ty>
std::size_t inlvector<Ty>::capacity() const
{
    return m_buffer.size();
}

template <class Ty>
void inlvector<Ty>::push_back(const Ty &value)
{
    grow(std::size_t{1});
    construct(m_size, value);
    ++m_size;
}

template <class Ty>
void inlvector<Ty>::pop_back()
{
    assert(m_size > 0); // maybe throw instead

    destruct(m_size - 1);
    --m_size;
}

template <class Ty>
void inlvector<Ty>::grow(std::size_t amount)
{
    if (m_buffer.size() - m_size < amount)
        reallocate(m_size + amount);
}

template <class Ty>
void inlvector<Ty>::reallocate(std::size_t min_size)
{
    assert(min_size > m_size);

    auto capacity = std::max(min_size, m_buffer.size() + std::max(m_buffer.size() / std::size_t{2}, std::size_t{1})); // growth factor of 1.5ish

    auto buffer = memory_block<Ty>(capacity);
    std::uninitialized_copy_n(m_buffer.data(), m_size, buffer.data()); // copy each element from old buffer to new buffer w/ placement new

    destruct_all(); // clean up the old buffer (call destructors on each of the old elements)

    m_buffer = std::move(buffer); // take ownership of the new buffer
}

template <class Ty>
void inlvector<Ty>::construct(std::size_t index, const Ty &value)
{
    new (m_buffer.at(index)) Ty(value); // placement new w/ copy constructor
}

template <class Ty>
void inlvector<Ty>::destruct(std::size_t index)
{
    assert(index < m_size);

    m_buffer.at(index)->~Ty(); // explictly call destructor (because we used placement new)
}

template <class Ty>
void inlvector<Ty>::destruct_all()
{
    for (auto index = std::size_t{0}; index != m_size; ++index)
        destruct(index);
}