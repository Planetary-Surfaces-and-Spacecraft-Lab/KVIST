#include <vector>
#include <exception>
// #include "ian_numerics_library.h"

// Ian Numerics Library
// Based on https://codereview.stackexchange.com/questions/211241/simple-re-implementation-of-stdvector
// The intention behind this vector is that it allows a scalar function to
// initialize a vector
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

    void swap(memory_block &other)
    {
        using std::swap;
        swap(m_size, other.m_size);
        swap(m_buffer, other.m_buffer);
    }

    memory_block &operator=(memory_block const &other)
    {
        auto temp = memory_block(other);
        swap(temp);
        return *this;
    }

    memory_block &operator=(memory_block &&other)
    {
        auto temp = memory_block(std::move(other));
        swap(temp);
        return *this;
    }

    std::size_t size() const
    {
        return m_size;
    }

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
    inlvector();
    explicit inlvector(std::size_t count);
    inlvector(std::size_t count, const Ty &val);

    inlvector(const inlvector &other);
    inlvector(inlvector &&other);

    ~inlvector();

    void swap(inlvector &other);

    inlvector &operator=(const inlvector &other);
    inlvector &operator=(inlvector &&other);

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


// Implemented in math_util.cpp
inlvector<double> allocate_linspace(double min, double max, const int n);

// Implemented in math_util.cpp
std::vector<double> vectorized_sqrt(std::vector<double> input);
