#include "math_util.h"
#include <math.h>
#include <chrono>
#include <iostream>
#include <numeric>

inlvector<double> allocate_linspace(double min, double max, const int n)
{
    auto t1 = std::chrono::high_resolution_clock::now();
    inlvector<double> ret(n);
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << ret[5] << std::endl;
    std::chrono::duration<double, std::micro> ms_double = t2 - t1;
    std::cout << "allocationg time = " << ms_double.count() << "micros\n";

    t1 = std::chrono::high_resolution_clock::now();
    const double dx = (max - min) / (double(n) - 1.0);
    t2 = std::chrono::high_resolution_clock::now();
    ms_double = t2 - t1;
    std::cout << "scalar compute time = " << ms_double.count() << "micros\n";

    t1 = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < n; ++i)
    {
        ret[i] = min + (i * dx);
    }
    t2 = std::chrono::high_resolution_clock::now();
    ms_double = t2 - t1;
    std::cout << "computing time = " << ms_double.count() << "micros\n";
    return ret;
}

std::vector<double> vectorized_sqrt(std::vector<double> input)
{
    std::vector<double> output(input.size());
    int i = 0;

    for (int i = 0; i < input.size(); ++i)
    {
        output[i] = sqrt(input[i]);
    }

    return output;
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
inlvector<Ty> &inlvector<Ty>::operator=(inlvector &&other)
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