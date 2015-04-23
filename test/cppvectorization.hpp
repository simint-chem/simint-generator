#ifndef CPPVECTORIZATION_H
#define CPPVECTORIZATION_H

#include "vectorization.h"

template <typename T>
struct AlignedAllocator
{
    typedef T value_type;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;

    typedef T * pointer;
    typedef const T * const_pointer;

    typedef T & reference;
    typedef const T & const_reference;

    AlignedAllocator() throw () { }

    template <typename T2>
    AlignedAllocator(const AlignedAllocator<T2> &) throw () { }

    ~AlignedAllocator() throw () { }

    pointer address(reference r)
    {
        return &r;
    }

    const_pointer address(const_reference r) const
    { 
        return &r;
    }

    pointer allocate(size_type n)
    {
        return (pointer)ALLOC(n*sizeof(value_type));
    }

    void deallocate(pointer p, size_type)
    {
        FREE(p);
    }

    void construct(pointer p, const value_type & v)
    {
        new (p) value_type(v);
    }

    void destroy(pointer p)
    {
        p->~value_type();
    }

    size_type max_size() const throw ()
    {
        return size_type(-1) / sizeof(value_type);
    }

    bool operator!=(const AlignedAllocator<T> & other) const
    {
        return false;
    }

    bool operator==(const AlignedAllocator<T> & other) const 
    {
        return true;
    }

    template <typename T2>
    struct rebind
    {
        typedef AlignedAllocator<T2> other;
    };
};

#endif
