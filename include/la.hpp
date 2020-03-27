#pragma once

#include <iostream>
#include <compare>
#include <boost/multiprecision/cpp_int.hpp>

namespace la {
    template<typename integer_type>
    class rational_t
    {
    public:
        using int_t = integer_type;

    private:
        using self_t = rational_t<integer_type>;

        int_t _num;
        int_t _den;

    public:

        rational_t(const int_t& num, const int_t& den);
        rational_t(int_t&& num, int_t&& den);
        rational_t(const int_t& value);
        rational_t(int_t&& value);

        int_t& numerator();
        const int_t& numerator() const;

        int_t& denominator();
        const int_t& denominator() const;

        int sign() const;

        void simplify();

        self_t& operator =(const self_t& other);
        self_t& operator =(self_t&& other);

        bool operator ==(const self_t& other);
        bool operator !=(const self_t& other);
        std::strong_ordering operator <=>(const self_t& other);

        self_t& operator ++();
        self_t& operator --();
        self_t operator ++(int _);
        self_t operator --(int _);

        self_t& operator +=(const self_t& other);
        self_t& operator -=(const self_t& other);
        self_t& operator *=(const self_t& other);
        self_t& operator /=(const self_t& other);

        friend int_t floor(const self_t& value);
        friend int_t ceil(const self_t& value);
    };
}
