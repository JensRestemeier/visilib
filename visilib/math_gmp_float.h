#pragma once
#ifndef _MATH_GMP_FLOAT_H_
#define _MATH_GMP_FLOAT_H_

#include <gmp.h>
struct GmpFloat;

inline GmpFloat sqrt(const GmpFloat& val);
inline std::ostream& operator<< (std::ostream& stream, const GmpFloat& val);
inline GmpFloat operator+(const GmpFloat& lhs, const GmpFloat& rhs);
inline GmpFloat operator-(const GmpFloat& lhs, const GmpFloat& rhs);
inline GmpFloat operator/(const GmpFloat& lhs, const GmpFloat& rhs);
inline GmpFloat operator*(const GmpFloat& lhs, const GmpFloat& rhs);
inline double to_double(const GmpFloat& rhs);

struct GmpFloat {
    mpf_t v;
    void init()
    {
        mpf_init2(v, 256);
    }
    inline void copy(const GmpFloat& rhs)
    {
        mpf_set(v, rhs.v);
    }
    GmpFloat() {
        init();
        mpf_set_d(v, 0);
    }
    GmpFloat(double d) {
        init();
        mpf_set_d(v, d);
    }
    GmpFloat(const GmpFloat& rhs)
    {
        init();
        copy(rhs);
    }
    const GmpFloat& operator=(const GmpFloat& rhs)
    {
        if (this != &rhs)
        {
            copy(rhs);
        }
        return *this;
    }
    inline GmpFloat operator-() const
    {
        GmpFloat tmp;
        mpf_neg(tmp.v, v);
        return tmp;
    }
    inline const GmpFloat& operator*=(const GmpFloat& rhs)
    {
        *this = *this * rhs;
        return *this;
    }
    inline const GmpFloat& operator/=(const GmpFloat& rhs)
    {
        *this = *this / rhs;
        return *this;
    }
    inline const GmpFloat& operator-=(const GmpFloat& rhs)
    {
        *this = *this - rhs;
        return *this;
    }
    inline const GmpFloat& operator+=(const GmpFloat& rhs)
    {
        *this = *this + rhs;
        return *this;
    }
    inline bool operator<(const GmpFloat& rhs) const
    {
        return mpf_cmp(v, rhs.v) < 0;
    }
    inline bool operator>(const GmpFloat& rhs) const
    {
        return mpf_cmp(v, rhs.v) > 0;
    }
    inline bool operator<=(const GmpFloat& rhs) const
    {
        return mpf_cmp(v, rhs.v) <= 0;
    }
    inline bool operator>=(const GmpFloat& rhs) const
    {
        return mpf_cmp(v, rhs.v) >= 0;
    }
};

inline GmpFloat sqrt(const GmpFloat& val)
{
    GmpFloat tmp;
    mpf_sqrt(tmp.v, val.v);
    return tmp;
}
inline std::ostream& operator<< (std::ostream& stream, const GmpFloat& val)
{
    return stream << mpf_get_d(val.v);
}
inline GmpFloat operator+(const GmpFloat& lhs, const GmpFloat& rhs)
{
    GmpFloat tmp;
    mpf_add(tmp.v, lhs.v, rhs.v);
    return tmp;
}
inline GmpFloat operator-(const GmpFloat& lhs, const GmpFloat& rhs)
{
    GmpFloat tmp;
    mpf_sub(tmp.v, lhs.v, rhs.v);
    return tmp;
}
inline GmpFloat operator/(const GmpFloat& lhs, const GmpFloat& rhs)
{
    GmpFloat tmp;
    mpf_div(tmp.v, lhs.v, rhs.v);
    return tmp;
}
inline GmpFloat operator*(const GmpFloat& lhs, const GmpFloat& rhs)
{
    GmpFloat tmp;
    mpf_mul(tmp.v, lhs.v, rhs.v);
    return tmp;
}
inline double to_double(const GmpFloat& rhs)
{
    return mpf_get_d(rhs.v);
}

#endif _MATH_GMP_FLOAT_H_
