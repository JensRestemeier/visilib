#pragma once
#ifndef _MATH_REAL_EXPR_H_
#define _MATH_REAL_EXPR_H_

#if 0
#include <Expr.h>
#include <ExprUtil.h>
#else
//#include <CORE/Expr.h>
#define CORE_LEVEL 4
#include <CORE.h>

typedef CORE::Expr Expr;
typedef CORE::BigFloat BigFloat;
#endif


struct RealExpr;

inline std::ostream& operator<< (std::ostream& stream, const RealExpr& val);
inline RealExpr operator+(const RealExpr& lhs, const RealExpr& rhs);
inline RealExpr operator-(const RealExpr& lhs, const RealExpr& rhs);
inline RealExpr operator/(const RealExpr& lhs, const RealExpr& rhs);
inline RealExpr operator*(const RealExpr& lhs, const RealExpr& rhs);

struct RealExpr {
    Expr expr;
    void init()
    {
    }
    inline void copy(const RealExpr& rhs)
    {
        expr= rhs.expr;
    }
    RealExpr() {
        init();
    }
    RealExpr(double d) {
        init();
        expr = BigFloat(d);
    }
    RealExpr(const RealExpr& rhs)
    {
        init();
        copy(rhs);
    }
    RealExpr(const Expr& _expr)
    {
        init();
        expr = _expr;
    }
    const RealExpr& operator=(const RealExpr& rhs)
    {
        if (this != &rhs)
        {
            copy(rhs);
        }
        return *this;
    }
    inline RealExpr operator-() const
    {
        return RealExpr(-expr);
    }
    inline const RealExpr& operator*=(const RealExpr& rhs)
    {
        Expr tmp(expr* rhs.expr);
        expr=(tmp);
        return *this;
    }
    inline const RealExpr& operator/=(const RealExpr& rhs)
    {
        Expr tmp(expr / rhs.expr);
        expr=(tmp);
        return *this;
    }
    inline const RealExpr& operator-=(const RealExpr& rhs)
    {
        Expr tmp(expr - rhs.expr);
        expr=(tmp);
        return *this;
    }
    inline const RealExpr& operator+=(const RealExpr& rhs)
    {
        Expr tmp(expr + rhs.expr);
        expr=(tmp);
        return *this;
    }
    inline bool operator<(const RealExpr& rhs) const
    {
        return expr < rhs.expr;
    }
    inline bool operator>(const RealExpr& rhs) const
    {
        return expr > rhs.expr;
    }
    inline bool operator<=(const RealExpr& rhs) const
    {
        return expr <= rhs.expr;
    }
    inline bool operator>=(const RealExpr& rhs) const
    {
        return expr >= rhs.expr;
    }

    inline bool isZero() const {
        return expr == Expr(0);
    }

    RealExpr abs() const 
    {
#if true
        //return RealExpr(expr * Real(expr < 0 ? -1 : 1));
        //return RealExpr((expr * expr).sqrt());
        if (expr < 0)
        {
            return RealExpr( - expr);
        }
        else
        {
            return *this;
        }
#else
        return RealExpr(expr.abs());
#endif
    }
    static RealExpr tolerance() {
        return RealExpr(1e-20);
    }
};

inline RealExpr sqrt(const RealExpr& x)
{
    RealExpr tmp(x);
    return RealExpr(sqrt(tmp.expr));
}
inline std::ostream& operator<< (std::ostream& stream, const RealExpr& val)
{
    return stream << val.expr;
}
inline RealExpr operator+(const RealExpr& lhs, const RealExpr& rhs)
{
    return RealExpr(lhs.expr + rhs.expr);
}
inline RealExpr operator-(const RealExpr& lhs, const RealExpr& rhs)
{
    return RealExpr(lhs.expr - rhs.expr);
}
inline RealExpr operator/(const RealExpr& lhs, const RealExpr& rhs)
{
    return RealExpr(lhs.expr / rhs.expr);
}
inline RealExpr operator*(const RealExpr& lhs, const RealExpr& rhs)
{
    return RealExpr(lhs.expr * rhs.expr);
}
inline double to_double(const RealExpr& rhs)
{
#if 0
    RealExpr tmp(rhs);
    return tmp.expr.approx().to_double();
#else
    RealExpr tmp(rhs);
    return tmp.expr.approx().get_d();
#endif
}

#endif // _MATH_REAL_EXPR_H_
