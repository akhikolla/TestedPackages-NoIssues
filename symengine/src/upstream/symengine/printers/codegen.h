#ifndef SYMENGINE_CODEGEN_H
#define SYMENGINE_CODEGEN_H

#include <symengine/visitor.h>
#include <symengine/printers/strprinter.h>
#include <symengine/symengine_exception.h>

namespace SymEngine
{

class CodePrinter : public BaseVisitor<CodePrinter, StrPrinter>
{
public:
    using StrPrinter::str_;
    using StrPrinter::apply;
    using StrPrinter::bvisit;
    void bvisit(const Basic &x);
    void bvisit(const Complex &x);
    void bvisit(const Interval &x);
    void bvisit(const Contains &x);
    void bvisit(const Piecewise &x);
    void bvisit(const Rational &x);
    void bvisit(const EmptySet &x);
    void bvisit(const FiniteSet &x);
    void bvisit(const UniversalSet &x);
    void bvisit(const Abs &x);
    void bvisit(const Ceiling &x);
    void bvisit(const Truncate &x);
    void bvisit(const Max &x);
    void bvisit(const Min &x);
    void bvisit(const Constant &x);
    void bvisit(const NaN &x);
    void bvisit(const Equality &x);
    void bvisit(const Unequality &x);
    void bvisit(const LessThan &x);
    void bvisit(const StrictLessThan &x);
    void bvisit(const UnivariateSeries &x);
    void bvisit(const Derivative &x);
    void bvisit(const Subs &x);
    void bvisit(const GaloisField &x);
};

class C89CodePrinter : public BaseVisitor<C89CodePrinter, CodePrinter>
{
public:
    using CodePrinter::str_;
    using CodePrinter::apply;
    using CodePrinter::bvisit;
    void bvisit(const Infty &x);
    void _print_pow(std::ostringstream &o, const RCP<const Basic> &a,
                    const RCP<const Basic> &b);
};

class C99CodePrinter : public BaseVisitor<C99CodePrinter, C89CodePrinter>
{
public:
    using C89CodePrinter::str_;
    using C89CodePrinter::apply;
    using C89CodePrinter::bvisit;
    void bvisit(const Infty &x);
    void _print_pow(std::ostringstream &o, const RCP<const Basic> &a,
                    const RCP<const Basic> &b);
    void bvisit(const Gamma &x);
    void bvisit(const LogGamma &x);
};

class JSCodePrinter : public BaseVisitor<JSCodePrinter, CodePrinter>
{
public:
    using CodePrinter::str_;
    using CodePrinter::apply;
    using CodePrinter::bvisit;
    void bvisit(const Constant &x);
    void _print_pow(std::ostringstream &o, const RCP<const Basic> &a,
                    const RCP<const Basic> &b);
    void bvisit(const Abs &x);
    void bvisit(const Sin &x);
    void bvisit(const Cos &x);
    void bvisit(const Max &x);
    void bvisit(const Min &x);
};
}

#endif // SYMENGINE_CODEGEN_H
