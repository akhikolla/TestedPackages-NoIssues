# ifndef HESIM_MATH_COMPOSITE_H
# define HESIM_MATH_COMPOSITE_H
#include <RcppArmadillo.h>

namespace hesim {

/** @ingroup math 
 * Functions and classes for numerical computing including integration 
 * and root finding.
 */
namespace math {

namespace detail {

/***************************************************************************//** 
 * Composite rule given a fixed sample.
 * Integrate a function with values y at points x using a composite rule. 
 * @param x_first, x_last Iterators defining the values of the points x.  
 * @param y_first Initial value of y.
 * @param rule A function denoting the integration rule.
 * @return The integral of the function approximated by the integration rule.
 ******************************************************************************/ 
template <typename Func, typename InputIt>
double composite(InputIt x_first, InputIt x_last, InputIt y_first, Func rule){
  double area = 0;
  InputIt y_it = y_first;
  for (InputIt x_it = x_first; x_it != (x_last - 1); ++x_it){
    area += rule(*x_it, *(x_it + 1), *y_it, *(y_it + 1));
    ++y_it;
  }
  return area;
}

} // end namespace detail

/***************************************************************************//** 
 * Composite trapezoid rule given a fixed sample.
 * Integrate a function with values y at points x using a composite trapezoid rule. 
 * @param x_first, x_last Iterators defining the values of the points x.  
 * @param y_first Initial value of y.
 * @return The integral of the function.
 ******************************************************************************/ 
template <typename InputIt>
inline double trapz(InputIt x_first, InputIt x_last, InputIt y_first){
  auto fun = [](double a, double b, double y_a, double y_b){
    return (b - a)/2 * (y_a + y_b);
  };
  return detail::composite(x_first, x_last, y_first, fun);
}

/***************************************************************************//** 
 * Left Riemann sum given a fixed sample.
 * Integrate a function with values y at points x by using approximating a 
 * function with points at the left point of intervals.
 * @param x_first, x_last Iterators defining the values of the points x.  
 * @param y_first Initial value of y.
 * @return The integral of the function.
 ******************************************************************************/ 
template <typename InputIt>
inline double riemann_left(InputIt x_first, InputIt x_last, InputIt y_first){
  auto fun = [](double a, double b, double y_a, double y_b){
    return (b - a) * y_a;
  };
  return detail::composite(x_first, x_last, y_first, fun);
}

/***************************************************************************//** 
 * Right Riemann sum given a fixed sample.
 * Integrate a function with values y at points x by using approximating a
 * function with points at the right point of intervals.
 * @param x_first, x_last Iterators defining the values of the points x.  
 * @param y_first Initial value of y.
 * @return The integral of the function.
 ******************************************************************************/ 
template <typename InputIt>
inline double riemann_right(InputIt x_first, InputIt x_last, InputIt y_first){
  auto fun = [](double a, double b, double y_a, double y_b){
    return (b - a) * y_b;
  };
  return detail::composite(x_first, x_last, y_first, fun);
}

} // end namespace math

} // end namespace hesim

# endif
