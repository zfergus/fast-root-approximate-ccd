#pragma once

#include "cubic.hpp"

namespace ccd {

/// @brief Compute the sign of a number.
/// @tparam T Type of the number.
/// @param x Number to compute the sign of.
/// @return -1 if the number is negative, 1 if the number is positive, 0 otherwise.
template <typename T> inline int sgn(T x) { return (T(0) < x) - (x < T(0)); }

/// @brief Compute the roots of a quadratic equation.
///
/// The roots are sorted in ascending order. Undefined behavior if the equation
/// has no real roots. If the equation has a double root, the root is repeated.
///
/// @param a Coefficient of the quadratic term.
/// @param b Coefficient of the linear term.
/// @param c Constant term.
/// @return Roots of the quadratic equation.
std::array<double, 2>
solve_quadratic_equation(const double a, const double b, const double c);

/// @brief Perform a modified Newton-Raphson root finding algorithm to find the roots of a cubic equation.
/// @param f Cubic equation to find the roots of.
/// @param x0 Initial guess for the root.
/// @param tolerance Tolerance for the root finding algorithm.
/// @return Root of the cubic equation.
double newton_raphson(
    const CubicEquation& f, const double x0, const double tolerance = 1e-6);

/// @brief Perform a Newton-Raphson root finding algorithm to find the roots of a cubic equation.
/// @param f Cubic equation to find the roots of.
/// @param x0 Initial guess for the root.
/// @param tolerance Tolerance for the root finding algorithm.
/// @return Root of the cubic equation.
double modified_newton_raphson(
    const CubicEquation& f,
    const double x0,
    const double locally_min_gradient,
    const double tolerance = 1e-6,
    const int max_iter = 100);

} // namespace ccd
