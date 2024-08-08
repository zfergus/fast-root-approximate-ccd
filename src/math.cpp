#include "math.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

namespace ccd {

std::array<double, 2>
solve_quadratic_equation(const double a, const double b, const double c)
{
    const double delta = b * b - 4 * a * c;
    assert(delta >= -std::numeric_limits<double>::epsilon());

    double tmp = b;
    if (delta > std::numeric_limits<double>::epsilon()) {
        tmp += sgn(b) * std::sqrt(delta);
    }

    std::array<double, 2> roots = { { -2 * c / tmp, -tmp / (2 * a) } };
    if (roots[0] > roots[1])
        std::swap(roots[0], roots[1]);
    return roots;
}

double
newton_raphson(const CubicEquation& f, const double x0, const double tolerance)
{
    double prev_x, x = x0;
    do {
        prev_x = x;
        x -= std::clamp(f(x) / f.derivative(x), -1.0, 1.0);
    } while (std::abs(x - prev_x) > tolerance);
    return x;
}

double modified_newton_raphson(
    const CubicEquation& f,
    const double x0,
    const double locally_min_gradient,
    const double tolerance,
    const int max_iter)
{
    double prev_x, x = x0;
    int iter = 0;
    do {
        prev_x = x;
        x -= std::clamp(f(x) / locally_min_gradient, -1.0, 1.0);
    } while (std::abs(x - prev_x) > tolerance && ++iter < max_iter);
    return x;
}

} // namespace ccd
