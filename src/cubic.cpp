#include "cubic.hpp"
#include "math.hpp"

#include <algorithm>
#include <cassert>

namespace ccd {

std::array<double, 2> CubicEquation::extrema() const
{
    return solve_quadratic_equation(3 * a, 2 * b, c);
}

std::vector<std::array<double, 2>> CubicEquation::monotonic_intervals() const
{
    constexpr double x_start = 0, x_end = 1;
    auto [x1, x2] = extrema();
    assert(x1 <= x2);
    x1 = std::clamp(x1, x_start, x_end);
    x2 = std::clamp(x2, x_start, x_end);

    const double sgn_f_x1 = sgn((*this)(x1));
    const double sgn_f_x2 = sgn((*this)(x2));

    std::vector<std::array<double, 2>> monotonic_intervals;
    if (x_start < x1 && sgn((*this)(x_start)) != sgn_f_x1) { // (x_start, x1]
        monotonic_intervals.push_back({ { x_start, x1 } });
    }
    if (x1 < x2 && sgn_f_x1 != sgn_f_x2) { // (x1, x2]
        monotonic_intervals.push_back({ { x1, x2 } });
    }
    if (x2 < x_end && sgn_f_x2 != sgn((*this)(x_end))) { // (x2, x_end]
        monotonic_intervals.push_back({ { x2, x_end } });
    }

    return monotonic_intervals;
}

bool CubicEquation::is_nearly_quadratic(const double tol) const
{
    return std::abs(a) < tol || std::abs(a / (b != 0 ? b : 1)) < tol;
}

bool CubicEquation::is_nearly_linear(const double tol) const
{
    return is_nearly_quadratic(tol)
        && (std::abs(b) < tol || std::abs(b / (c != 0 ? c : 1)) < tol);
}

bool CubicEquation::is_nearly_constant(const double tol) const
{
    return is_nearly_linear(tol)
        && (std::abs(c) < tol || std::abs(c / (d != 0 ? d : 1)) < tol);
}

std::ostream& operator<<(std::ostream& os, const CubicEquation& f)
{
    return (os << f.a << "x^3 + " << f.b << "x^2 + " << f.c << "x + " << f.d);
}

} // namespace ccd