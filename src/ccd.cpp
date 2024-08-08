#include "ccd.hpp"
#include "autogen.hpp"
#include "geometry.hpp"
#include "math.hpp"

namespace ccd {

bool point_triangle_ccd(
    const Eigen::Ref<const Eigen::Vector3d>& p_t0,
    const Eigen::Ref<const Eigen::Vector3d>& t0_t0,
    const Eigen::Ref<const Eigen::Vector3d>& t1_t0,
    const Eigen::Ref<const Eigen::Vector3d>& t2_t0,
    const Eigen::Ref<const Eigen::Vector3d>& p_t1,
    const Eigen::Ref<const Eigen::Vector3d>& t0_t1,
    const Eigen::Ref<const Eigen::Vector3d>& t1_t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2_t1,
    double& toi)
{
    return fast_approximate_root_ccd(
        autogen::point_triangle_ccd_equation(
            p_t0.x(), p_t0.y(), p_t0.z(), t0_t0.x(), t0_t0.y(), t0_t0.z(),
            t1_t0.x(), t1_t0.y(), t1_t0.z(), t2_t0.x(), t2_t0.y(), t2_t0.z(),
            p_t1.x(), p_t1.y(), p_t1.z(), t0_t1.x(), t0_t1.y(), t0_t1.z(),
            t1_t1.x(), t1_t1.y(), t1_t1.z(), t2_t1.x(), t2_t1.y(), t2_t1.z()),
        [&](const double toi) {
            return is_point_inside_triangle(
                (p_t1 - p_t0) * toi + p_t0, (t0_t1 - t0_t0) * toi + t0_t0,
                (t1_t1 - t1_t0) * toi + t1_t0, (t2_t1 - t2_t0) * toi + t2_t0);
        },
        toi);
}

bool edge_edge_ccd(
    const Eigen::Ref<const Eigen::Vector3d>& ea0_t0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1_t0,
    const Eigen::Ref<const Eigen::Vector3d>& eb0_t0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1_t0,
    const Eigen::Ref<const Eigen::Vector3d>& ea0_t1,
    const Eigen::Ref<const Eigen::Vector3d>& ea1_t1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0_t1,
    const Eigen::Ref<const Eigen::Vector3d>& eb1_t1,
    double& toi)
{
    return fast_approximate_root_ccd(
        autogen::edge_edge_ccd_equation(
            ea0_t0.x(), ea0_t0.y(), ea0_t0.z(), ea1_t0.x(), ea1_t0.y(),
            ea1_t0.z(), eb0_t0.x(), eb0_t0.y(), eb0_t0.z(), eb1_t0.x(),
            eb1_t0.y(), eb1_t0.z(), ea0_t1.x(), ea0_t1.y(), ea0_t1.z(),
            ea1_t1.x(), ea1_t1.y(), ea1_t1.z(), eb0_t1.x(), eb0_t1.y(),
            eb0_t1.z(), eb1_t1.x(), eb1_t1.y(), eb1_t1.z()),
        [&](const double toi) {
            return are_edges_intersecting(
                (ea0_t1 - ea0_t0) * toi + ea0_t0,
                (ea1_t1 - ea1_t0) * toi + ea1_t0,
                (eb0_t1 - eb0_t0) * toi + eb0_t0,
                (eb1_t1 - eb1_t0) * toi + eb1_t0);
        },
        toi);
}

bool fast_approximate_root_ccd(
    CubicEquation d,
    const std::function<bool(const double)>& is_inside,
    double& toi)
{
    const double d0 = d(0);
    if (d0 == 0 && is_inside(0)) {
        toi = 0;
        return true;
    } else if (d0 < 0) {
        d *= -1;
    }

    if (d.is_nearly_constant()) {
        toi = 0;
        return d.d == 0;
    } else if (d.is_nearly_linear()) {
        toi = -d.d / d.c;
        return toi >= 0 && toi <= 1 && is_inside(toi);
    } else if (d.is_nearly_quadratic()) {
        if (d.c * d.c - 4 * d.b * d.d < 0) // no real roots
            return false;
        for (const double root : solve_quadratic_equation(d.b, d.c, d.d)) {
            if (0 <= root && root <= 1 && is_inside(root)) {
                toi = root;
                return true;
            }
        }
        return false; // no roots in [0, 1]
    } else if (4 * d.b * d.b - 12 * d.a * d.c < 0) {
        return false;
    }

    // for (const auto& [a, b, min_gradient] : d.monotonic_intervals()) {
    //     toi = newton_raphson(d, a);
    //     if (0 <= toi && toi <= 1 && is_inside(toi)) {
    //         return true;
    //     }
    // }
    // return false;

    const auto root_interval = determine_cubic_root_interval(d);
    if (!root_interval || root_interval->a > root_interval->b) {
        return false;
    }

    toi = modified_newton_raphson(
        d, root_interval->a, root_interval->min_gradient);

    return toi >= 0 && toi <= 1 && is_inside(toi);
}

std::optional<RootInterval>
determine_cubic_root_interval(const CubicEquation& d)
{
    constexpr double t0 = 0, t1 = 1;
    assert(d(t0) >= 0);

    const auto [tm0, tm1] = d.extrema();
    assert(tm0 <= tm1);

    assert(d.a != 0);
    if (d.a > 0) {
        // Case 1 (Interval II)
        const double t_min = tm1, t_max = tm0;
        if (d(t_min) > 0 || t0 >= t_min) {
            // assert(d(t0) > 0);
            return std::nullopt;
        }
        return RootInterval(t_max, t_min, d.derivative((t_min + t_max) / 2));
    }

    const double t_min = tm0, t_max = tm1;
    if (d(t_min) >= 0) {
        // Case 2 (Interval III or may be an invalid interval if t_max > t1)
        return RootInterval(t_max, t1, d.derivative(t1));
    } else {
        // Case 3
        return (t0 < t_min)
            ? RootInterval(t0, t_min, d.derivative(t0))  // Interval I
            : RootInterval(t_max, t1, d.derivative(t1)); // Interval III
    }
}

} // namespace ccd