#include <catch2/catch_all.hpp>

#include <ccd.hpp>

using namespace ccd;

static const double EPSILON = std::numeric_limits<float>::epsilon();

TEST_CASE("Point-Triangle CCD", "[ccd][point-triangle]")
{
    // point
    double v0z = GENERATE(0.0, -1.0);
    Eigen::Vector3d v0(0, 1, v0z);
    // triangle = (v1, v2, v3)
    Eigen::Vector3d v1(-1, 0, 1);
    Eigen::Vector3d v2(1, 0, 1);
    Eigen::Vector3d v3(0, 0, -1);

    // displacements
    double u0y =
        -GENERATE(-1.0, 0.0, 0.5 - EPSILON, 0.5, 0.5 + EPSILON, 1.0, 2.0);
    double u0z = GENERATE(-EPSILON, 0.0, EPSILON);
    Eigen::Vector3d u0(0, u0y, u0z);
    double u1y =
        GENERATE(-1.0, 0.0, 0.5 - EPSILON, 0.5, 0.5 + EPSILON, 1.0, 2.0);
    Eigen::Vector3d u1(0, u1y, 0);

    bool is_collision_expected = ((-u0y + u1y >= 1) && (v0z + u0z >= v3.z()));

    double toi;
    bool is_colliding = point_triangle_ccd(
        v0, v1, v2, v3, v0 + u0, v1 + u1, v2 + u1, v3 + u1, toi);

    CAPTURE(v0z, u0y, u0z, u1y, EPSILON, toi);
    REQUIRE(is_colliding == is_collision_expected);
}

TEST_CASE("Zhongshi test case", "[ccd][point-triangle]")
{
    double qy = GENERATE(-EPSILON, 0, EPSILON);

    Eigen::Vector3d q;
    q << 0, qy, 0;

    Eigen::Vector3d b0;
    b0 << 0, 0, 0;
    Eigen::Vector3d b1;
    b1 << 0, 1, 0;
    Eigen::Vector3d b2;
    b2 << 1, 0, 0;

    Eigen::Vector3d t0;
    t0 << 0, 0, 1;
    Eigen::Vector3d t1;
    t1 << 0, 1, 1;
    Eigen::Vector3d t2;
    t2 << 1, 0, 1;

    Eigen::Vector3d q1;
    q1 << 0, qy, 0;

    bool is_collision_expected = q.y() >= 0;

    double toi;
    bool is_colliding = point_triangle_ccd(q, b0, b1, b2, q1, t0, t1, t2, toi);

    CAPTURE(qy);
    CHECK(is_colliding == is_collision_expected);
}

TEST_CASE("Bolun test case", "[ccd][point-triangle]")
{
    Eigen::Vector3d x0(0.1, 0.1, 0.1), x1(0, 0, 1), x2(1, 0, 1), x3(0, 1, 1),
        x0b(0.1, 0.1, 0.1), x1b(0, 0, 0), x2b(0, 1, 0), x3b(1, 0, 0);

    bool is_collision_expected = true;

    double toi;
    bool is_colliding =
        point_triangle_ccd(x0, x1, x2, x3, x0b, x1b, x2b, x3b, toi);

    CHECK(is_colliding == is_collision_expected);
}

TEST_CASE("No Zero ToI CCD", "[ccd][point-triangle]")
{
    Eigen::Vector3d p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1;
    p_t0 << 0.0133653, 0.100651, -0.0215935;
    t0_t0 << 0.0100485, 0.0950896, -0.0171013;
    t1_t0 << 0.0130388, 0.100666, -0.0218112;
    t2_t0 << 0.015413, 0.100554, -0.0202265;
    p_t1 << 0.0133652999767858, 0.099670000268615, -0.0215934999996444;
    t0_t1 << 0.0100484999799995, 0.0941086002577558, -0.0171012999972189;
    t1_t1 << 0.0130387999724314, 0.0996850002629403, -0.0218111999936902;
    t2_t1 << 0.0154129999740718, 0.0995730002646605, -0.020226499996014;

    double toi;
    bool is_impacting = point_triangle_ccd(
        p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1, toi);

    CAPTURE(toi);
    CHECK(!is_impacting);
}

TEST_CASE(
    "Point-triangle false-negative", "[ccd][point-triangle][false-negative]")
{
    Eigen::Vector3d p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1;

    SECTION("Case 1")
    {
        // clang-format off
        p_t0  << -4375459 / 1048576.0, 13354865 / 33554432.0, -11393047 /  8388608.0;
        t0_t0 << -8363723 / 2097152.0, 16021805 / 33554432.0, -15346925 / 16777216.0;
        t1_t0 << -4334121 / 1048576.0,  3618099 /  8388608.0, -12072633 /  8388608.0;
        t2_t0 << -4482915 / 1048576.0,  3521177 /  8388608.0, -10864455 /  8388608.0;
        p_t1  << -4219061 / 1048576.0, 10201537 / 16777216.0,  -3974493 /  4194304.0;
        t0_t1 << -1051037 /  262144.0,  9399217 / 16777216.0,  -6996963 /  8388608.0;
        t1_t1 << -8539773 / 2097152.0,  4662053 /  8388608.0, -10288603 /  8388608.0;
        t2_t1 << -8621451 / 2097152.0,  9321975 / 16777216.0,   -151795 /   131072.0;
        // clang-format on
    }
    SECTION("Case 2")
    {
        // clang-format off
        p_t0  <<  -5066329 / 1048576.0,  3747703 /  8388608.0,  -9110393 / 4194304.0;
        t0_t0 <<  -5097349 / 1048576.0, 15565867 / 33554432.0,  -8419809 / 4194304.0;
        t1_t0 << -10108719 / 2097152.0,  4441449 /  8388608.0,  -2398085 / 1048576.0;
        t2_t0 << -10030265 / 2097152.0,  5321347 /  8388608.0,  -4791705 / 2097152.0;
        p_t1  <<  -9984255 / 2097152.0,  8516467 / 16777216.0,  -7675421 / 4194304.0;
        t0_t1 <<  -5031053 / 1048576.0,  7841285 / 16777216.0, -14946117 / 8388608.0;
        t1_t1 <<  -4854525 / 1048576.0,  4763077 /  8388608.0,   -269775 /  131072.0;
        t2_t1 <<  -9736347 / 2097152.0, 10062471 / 16777216.0,  -4408415 / 2097152.0;
        // clang-format on
    }

    double toi;
    bool is_impacting = point_triangle_ccd(
        p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1, toi);

    CAPTURE(toi);
    CHECK(is_impacting);
}