#include <catch2/catch_all.hpp>

#include <ccd.hpp>
#include <autogen.hpp>

using namespace ccd;

static const double EPSILON = std::numeric_limits<float>::epsilon();

TEST_CASE("Edge-Edge CCD", "[ccd][3D][edge-edge][!mayfail]")
{
    Eigen::Vector3d ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1,
        eb1_t1;
    bool is_collision_expected;
    bool conservative_check = true;

    std::string name;
    SECTION("General")
    {
        double uy = GENERATE(-1.0, 0.0, 1 - EPSILON, 1.0, 1 + EPSILON, 2.0);
        double e1x = GENERATE(
            -1 - EPSILON, -1, -1 + EPSILON, -0.5, 0, 0.5, 1 - EPSILON, 1,
            1 + EPSILON);

        ea0_t0 << -1, -1, 0;
        ea1_t0 << 1, -1, 0;
        eb0_t0 << e1x, 1, -1;
        eb1_t0 << e1x, 1, 1;

        Eigen::Vector3d u0, u1;
        SECTION("moving")
        {
            name = "General moving case";
            u0 << 0, uy, 0;
            u1 << 0, -uy, 0;
            is_collision_expected = uy >= 1.0 && e1x >= -1 && e1x <= 1;
        }
        SECTION("fixed")
        {
            name = "General fixed case";
            u0 << 0, 2 * uy, 0;
            u1.setZero();
            is_collision_expected = uy >= 2.0 && e1x >= -1 && e1x <= 1;
        }

        ea0_t1 = ea0_t0 + u0;
        ea1_t1 = ea1_t0 + u0;
        eb0_t1 = eb0_t0 + u1;
        eb1_t1 = eb1_t0 + u1;
    }
    SECTION("Double root test case 1")
    {
        name = "Double root test case 1";
        ea0_t0 << -3.0022200, 0.2362580, 0.0165247;
        ea1_t0 << -3.2347850, 0.8312380, -0.1151003;
        eb0_t0 << -3.0319900, 0.3148750, 0.0000000;
        eb1_t0 << -2.8548800, 0.0900349, 0.0000000;
        ea0_t1 << -2.8995600, 0.0345838, 0.0638580;
        ea1_t1 << -3.1716930, 0.6104858, -0.0713340;
        eb0_t1 = eb0_t0;
        eb1_t1 = eb1_t0;

        is_collision_expected = true;
    }
    SECTION("Double root test case 2")
    {
        name = "Double root test case 2";
        ea0_t0 << 0, 0, 1;
        ea1_t0 << 0, 1, 1;
        eb0_t0 << 0.1, 0.2, 2;
        eb1_t0 << 0.1, 0.2, -1;
        ea0_t1 << 1, 1, 0;
        ea1_t1 << 0, 0, 0;
        eb0_t1 = eb0_t0;
        eb1_t1 = eb1_t0;

        const double t = GENERATE(0.5, 0.8, 0.88, 0.9, 1.0);
        ea0_t1 = (ea0_t1 - ea0_t0) * t + ea0_t0;
        ea1_t1 = (ea1_t1 - ea1_t0) * t + ea1_t0;

        bool is_collision_expected = true;
    }
    SECTION("Slow Case 1")
    {
        name = "Slow Case 1";
        ea0_t0 << 1, 0.50803125, 2.10835646075301e-18;
        ea1_t0 << -2.38233935445388e-18, 0.50803125, 1;
        eb0_t0 << -4.99999999958867e-07, 0.5, 0;
        eb1_t0 << -4.99999999958867e-07, 0.5, 1;
        ea0_t1 << 1, 0.47124375, 4.11078309465837e-18;
        ea1_t1 << -2.8526707189104e-18, 0.47124375, 1;
        eb0_t1 << -4.99999999958867e-07, 0.5, 0;
        eb1_t1 << -4.99999999958867e-07, 0.5, 1;

        is_collision_expected = true;
        // conservative_check = false;
    }
    SECTION("Slow Case 2")
    {
        name = "Slow Case 2";
        ea0_t0 << 1.00002232466453, 0.500004786049044, -2.06727783590977e-05;
        ea1_t0 << 1.64687846177844e-05, 0.499996645067319, 1.63939999009028e-05;
        eb0_t0 << 1, 0.5, 0;
        eb1_t0 << 0, 0.5, 0;
        ea0_t1 << 1.00294282700155, 0.498652627047143, 0.003626320742036;
        ea1_t1 << -0.00219276550735626, 0.500871179186644, -0.00315828804921928;
        eb0_t1 << 1, 0.5, 0;
        eb1_t1 << 0, 0.5, 0;

        is_collision_expected = true;
        // conservative_check = false;
    }
    INFO(name);

    CAPTURE(autogen::edge_edge_ccd_equation(
        ea0_t0.x(), ea0_t0.y(), ea0_t0.z(), ea1_t0.x(), ea1_t0.y(), ea1_t0.z(),
        eb0_t0.x(), eb0_t0.y(), eb0_t0.z(), eb1_t0.x(), eb1_t0.y(), eb1_t0.z(),
        ea0_t1.x(), ea0_t1.y(), ea0_t1.z(), ea1_t1.x(), ea1_t1.y(), ea1_t1.z(),
        eb0_t1.x(), eb0_t1.y(), eb0_t1.z(), eb1_t1.x(), eb1_t1.y(),
        eb1_t1.z()));

    double toi;
    bool is_colliding = edge_edge_ccd(
        ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1, toi);
    if (conservative_check) {
        CHECK((is_colliding || !is_collision_expected));
    } else {
        CHECK(is_colliding == is_collision_expected);
    }
}
