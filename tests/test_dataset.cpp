#include <catch2/catch_all.hpp>

#include <ccd.hpp>

#include <ccd_io/read_ccd_queries.hpp>

#include <filesystem>

enum Case { POINT_TRIANGLE, EDGE_EDGE };

TEST_CASE("Dataset", "[ccd][point-triangle][edge-edge][dataset][!mayfail]")
{
    using namespace ccd;
    namespace fs = std::filesystem;

    const std::string folder_name = GENERATE(
        "erleben-sliding-spike", "erleben-spike-wedge", "erleben-sliding-wedge",
        "erleben-wedge-crack", "erleben-spike-crack", "erleben-wedges",
        "erleben-cube-cliff-edges", "erleben-spike-hole",
        "erleben-cube-internal-edges", "erleben-spikes", "unit-tests", "chain",
        "cow-heads", "golf-ball", "mat-twist");
    const Case test_case = GENERATE(POINT_TRIANGLE, EDGE_EDGE);

    const fs::path path = fs::path(CCD_IO_SAMPLE_QUERIES_DIR) / folder_name
        / (test_case == POINT_TRIANGLE ? "vertex-face" : "edge-edge");
    REQUIRE(fs::exists(path));
    REQUIRE(fs::is_directory(path));

    for (const auto& f : fs::directory_iterator(path)) {
        if (f.path().extension() != ".csv") {
            continue;
        }

        CAPTURE(f.path().string());
        const std::vector<ccd_io::CCDQuery> queries =
            ccd_io::read_ccd_queries(f.path().string());

        for (const auto& [vertices, is_collision_expected] : queries) {
            Eigen::Map<const Eigen::Matrix<double, 8, 3>> V(&vertices[0][0]);

            bool result;
            double toi;
            switch (test_case) {
            case POINT_TRIANGLE:
                result = point_triangle_ccd(
                    V.row(0), V.row(1), V.row(2), V.row(3), V.row(4), V.row(5),
                    V.row(6), V.row(7), toi);
                break;
            case EDGE_EDGE:
                result = edge_edge_ccd(
                    V.row(0), V.row(1), V.row(2), V.row(3), V.row(4), V.row(5),
                    V.row(6), V.row(7), toi);
            }

            CHECK((result || !is_collision_expected)); // conservative check
        }
    }
}