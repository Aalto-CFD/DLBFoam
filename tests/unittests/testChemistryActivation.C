#include "catch.hpp"
#include "ChemistryActivation.H"

using namespace Foam;

TEST_CASE("startTime/duration without repeat")
{
    const scalar start = 10.0;
    const scalar duration = 5.0;
    const scalar repeat = 0;

    CHECK(!chemistryActive(start, duration, repeat, 0.0));
    CHECK(!chemistryActive(start, duration, repeat, 9.999));
    CHECK(chemistryActive(start, duration, repeat, 10.0));
    CHECK(chemistryActive(start, duration, repeat, 12.5));
    CHECK(chemistryActive(start, duration, repeat, 15.0));
    CHECK(!chemistryActive(start, duration, repeat, 15.0001));
}

TEST_CASE("startTime/duration with repeat (non-wrapping)")
{
    const scalar start = 2.0;
    const scalar duration = 3.0;
    const scalar repeat = 10;

    CHECK(chemistryActive(start, duration, repeat, 2.0));
    CHECK(chemistryActive(start, duration, repeat, 4.0));
    CHECK(chemistryActive(start, duration, repeat, 5.0));
    CHECK(!chemistryActive(start, duration, repeat, 6.0));
    CHECK(chemistryActive(start, duration, repeat, 12.0));
    CHECK(!chemistryActive(start, duration, repeat, 16.0));
}

TEST_CASE("startTime/duration with repeat (wrapped interval)")
{
    const scalar start = 8.0;
    const scalar duration = 5.0;
    const scalar repeat = 10;

    CHECK(!chemistryActive(start, duration, repeat, 7.9));
    CHECK(chemistryActive(start, duration, repeat, 8.0));
    CHECK(chemistryActive(start, duration, repeat, 9.5));
    CHECK(chemistryActive(start, duration, repeat, 10.0)); // 10 % 10 == 0
    CHECK(!chemistryActive(start, duration, repeat, 0.0));
    CHECK(!chemistryActive(start, duration, repeat, 3.0));
}

TEST_CASE("start/duration sentinel (undefined via -great, great) are treated as always active")
{
    const scalar repeat = 0;

    // start undefined (sentinel -great)
    CHECK(chemistryActive(-great, 1.0, repeat, 1.234));

    // duration undefined (sentinel great)
    CHECK(chemistryActive(1.0, great, repeat, 1.234));

    // both undefined
    CHECK(chemistryActive(-great, great, repeat, 1000.0));

}

TEST_CASE("startTime/duration with repeat (startTime effect)")
{
    const scalar start = 700;
    const scalar duration = 50;
    const scalar repeat = 720;

    CHECK(!chemistryActive(start, duration, repeat, 0.0));
    CHECK(!chemistryActive(start, duration, repeat, 50.0));
    CHECK(!chemistryActive(start, duration, repeat, 699.999));
    CHECK(chemistryActive(start, duration, repeat, 700.0));
    CHECK(chemistryActive(start, duration, repeat, 749.999));
    CHECK(!chemistryActive(start, duration, repeat, 750.001));
    CHECK(!chemistryActive(start, duration, repeat, 1419.999));
    CHECK(chemistryActive(start, duration, repeat, 1420.0));
    CHECK(chemistryActive(start, duration, repeat, 1420.001));
    CHECK(chemistryActive(start, duration, repeat, 1469.999));
    CHECK(!chemistryActive(start, duration, repeat, 1470.001));


}