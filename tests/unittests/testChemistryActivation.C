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

    CHECK(chemistryActive(start, duration, repeat, timeValueFromUserTime(2.0, repeat)));
    CHECK(chemistryActive(start, duration, repeat, timeValueFromUserTime(4.0, repeat)));
    CHECK(chemistryActive(start, duration, repeat, timeValueFromUserTime(5.0, repeat)));
    CHECK(!chemistryActive(start, duration, repeat, timeValueFromUserTime(6.0, repeat)));
    CHECK(chemistryActive(start, duration, repeat, timeValueFromUserTime(12.0, repeat)));
    CHECK(!chemistryActive(start, duration, repeat, timeValueFromUserTime(16.0, repeat)));
}

TEST_CASE("startTime/duration with repeat (wrapped interval)")
{
    const scalar start = 8.0;
    const scalar duration = 5.0;
    const scalar repeat = 10;

    CHECK(!chemistryActive(start, duration, repeat, timeValueFromUserTime(7.9, repeat)));
    CHECK(chemistryActive(start, duration, repeat, timeValueFromUserTime(8.0, repeat)));
    CHECK(chemistryActive(start, duration, repeat, timeValueFromUserTime(9.5, repeat)));
    CHECK(chemistryActive(start, duration, repeat, timeValueFromUserTime(10.0, repeat))); // 10 % 10 == 0
    CHECK(chemistryActive(start, duration, repeat, timeValueFromUserTime(0.0, repeat)));
    CHECK(chemistryActive(start, duration, repeat, timeValueFromUserTime(3.0, repeat)));
    CHECK(!chemistryActive(start, duration, repeat, timeValueFromUserTime(3.1, repeat)));
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

    // ensure folding still works for repeat when sentinel is not present
    const scalar rep = 10;
    // extreme start sentinel combined with repeat should still be treated as always active
    CHECK(chemistryActive(great, 3.0, rep, timeValueFromUserTime(5.0, rep)));
}