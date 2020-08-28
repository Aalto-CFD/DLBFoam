#pragma once

#include "result.H"

#include <vector>

struct Runner
{

    template <class BenchmarkType>
    static Result run(BenchmarkType benchmark, size_t n_times)
    {

        std::vector<double> times;
        clockTime           time;
        for(size_t i = 0; i < n_times; ++i)
        {
            time.timeIncrement();
            benchmark.run();
            times.push_back(time.timeIncrement());
            Pstream::waitRequests(); // not sure if this is a barrier but it
                                     // would be good to have one here...
        }
        // Remove the first logged instance to treat it as a dry run
        times.erase(times.begin());
        return Result(benchmark.get_info(), times);
    }
};
