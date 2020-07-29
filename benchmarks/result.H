#pragma once

#include <vector>
#include <string>


struct Result{

    Result(const std::string& name_, const std::vector<double>& times)
    : name(name_), samples(times.size()), average(calc_avg(times)), standard_dev(calc_std(times))
    {}

    double calc_avg(const std::vector<double>& times) const{
        return std::accumulate(times.begin(), times.end(), 0.0) / times.size();
    }

    double calc_std(const std::vector<double>& times) const {
        auto mean = calc_avg(times);
        double sum = 0.0;
        for (const auto& time : times) {
            sum += (mean - time) * (mean - time);
        }
        return sum / times.size(); //!

    }


    std::string to_string() const {
        std::string ret;
        ret += name;
        ret += std::string(" Samples :") + std::to_string(samples);  
        ret += std::string(" Mean runtime: ") + std::to_string(average);
        ret += std::string(" Standard deviation: ") + std::to_string(standard_dev);
        return ret;       
    }

    std::string name;
    size_t samples;
    double average, standard_dev;

};