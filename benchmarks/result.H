#pragma once

#include <vector>
#include <string>
#include <iomanip>
#include "benchmark_info.H"

struct Result{

    Result(const BenchmarkInfo& info_, const std::vector<double>& times)
    : info(info_), samples(times.size()), average(calc_avg(times)), standard_dev(calc_std(times))
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


    static std::string get_header() {

        std::ostringstream ss;

        ss << std::setw(20) << "Model";
        ss << std::setw(20) << "Function";
        ss << std::setw(20) << "Balancer";
        ss << std::setw(20) << "Init. Condition";
        ss << std::setw(10) << "N samples";
        ss << std::setw(10) << "Mean";
        ss << std::setw(10) << "Stdev";
        return ss.str();
    }

    static std::string get_header_csv() {
    

        std::ostringstream ss;
        ss << "Model," << "Function," << "Balancer," << "Init. Condition," << "N samples," << "Mean," << "Stdev"; 
        return ss.str();

    }

    std::string to_string() const{

        std::ostringstream ss;
        ss << std::setw(20) << info.chem_model_name;
        ss << std::setw(20) << info.benchmarked_function;
        ss << std::setw(20) << info.balancer_name;
        ss << std::setw(20) << info.initial_condition_name;
        ss << std::setw(10) << samples;
        ss << std::setw(10) << average;
        ss << std::setw(10) << standard_dev;

        return ss.str();


    }

    std::string to_csv() const{

        std::ostringstream ss;
        ss << info.chem_model_name; 
        ss << ",";
        ss << info.benchmarked_function;
        ss << ",";
        ss << info.balancer_name;
        ss << ",";
        ss << info.initial_condition_name;
        ss << ",";
        ss << samples;
        ss << ",";
        ss << average;
        ss << ",";
        ss << standard_dev;
        
        return ss.str();


    }


    BenchmarkInfo info;
    size_t samples;
    double average, standard_dev;

};
/*
std::string chem_model_name;
    std::string benchmarked_function;
    std::string balancer_name;
    std::string initial_condition_name;
*/
Foam::Ostream& operator<<(Foam::Ostream& os, const Result& r) {

    os << r.to_string();
    return os;

}