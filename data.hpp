#ifndef DATA_HPP
#define DATA_HPP

#include <vector>

class event_data
{
    public:
        // A T-indexed array of fluctuations
        std::vector<double> avg_E;
        std::vector<double> E_fluc;
        // A T-indexed array of magnetizations
        std::vector<double> avg_M;
        std::vector<double> M_fluc;

        // Heat capacity
        std::vector<double> C_v;

        // Entropy
        std::vector<double> S;

        // magnetization
        std::vector<double> M_sweeps;
        std::vector<double> E_sweeps;

        std::vector<double> T;

        std::vector<std::vector<double>> cf;
        std::vector<unsigned> cf_len;
};

#endif
