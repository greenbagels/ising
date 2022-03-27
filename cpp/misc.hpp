#ifndef MISC_HPP
#define MISC_HPP

#include <vector>

namespace ising
{
    enum MC_BACKENDS
    {
        BACKEND_SINGLE_FLIP = 0,
        BACKEND_WOLFF = 1,
        NUM_BACKENDS = 2
    };

    struct event_data
    {
            // A T-indexed array of fluctuations
            std::vector<double> avg_E;
            std::vector<double> E_fluc;
            // A T-indexed array of magnetizations
            std::vector<double> avg_M;
            std::vector<double> M_fluc;

            // magnetization
            std::vector<double> M_sweeps;
            std::vector<double> E_sweeps;

            std::vector<double> T;

            std::vector<std::vector<double>> cf;
            std::vector<unsigned> cf_len;
    };
}

#endif
