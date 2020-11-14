#ifndef DATA_HPP
#define DATA_HPP

#include <vector>

class event_data
{
    public:
        // A T-indexed array of fluctuations
        std::vector<double> E;
        std::vector<double> E_fluc;
        // A T-indexed array of magnetizations
        std::vector<double> M;
        std::vector<double> M_fluc;

        std::vector<double> T;
};

#endif
