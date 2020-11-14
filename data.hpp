#ifndef DATA_HPP
#define DATA_HPP

#include <vector>

class event_data
{
    public:
        // Energy
        std::vector<double> E;
        std::vector<double> equilibrium_E;
        // Magnetization
        std::vector<double> m;
        std::vector<double> equilibrium_m;

        std::vector<double> T;
};

#endif
