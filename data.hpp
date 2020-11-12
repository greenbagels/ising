#ifndef DATA_HPP
#define DATA_HPP

#include <vector>

class event_data
{
    public:
        // Energy
        std::vector<double> E;
        double equilibrium_E;
        // Magnetization
        std::vector<double> m;
        double equilibrium_m;
};

#endif
