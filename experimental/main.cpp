
#include <cstdlib>
#include <iostream>
#include <string>
#include "wolff.hpp"

int main(int argc, char *argv[])
{

    int sweeps = atoi(argv[1]);
    int width = 400;
    int scale = 4;
    double temp = 4.;
    std::cerr << "Initializing sim with parameters:\n" <<
        "    sweeps = " << sweeps << "\n    scale = " << scale
        << "\n     temp = " << temp << "\n";

    wolff sim(sweeps, width, scale, temp); 
    sim.print_to_file("initial_state.png");
    sim.run();
    sim.print_to_file("T_" + std::to_string(temp) + "_final_t_"
            + std::to_string(sweeps) + ".png");
    return 0;
}
