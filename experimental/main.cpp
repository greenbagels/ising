
#include <cstdlib>
#include <iostream>
#include <string>
#include "wolff.hpp"

int main(int argc, char *argv[])
{

    int sweeps = atoi(argv[1]);
    int width = 400;
    int scale = 4;
    double temp = 2.27;
    std::cerr << "Initializing sim with parameters:\n" <<
        "    sweeps = " << sweeps << "\n    scale = " << scale
        << "\n     temp = " << temp << "\n";

    wolff sim(sweeps, width, scale, temp); 
    sim.print_to_file("initial_state.png");
    for (auto i = 1; i <= 100; i++)
    {
        sim.run();
        sim.print_to_file("final_t_" + std::to_string(i*sweeps) + ".png");
    }
    return 0;
}
