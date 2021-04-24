
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fmt/core.h>
#include <iostream>
#include <string>
#include "wolff.hpp"

int main(int argc, char *argv[])
{
    // argv: {executable, sweeps, steps, width, temp}
    if (argc < 5)
    {
        std::cerr << "Usage: wolff sweeps steps width temp" << std::endl;
        return 1;
    }
    // Get configuration specified at runtime
    int sweeps = std::atoi(argv[1]);
    int steps = std::atoi(argv[2]);
    int width = std::atoi(argv[3]);
    double temp = std::atof(argv[4]);
    // Scale up our image 4x for clarity
    int scale = 4;

    std::cerr << "Initializing sim with parameters:\n" <<
        "    sweeps = " << sweeps << "\n    scale = " << scale
        << "\n     temp = " << temp << "\n";

    // Initialize the sim
    wolff sim(width, scale, temp);
    // Let's make a directory for our output data
    // Example: "output/400/2.4/<imgs>"
    std::string output_dir = "output/" + std::to_string(width) + "/" + std::to_string(temp);
    std::filesystem::create_directories(output_dir);
    auto format_width = std::strlen(argv[1]);
    sim.print_to_file(output_dir + fmt::format("/{:0{}d}", 0, format_width) + ".png");
    for (auto i = 0; i < steps; i++)
    {
        sim.run(sweeps / steps);
        sim.print_to_file(output_dir + fmt::format("/{:0{}d}", sweeps/steps * (i+1), format_width) + ".png");
    }
    return 0;
}
