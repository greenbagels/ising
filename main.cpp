// C++ Standard Libraries
#include <iostream>
#include <stdexcept>
#include <string>

// C Standard Libraries
#include <cmath>

// System Libraries
#include <omp.h>
#include <png++/png.hpp>
#include <boost/program_options.hpp>

// Project Libraries
#include "ising.hpp"

int main(int argc, char* argv[])
{
    /*
     * argv: {"./ising size iters mode
     *
     * size: width of grid
     * iters: number of sweeps
     * mode: which operation to run
     *   - 1: print only final image
     *   - 2: print
     *
     */
    namespace po = boost::program_options;
    po::options_description desc("Supported options");

    desc.add_options()
        ("help", "Print help info")
        ("benchmark", "Perform single-iteration benchmark")
        ("cli", "Print images to the terminal")
        ("width", po::value<std::size_t>()->default_value(16), "Set grid width and height")
        ("nimg", po::value<std::size_t>()->default_value(1), "Set number of images to save")
        ("scale", po::value<unsigned>()->default_value(1), "Set post-process image scaling")
        ("sweeps", po::value<std::size_t>()->default_value(16*16*100), "Set minimum spin-flip iteration count")
        ("temp", po::value<double>()->default_value(2.5), "Set lattice temperature")
        ("fstr", po::value<double>()->default_value(0.0), "Set external field strength")
        ("neighbors", po::value<unsigned>()->default_value(4), "Set number of cell neighbors")
        ("backend", po::value<std::string>()->default_value(std::string("cpu-serial")), "Set parallelism backend");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
        std::cout << desc << "\n";
        return 1;
    }

    if (vm["backend"].as<std::string>() != "cpu-serial" &&
        vm["backend"].as<std::string>() != "cpu-parallel" &&
        vm["backend"].as<std::string>() != "gpu")
    {
        std::cerr << "Invalid backend specified. \
            Valid options are 'cpu-serial', 'cpu-parallel', and 'gpu'\n";
        return 1;
    }

    if (vm["temp"].as<double>() < 0)
    {
        std::cerr << "Temperature cannot be negative!\n";
        return 1;
    }

    if (vm["nimg"].as<std::size_t>() > vm["sweeps"].as<std::size_t>())
    {
        std::cerr << "Cannot save more images than sweeps!\n";
        return 1;
    }

    ising sim(vm["sweeps"].as<std::size_t>(),
              vm["width"].as<std::size_t>(),
              vm["neighbors"].as<unsigned>(),
              vm["nimg"].as<std::size_t>(),
              vm["scale"].as<unsigned>(),
              vm["temp"].as<double>(),
              vm["fstr"].as<double>(),
              vm["backend"].as<std::string>());

    if (vm.count("cli"))
    {
        sim.set_display_mode(0);
    }
    else
    {
        sim.set_display_mode(1);
    }

    if (vm.count("benchmark"))
    {
        sim.set_benchmark_mode(1);
    }
    else
    {
        sim.set_benchmark_mode(0);
    }

    sim.run();

    return 0;
}
