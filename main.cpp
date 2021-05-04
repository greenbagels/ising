/*! @file main.cpp
 *  @brief Sample test program for the ising::simulation class
 *  @author Sameed Pervaiz (pervaiz.8@osu.edu)
 *  @copyright GPLv3
 *  @date 2021-05-04
 */

// For time benchmarking
#include <chrono>

// C Standard Library headers
#include <cmath>
#include <cstdlib>
#include <cstring>

// C++ Standard Library headers
#include <filesystem>
/* Not supported by many compilers, but essentially the same
 * as libfmt.
#include <format>
 */
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <thread>

// Project headers
#include "sim.hpp"
#include "misc.hpp"

// Third-party headers
// Commonly provided by 'libfmt' packages (at least, on Debian). Check out
// their GitHub, or your system's packaging documentation for info on obtaining
// a copy.
#include <fmt/core.h>
#include <boost/program_options.hpp>

void do_block_calc(double Tmin, double dT, int count, double H, int width, bool randomize,
        double *U_avg, double *U_fluc, double *M_avg, double *M_fluc, double *correl_fun, int backend);

int main(int argc, char *argv[])
{
    namespace po = boost::program_options;
    po::options_description desc("Supported options");

    desc.add_options()
        ("help", "Print help info")
        ("width", po::value<int>()->default_value(16), "Set grid width and height")
        ("scale", po::value<int>()->default_value(1), "Set post-process image scaling")
        ("tmin", po::value<double>()->default_value(0.), "Set starting lattice temperature")
        ("tmax", po::value<double>()->default_value(5.), "Set final lattice temperature")
        ("dt", po::value<double>()->default_value(0.1), "Set lattice temperature step")
        ("field", po::value<double>()->default_value(0.0), "Set external field strength")
        ("backend", po::value<int>()->default_value(ising::BACKEND_WOLFF), "Set equilbration backend")
        ("threads", po::value<unsigned>()->default_value(1u), "Set number of threads for temp benchmarking")
        ("output-dir", po::value<std::string>()->default_value("output"), "Set output directory")
        ("randomize", po::bool_switch(), "Use random initial grid");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
        std::cout << desc << "\n";
        return 1;
    }

    auto width = vm["width"].as<int>();
    auto scale = vm["scale"].as<int>();
    auto tmin = vm["tmin"].as<double>();
    auto tmax = vm["tmax"].as<double>();
    auto dt = vm["dt"].as<double>();
    auto field = vm["field"].as<double>();
    auto backend = vm["backend"].as<int>();
    auto N = static_cast<unsigned>((tmax - tmin) / dt + 1);
    auto nthreads = std::min(vm["threads"].as<unsigned>(), N);
    auto output_dir = vm["output-dir"].as<std::string>();
    auto randomize = vm["randomize"].as<bool>();


    std::cout << "Initializing sim with parameters:"
        << "\n\twidth = " << width
        << "\n\tscale = " << scale
        << "\n\ttmin = " << tmin
        << "\n\ttmax = " << tmax
        << "\n\tdt = " << dt
        << "\n\tfield = " << field
        << "\n\tbackend = " << (backend ? "Wolff" : "Single Flip")
        << "\n\tthreads = " << nthreads
        << "\n\toutput-dir = " << output_dir
        << "\n\trandomize = " << randomize
        << std::endl;

    // Let's make a directory for our output data
    // Example: "output/<width>/<temp>/<field>/<imgs>"
    output_dir += "/" + std::to_string(width) + "/";
    std::filesystem::create_directories(output_dir);
    // auto format_width = std::strlen(argv[1]);
    // Automatically managed array for time benchmarking
    // std::unique_ptr<double[]> times(new double[sweeps]);
    // Initialize the sim... if we weren't benchmarking
    //wolff sim(width, scale, temp, field, randomize);

    // This records our variables as a function of temperature
    std::ofstream output_file(output_dir + "temp_vars.dat");
    output_file << "#T    U    U_fluc    M    M_fluc\n";

    std::unique_ptr<double[]> U_avg (new double[N]());
    std::unique_ptr<double[]> U_fluc(new double[N]());
    std::unique_ptr<double[]> M_avg (new double[N]());
    std::unique_ptr<double[]> M_fluc(new double[N]());

    // We need a correlation function for each temperature. But for a given
    // temperature, the CF is a function of x in [1, width/2]
    std::unique_ptr<double[]> correl_fun(new double[N * (width / 2)]());

    std::unique_ptr<std::thread[]> threads(new std::thread[nthreads]);

    int offset = 0;
    for (auto i = 0; i < nthreads; i++)
    {
        int blocksize = N / nthreads;
        // There could be a remainder left; we could allot it all to the last
        // thread (the lazy way), or, to make the work-items more evenly
        // distributed, we add an extra split it across the first (remainder)
        // items.
        if (i < N % nthreads)
        {
            blocksize++;
        }
        std::cerr << "Thread " << i << " handles indices in range [" << offset << "," << offset + blocksize - 1 << "]\n";

        std::cout << "Spawning thread " << i << std::endl;
        threads[i] = std::thread(do_block_calc, tmin + dt * offset, dt, blocksize, field, width,
                randomize, &U_avg[offset], &U_fluc[offset], &M_avg[offset], &M_fluc[offset], &correl_fun[offset * (width/2)],
                backend);
        offset += blocksize;
    }

    // Thread join loop, to synchronize
    for (auto i = 0; i < nthreads; i++)
    {
        threads[i].join();
        std::cout << "Thread " << i << " joined!"<< std::endl;
    }

    for (auto i = 0; i < N; i++)
    {
        output_file << tmin + i*dt << " " << U_avg[i] << " " << U_fluc[i]
            << " " << M_avg[i] << " " << M_fluc[i] << std::endl;
        std::ofstream cf_ofile(output_dir + "cf_" + std::to_string(tmin+i*dt) +".dat");

        for (auto j = 0; j < width/2; j++)
        {
            cf_ofile << j << " " << correl_fun[i * (width/2) + j] << std::endl;
        }
    }

    return 0;
}

// TODO: rename
void do_block_calc(double Tmin, double dT, int count, double H, int width,
        bool randomize, double *U_avg, double *U_fluc, double *M_avg, double *M_fluc, double *correl_fun, int backend)
{
    // spawning many ephemeral threads is inefficient, so we use persistent
    // pthreads to handle our data.
    //
    // let's say we're calculating N different temperatures. We want to evolve
    // an ensemble of N different configurations through time, measuring the
    // dynamical variables of each configuration through time.
    // Now, make a sim object for this block
    ising::simulation sim(width, Tmin, H, randomize, backend);
    // Time for pre-sweeps! We should probable scale this with the temperature,
    // but for now, we just do a large amount for each starting temperature, and
    // hope the scaling works out well by reusing the grid
    sim.run(20000);

    for (auto i = 0; i < count; i++)
    {
        auto T = Tmin + i * dT;
        /* This code is to check for bad initialization of our arrays
        if (U_fluc[i] || U_avg[i] || M_fluc[i] || M_fluc[i])
            std::terminate();
        */
        sim.set_T(T);
        // Now, equilibrate 500 times to adjust from the old temperature
        sim.run(500);

        // Now, start the averaging loop!
        std::vector<double> M, U;
        for (auto n = 0; n < 50000; n++)
        {
            // We want to record data while minimizing autocorrelations. To do this, we
            // iterate until at least width*width spins have flipped; we cannot stop in
            // the middle of an iteration without violating detailed balance, so we just
            // aim to overshoot this target as little as possible.
            for (int flips = 0; flips < width * width; flips += sim.iterate());

            // Now, push back the variables of interest
            M.push_back(sim.get_M());
            U.push_back(sim.get_U());
            /*
            if (std::fabs(sim.get_U()) > 2*width*width)
            {
                std::terminate();
            }
            */

            // The average spin is just the Magnetization / (width * width), so
            // we calculate <spin1 * spin2> - <spin>, where spin1 and spin2 are
            // a distance r apart:

            // Consider loop re-odering if this is too slow!
            for (auto r = 0; r < width / 2; r++)
            {
                correl_fun[i * (width/2) + r] += sim.calc_cf(r);
            }
        }

        for (auto r = 0; r < width / 2; r++)
        {
            correl_fun[i * (width/2) + r] /= 50000;
        }

        for (auto j = 0; j < M.size(); j++)
        {
            // TODO: if catastrophic cancellation is an issue, implement
            // Kahan summation.
            auto fixedm = M[j];
            if (T < 2.27)
                fixedm = std::fabs(fixedm);
            M_avg[i] += fixedm;
            U_avg[i] += U[j];
        }

        M_avg[i] /= M.size();
        U_avg[i] /= U.size();

        // The fluctuations are just the same thing as the variance, except
        // we need to remember that E[X^2]-E[X]^2 is not a good idea with floats.
        // So like with the single-spin-flip case, we calculate the average of
        // the residuals (using Kahan summation if roundoff becomes problematic)

        for (auto j = 0; j < M.size(); j++)
        {
            auto fixedm = M[j];
            if (T < 2.27)
                fixedm = std::fabs(fixedm);
            M_fluc[i] += (fixedm - M_avg[i]) * (fixedm - M_avg[i]);
            U_fluc[i] += (U[j] - U_avg[i]) * (U[j] - U_avg[i]);
        }
        M_fluc[i] /= (M.size() - 1);
        U_fluc[i] /= (U.size() - 1);
    }
}
