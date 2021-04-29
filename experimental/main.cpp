/*! @file main.cpp
 *  @brief Sample test program for the wolff class
 *  @author Sameed Pervaiz (pervaiz.8@osu.edu)
 *  @copyright GPLv3
 *  @date 2021-04-21
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
#include "wolff.hpp"

// Third-party headers
// Commonly provided by 'libfmt' packages (at least, on Debian). Check out
// their GitHub, or your system's packaging documentation for info on obtaining
// a copy.
#include <fmt/core.h>
#include <omp.h>

void do_block_calc(double Tmin, double dT, int count, double H, int width, bool randomize,
        double *U_avg, double *U_fluc, double *M_avg, double *M_fluc);

int main(int argc, char *argv[])
{
    // argv: {executable, width, Tmin, Tmax, dT, field, scale}

    if (argc < 7)
    {
        std::cerr << "Usage: wolff width tmin tmax dt field scale" << std::endl;
        return 1;
    }
    // Get configuration specified at runtime
    int width = std::atoi(argv[1]);
    double tmin = std::atof(argv[2]);
    double tmax = std::atof(argv[3]);
    double dt = std::atof(argv[4]);
    double field = std::atof(argv[5]);
    int scale = std::atoi(argv[6]);

    // Should we use random ICs, or a predetermined seed for testing?
    bool randomize = true;

    std::cout << "Initializing sim with parameters:"
        << "\n\twidth = " << width
        << "\n\ttmin = " << tmin
        << "\n\ttmax = " << tmax
        << "\n\tdt = " << dt
        << "\n\tfield = " << field
        << "\n\tscale = " << scale
        << std::endl;

    // Let's make a directory for our output data
    // Example: "output/<width>/<temp>/<field>/<imgs>"
    // std::string output_dir = "output/" + std::to_string(width) + "/" + std::to_string(temp);
    // std::filesystem::create_directories(output_dir);
    // auto format_width = std::strlen(argv[1]);
    // Automatically managed array for time benchmarking
    // std::unique_ptr<double[]> times(new double[sweeps]);
    // Initialize the sim... if we weren't benchmarking
    //wolff sim(width, scale, temp, field, randomize);

    // This records our variables as a function of temperature
    std::ofstream output_file("temp_data_long_run.dat");
    output_file << "#T    U    U_fluc    M    M_fluc\n";

    auto nthreads = std::thread::hardware_concurrency();
    auto N = static_cast<int>((tmax - tmin)/dt + 1);

    std::unique_ptr<double[]> U_avg (new double[N]());
    std::unique_ptr<double[]> U_fluc(new double[N]());
    std::unique_ptr<double[]> M_avg (new double[N]());
    std::unique_ptr<double[]> M_fluc(new double[N]());

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

        std::cout << "Spawning thread " << i << std::endl;
        threads[i] = std::thread(do_block_calc, tmin + dt * offset, dt, blocksize, field, width,
                randomize, &U_avg[offset], &U_fluc[offset], &M_avg[offset], &M_fluc[offset]);
        offset += blocksize;
    }

    // Thread join loop, to synchronize
    for (auto i = 0; i < nthreads; i++)
    {
        threads[i].join();
        std::cout << "Thread" << i << " joined!"<< std::endl;
    }

    for (auto i = 0; i < N; i++)
    {
        output_file << tmin + i*dt << " " << U_avg[i] << " " << U_fluc[i]
            << " " << M_avg[i] << " " << M_fluc[i] << std::endl;
    }
    // Old loop for benchmarking, etc.
    // TODO: move this to its own function.

    /*
    sim.print_to_file(output_dir + fmt::format("/{:0{}d}", 0, format_width) + ".png");
    for (auto i = 0; i < steps; i++)
    {
        // We keep an interior loop instead of just calling sim.run(sweeps/steps)
        // to allow later changes like adding std::chrono benchmarking (which we
        // did in the past) or other inter-flip steps
        for (auto j = 0; j < sweeps / steps; j++)
        {
            sim.iterate();
        }
        sim.print_to_file(output_dir + fmt::format("/{:0{}d}", sweeps/steps * (i+1), format_width) + ".png");
    }
    */

    return 0;
}

// TODO: rename
void do_block_calc(double Tmin, double dT, int count, double H, int width,
        bool randomize, double *U_avg, double *U_fluc, double *M_avg, double *M_fluc)
{
    // spawning many ephemeral threads is inefficient, so we use persistent
    // pthreads to handle our data.
    //
    // let's say we're calculating N different temperatures. We want to evolve
    // an ensemble of N different configurations through time, measuring the
    // dynamical variables of each configuration through time.

    // Now, make a sim object for this block
    wolff sim(width, Tmin, H, randomize);
    // Time for pre-sweeps! We should probable scale this with the temperature,
    // but for now, we just do a large amount for each starting temperature, and
    // hope the scaling works out well by reusing the grid
    sim.run(50000);

    for (auto T = Tmin; T < Tmin + count * dT; T += dT)
    {
        int idx = static_cast<int>((T - Tmin) / dT);
        sim.set_T(T);
        // Now, equilibrate 5000 times to adjust from the old temperature
        sim.run(5000);

        // Now, start the averaging loop!
        std::vector<double> M, U;
        for (auto n = 0; n < 5000; n++)
        {
            // We want to record data while minimizing autocorrelations. To do this, we
            // iterate until at least width*width spins have flipped; we cannot stop in
            // the middle of an iteration without violating detailed balance, so we just
            // aim to overshoot this target as little as possible.
            for (int flips = 0; flips < width * width; flips += sim.iterate());

            // Now, push back the variables of interest
            M.push_back(sim.get_M());
            U.push_back(sim.get_U());
        }

        for (auto i = 0; i < M.size(); i++)
        {
            // TODO: if catastrophic cancellation is an issue, implement
            // Kahan summation.
            M_avg[idx] += M[i] / M.size();
            M_fluc[idx] += M[i] * M[i] / M.size();
            U_avg[idx] += U[i] / U.size();
            U_fluc[idx] += U[i] * U[i] / U.size();
        }

        M_fluc[idx] -= M_avg[idx] * M_avg[idx];
        U_fluc[idx] -= U_avg[idx] * U_avg[idx];
    }
}
