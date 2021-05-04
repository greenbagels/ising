/*! @file sim.hpp
 *  @brief Ising model simulation class interface
 *  @author Sameed Perviaz
 *  @copyright GPLv3
 *  @date 2021-05-04
 */

#ifndef SIM_HPP
#define SIM_HPP

#include <string>
#include <random>
#include <functional>
#include <memory>

namespace ising
{
    class simulation
    {
        public:
            /*! Construct the class and intiialize class invariants */
            simulation(int w, double temp, double field,
                    bool randomize, int backend);

            /*! Performs one iteration of the Wolff algorithm, preserving detailed balance */
            int iterate();
            /*! Convenience wrapper for performing multiple iterations sequentially */
            void run(int sweeps);
            /*! Uses libpng to print a grayscale representation of the current grid state */
            void print_to_file(std::string filename, int scale);

            /*! Set the simulation temperature */
            void set_T(double new_temp);

            /*! Sets the grid state to a random configuration of 1s or 0s */
            void randomize_grid();

            /*! Calculates the equal-time correlation function at separation r */
            double calc_cf(int r);

            /*! Exposes the current temperature */
            double get_T() const;
            /*! Exposes the total magnetization */
            double get_M() const;
            /*! Exposes the total energy */
            double get_U() const;

        private:
            /*! Convenience wrapper for <random> uniform int distribution */
            int rand_int(int a, int b);
            /*! Convenience wrapper for <random> uniform real distribution */
            double rand_dbl(double a, double b);

            // Total calculations (performs n^2 calculations)
            /*! Calculates the *total* grid magnetization */
            void calculate_M();
            /*! Calculates the *total* grid interaction energy, including the field */
            void calculate_U();

            /*! Gets a pointer to the underlying spin element, stored row-major in grid */
            int *get_spin(int x, int y);

            /*! Probabilistically attempts to flip a single spin */
            int single_spin_flip();
            /*! Probabilistically forms a cluster and flips the built cluster */
            int wolff_cluster_flip();

            // Grid variables
            std::unique_ptr<int[]> grid;
            int width;

            // Determines whether to reset the grid bet
            bool sanitize_every_run;

            // Dynamics
            double T;
            double H;
            // Magnetization, aka total spin
            double M;
            // Total energy (incl. field energy)
            double U;

            // Algorithm backend (i.e. single-flip metropolis, etc.)
            int algo_backend;

            // Command pattern for the backend
            std::vector<std::function<int()>> spin_algos =
            {
                [this](){return this->single_spin_flip();},
                [this](){return this->wolff_cluster_flip();}
            };
            // Utilities
            std::mt19937 engine;
    };
}
#endif
