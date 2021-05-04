/*! @file wolff.hpp
 *  @brief wolff cluster flip class header
 *  @author Sameed Perviaz
 *  @copyright GPLv3
 *  @date 2021-04-21
 */

#ifndef WOLFF_HPP
#define WOLFF_HPP

#include <string>
#include <random>

class wolff
{
    public:
        /*! Construct the class and intiialize class invariants */
        wolff(int w, double temp, double field, bool randomize);
        /*! We need this explicitly to free our grid
         * TODO: wrap pointer in std::unique_ptr so that we can rely on default
         * dtors
         */
        ~wolff();

        /*! Performs one iteration of the Wolff algorithm, preserving detailed balance */
        int iterate();
        /*! Convenience wrapper for performing multiple iterations sequentially */
        void run(int sweeps);
        /*! Uses libpng to print a grayscale representation of the current grid state */
        void print_to_file(std::string filename, int scale);

        /*! Set the simulation temperature */
        void set_T(double new_temp);

        // Total calculations (performs n^2 calculations)
        /*! Calculates the *total* grid magnetization */
        void calculate_M();
        /*! Calculates the *total* grid interaction energy, including the field */
        void calculate_U();

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

        /*! Sets the grid state to a random configuration of 1s or 0s */
        void randomize_grid();

        /*! Gets a pointer to the underlying spin element, stored row-major in grid */
        int *get_spin(int x, int y);

        // Grid variables
        int* grid;
        int width;

        // Dynamics
        double T;
        double H;
        // Magnetization, aka total spin
        double M;
        // Total energy (incl. field energy)
        double U;

        // Utilities
        std::mt19937 engine;
};

#endif
