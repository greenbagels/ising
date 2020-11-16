#ifndef ISING_HPP
#define ISING_HPP

#include <memory>
#include <string>

#include "data.hpp"

struct neighbors
{
    int top;
    int bottom;
    int left;
    int right;
};

class ising
{
    public:
        ising(std::size_t sweeps, std::size_t width, unsigned num_neighbors,
                std::size_t nimg, unsigned scale, double temp, double fstr,
                std::string backend);
        double calc_deltaU(unsigned i, unsigned j);
        void set_display_mode(unsigned char mode);
        void set_benchmark_mode(unsigned char mode);
        void run();

    private:
        void initialize_spins();
        void flip_spin(std::size_t i, std::size_t j);
        int get_spin(std::size_t i, std::size_t j) const;
        neighbors get_neighbors(unsigned i, unsigned j);
        double calc_totalU();
        double calc_totalM();
        void save_png_snapshot(const char* fname);
        void print_snapshot();

        std::size_t sweeps;
        std::size_t width;
        unsigned num_neighbors;
        unsigned dim;
        unsigned scale;
        std::size_t nimg;
        double temp;
        double field_strength;
        char backend;

        double total_U;
        double total_M;

        std::unique_ptr<char[]> spins;

        unsigned char display_mode;
        unsigned char benchmark_mode;

        event_data data;
};

#endif
