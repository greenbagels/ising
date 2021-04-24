#ifndef WOLFF_HPP
#define WOLFF_HPP

#include <string>
#include <random>

class wolff
{
    public:
        wolff(int w, int scale, double temp);
        ~wolff();

        void iterate();
        void run(int sweeps);
        void print_to_file(std::string filename);
    private:
        void randomize_grid();
        int rand_int(int a, int b);
        double rand_dbl(double a, double b);
        int *get_spin(int x, int y);

        // Grid variables
        int* grid;
        int width;
        int scale;

        // Dynamics
        double T;
        // double H;

        // Utilities
        std::mt19937 engine;
};

#endif
