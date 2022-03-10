
#include "minimizer.hpp"

#include <cmath>
#include <cstdio>

#include <array>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlinear.h>

// Lazy
using arr = std::vector<double>;
using arrref = std::vector<double>&;

void read_data(std::string fname, arrref t1, arrref y1, arrref dy1);
int exp_f (const gsl_vector *x, void *data, gsl_vector *f);
int exp_df (const gsl_vector *x, void *data, gsl_matrix *J);
int exp_fvv (const gsl_vector *x, const gsl_vector *v, void *data, gsl_vector *fvv);

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage: ./executable filename\n";
        return 1;
    }

    arr t1,y1,dy1;
    std::string input_fname(argv[1]);
    double T;
    // The input is of the form cf_%lf.dat
    std::sscanf(input_fname.c_str(), "cf_%lf.dat", &T);

    // Read chan1, chan2 data
    read_data(input_fname, t1, y1, dy1);

    std::vector<double> init_vals = {{y1[0], -5.}};

    // Fit the first set of parameters, aka CH_A, the input signal
    minimizer::nonlinear_ls<std::vector<double>> fitter(t1, y1, dy1, exp_f, exp_df, exp_fvv, init_vals); 
    fitter.fit();
    fitter.print_results();

    // We're going to loop over the files in bash, so just append for now
    std::ofstream cl_file("corr_len.dat", std::ios::app);

    cl_file << T << " " << 1. / fitter.parameter(1) << std::endl;
    return 0;
}

void read_data(std::string fname, arrref t1, arrref y1, arrref dy1)
{
    std::size_t line_count = 0u;
    std::ifstream infile(fname);
    if (!infile.is_open())
    {
        throw std::runtime_error("Error opening input file!\n");
    }
    for (std::string temp; std::getline(infile, temp);)
    {
        line_count++;
    }

    // Clears EOF bit as of C++11
    infile.clear();
    infile.seekg(0);

    for (auto i = 0u; i < line_count; i++)
    {
        double temp;
        infile >> temp;
        t1.push_back(temp);
        infile >> temp;
        y1.push_back(temp);
        infile >> temp;
        dy1.push_back(1/(temp*temp));
    }
}

int exp_f (const gsl_vector *x, void *data, gsl_vector *f)
{
    // y = A * exp(B*x);

    minimizer::user_data *sdata = reinterpret_cast<minimizer::user_data*>(data);

    std::size_t n = sdata->n;
    double *t = sdata->t;
    double *y = sdata->y;

    double A = gsl_vector_get(x, 0);
    double B = gsl_vector_get(x, 1);

    for (auto i = 0u; i < n; i++)
    {
        double yi = A * std::exp(B * t[i]);
        gsl_vector_set(f, i, yi - y[i]);
    }

    return GSL_SUCCESS;
}

int exp_df (const gsl_vector *x, void *data, gsl_matrix *J)
{
    // y = A * exp(B*x)

    // dy/dA = exp(B*x)
    // dy/dB = A * x * exp(B*x)

    minimizer::user_data *sdata = reinterpret_cast<minimizer::user_data*>(data);

    std::size_t n = sdata->n;
    double *t = sdata->t;
    // double *y = sdata->y;

    double A = gsl_vector_get(x, 0);
    double B = gsl_vector_get(x, 1);

    for (auto i = 0u; i < n; i++)
    {
        gsl_matrix_set (J, i, 0, std::exp(B * t[i]));
        gsl_matrix_set (J, i, 1, A * t[i] * std::exp(B * t[i]));
    }

    return GSL_SUCCESS;
}

int exp_fvv (const gsl_vector *x, const gsl_vector *v, void *data, gsl_vector *fvv)
{
    // y = A * exp(B*x)

    // First partials
    // dy/dA = exp(B*x))
    // dy/dB = A * x * exp(B*x)

    // Second non-mixed partials
    // d^2y/dA^2 = 0
    // d^2y/dB^2 = A * x^2 * exp(B*x)
    
    // Mixed partials (equal by clairaut, so double them)
    // d^2y/dAdB = x * exp(B*x)
    
    minimizer::user_data *sdata = reinterpret_cast<minimizer::user_data*>(data);

    auto n = sdata->n;
    auto t = sdata->t;
    // auto y = sdata->y;

    double A = gsl_vector_get(x, 0);
    double B = gsl_vector_get(x, 1);

    for (auto i = 0u; i < n; i++)
    {
        auto va = gsl_vector_get(v, 0);
        auto vb = gsl_vector_get(v, 1);

        auto fab = t[i] * std::exp(B*t[i]);
        auto fbb = A * t[i] * t[i] * std::exp(B*t[i]);

        auto fvvi = vb * vb * fbb + 2. * (va * vb * fab);

        gsl_vector_set(fvv, i, fvvi);
    }
    return GSL_SUCCESS;
}

