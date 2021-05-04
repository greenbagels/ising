# ising++

A small C++ library to perform Ising model simulations of a ferromagnet.

# Usage

The sample test program contained in `main.cpp` is configured by default to
compute spin-spin correlation functions as a function of temperature, as well
as measure average equilbrium energies and magnetizations (with their variances)
for a sane set of default options. To do this, just build and run:

```
make
./ising
```

But if you want, you can tailor the program to run large sims, sims of different
dimensions, saving image snapshots, and more. To see a full list of options, pass
the `--help` option:

```
./ising --help
```

The test program (with its Makefile) requires C++17 for `std::filesystem` goodies.
The library itself should only require C++11 or so, but I haven't double-checked
for sure (so YMMV on older systems!).

# Documentation

Documentation is a little short, and is only found in the code for now.

# License

All components are licensed under GPLv3 unless otherwise noted.
