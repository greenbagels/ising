#! /usr/bin/env gnuplot

file1 = "output_6.dat"
file2 = "output_8.dat"
file3 = "output_10.dat"

set title "Ising Magnet, h=0"
set term pngcairo enhanced

set xlabel "k_B T/eps"
set ylabel "E/eps"
set output "1.png"
plot file1 using 1:2 title "n=6", \
     file2 using 1:2 title "n=8", \
     file3 using 1:2 title "n=10"

set xlabel "k_B T/eps"
set ylabel "k_B C_v/(eps^2)"
set output "2.png"
set yrange [0:100]
plot file1 using 1:($1 == 0 ? NaN : $3/($1*$1)) title "n=6", \
     file2 using 1:($1 == 0 ? NaN : $3/($1*$1)) title "n=8", \
     file3 using 1:($1 == 0 ? NaN : $3/($1*$1)) title "n=10"
unset yrange

set xlabel "k_B T/eps"
set ylabel "m"
set output "3.png"
plot file1 using 1:4 title "n=6", \
     file2 using 1:4 title "n=8", \
     file3 using 1:4 title "n=10"

set xlabel "k_B T/eps"
set ylabel "|m|"
set output "4.png"
plot file1 using 1:(abs($4)) title "n=6", \
     file2 using 1:(abs($4)) title "n=8", \
     file3 using 1:(abs($4)) title "n=10"

set xlabel "k_B T/eps"
set ylabel "Susceptibility"
set output "5.png"
plot file1 using 1:($1 == 0 ? NaN : $5/$1) title "n=6", \
     file2 using 1:($1 == 0 ? NaN : $5/$1) title "n=8", \
     file3 using 1:($1 == 0 ? NaN : $5/$1) title "n=10"
