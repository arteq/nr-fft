#!/usr/bin/gnuplot

set title "Widma mocy"
plot 'sin1_fft.dat' w l t "Sin #1", 'sin2_fft.dat' w l t "Sin #2", 'gauss1_fft.dat' w l t "Gauss #1", 'gauss2_fft.dat' w l t "Gauss #2"

set term png enh size 800,600
set output 'plot.png'
plot 'lorentz_fft.dat' w l t "Lorentz"
