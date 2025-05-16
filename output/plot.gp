reset
set terminal pngcairo size 1000,2000 enhanced font "Arial,12"
set output 'history.png'

set multiplot layout 3,1 title "Tr(A_μ A^μ) Analysis" font ",11"

# 1. Scatter plot: Im vs Re
set xlabel "Re[Tr(A_μ A^μ)]" font ",11"
set ylabel "Im[Tr(A_μ A^μ)]" font ",11"
set grid
set key off
plot '../output/virial.dat' using 1:2 with points lc rgb "blue"

# 2. Real part history
set xlabel "HMC time" font ",11"
set ylabel "Re[Tr(A_μ A^μ)]" font ",11"
plot '../output/virial.dat' using 0:1 with lines lw 1 lc rgb "red" title "Re"

# 3. Imaginary part history
set xlabel "HMC time" font ",11"
set ylabel "Im[Tr(A_μ A^μ)]" font ",11"
plot '../output/virial.dat' using 0:2 with lines lw 1 lc rgb "orange" title "Im"

unset multiplot
