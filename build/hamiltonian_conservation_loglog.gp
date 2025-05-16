set terminal pngcairo size 1200,800 enhanced font 'Helvetica,14'
set output 'hamiltonian_conservation_loglog.png'
set logscale x
set logscale y
set xlabel 'Number of Steps'
set ylabel 'Hamiltonian |ΔH|'
set title 'RATTLE Hamiltonian Conservation'
set grid
plot 'hamiltonian_conservation.dat' using 1:2 with linespoints lw 2 title 'ΔH'
unset output
