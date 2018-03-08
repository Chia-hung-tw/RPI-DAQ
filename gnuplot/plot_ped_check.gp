data = sprintf("gnuplot/ped_check.dat")

set key font ",14"
set key bottom
set pointsize 0.7

set ylabel "ADC"
set ylabel font ",12" offset char -1,0
set xlabel "channel"
set xlabel font ",12" offset char 0,0.5
set ytics font ",12"
set xtics font ",12"
set yrange[100:]

set terminal png size 1080,810

set output sprintf('Plot_out/ped_check.png')
plot data using 1:2 with points pt 7 lc rgb "red"  title "HighGain",\
     data using 1:3 with points pt 21 lc rgb "blue" title "LowGain"