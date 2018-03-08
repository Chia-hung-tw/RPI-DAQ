
data = sprintf("gnuplot/dac_adc.dat")

#set key top left
set key bot right
set pointsize 0.7
set ylabel "ADC"
set ylabel font ",12" offset char -1,0
set xlabel "DAC"
set xlabel font ",12" offset char 0,0.5
set ytics font ",12"
set xtics font ",12"
set yrange[200:]
#set yrange[0:200]



plot data using 2:3 with points pt 21 lc rgb "red" title "chip0 HG",\
     data using 2:4 with points pt 7 lc rgb "red" title "chip0 LG",\
     data using 2:5 with points pt 21 lc rgb "yellow" title "chip1 HG",\
     data using 2:6 with points pt 7 lc rgb "yellow" title "chip1 LG",\
     data using 2:7 with points pt 21 lc rgb "green" title "chip2 HG",\
     data using 2:8 with points pt 7 lc rgb "green" title "chip2 LG",\
     data using 2:9 with points pt 21 lc rgb "blue"  title "chip3 HG",\
     data using 2:10 with points pt 7 lc rgb "blue" title "chip3 LG"

#plot data using 2:11 with points pt 21 lc rgb "red" title "chip0 TOT",\
#     data using 2:12 with points pt 21 lc rgb "yellow" title "chip1 TOT",\
#     data using 2:13 with points pt 21 lc rgb "green" title "chip2 TOT",\
#     data using 2:14 with points pt 21 lc rgb "blue" title "chip3 TOT"