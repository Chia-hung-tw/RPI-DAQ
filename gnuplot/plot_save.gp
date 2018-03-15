fileN = "Module_1_RUN_080318_0825"
plotset = "LGSCA0"

set style fill solid noborder
set style circle radius 0.45
set terminal png size 1080,810

# start value for H
h1 = 1/360.0
# end value for H
h2 = 250/360.0
# creating the palette (its the same as the above defined one)
set palette model HSV functions (1-gray)*(h2-h1)+h1,1,0.68
set termoption enhanced
set ylabel "mm"
set ylabel font ",12" offset char 1,0
set xlabel "mm"
set xlabel font ",12" offset char 0,0.5
set ytics font ",12"
set ytics format "{/:Bold %.0s}"
set xtics font ",12"
set xtics format "{/:Bold %.0s}"
set cbtics font ",12"
set cbtics format "{/:Bold %.0f}"
set cblabel "avg (ADC)" offset char 1.5,0
set cblabel font ",16" 
set title   font ",16" noenhanced
unset key

data = sprintf("gnuplot/%s.dat",plotset)
frame = "gnuplot/gnu_frame.dat"
Title = sprintf("%s %s",fileN,plotset)

Title_avg = sprintf("%s avg",Title)
Title_noise = sprintf("%s noise",Title)
Title_chinfo = sprintf("%s CHinfo",Title)

Save_avg = sprintf("%s_%s_avg.png",fileN,plotset)
Save_noise = sprintf("%s_%s_noise.png",fileN,plotset)
Save_chinfo = sprintf("%s_%s_CHinfo.png",fileN,plotset)


set title Title_avg
set output sprintf('Plot_out/%s',Save_avg)
plot data using 1:2:3 with circles lc palette, \
     data using 1:2:(sprintf("%i",$3)) with labels font ",12", \
     frame with lines lw 2 lc "black"

set cblabel "noise (ADC)" offset char 1.5,0

set title Title_noise
set output sprintf('Plot_out/%s',Save_noise)
plot data using 1:2:4 with circles lc palette, \
     data using 1:2:4 with labels font ",10", \
     frame with lines lw 2 lc "black"

set title Title_chinfo
set output sprintf('Plot_out/%s',Save_chinfo)
plot data using 1:2:4 with circles lc palette, \
     data using 1:2:(sprintf("(%i, %i)", $5, $6)) with labels font ",10" ,\
     frame with lines lw 2 lc "black"