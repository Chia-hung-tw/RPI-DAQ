echo "./take_pedestal"
./take_pedestal
echo "./decode_ana"
./decode_ana
echo "gnuplot gnuplot/plot_ped_check.gp"
gnuplot gnuplot/plot_ped_check.gp

echo "reading information for gnuplot..."
linec=0
echo ""
while IFS='' read -r line || [[ -n "$line" ]]; do
    let "linec = linec+1"
    if [ $linec == 1 ]; then
	echo "Run information: $line"
	file=$line
    fi
done < "./gnuplot/plotconfig.txt"

mv Plot_out/ped_check.png Plot_out/${file}_ped_check.png
eog -w Plot_out/${file}_ped_check.png &
