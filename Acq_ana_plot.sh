argc=$#
if [ $argc != 0 ]; then
    echo "Input rawfile mode"
    filecount=0
    while [ $filecount != $argc ]; do
	let "filecount = filecount + 1"
	tolog=$1
	sed -i '1c '$tolog'' runinfo.log
	shift	
	echo "./decode_ana"
	./decode_ana
	echo "reading information for gnuplot..."
	linec=0
	echo ""
	while IFS='' read -r line || [[ -n "$line" ]]; do
	    let "linec = linec+1"
	    if [ $linec == 1 ]; then
		echo "Run information: $line"
		file=$line
	    elif [ $linec == 2 ]; then
		config1=$line
		if [ $config1 == "CC" ]; then
		    echo "Plotting Connected Channel ..."
		elif [ $config1 == "NC" ]; then
		    echo "Plotting Non-Connected Channel ..."
		else
		    echo "wrong plotting config! config line"$linec "is wired!"
		fi
	    elif [ $linec == 3 ]; then
		config2=$line
		if [ $config2 == "HG" ]; then
		    echo "Plotting High Gain ..."
		elif [ $config2 == "LG" ]; then
		    echo "Plotting Low Gain  ..."
		elif [ $config2 == "HL" ]; then
		    echo "Plotting both High Gain and Low Gain  ..."
		else
		    echo "wrong plotting config! config line"$linec "is wired!"
		fi
	    elif [ $linec == 4 ]; then
		config3=$line
		echo "plotting SCA"$config3
	    fi
	done < "./gnuplot/plotconfig.txt"

	if [ $config2 != "HL" ]; then
	    plotset=$config2"SCA"$config3
	    sed -i '1c fileN = "'$file'"' gnuplot/plot_save.gp
	    sed -i '2c plotset = "'$plotset'"' gnuplot/plot_save.gp
	    echo "saving plots..."
	    gnuplot gnuplot/plot_save.gp
	    echo "Output plot will be: ./Plot_out/"$file"_"$plotset"_*.png"
	else
	    echo "plotting HG first ..."
 	    plotset="HGSCA"$config3
	    sed -i '1c fileN = "'$file'"' gnuplot/plot_save.gp
	    sed -i '2c plotset = "'$plotset'"' gnuplot/plot_save.gp
	    echo "saving plots..."
	    gnuplot gnuplot/plot_save.gp
	    echo "Output plot will be: ./Plot_out/"$file"_"$plotset"_*.png"

	    echo ""
	    echo "plotting LG ... "
	    plotset="LGSCA"$config3
	    sed -i '1c fileN = "'$file'"' gnuplot/plot_save.gp
	    sed -i '2c plotset = "'$plotset'"' gnuplot/plot_save.gp
	    echo "saving plots..."
	    gnuplot gnuplot/plot_save.gp
	    echo "Output plot will be: ./Plot_out/"$file"_"$plotset"_*.png"
	fi
    done

else
    echo "./acquisition"
    ./acquisition
    echo "./decode_ana"
    ./decode_ana
    echo "reading information for gnuplot..."
    linec=0
    echo ""
    while IFS='' read -r line || [[ -n "$line" ]]; do
	let "linec = linec+1"
	if [ $linec == 1 ]; then
	    echo "Run information: $line"
	    file=$line
	elif [ $linec == 2 ]; then
	    config1=$line
	    if [ $config1 == "CC" ]; then
		echo "Plotting Connected Channel ..."
	    elif [ $config1 == "NC" ]; then
		echo "Plotting Non-Connected Channel ..."
	    else
		echo "wrong plotting config! config line"$linec "is wired!"
	    fi
	elif [ $linec == 3 ]; then
	    config2=$line
	    if [ $config2 == "HG" ]; then
		echo "Plotting High Gain ..."
	    elif [ $config2 == "LG" ]; then
		echo "Plotting Low Gain  ..."
	    elif [ $config2 == "HL" ]; then
		echo "Plotting both High Gain and Low Gain  ..."
	    else
		echo "wrong plotting config! config line"$linec "is wired!"
	    fi
	elif [ $linec == 4 ]; then
	    config3=$line
	    echo "plotting SCA"$config3
	fi
    done < "./gnuplot/plotconfig.txt"

    if [ $config2 != "HL" ]; then
	plotset=$config2"SCA"$config3
	#sed -i '1c fileN = "'$file'"' gnuplot/plot_show.gp
	#sed -i '2c plotset = "'$plotset'"' gnuplot/plot_show.gp
	#gnuplot gnuplot/plot_show.gp

	echo "saving plots..."
	gnuplot gnuplot/plot_save.gp
	echo "Output plot will be: ./Plot_out/"$file"_"$plotset"_*.png"
	eog -w './Plot_out/'$file'_'$plotset'_'*'.png' &
    else
	echo "plotting HG first ..."
 	plotset="HGSCA"$config3
	#sed -i '1c fileN = "'$file'"' gnuplot/plot_show.gp
	#sed -i '2c plotset = "'$plotset'"' gnuplot/plot_show.gp
	#gnuplot gnuplot/plot_show.gp

	echo "saving plots..."
	sed -i '1c fileN = "'$file'"' gnuplot/plot_save.gp
	sed -i '2c plotset = "'$plotset'"' gnuplot/plot_save.gp	
	gnuplot gnuplot/plot_save.gp
	echo "Output plot will be: ./Plot_out/"$file"_"$plotset"_*.png"

	echo ""
	echo "plotting LG ... "
	plotset="LGSCA"$config3
	#sed -i '1c fileN = "'$file'"' gnuplot/plot_show.gp
	#sed -i '2c plotset = "'$plotset'"' gnuplot/plot_show.gp
	#gnuplot gnuplot/plot_show.gp

	echo "saving plots..."
	sed -i '1c fileN = "'$file'"' gnuplot/plot_save.gp
	sed -i '2c plotset = "'$plotset'"' gnuplot/plot_save.gp
	gnuplot gnuplot/plot_save.gp
	echo "Output plot will be: ./Plot_out/"$file"_"$plotset"_*.png"
	eog -w './Plot_out/'$file'_'*'.png' &
    fi
    echo "done!"    
fi
