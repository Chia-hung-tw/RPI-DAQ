all: acquisition decode_ana take_pedestal

####################################

acquisition: main.o gpiolib.o
	gcc main.o gpiolib.o bcm2835.o -lm -o acquisition

main.o: main.c
	gcc -c -I ./RPi_software/bcm2835-1.52/src ./RPi_software/bcm2835-1.52/src/bcm2835.c main.c

gpiolib.o: gpiolib.c
	gcc -c -I ./RPi_software/bcm2835-1.52/src ./RPi_software/bcm2835-1.52/src/bcm2835.c gpiolib.c

####################################

decode_ana: decodefromraw.o CHMapping.o
	g++ decodefromraw.o CHMapping.o -o decode_ana

decodefromraw.o: ./quick_analysis/decode_main.cc ./quick_analysis/CHMapping.h
	g++ -g -c ./quick_analysis/decode_main.cc -o decodefromraw.o

CHMapping.o: ./quick_analysis/CHMapping.cc ./quick_analysis/CHMapping.h
	g++ -g -c ./quick_analysis/CHMapping.cc

#####################################

take_pedestal: main_ped.o gpiolib.o
	gcc main_ped.o gpiolib.o bcm2835.o -lm -o take_pedestal

main_ped.o: main_ped.c
	gcc -c -I ./RPi_software/bcm2835-1.52/src ./RPi_software/bcm2835-1.52/src/bcm2835.c main_ped.c

#####################################

clean:
	rm -f acquisition decode_ana take_pedestal
	rm -f *.o
