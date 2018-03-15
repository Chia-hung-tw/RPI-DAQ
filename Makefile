CPPFLAG=-g -std=c++11
all: decode_ana
####################################

decode_ana: decodefromraw.o CHMapping.o
	g++ $(CPPFLAG) decodefromraw.o CHMapping.o -o decode_ana

decodefromraw.o: ./quick_analysis/decode_main.cc ./quick_analysis/CHMapping.h
	g++ -c $(CPPFLAG) ./quick_analysis/decode_main.cc -o decodefromraw.o

CHMapping.o: ./quick_analysis/CHMapping.cc ./quick_analysis/CHMapping.h
	g++ -c $(CPPFLAG) ./quick_analysis/CHMapping.cc

#####################################


clean:
	rm -f decode_ana
	rm -f *.o
