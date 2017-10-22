CXX=g++
CFLAGS=-g
LFLAGS=-lm

all: fft.o
	$(CXX) $(LFLAGS) -o fft fft.o

fft.o:
	$(CXX) $(CFLAGS) -c fft.cpp

clean:
	rm -f fft.o fft
	rm -f sin1_org.dat sin1_fft.dat
	rm -f sin2_org.dat sin2_fft.dat
	rm -f gauss1_org.dat gauss1_fft.dat
	rm -f gauss2_org.dat gauss2_fft.dat
	rm -f lorentz_fft.dat
