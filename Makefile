CCX = g++
CCXFLAGS = -O3  -funroll-loops -march=native -std=c++11 -pthread  -fopenmp  -I. -I./include 
INCLUDEPATH = -I/home/han/MKHE/code/dependecy/NTL/ntl-11.5.1-exce-install/include -I. -I./include 
LIBPATH= -L/home/han/MKHE/code/dependecy/NTL/ntl-11.5.1-exce-install/lib 
DEPS = -lntl -lgmp -lfftw3 -lm

all: clean test
	 
# clean:
# 	$(RM) test test.o lwehe.o ntruhe.o fft.o sampler.o keygen.o libfinal.a libmkhe.a

# clean:
# 	$(RM) test test.o lwehe.o ntruhe.o fft.o sampler.o keygen.o libfinal.a poly.o libmkhe.a MKHEparams.o MKHEkeygen.o MKLwe.o MKHEscheme.o

clean:
	$(RM) test libmkhe.a MKHEscheme.o MKHEkeygen.o MKHEparams.o poly.o  libmkhe.a 
test: FINAL.h libfinal.a libmkhe.a 
	$(CCX) $(CCXFLAGS)  -o test test.cpp  libmkhe.a libfinal.a $(DEPS)



libfinal.a: include/params.h   fft.o sampler.o
	$(AR) -q libfinal.a   fft.o sampler.o


fft.o: include/fft.h
	$(CCX) $(CCXFLAGS) -c src/fft.cpp

sampler.o: include/sampler.h include/params.h src/sampler.cpp
	$(CCX) $(CCXFLAGS) -c src/sampler.cpp


libmkhe.a: MKHEscheme.o NTRUhe.o RLWEhe.o MKLwe.o MKHEkeygen.o MKHEparams.o poly.o  libfinal.a
	$(AR) -q libmkhe.a  MKHEscheme.o NTRUhe.o RLWEhe.o MKLwe.o MKHEkeygen.o MKHEparams.o poly.o

poly.o: include/poly.h libfinal.a
	$(CCX) $(CCXFLAGS)  -c src/poly.cpp
	
MKHEparams.o: include/MKHEparams.h  libfinal.a 
	$(CCX) $(CCXFLAGS) -c src/MKHEparams.cpp

NTRUhe.o: include/NTRUhe.h MKHEparams.o poly.o
	$(CCX) $(CCXFLAGS) -c src/NTRUhe.cpp

MKLwe.o: include/MKLwe.h poly.o 
	$(CCX) $(CCXFLAGS) -c src/MKLwe.cpp

MKHEkeygen.o: include/MKHEkeygen.h NTRUhe.o MKHEparams.o poly.o
	$(CCX) $(CCXFLAGS) -c src/MKHEkeygen.cpp


RLWEhe.o: include/RLWEhe.h MKHEkeygen.o
	$(CCX) $(CCXFLAGS) -c src/RLWEhe.cpp

MKHEscheme.o: include/MKHEscheme.h NTRUhe.o MKLwe.o RLWEhe.o
	$(CCX) $(CCXFLAGS) -c src/MKHEscheme.cpp

