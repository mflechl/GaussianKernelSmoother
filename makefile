CXX=g++  -std=c++11
INCDIR=.
ROOTINC=$(shell root-config --incdir)
ROOTLIB=$(shell root-config --libs)

all: GaussianKernelSmoother.o test4

GaussianKernelSmoother.o: GaussianKernelSmoother.C
	$(CXX) -I$(INCDIR) -I$(ROOTINC) $(ROOTLIB) -c GaussianKernelSmoother.C

test4: test4.C
	$(CXX) -I$(INCDIR) -I$(ROOTINC) $(ROOTLIB) -o $@ GaussianKernelSmoother.o test4.C

clean:
	rm *.o test4
