# Define these macros to match with your project's filenames
TARGET = RatesMC2
OBJS = RatesMC2.o

#for icc
#LIBS = -lgsl -lgslcblas -lm -lstdc++
#CPP = icc
#CPPFLAGS = -O3 -parallel

#for gcc
LIBS = -lgsl -lgslcblas -lm -lstdc++
CPP = gcc
CPPFLAGS = -O3 -Wall


.SUFFIXES: .cpp

$(TARGET): $(OBJS)
	$(CPP) $(CPPFLAGS) -o $@ $(OBJS) $(LIBS)

# `make clean' removes extra files that crop up during development
clean:
	-rm $(OBJS) 
	-rm *~

.cpp.o:
	$(CPP) $(CPPFLAGS) -c $< 

docs:
	cd doc && ${MAKE}
