CC=clang-omp++
CFLAGS=-fopenmp -Imodels -C src
LFLAGS=-liomp
LDFLAGS=-L/usr/local/lib -L/Users/Lucy/libomp_oss/exports/mac_32e/lib.thin -lgsl

BASEOBJ = data.o likelihood.o main.o MCMC.o model.o parameter.o particle.o pfilter.o trajectory.o
SIROFFSPRINGOBJ = SIR_offspring_distribution.o

SIRoffspring : $(BASEOBJ) $(SIROFFSPRINGOBJ)
	$(CC) $(CFLAGS) $(LFLAGS) $(LDFLAGS) $(BASEOBJ) $(SIROFFSPRINGOBJ) -o SIRoffspring

data.o : data.cpp data.h
	$(CC) $(CFLAGS) data.cpp -o data.o

likelihood.o : likelihood.h likelihood.cpp parameter.h
	$(CC) $(CFLAGS) likelihood.cpp

MCMC.o : MCMC.h MCMC.cpp pfilter.h model.h parameter.h trajectory.h data.h likelihood.h particle.h
	$(CC) $(CFLAGS) MCMC.cpp

model.o : model.h model.cpp trajectory.h parameter.h
	$(CC) $(CFLAGS) model.cpp

parameter.o : parameter.h parameter.cpp
	$(CC) $(CFLAGS) parameter.cpp

particle.o : particle.h particle.cpp trajectory.h
	$(CC) $(CFLAGS) particle.cpp

pfilter.o : pfilter.h pfilter.cpp model.h parameter.h trajectory.h data.h likelihood.h particle.h
	$(CC) $(CFLAGS) trajectory.cpp

trajectory.o : trajectory.h trajectory.cpp
	$(CC) $(CFLAGS) trajectory.cpp

SIR_offspring_distribution.o : model.h SIR_offspring_distribution.cpp
	$(CC) $(CFLAGS) SIR_offspring_distribution.cpp

main.o : main.cpp MCMC.h pfilter.h trajectory.h model.h parameter.h data.h likelihood.h particle.h
	$(CC) $(CFLAGS) main.cpp

clean : rm *.o *~ SIRoffspring

tar : tar cfv *.h *.cpp models/SIR_offspring_distribution.cpp
