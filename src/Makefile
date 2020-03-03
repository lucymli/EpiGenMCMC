CC=g++
CFLAGS=-fopenmp -std=c++11
LDFLAGS=-L/usr/local/lib
LGSLFLAGS=-lgsl -lgslcblas -lm

BASEOBJ = data.o likelihood.o main.o MCMC.o model.o parameter.o particle.o pfilter.o trajectory.o
COMMAND = $(CC) $(CFLAGS) $(LDFLAGS) $(BASEOBJ)

SIROFFSPRINGOBJ = SIR_offspring_distribution.o
SIRoffspring : $(BASEOBJ) $(SIROFFSPRINGOBJ)
	$(COMMAND) $(SIROFFSPRINGOBJ) -o SIRoffspring $(LGSLFLAGS)

BRANCHINGEBOLALIBERIAOBJ = branching_process_EbolaLiberia.o
EbolaLiberiaBranching : $(BASEOBJ) $(BRANCHINGEBOLALIBERIAOBJ)
	$(COMMAND) $(BRANCHINGEBOLALIBERIAOBJ) -o EbolaLiberiaBranching $(LGSLFLAGS)

BRANCHINGPOLIOPAKENVOBJ = branching_process_PakEnv.o
PolioPakEnvBranching : $(BASEOBJ) $(BRANCHINGPOLIOPAKENVOBJ)
	$(COMMAND) $(BRANCHINGPOLIOPAKENVOBJ) -o PolioPakEnvBranching $(LGSLFLAGS)

BRANCHINGCOVID19OBJ = branching_process_covid19.o
covid19branching : $(BASEOBJ) $(BRANCHINGCOVID19OBJ)
	$(COMMAND) $(BRANCHINGCOVID19OBJ) -o covid19branching $(LGSLFLAGS)

all : SIRoffspring PolioPakEnvBranching covid19branching

MODELPATH=models

data.o : data.cpp data.h
	$(CC) $(CFLAGS) -c data.cpp

likelihood.o : likelihood.h likelihood.cpp parameter.h
	$(CC) $(CFLAGS) -c likelihood.cpp

MCMC.o : MCMC.h MCMC.cpp pfilter.h model.h parameter.h trajectory.h data.h likelihood.h particle.h
	$(CC) $(CFLAGS) -c MCMC.cpp

model.o : model.h model.cpp trajectory.h parameter.h
	$(CC) $(CFLAGS) -c model.cpp

parameter.o : parameter.h parameter.cpp
	$(CC) $(CFLAGS) -c parameter.cpp

particle.o : particle.h particle.cpp trajectory.h
	$(CC) $(CFLAGS) -c particle.cpp

pfilter.o : pfilter.h pfilter.cpp model.h parameter.h trajectory.h data.h likelihood.h particle.h
	$(CC) $(CFLAGS) -c pfilter.cpp

trajectory.o : trajectory.h trajectory.cpp
	$(CC) $(CFLAGS) -c trajectory.cpp

main.o : main.cpp MCMC.h pfilter.h trajectory.h model.h parameter.h data.h likelihood.h particle.h
	$(CC) $(CFLAGS) -c main.cpp

$(SIROFFSPRINGOBJ) : model.h parameter.h trajectory.h $(MODELPATH)/SIR_offspring_distribution.cpp
	$(CC) $(CFLAGS) -c $(MODELPATH)/SIR_offspring_distribution.cpp

$(BRANCHINGEBOLALIBERIAOBJ) : model.h parameter.h trajectory.h $(MODELPATH)/branching_process_EbolaLiberia.cpp
	$(CC) $(CFLAGS) -c $(MODELPATH)/branching_process_EbolaLiberia.cpp

$(BRANCHINGPOLIOPAKENVOBJ) : model.h parameter.h trajectory.h $(MODELPATH)/branching_process_PakEnv.cpp
	$(CC) $(CFLAGS) -c $(MODELPATH)/branching_process_PakEnv.cpp

$(BRANCHINGCOVID19OBJ) : model.h parameter.h trajectory.h $(MODELPATH)/branching_process_covid19.cpp
	$(CC) $(CFLAGS) -c $(MODELPATH)/branching_process_covid19.cpp

clean :
	rm *.o
	find . -maxdepth 1 -type f -name '*~' -delete

tar : tar cfv *.h *.cpp models/SIR_offspring_distribution.cpp
