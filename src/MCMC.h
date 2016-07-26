//
//  MCMC.h
//  EpiGenMCMC
//
//  Created by Lucy Li on 08/04/2016.
//  Copyright (c) 2016 Lucy Li, Imperial College London. All rights reserved.
//

#ifndef EpiGenMCMC_MCMC_h
#define EpiGenMCMC_MCMC_h

#include <chrono>
#include <ctime>
//#include "parameter.h"
#include "pfilter.h"


namespace EpiGenMCMC_MCMC {
    void initialise_logfile(std::string, Parameter, MCMCoptions, unsigned int, unsigned int);
    void mcmc_log(std::string, Parameter, int, double, double);
    template <class Type> void reorder(std::vector<Type>&);
    void lhs(Parameter, std::vector <double> &, std::vector <double> &, int, int, bool);
    int mcmc (Parameter&, MCMCoptions&, Trajectory&, TimeSeriesData&, TreeData&, MultiTreeData&);
}


#endif
