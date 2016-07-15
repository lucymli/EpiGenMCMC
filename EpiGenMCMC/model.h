//
//  model.h
//  EpiGenMCMC
//
//  Created by Lucy Li on 13/04/2016.
//  Copyright (c) 2016 Lucy Li, Imperial College London. All rights reserved.
//

#ifndef EpiGenMCMC_model_h
#define EpiGenMCMC_model_h

#include "parameter.h"
#include "trajectory.h"


class Model{
    bool use_deterministic;
public:
    Model();
    void set_deterministic(bool);
    void simulate(std::vector<double> &, std::vector<std::string>&, Trajectory *, int, int, double, int, gsl_rng *);
};



#endif
