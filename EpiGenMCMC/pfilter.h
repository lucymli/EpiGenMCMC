//
//  pfilter.h
//  EpiGenMCMC
//
//  Created by Lucy Li on 08/04/2016.
//  Copyright (c) 2016 Lucy Li, Imperial College London. All rights reserved.
//

#ifndef EpiGenMCMC_pfilter_h
#define EpiGenMCMC_pfilter_h

#include <omp.h>
#include "model.h"
#include "data.h"
#include "likelihood.h"
#include "particle.h"



namespace EpiGenPfilter {
    double pfilter(Model & sim_model, Parameter & model_params, MCMCoptions & options, Particle &particles, Trajectory & output_traj, TimeSeriesData &epi_data, TreeData &tree_data, MultiTreeData &multitree_data);
}

#endif
