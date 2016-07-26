//
//  main.cpp
//  EpiGenMCMC
//
//  Created by Lucy Li on 07/04/2016.
//  Copyright (c) 2016 Lucy Li, Imperial College London. All rights reserved.
//


#include <iostream>
#include "MCMC.h"

int main(int argc, const char * argv[]) {
    Parameter model_params(argv[1]);
    std::cout << "Read in a total of " << model_params.get_total_params() << " parameters." << std::endl;
    MCMCoptions options(argv[2]);
    std::cout << "Finished reading in MCMC options" << std::endl;
    int total_data_dt = 0;
    double sim_dt = 1.0;
    int num_groups = 1;
    TimeSeriesData epi_data;
    TreeData tree_data;
    MultiTreeData multitree_data;
    if (options.which_likelihood < 2) {
        epi_data = *new TimeSeriesData(argv[4]);
        total_data_dt = epi_data.get_T();
        num_groups = epi_data.get_num_group();
        options.num_groups = num_groups;
        std::cout << "Finished reading in epi data. Total time steps: " << total_data_dt << std::endl;
        sim_dt = epi_data.get_dt();
    }
    if (options.which_likelihood!=1) {
        int i = 4;
        if (options.which_likelihood==0) i = 5;
        if (options.num_trees > 1) {
            multitree_data = *new MultiTreeData(argv[i]);
            total_data_dt = multitree_data.get_T();
            sim_dt = multitree_data.get_dt();
        }
        else {
            tree_data = *new TreeData(argv[i]);
            total_data_dt = tree_data.get_T();
            sim_dt = tree_data.get_dt();
        }
        std::cout << "Finished reading in tree data. Total time steps: " << total_data_dt << std::endl;
    }
    options.total_dt = total_data_dt;
    options.sim_dt = sim_dt;
    Trajectory init_traj(total_data_dt, num_groups);
    if (options.save_traj) init_traj.initialise_file(options.traj_filename, options.pfilter_every);
    init_traj.initialise_states(argv[3]);
    int mcmc_out = EpiGenMCMC_MCMC::mcmc (model_params, options, init_traj, epi_data, tree_data, multitree_data);
    std::cout << "MCMC complete" << std::endl;
    return mcmc_out;
}
