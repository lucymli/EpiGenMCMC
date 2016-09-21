//
//  main.cpp
//  EpiGenMCMC
//
//  Created by Lucy Li on 07/04/2016.
//  Copyright (c) 2016 Lucy Li, Imperial College London. All rights reserved.
//

#include <omp.h>
#include <iostream>
#include "model.h"

// Argument order:
//   1. Parameters
//   2. Replicates
//   3. Total dt
//   4. Size of simulation time step
//   5. Sum every time steps
//   6. Number of groups
//   7. Seed
//   8. Number of threads
//   9. Filename of inital states
//   10. Filename of output
int main(int argc, const char * argv[]) {
    // Read in parameter values etc.
    Parameter model_params(argv[1]);
    std::vector <double> values = model_params.get_values_vector();
    std::vector <std::string> param_names = model_params.get_names_vector();
    std::cout << "Read in a total of " << model_params.get_total_params() << " parameters." << std::endl;
    int replicates = std::stoi(argv[2]);
    int total_dt = std::stoi(argv[3]);
    double dt_size = std::stod(argv[4]);
    int sum_every = std::stoi(argv[5]);
    int num_groups = std::stoi(argv[6]);
    int seed_num = std::stoi(argv[7]);
    int num_threads = std::stoi(argv[8]);
    std::string traj_input = argv[9];
    std::string traj_output = argv[10];
    // Initialise trajectories
    Trajectory init_traj(total_dt, num_groups);
    init_traj.resize(total_dt, num_groups);
    init_traj.initialise_states(traj_input);
    init_traj.initialise_file(traj_output, sum_every);
    Model sim_model;
    const gsl_rng_type * rng_t = gsl_rng_rand;
    gsl_rng **rng;
    rng = (gsl_rng **) malloc(num_threads * sizeof(gsl_rng *));
    for (int i = 0; i != num_threads; ++i) {
        rng[i] = gsl_rng_alloc(rng_t);
        gsl_rng_set(rng[i], seed_num+i);
    }
    // A single simulation
    if (replicates==1) {
        sim_model.simulate(values, param_names, &init_traj, 0, total_dt, dt_size, seed_num+1, rng[0]);
        init_traj.print_to_file(traj_output, sum_every, true);
    } else {
        // Parallelised simulations
        std::vector <Trajectory *> trajectories;
        for (int i=0; i<replicates; ++i) {
            trajectories.push_back(new Trajectory(init_traj));
        }
#pragma omp parallel for schedule(static,1)
        for (int tn=0; tn<num_threads; ++tn) {
            Model sim_model_tn;
            for (int i=tn; i<replicates; ++i) {
                sim_model_tn.simulate(values, param_names, trajectories[i], 0, total_dt, dt_size, seed_num+i, rng[tn]);
                std::cout << i << "\t";
            }
        }
        for (int i=0; i<replicates; ++i) {
            trajectories[i]->print_to_file(i, traj_output, sum_every, true);
        }
    }
    std::cout << "Simulation complete" << std::endl;
    return 0;
}
