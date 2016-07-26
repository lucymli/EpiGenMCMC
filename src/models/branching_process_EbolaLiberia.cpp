//
//  branching_process_EbolaLiberia.cpp
//  EpiGenMCMC
//
//  Created by Lucy Li on 24/06/2016.
//  Copyright (c) 2016 Lucy Li, Imperial College London. All rights reserved.
//

#include <stdio.h>
#include "../model.h"


void Model::simulate(std::vector<double> & model_params, std::vector<std::string> & param_names, Trajectory * traj, int start_dt, int end_dt, double step_size, int total_dt, gsl_rng * rng) {
    if ((traj->get_state(0)+traj->get_state(1)) < 1.0) {
        return;
    }
    
    // /* For slightly faster implementation, call parameters by index
    double R0_0 = model_params[0];
    double R0_1 = model_params[1];
    double R0_T0 = model_params[2];
    double R0_now = 0.0;
    double k = model_params[3];
    double rateE2I = model_params[4];
    double rateI2R = model_params[5];
    double Re = 0.0;
    double recoveries=0.0;
    double new_infections=0.0;
    double durI = 0.0;
    if (start_dt < step_size) {  // Set initial number of infected
        traj->resize_recoveries(total_dt);
        int init_inf = (int)round(model_params[6]);
        traj->set_state(init_inf, 0);
        for (int i=0; i!=init_inf; ++i) {
            durI = gsl_ran_exponential(rng, 1.0/rateE2I)+gsl_ran_exponential(rng, 1.0/rateI2R);
            durI = start_dt + (int)round(durI/step_size);
            traj->add_recovery_time(durI);
        }
    }
    // */
    for (int t=start_dt; t<end_dt; ++t) {
        //
        // Transitions
        //
        // Recoveries: I --> R
        recoveries = traj->num_recover_at(t-start_dt);
        if (recoveries > 0) {
            traj->set_traj(recoveries, t-start_dt);
            if (t*step_size<(R0_T0)) R0_now = R0_0;
            else R0_now = R0_1;
            new_infections = gsl_ran_negative_binomial(rng, k/(k+R0_now), k*recoveries);
            if (new_infections > 0) {
                for (int i=0; i!=new_infections; ++i) {
                    durI = gsl_ran_exponential(rng, 1.0/rateE2I)+gsl_ran_exponential(rng, 1.0/rateI2R);
                    durI = t + (int)round(durI/step_size);
                    traj->add_recovery_time(durI-start_dt);
                }
            }
            traj->set_state(traj->get_state(0)+new_infections-recoveries, 0);
        }
        double num_infected = traj->get_state(0);
        double curr_coal_rate = 1.0/num_infected;
        // Record 1/N for coalescent rate calculation
        if (num_infected > 0.0) {
            traj->set_traj2(curr_coal_rate*Re/((1.0/rateI2R)+(1.0/rateE2I))*(1.0+1.0/k), t-start_dt);
        }
        else {
            traj->set_traj2(0.0, t-start_dt);
            break;
        }
    }
    traj->delete_recoveries_before(end_dt-start_dt);
}