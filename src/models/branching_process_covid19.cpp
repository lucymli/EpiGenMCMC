//
//  branching_process_covid19.cpp
//  EpiGenMCMC
//
//  Created by Lucy Li on Feb 28, 2020
//

#include <stdio.h>
#include "../model.h"


void Model::set_custom_prob(double alpha, double scale) {
    custom_prob.resize(0);
    double w;
    for (int i=0; i<5000; ++i) {
        w = gsl_ran_gamma(model_rng, alpha, scale);
        custom_prob.push_back(w);
    }
}

void Model::simulate(std::vector<double> & model_params, std::vector<std::string> & param_names, Trajectory * traj, int start_dt, int end_dt, double step_size, int total_dt, gsl_rng * rng) {
    if ((traj->get_state(0)+traj->get_state(1)) < 1.0) {
        return;
    }
    // /* For slightly faster implementation, call parameters by index
    int R0_array_size = 7;
    double R0_array[R0_array_size];
    std::copy(model_params.begin(), model_params.begin()+R0_array_size, R0_array);
    int R0_change_times[R0_array_size-1];
    std::copy(model_params.begin()+8, model_params.begin()+8+R0_array_size-1, R0_change_times);
    double cv = model_params[7];
    double alpha = model_params[27];
    double scale = model_params[28];
    double Re = R0_array[0];
    double recoveries=0.0;
    double new_infections=0.0;
    double durI = 0.0;
    double curr_coal_rate = 0.0;
    double k;
    if (custom_prob.size()<1000) {
        model_rng = rng;
        set_custom_prob(alpha, scale);
    }
    int ran_num = (int)gsl_ran_flat(rng, 0, custom_prob.size());
    if (start_dt < step_size) {  // Set initial number of infected
        traj->resize_recoveries(total_dt+1);
        int init_inf = (int)round(model_params[29]);
        traj->set_state(init_inf, 0);
        for (int i=0; i!=init_inf; ++i) {
            durI = custom_prob[ran_num];
            ++ran_num;
            if (ran_num == custom_prob.size()) ran_num=0;
            durI = (int)(durI/step_size);
            traj->add_recovery_time(durI);
        }
    }
    // */
    double num_infected = traj->get_state(0);
    for (int t=start_dt; t<end_dt; ++t) {
        //
        // Transitions
        //
        Re = R0_array[0];
        for (int R0_i=0; R0_i<(R0_array_size-1); R0_i++) {
            if (t >= R0_change_times[R0_i]) {
                Re *= R0_array[R0_i+1];
            }
        }
        k = Re/(cv-1);
        // Recoveries: I --> R
        recoveries = traj->num_recover_at(t-start_dt);
         if (recoveries > 5000) { // Assume that if more than 3000 people recover within a day, then the epidemic is too large.
             num_infected = 0.0;
             recoveries = 0.0;
         }
        if (recoveries > 0) {
//          traj->set_traj(0, recoveries, t-start_dt);
            new_infections = gsl_ran_negative_binomial(rng, k/(k+Re), k*recoveries);
            if (new_infections > 0) {
                for (int i=0; i!=new_infections; ++i) {
                    durI = custom_prob[ran_num];
                    ++ran_num;
                    if (ran_num == custom_prob.size()) ran_num=0;
                    durI = (int)(durI/step_size);
                    traj->add_recovery_time(durI+t-start_dt);
                }
            }
            //            traj->set_state(traj->get_state(0)+new_infections-recoveries, 0);
            num_infected += new_infections - recoveries;
        }
        // Record 1/N for coalescent rate calculation
        if (num_infected > 0.0) {
            traj->set_traj(0, new_infections, t-start_dt);
            traj->set_traj(1, num_infected, t-start_dt);
            traj->set_traj(2, Re/(alpha*scale)*(1.0+1.0/k), t-start_dt);
        }
        else {
            traj->set_traj(0, 0.0, t-start_dt);
            traj->set_traj(1, 0.0, t-start_dt);
            traj->set_traj(2, 0.0, t-start_dt);
            break;
        }
    }
    traj->set_state(num_infected, 0);
    traj->delete_recoveries_before(end_dt-start_dt);
}

