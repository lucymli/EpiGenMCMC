//
//  branching_process_PakEnv_FATA_KP.cpp
//  EpiGenMCMC
//
//  Created by Lucy Li on 14/06/2016.
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
    double R0_1 = model_params[1]*R0_0;
    double R0_2 = model_params[2]*R0_1;
    double R0_3 = model_params[3]*R0_2;
    double R0_4 = model_params[4]*R0_3;
    double R0_5 = model_params[5]*R0_4;
    double R0_T0 = model_params[6];
    double R0_T1 = model_params[7];
    double R0_T2 = model_params[8];
    double R0_T3 = model_params[9];
    double R0_T4 = model_params[10];
    double R0_now = 0.0;
    double amplitude = model_params[11];
    double phase = model_params[12];
    double k = model_params[13];
    double rateE2I = model_params[14];
    double rateI2R = model_params[15];
    double Re = 0.0;
    double recoveries=0.0;
    double new_infections=0.0;
    double durI = 0.0;
    double ran_unif_num1 = 0.1;
    double ran_unif_num2 = 0.1;
    if (start_dt < step_size) {  // Set initial number of infected
        traj->resize_recoveries(total_dt);
        int init_inf = (int)round(model_params[17]);
        traj->set_state(init_inf, 0);
        for (int i=0; i!=init_inf; ++i) {
            ran_unif_num1 = gsl_ran_flat(rng, 0.0000001, 1);
            ran_unif_num2 = gsl_ran_flat(rng, 0.0000001, 1);
            durI = -log(ran_unif_num1)/rateE2I-log(ran_unif_num2)/rateI2R;
            //            durI = gsl_ran_exponential(rng, 1.0/rateE2I)+gsl_ran_exponential(rng, 1.0/rateI2R);
            durI = (int)(durI/step_size);
            traj->add_recovery_time(durI);
        }
    }
    double a = rateE2I*rateI2R/(rateI2R-rateE2I)/200.0;
    // */
    double num_infected = traj->get_state(0);
    for (int t=start_dt; t<end_dt; ++t) {
        //
        // Transitions
        //
        // Recoveries: I --> R
        recoveries = traj->num_recover_at(t-start_dt);
        if (recoveries > 3000) { // Assume that if more than 3000 people recover within a tenth of a day, then the epidemic is too large.
            num_infected = 0.0;
        }
        else if (recoveries > 0) {
            traj->set_traj(recoveries, t-start_dt);
            if (t<(R0_T0)) R0_now = R0_0;
            else if (t < R0_T1) R0_now = R0_1;
            else if (t < R0_T2) R0_now = R0_2;
            else if (t < R0_T3) R0_now = R0_3;
            else if (t < R0_T4) R0_now = R0_4;
            else R0_now = R0_5;
            double time_of_year = ((double)t*step_size)+phase;
            Re = R0_now * (1.0 + amplitude * cos(2.0*M_PI*time_of_year));
            new_infections = gsl_ran_negative_binomial(rng, k/(k+Re), k*recoveries);
            if (new_infections > 100) {
                double infection_counter = 0;
                for (int i=0; i!=total_dt; ++i) {
                    double recover_after = std::max(1.0, round(new_infections*a*(exp(-rateE2I*(i+1.0)*step_size)-exp(-rateI2R*(i+1.0)*step_size))));
                    infection_counter += recover_after;
                    if (infection_counter > new_infections) recover_after -= infection_counter-new_infections;
                    traj->add_recovery_time(i+t-start_dt, (int)recover_after);
                    if (infection_counter >= new_infections) break;
                }
            }
            else if (new_infections > 0) {
                for (int i=0; i!=new_infections; ++i) {
                    ran_unif_num1 = gsl_ran_flat(rng, 0.0000001, 1);
                    ran_unif_num2 = gsl_ran_flat(rng, 0.0000001, 1);
                    durI = -log(ran_unif_num1)/rateE2I - log(ran_unif_num2)/rateI2R;
                    //                    durI = gsl_ran_exponential(rng, 1.0/rateE2I)+gsl_ran_exponential(rng, 1.0/rateI2R);
                    durI = (int)(durI/step_size);
                    traj->add_recovery_time(durI+t-start_dt);
                }
            }
            //            traj->set_state(traj->get_state(0)+new_infections-recoveries, 0);
            num_infected += new_infections - recoveries;
        }
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
    traj->set_state(num_infected, 0);
    traj->delete_recoveries_before(end_dt-start_dt);
}