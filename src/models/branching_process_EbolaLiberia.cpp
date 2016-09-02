//
//  branching_process_EbolaLiberia.cpp
//  EpiGenMCMC
//
//  Created by Lucy Li on 24/06/2016.
//  Copyright (c) 2016 Lucy Li, Imperial College London. All rights reserved.
//

#include <stdio.h>
#include "../model.h"



void Model::set_custom_prob(double alpha, double scale) {
    custom_prob.resize(0);
    bool stop = false;
    double w = 0.0;
    while ((!stop)&(custom_prob.size()<400)) {
        w = gsl_cdf_gamma_P(custom_prob.size(), alpha, scale);
        custom_prob.push_back(w);
        stop = w > 0.99999;
    }
    for (int i=1; i!=custom_prob.size(); ++i) custom_prob[i-1] = custom_prob[i]-custom_prob[i-1];
    custom_prob.pop_back();
}


void Model::simulate(std::vector<double> & model_params, std::vector<std::string> & param_names, Trajectory * traj, int start_dt, int end_dt, double step_size, int total_dt, gsl_rng * rng) {
    if ((traj->get_state(0)+traj->get_state(1)) < 1.0) {
        return;
    }
    double R0_0=1.0;
    double R0_1=1.0;
    double R0_2=1.0;
    double R0_3=1.0;
    double R0_4=1.0;
    double R0_5=1.0;
    double R0_T0=0.0;
    double R0_T1=100000.0;
    double R0_T2=100000.0;
    double R0_T3=100000.0;
    double R0_T4=100000.0;
    double k=1.0;
    double alpha=1.0;
    double scale=1.0;
    // /* For slightly faster implementation, call parameters by index
    for (int i=0; i!=param_names.size(); ++i) {
        if (param_names[i]=="R0_0") R0_0 = model_params[i];
        if (param_names[i]=="R0_1") R0_1 = model_params[i];
        if (param_names[i]=="R0_2") R0_2 = model_params[i];
        if (param_names[i]=="R0_3") R0_0 = model_params[i];
        if (param_names[i]=="R0_4") R0_0 = model_params[i];
        if (param_names[i]=="R0_5") R0_0 = model_params[i];
        if (param_names[i]=="R0_T0") R0_T0 = model_params[i];
        if (param_names[i]=="R0_T1") R0_T1 = model_params[i];
        if (param_names[i]=="R0_T2") R0_T2 = model_params[i];
        if (param_names[i]=="R0_T3") R0_T3 = model_params[i];
        if (param_names[i]=="R0_T4") R0_T4 = model_params[i];
        if (param_names[i]=="k") k = model_params[i];
        if (param_names[i]=="alpha") alpha = model_params[i];
        if (param_names[i]=="scale") scale = model_params[i];
    }
    double R0_now = 0.0;
    double recoveries=0.0;
    double new_infections=0.0;
    double durI = 0.0;
    double ran_unif_num1, ran_unif_num2;
    if (custom_prob.size()==0) set_custom_prob(alpha, scale);
    if (start_dt < step_size) {  // Set initial number of infected
        traj->resize_recoveries(total_dt);
        int init_inf = (int)round(model_params[14]);
        traj->set_state(init_inf, 0);
        for (int i=0; i!=init_inf; ++i) {
//            durI = gsl_ran_gamma(rng, alpha, scale);
            ran_unif_num1 = gsl_ran_flat(rng, 0.0000001, 1);
            ran_unif_num2 = gsl_ran_flat(rng, 0.0000001, 1);
            durI = -log(ran_unif_num1)*scale-log(ran_unif_num2)*scale;
            durI = (int)(durI/365.0/step_size);
            traj->add_recovery_time(durI);
        }
    }
    double num_infected = traj->get_state(0);
    for (int t=start_dt; t<end_dt; ++t) {
        //
        // Transitions
        //
        // Recoveries: I --> R
        recoveries = traj->num_recover_at(t-start_dt);
        if (recoveries > 100000.0) { // If the epidemic is too large, set num_infected to 0, so that likelihood is 0.
            num_infected = 0.0;
            traj->set_traj(0.0, t-start_dt);
        }
        else if (recoveries > 0) {
            traj->set_traj(recoveries, t-start_dt);
            if (t*step_size<(R0_T0)) R0_now = R0_0;
            else if (t*step_size<(R0_T1)) R0_now = R0_1;
            else if (t*step_size<(R0_T2)) R0_now = R0_2;
            else if (t*step_size<(R0_T3)) R0_now = R0_3;
            else if (t*step_size<(R0_T4)) R0_now = R0_4;
            else {
                R0_now = R0_5;
            }
            new_infections = gsl_ran_negative_binomial(rng, k/(k+R0_now), k*recoveries);
            if (new_infections > 1000) {
                std::vector <unsigned int> a (106, 0);
                gsl_ran_multinomial(rng, 106, (unsigned int)new_infections, &custom_prob[0], &a[0]);
                for (int i=0; i!=a.size(); ++i) {
//                    durI = gsl_ran_gamma(rng, alpha, scale);

//                    ran_unif_num1 = gsl_rng_uniform_pos(rng);
//                    ran_unif_num2 = gsl_rng_uniform_pos(rng);
//                    durI = -log(ran_unif_num1)*scale-log(ran_unif_num2)*scale;
//                    durI = (int)(durI/365.0/step_size);
                    
//                    durI = durI_vec[gsl_rng_uniform_int(rng, 1000)];
//                    traj->add_recovery_time(durI+t-start_dt);
                    traj->add_recovery_time(i+t-start_dt, a[i]);
                }
            }
            else if (new_infections > 0) {
                for (int i=0; i!=new_infections; ++i) {
                    ran_unif_num1 = gsl_rng_uniform_pos(rng);
                    ran_unif_num2 = gsl_rng_uniform_pos(rng);
                    durI = -log(ran_unif_num1)*scale-log(ran_unif_num2)*scale;
                    durI = (int)(durI/365.0/step_size);
                    traj->add_recovery_time(durI+t-start_dt);
                }
            }
            num_infected += new_infections - recoveries;
        }
        double curr_coal_rate = 1.0/num_infected;
        // Record 1/N for coalescent rate calculation
        if (num_infected > 0.0) {
            traj->set_traj2(curr_coal_rate*R0_now/(alpha*scale)*(1.0+1.0/k), t-start_dt);
        }
        else {
            traj->set_traj2(0.0, t-start_dt);
            break;
        }
    }
    traj->set_state(num_infected, 0);
    traj->delete_recoveries_before(end_dt-start_dt);
}