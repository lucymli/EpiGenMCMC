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
    /* If the order of parameters is not known beforehand
     double rateI2R=1.0;
     double k=1.0;
     double Rt=-1;
     double Beta=-1;
     for (int i=0; i!=param_names.size(); ++i) {
     if (param_names[i] == "rateI2R") rateI2R = model_params[i];
     if (param_names[i] == "k") k = model_params[i];
     if (param_names[i] == "R0") {
     Rt = model_params[i] / traj->get_init_state(0) * traj->get_state(0);
     }
     if (param_names[i] == "beta") {
     Beta = model_params[i];
     }
     }
     if (Beta < 0) Beta = Rt * rateI2R / traj->get_state(0);
     if (Rt < 0) Rt = Beta / rateI2R * traj->get_state(0);
     */
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
    double rate_incub = model_params[16];
    double Re = 0.0;
    double S2E = 0.0;
    double E2I = 0.0;
    double I2R = 0.0;
    double total_infectious=0.0;
    double new_infections=0.0;
    double prob_recovery = rateI2R*step_size;
    double prob_infectious = rateE2I*step_size;
    double prob_incub = rate_incub*step_size;
    double infected = 0.0;
    double latent = 0.0;
    int num_incub_states = 7;
    std::vector < double> incub_trans (num_incub_states, 0.0);
    if (start_dt < step_size) {  // Set initial number of infected
        traj->set_state(round(model_params[10]), 1);
    }
    // */
    for (int t=start_dt; t<end_dt; ++t) {
        if (t*step_size<(R0_T0)) R0_now = R0_0;
        else if (t*step_size < R0_T1) R0_now = R0_1;
        else if (t*step_size < R0_T2) R0_now = R0_2;
        else if (t*step_size < R0_T3) R0_now = R0_3;
        else if (t*step_size < R0_T4) R0_now = R0_4;
        else R0_now = R0_5;
        double time_of_year = ((double)t*step_size)+phase;
        Re = R0_now * (1.0 + amplitude * cos(2.0*M_PI*time_of_year));
        //
        // Transitions
        //
        // Recoveries: I --> R
        if (traj->get_state(1) > 0) {
            infected = traj->get_state(1);
            if (prob_recovery >= 1.0) I2R = infected;
            else if (use_deterministic) {
                I2R = infected*prob_recovery;
            }
            else {
                if (false){//infected*prob_recovery > 50.0) {
//                    std::normal_distribution<double>distribution(infected*prob_recovery, sqrt(infected*prob_recovery*(1.0-prob_recovery)));
                    I2R = round(infected*prob_recovery);//round(distribution (rng));
                }
                else {
//                    std::binomial_distribution<int> distribution(infected, prob_recovery);
//                    I2R = distribution(rng);
                    I2R = gsl_ran_binomial(rng, prob_recovery, infected);
                }
            }
        }
        // Becoming infectious: E --> I
        if (traj->get_state(0) > 0) {
            latent = traj-> get_state(0);
            if (false){//latent * prob_infectious > 50.0) {
//                std::normal_distribution<double>distribution(latent*prob_infectious, sqrt(latent*prob_infectious*(1.0-prob_infectious)));
//                E2I = round(distribution (rng));
                E2I = round(latent*prob_infectious);
            }
            else {
//                std::binomial_distribution<int> distribution (latent, prob_infectious);
//                E2I = distribution (rng);
                E2I = gsl_ran_binomial(rng, prob_infectious, latent);
            }
        }
        // New infections: S --> I
        if (I2R > 0.0) {
            total_infectious += I2R;
            traj->set_traj(I2R, t-start_dt);
            if (use_deterministic) {
                S2E = Re * I2R;
            }
            else {
                // Draw from the negative binomial distribution (gamma-poisson mixture) to determine
                // number of secondary infections
                if (false){//Re*I2R > 50.0) {
                    S2E = round(Re*I2R);
                }
                else {
//                    double alpha = k*I2R;
//                    double scale = Re / k;
//                    double gamma_number;
//                    std::gamma_distribution<double> gam (alpha, scale);
//                    gamma_number = gam(rng);
//                    std::poisson_distribution<int> pois (gamma_number);
//                    S2E = pois(rng);
                    S2E=gsl_ran_negative_binomial(rng, k/(k+Re), k*I2R);
                }
            }
            if (S2E > 0) {
                new_infections += S2E;
//                traj->set_traj(S2E, t-start_dt);
            }
            traj->set_state(traj->get_state(0)+S2E-E2I, 0); //Infectious
            traj->set_state(traj->get_state(1)+E2I-I2R, 1); //Infectious
        }
        // Incubation period
//        double incub_num;
//        for (int incub_i = 0; incub_i!=num_incub_states; ++incub_i) {
//            incub_num = traj->get_state(incub_i+2);
//            std::binomial_distribution <int> distribution(incub_num, prob_incub);
//            incub_trans[incub_i] = distribution(rng);
//            if (incub_i == 0) {
//                traj->set_state(incub_num-incub_trans[incub_i]+S2E, incub_i+2);
//            }
//            else {
//                traj->set_state(incub_num-incub_trans[incub_i]+incub_trans[incub_i-1], incub_i+2);
//            }
//            traj->set_traj(incub_trans[num_incub_states-1], t-start_dt);
//        }
        double num_infected = traj->get_state(1);
        double curr_coal_rate = 1.0/num_infected;
        // Record 1/N for coalescent rate calculation
        if (num_infected > 0.0) {
            traj->set_traj2(curr_coal_rate*Re/((1.0/rateI2R)+(1.0/rateE2I))*(1.0+1.0/k), t-start_dt);
        }
        else {
            traj->set_traj2(0.0, t-start_dt);
            break;
        }
        S2E=0.0;
        E2I=0.0;
        I2R=0.0;
    }

//    // Calculate coalescent rate as 1/I(t) / Tg * Rt * (1 + 1/k)
//    for (int t=start_dt; t<end_dt; ++t) {
//        double coal_rate = traj->get_traj2(t-start_dt)*Re/((1.0/rateI2R)+(1.0/rateE2I))*(1.0+1.0/k);
//        traj->set_traj2(coal_rate, t-start_dt);
//    }
}