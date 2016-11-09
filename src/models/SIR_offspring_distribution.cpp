//
//  SIR_offspring_distribution.cpp
//  EpiGenMCMC
//
//  Created by Lucy Li on 13/04/2016.
//  Copyright (c) 2016 Lucy Li, Imperial College London. All rights reserved.
//

#include "../model.h"

void Model::simulate(std::vector<double> & model_params, std::vector<std::string> & param_names, Trajectory * traj, int start_dt, int end_dt, double step_size, int total_dt,  gsl_rng * rng) {
    if (traj->get_state(1) < 1.0) {
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
    double rateI2R=model_params[2];
    double k=model_params[1];
//    double Rt=model_params[0];
//    double Beta = Rt * rateI2R / traj->get_state(0);
    double Beta = model_params[0];
    double Rt = Beta/rateI2R*traj->get_state(0);
    double S2I = 0.0;
    double I2R = 0.0;
    double total_infectious=0.0;
    double new_infections=0.0;
    double p = rateI2R*step_size;
    // */
    for (int t=start_dt; t<end_dt; ++t) {
        //
        // Transitions
        //
        // Recoveries: I --> R
        if (traj->get_state(1) > 0) {
            double infected = traj->get_state(1);
            if (p >= 1.0) I2R = infected;
            else if (use_deterministic) {
                I2R = infected*p;
            }
            else {
                I2R = gsl_ran_binomial(rng, p, infected);
            }
        }
        // Infections: S --> I
        if (I2R > 0.0) {
            total_infectious += I2R;
            traj->set_state(traj->get_state(2)+I2R, 2); //Recovered
            traj->set_traj(1, I2R, t-start_dt); // People are sampled, i.e. appear in the time-series at the time of recovery
            double currS = traj->get_state(0);
            if (currS > 0) {
                if (use_deterministic) {
                    S2I = Beta*currS*traj->get_state(1)*step_size;
                }
                else {
                    // Current reproductive number
                    Rt = Beta / rateI2R * currS;
                    // Draw from the negative binomial distribution (gamma-poisson mixture) to determine
                    // number of secondary infections
                    S2I = gsl_ran_negative_binomial(rng, k/(k+Rt), k*I2R);
                }
                if (S2I > 0) {
                    S2I = std::min(currS, S2I);
                    new_infections += S2I;
                    traj->set_state(currS-S2I, 0); // Susceptible
                }
            }
            traj->set_state(traj->get_state(1)+S2I-I2R, 1); //Infectious
        }
        double num_infected = traj->get_state(1);
        // Record 1/N for coalescent rate calculation
        if (num_infected > 0.0) {
            traj->set_traj(2, num_infected, t-start_dt);
            traj->set_traj(3, S2I/I2R/Tg*(1.0+1.0/k), t-start_dt);
        }
        else {
            traj->set_traj2(0.0, t-start_dt);
            break;
        }
        S2I=0.0;
        I2R=0.0;
    }
}
