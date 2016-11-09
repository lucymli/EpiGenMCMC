//
//  SEIR_Ebola.cpp
//  EpiGenMCMC
//
//  Created by Lucy Li on 13/04/2016.
//  Copyright (c) 2016 Lucy Li, Imperial College London. All rights reserved.
//

#include "../model.h"

void Model::simulate(std::vector<double> & model_params, std::vector<std::string> & param_names, Trajectory * traj, int start_dt, int end_dt, double step_size, int total_dt,  gsl_rng * rng) {
    if ((traj->get_state(1)+traj->get_state(2)) < 1.0) {
        return;
    }
    double Beta = model_params[0];
    double k=model_params[1];
    double rateE2I=model_params[2];
    double rateI2R=model_params[3];
    double R0_reduce=model_params[6];
    int change_T=model_params[7]/step_size;
    double Rt = Beta/rateI2R*traj->get_state(0);
    double total_infectious=0.0;
    double new_infections=0.0;
    double divisions = 10.0;
    double p1 = rateE2I*step_size/divisions;
    double p2 = rateI2R*step_size/divisions;
    double S2E = 0.0;
    double E2I = 0.0;
    double I2R = 0.0;
    double currS;
    double latent;
    double infectious;
    double Tg = 1.0/rateE2I + 1.0/rateI2R;
    // */
    for (int t=start_dt; t<end_dt; ++t) {
        double sub_t_I2R = 0.0;
        for (int sub_t=0; sub_t<(int) divisions; ++sub_t) {
            //
            // Transitions
            //
            // Recoveries: I --> R
            currS = traj->get_state(0);
            latent = traj->get_state(1);
            infectious = traj->get_state(2);
            if (infectious > 0) {
                if (p2 >= 1.0) I2R = infectious;
                else if (use_deterministic) {
                    I2R = infectious*p2;
                }
                else {
                    I2R = gsl_ran_binomial(rng, p2, infectious); // Recoveries
                }
            }
            // Becoming infectious: E --> I
            if (latent > 0) {
                if (p1 >= 1.0) E2I = latent;
                else if (use_deterministic) {
                    E2I = latent * p1;
                }
                else {
                    E2I = gsl_ran_binomial(rng, p1, latent); // Becoming infectious
                }
            }
            traj->set_state(infectious-I2R+E2I, 2); // Infectious
            if (use_deterministic) {
                S2E = Beta * currS * step_size * infectious;
            } else {
                // Infections: S --> E
                if (I2R > 0.0) {
                    total_infectious += I2R;
                    sub_t_I2R += I2R;
                    //                traj->set_traj(I2R, t-start_dt); // People are sampled, i.e. appear in the time-series at the time of recovery
                    if (currS > 0) {
                        // Current reproductive number
                        Rt = Beta / rateI2R * currS;
                        if (change_T>=t) Rt *= R0_reduce;
                        // Draw from the negative binomial distribution (gamma-poisson mixture) to determine
                        // number of secondary infections
                        S2E = gsl_ran_negative_binomial(rng, k/(k+Rt), k*I2R);  // New infections
                        if (S2E > 0) {
                            S2E = std::min(currS, S2E);
                            new_infections += S2E;
                            traj->set_state(currS-S2E, 0); // Susceptible
                        }
                    }
                }
            }
            traj->set_state(latent+S2E-E2I, 1); // Latent
            S2E=0.0;
            E2I=0.0;
            I2R=0.0;
        }
        traj->set_traj(1, sub_t_I2R, t-start_dt);
        double N = traj->get_state(1)+traj->get_state(2);
        // Record 1/N for coalescent rate calculation
        if (N > 0.0) {
            traj->set_traj(2, N, t-start_dt);
            traj->set_traj(3, Rt/Tg*(1.0+1.0/k), t-start_dt);
        }
        else {
            traj->set_traj(2, 0.0, t-start_dt);
            traj->set_traj(3, 0.0, t-start_dt);
            break;
        }
    }
}
