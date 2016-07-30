//
//  SEIR_2age_open_offspring.cpp
//  EpiGenMCMC
//
//  Created by Lucy Li on 18/05/2016.
//  Copyright (c) 2016 Lucy Li, Imperial College London. All rights reserved.
//

#include "../model.h"

void Model::simulate(std::vector<double> & model_params, std::vector<std::string> & param_names, Trajectory * traj, int start_dt, int end_dt, double step_size, int total_dt, gsl_rng * rng) {
    if ((traj->get_state(1)+traj->get_state(2)+traj->get_state(6)+traj->get_state(5)) < 1.0) {
        return;
    }
    // states
    //  0-3: S1, E1, I1, R1
    //  4-7: S2, E2, I2, R2
    //  8-21: H1_1, H2_1, H3_1, H4_1, H5_1, H6_1, H7_1, H1_2...
    //
    // parameters
    //  0-1: beta1, beta2
    //  2-3: rateE2I, ratetransitions[2*age_i+
    //  4-5: seasonality amplitude, shift
    //  6: k
    //  7: infPerCase
    //  8: nat_hist_rate
    //  9: birth_rate
    //  10-11: death_rate1, 2
    //  12-13: N1, 2
    //  14-15: propS1, 2
    //  16: routine vaccination effectiveness
    //  17-18: sia_vaccination effectiveness, frequency
    
    /* If the order of parameters is not known beforehand
     double beta1 = -1;
     double beta2 = -1;
     double rateE2I = -1;
     double ratetransitions[2*age_i+ = -1;
     double seasonality = -1;
     double T_seasonality = -1;
     double k = -1;
     double infPerCase = -1;
     double nat_hist_rate = -1;
     double birth_rate = -1;
     double death_rate_1 = -1;
     double death_rate_2 = -1;
     double N1 = -1;
     double N2 = -1;
     double propS1 = -1;
     double propS2 = -1;
     double routine_vacc = -1;
     double sia_vacc = -1;
     double sia_every = -1;
     for (int i=0; i!=param_names.size(); ++i) {
        if (param_names[i] == "beta1") beta1 = model_params[i];
        if (param_names[i] == "beta2") beta2 = model_params[i];
        if (param_names[i] == "rateE2I") rateE2I = model_params[i];
        if (param_names[i] == "ratetransitions[2*age_i+") ratetransitions[2*age_i+ = model_params[i];
        if (param_names[i] == "seasonality") seasonality = model_params[i];
        if (param_names[i] == "T_seasonality") T_seasonality = model_params[i];
        if (param_names[i] == "k") k = model_params[i];
        if (param_names[i] == "infPerCase") infPerCase = model_params[i];
        if (param_names[i] == "nat_hist_rate") nat_hist_rate = model_params[i];
        if (param_names[i] == "birth_rate") birth_rate = model_params[i];
        if (param_names[i] == "death_rate_1") death_rate_1 = model_params[i];
        if (param_names[i] == "death_rate_2") death_rate_2 = model_params[i];
        if (param_names[i] == "N1") N1 = model_params[i];
        if (param_names[i] == "N2") N2 = model_params[i];
        if (param_names[i] == "propS1") propS1 = model_params[i];
        if (param_names[i] == "propS2") propS2 = model_params[i];
        if (param_names[i] == "routine_vacc") routine_vacc = model_params[i];
        if (param_names[i] == "sia_vacc") sia_vacc = model_params[i];
        if (param_names[i] == "sia_every") sia_every = model_params[i];
     }
     */
    // /* For slightly faster implementation, call parameters by index
//    double beta1 = model_params[0];
//    double beta2 = model_params[1];
//    double rateE2I = model_params[2];
//    double ratetransitions[2*age_i+ = model_params[3];
//    double seasonality = model_params[4];
//    double T_seasonality = model_params[5];
//    double k = model_params[6];
//    double nat_hist_rate = model_params[8];
//    double birth_rate = model_params[9];
//    double death_rate_1 = model_params[10];
//    double death_rate_2 = model_params[11];
////    double N1 = model_params[12];
////    double N2 = model_params[13];
////    double propS1 = model_params[14];
////    double propS2 = model_params[15];
//    double routine_vacc = model_params[16];
//    double sia_vacc = model_params[17];
//    double sia_every = model_params[18];
    
    double beta1 = model_params.get(0);
    double beta2 = model_params.get(1);
    double rateE2I = model_params.get(2);
    double rateI2R = model_params.get(3);
    double seasonality = model_params.get(4);
    double T_seasonality = model_params.get(5);
    double k = model_params.get(6);
    double nat_hist_rate = model_params.get(8);
    double birth_rate = model_params.get(9);
    double death_rate_1 = model_params.get(10);
    double death_rate_2 = model_params.get(11);
    //    double N1 = model_params.get(12);
    //    double N2 = model_params.get(13);
    //    double propS1 = model_params.get(14);
    //    double propS2 = model_params.get(15);
    double routine_vacc = model_params.get(16);
    double sia_vacc = model_params.get(17);
    double sia_every = model_params.get(18);
    
    int num_ages = 2;
    
//    double prop_in_1 = N1/(N1+N2);
//    double rate_mature = death_rate_2 * (1 - prop_in_1) / prop_in_1;
    std::vector <double> matrix;
    std::vector <double> matrix_prop(num_ages, 0.0);
    double curr_susc1 = traj->get_state(0);
    double curr_susc2 = traj->get_state(4);
    matrix.push_back(beta1*curr_susc1);
    matrix.push_back(beta1*beta2*curr_susc1);
    matrix.push_back(beta1*beta2*curr_susc2);
    matrix.push_back(beta1*beta2*curr_susc2);
    std::vector <double> Rts;
    Rts.push_back(1.0/ rateI2R * (matrix[0] + matrix[2]));
    Rts.push_back(1.0/ rateI2R* (matrix[1] + matrix[3]));
    
    
    
    double total_infectious=0.0;
    double new_infections=0.0;
    double p = rateI2R*step_size;
    // */
    std::vector<double> transitions (num_ages*3, 0.0);
    std::vector<double> currS (num_ages, 0.0);
    
    double num_infected = 1.0;
    
    int t = start_dt-1;
    while ((num_infected > 0.0) & (++t < end_dt)) {
        num_infected = 0.0;
        double total_pop_size = 0.0;
        for (int i=0; i!=4; ++i) {
            double state_size = traj->get_state(i);
            if (state_size*death_rate_1*step_size > 2) {
                std::binomial_distribution<int> distribution ((int)state_size, death_rate_1*step_size);
                double deaths1 = distribution(rng);
                traj->set_state(std::max(0.0, round(state_size-deaths1)), i);
            }
            total_pop_size += state_size;
        }
        for (int i=4; i!=8; ++i) {
            double state_size = traj->get_state(i);
            if (state_size*death_rate_2*step_size > 2) {
                std::binomial_distribution<int> distribution (state_size, death_rate_2*step_size);
                traj->set_state(std::max(0.0, round(state_size-distribution(rng))), i);
            }
            total_pop_size += state_size;
        }
//        std::binomial_distribution <int> births_distribution(total_pop_size, birth_rate);
        double num_births = round(total_pop_size * birth_rate * step_size);
        double num_routine_vacc = round(num_births * routine_vacc);
        traj->set_state(traj->get_state(0)+num_births-num_routine_vacc, 0);
        traj->set_state(traj->get_state(3)+num_routine_vacc, 3);
      
        
        //
        // Transitions
        //
        // Recoveries: I --> R
        for (int age_i=0; age_i!=num_ages; ++age_i) {
            double infected = traj->get_state(age_i*4+2);
            if (infected > 0) {
                if (p >= 1.0) transitions[2*age_i+age_i] = infected;
                else {
                    std::binomial_distribution<int> distribution(infected, p);
                    transitions[2*age_i+age_i] = std::min(infected, (double)distribution(rng));
                }
            }
        }
        for (int age_i=0; age_i!=num_ages; ++age_i) {
            currS[age_i] = traj->get_state(age_i*4+0);
        }
        // Infections: S --> E
        for (int age_i=0; age_i!=num_ages; ++age_i) {
            if (transitions[2*age_i+age_i] > 0.0) {
                total_infectious += transitions[2*age_i+age_i];
                traj->set_state(traj->get_state(age_i*4+3)-transitions[2*age_i+age_i], age_i*4+3); //Recovered
                // Current reproductive number
                // Draw from the negative binomial distribution (gamma-poisson mixture) to determine
                // number of secondary infections
                double alpha = k*transitions[2*age_i+age_i];
                double scale = Rts[age_i]/k*(1.0+seasonality*cos(2*M_PI*(double)t*step_size+T_seasonality));
                double gamma_number;
                std::gamma_distribution<double> gam (alpha, scale);
                gamma_number = gam(rng);
                std::poisson_distribution<int> pois (gamma_number);
                double inf = pois(rng);
                if (inf > 0) {
                    // Probability of infection in each age group by an infectious individual in age group age_i
                    for (int i=0; i!=num_ages; ++i) {
                        int index = i*num_ages+age_i;
                        matrix_prop[i] = matrix[index];
                    }
                    double A = std::accumulate(matrix_prop.begin(), matrix_prop.end(), 0.0);
                    for (int i=0; i!=matrix_prop.size(); ++i) {
                        matrix_prop[i] /= A;
                        transitions[i] += round(matrix_prop[i] * inf);
                    }
                }
            }
        }
        // Becoming infectious: E --> I
        for (int age_i=0; age_i!=num_ages; ++age_i) {
            double num_latent = traj->get_state(age_i*4+1);
            if (num_latent > 0) {
                std::binomial_distribution <int> distribution (num_latent, rateE2I*step_size);
                transitions[age_i+age_i] = distribution(rng);
                
            }
            if (transitions[age_i] > 0) {
                transitions[age_i] = std::min(currS[age_i], transitions[age_i]);
                new_infections += transitions[age_i];
                traj->set_state(currS[age_i]-transitions[age_i], age_i*4+0); // Susceptible
            }
            traj->set_state(num_latent+transitions[age_i]-transitions[age_i+age_i], age_i*4+1); // Exposed
            double previous_infectious = traj->get_state(age_i*4+2);
            double new_infectious = previous_infectious+transitions[age_i+age_i]-transitions[2*age_i+age_i];
            traj->set_state(new_infectious, age_i*4+2); //Infectious
            num_infected += new_infectious;
        }
        // Natural history compartments
        for (int age_i=0; age_i!=num_ages; ++age_i) {
            double curr_size;
            double transition1=0.0, transition2=0.0;
            for (int i=0; i!=7; ++i) {
                curr_size = traj->get_state(8+age_i*4+i);
                if (curr_size > 0) {
                    std::binomial_distribution <int> distribution (curr_size, nat_hist_rate*step_size);
                    transition1 = round(distribution(rng));
                }
                if (i==0) traj->set_state(curr_size+transitions[age_i]-transition1, 8+age_i*4+i);
                else if (i==6) {
                    traj->set_state(curr_size+transition2-transition1, 8+age_i*4+i);
                    traj->set_traj(transition1, t-start_dt, age_i);
                }
                else traj->set_state(curr_size+transition2-transition1, 8+age_i*4+i);
                transition2 = transition1;
            }
        }
        // SIA immunisation
        if ((t%(int)sia_every)==0) {
            double susceptible_children = traj->get_state(0);
            if (susceptible_children > 0) {
                double immunised = susceptible_children*sia_vacc;
                if (immunised < 100) {
                    std::binomial_distribution<int>distribution(susceptible_children, sia_vacc);
                    immunised = round(distribution(rng));
                }
                if (immunised > 0) {
                    traj->set_state(traj->get_state(0)-immunised, 0); // S
                    traj->set_state(traj->get_state(3)+immunised, 3); // R
                }
            }
        }
        // Record 1/N for coalescent rate calculation
        if (num_infected > 0.0) {
            traj->set_traj2(1.0/num_infected, t-start_dt);
        }
        else {
            traj->set_traj2(0.0, t-start_dt);
        }
        std::fill(transitions.begin(), transitions.end(), 0.0);
    }
    // Empirical reproductive number
    double Rt = 0.0;
    if (new_infections > 0.0) {
        Rt = new_infections/total_infectious;
    }
    // Calculate coalescent rate as 1/I(t) / Tg * Rt * (1 + 1/k)
    for (int t=start_dt; t<end_dt; ++t) {
        double coal_rate = traj->get_traj2(t-start_dt)*Rt*(rateE2I+rateI2R)*(1.0+1.0/k);
        traj->set_traj2(coal_rate, t-start_dt);
    }
    std::vector <double> ().swap(matrix);
    std::vector <double> ().swap(matrix_prop);
    std::vector <double> ().swap(Rts);
    std::vector <double> ().swap(transitions);
    std::vector <double> ().swap(currS);
}
