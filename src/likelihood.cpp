//
//  likelihood.cpp
//  EpiGenMCMC
//
//  Created by Lucy Li on 13/04/2016.
//  Copyright (c) 2016 Lucy Li, Imperial College London. All rights reserved.
//


#include "likelihood.h"

Likelihood::Likelihood() {}

double Likelihood::binomial_lik(double reporting_rate, double process, int size, int start, int end, int num_groups, bool return_log) {
    // Assuming no data is observed
    double loglik=0.0;
    for (int group_id=0; group_id!=num_groups; ++group_id) {
        double total_incidence = process;//0.0;
        //        for (int i=start; i!=end; ++i) {
        //            total_incidence += *(process+group_id*size+i-start);
        //        }
        loglik += log(1.0-reporting_rate) * (double)total_incidence;//dbinom(0.0, total_incidence, params[5], 1);
    }
    if (return_log) return(loglik);
    return(exp(loglik));
}

double Likelihood::binomial_lik(double reporting_rate, double process, std::vector<double>::iterator data, int size, int start, int end, int shift, int num_groups, bool return_log) {
    double loglik=0.0;
    for (int group_id=0; group_id!=num_groups; ++group_id) {
        double total_incidence = process;//0.0;
        double total_data = 0.0;
        for (int i=start; i!=end; ++i) {
            //            total_incidence += *(process+group_id*size+i-start);
//            double curr_inc = *(data+group_id*size+i-shift);
            double curr_inc = *(data+group_id*(size-shift)+i);
            total_data += curr_inc;
            if (total_data < 0) {
                if (curr_inc >= 0) {
                    total_data = curr_inc;
                }
            }
            
        }
        if (total_data < 0) {
            if (return_log) return(0.0);
            return(1.0);
        }
//        if (total_data < 1.0) {
//            if (return_log) return(0.0);
//            return(1.0);
//        }
//        double total_diff = abs(total_incidence-total_data/reporting_rate);
//        if (total_diff > (total_data/reporting_rate)) total_diff = total_data/reporting_rate;
//        if (return_log) return (log(1.0-total_diff/(total_data/reporting_rate)));
//        else return (1.0-total_diff/(total_data/reporting_rate));
        if (total_data > total_incidence) {
            if (return_log) return(-std::numeric_limits<double>::max());
            return(0.0);
        }
        if (total_data == 0.0) {
            loglik += log(1.0-reporting_rate) * (double)total_incidence;
        }
        else {
            loglik += log(gsl_ran_binomial_pdf(total_data, reporting_rate, total_incidence));
        }
    }
    if (return_log) return(loglik);
    return(exp(loglik));
}

double Likelihood::coalescent_lik(std::vector<double>::iterator sim_prev, std::vector<double>::iterator sim_coal_rate,
                                  std::vector<double>::iterator binomial, std::vector<double>::iterator intervals,
                                  std::vector<double>::iterator indices, int start, int end, int shift, bool return_log) {
    double weight = 0.0;
    for (int deltaT=start; deltaT!=end; ++deltaT) {
        // Loop over each simulation time step
        double NtoNe = *(sim_coal_rate+deltaT-start);
        double prev = *(sim_prev+deltaT-start);
        double coal_rate = NtoNe/prev;
        if (prev<1.0) coal_rate = 0.0;
        int first_index;
        if ((deltaT)==0) first_index = 0;
        else first_index = *(indices+deltaT-1)+1;
        int last_index = *(indices+deltaT);
        for (int event=first_index; event<=last_index; ++event) {
            // Loop over event during a simulation time step (either coalescence or sampling)
            double binom_coeff = *(binomial+event);
            if (binom_coeff > 0) {
                if ((prev <= 1.0) | ((prev*(prev-1.0)/2.0) < binom_coeff)) {
                    if (return_log) return(-std::numeric_limits<double>::max());
                    return (0.0);
                }
                double coal_rate_population = binom_coeff*coal_rate;
                double time_to_next_event = *(intervals+event);
                if (time_to_next_event < 0.0) {
                    // Coalescent ended interval
                    if (coal_rate == 0.0){
                        // If epidemic has died out before the most recent tip
                        if (return_log) return(-std::numeric_limits<double>::max());
                        return (0.0);
                    }
                    double time_to_coal = -*(intervals+event);
                    weight += log(coal_rate_population)-(coal_rate_population*time_to_coal);
                }
                else {
                    // Sampling ended interval, or the end of the simulation time step
                    double time_to_event = *(intervals+event);
                    weight += -(coal_rate_population*time_to_event);
                }
            }
        }
        if (!std::isfinite(weight)) {
            if (return_log) return (-std::numeric_limits<double>::max());
            return (0.0);
        }
    }
    if (return_log) return (weight);
    return (exp(weight));
}
