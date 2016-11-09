//
//  parameter.h
//  EpiGenMCMC
//
//  Created by Lucy Li on 08/04/2016.
//  Copyright (c) 2016 Lucy Li, Imperial College London. All rights reserved.
//

#ifndef EpiGenMCMC_parameter_h
#define EpiGenMCMC_parameter_h


#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <limits>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>


class Parameter {
    std::vector <double> parameter_values;
    std::vector <std::string> parameter_names;
    std::vector <bool> estimate;
    std::vector <std::string> transform;
    std::vector <std::string> prior;
    std::vector <double> prior_par_1;
    std::vector <double> prior_par_2;
    std::vector <double> prior_par_3;
    std::vector <std::string> proposal;
    std::vector <double> proposal_sd;
    std::vector <double> proposal_lower;
    std::vector <double> proposal_upper;
    std::vector <double> accepted;
    std::vector <double> rejected;
    std::vector <double> acceptance_rate;
    std::vector <int> params_to_estim;
    int total_params;
    int curr_param_to_estimate;
    double old_param_value;
    double optimal_acceptance;
    double lower_acceptance;
    double upper_acceptance;
    int adapt_every;
    int max_adapt_times;
    bool stop_adapting;
    double get_prior(double, int) const;
    double get_transform(double, std::string, bool) const;
public:
    Parameter();
    Parameter(const Parameter &);
    Parameter (std::string);
    bool is_estim(int) const;
    int get_curr_estim() const;
    void set_next_param();
    int get_total_params() const;
    int get_total_params_to_estim() const;
    int get_estim_index(int) const;
    double get(std::string) const;
    double get(int) const;
    double get_lower(int) const;
    double get_upper(int) const;
    void set(int, double);
    std::vector<double> get_values_vector() const;
    std::vector <std::string> get_names_vector() const;
    std::string getname(int) const;
    void reset();
    void accept();
    void reject();
    double get_acceptance() const;
    void stop_adapt();
    void start_adapt();
    void adapt();
    double propose(gsl_rng *);
    double get_prior_ratio() const;
    double get_prior_all() const;
    bool param_exists(std::string) const;
    void transform_param(int, bool);
    double get_lognormal_sd (double, double) const;
    double get_lognormal_mean (double, double) const;
};



struct MCMCoptions {
    int particles;
    int iterations;
    int log_every;
    int pfilter_every;
    double pfilter_threshold;
    int which_likelihood;
    int num_trees;
    std::string log_filename;
    std::string traj_filename;
    std::string model;
    int total_dt;
    double sim_dt;
    int num_groups;
    int seed;
    bool verbose;
    bool save_traj;
    bool use_lhs;
    int lhs_divides;
    int lhs_iterations;
    int num_threads;
    double heat_factor;
    int heat_length;
    double cool_rate;
    gsl_rng ** rng;
    MCMCoptions();
    MCMCoptions(const MCMCoptions&);
    MCMCoptions(std::string);
    
};




#endif
