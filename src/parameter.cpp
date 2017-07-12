//
//  parameter.cpp
//  EpiGenMCMC
//
//  Created by Lucy Li on 13/04/2016.
//  Copyright (c) 2016 Lucy Li, Imperial College London. All rights reserved.
//

#define MATHLIB_STANDALONE

#include <omp.h>
#include "parameter.h"


Parameter::Parameter() {}

Parameter::Parameter(const Parameter & copy_from) {
    parameter_values.clear();
    parameter_values.insert(parameter_values.begin(), copy_from.parameter_values.begin(), copy_from.parameter_values.end());
    parameter_names.clear();
    parameter_names.insert(parameter_names.begin(), copy_from.parameter_names.begin(), copy_from.parameter_names.end());
    estimate.clear();
    estimate.insert(estimate.begin(), copy_from.estimate.begin(), copy_from.estimate.end());
    transform.clear();
    transform.insert(transform.begin(), copy_from.transform.begin(), copy_from.transform.end());
    prior.clear();
    prior.insert(prior.begin(), copy_from.prior.begin(), copy_from.prior.end());
    prior_par_1.clear();
    prior_par_1.insert(prior_par_1.begin(), copy_from.prior_par_1.begin(), copy_from.prior_par_1.end());
    prior_par_2.clear();
    prior_par_2.insert(prior_par_2.begin(), copy_from.prior_par_2.begin(), copy_from.prior_par_2.end());
    prior_par_3.clear();
    prior_par_3.insert(prior_par_3.begin(), copy_from.prior_par_3.begin(), copy_from.prior_par_3.end());
    proposal.clear();
    proposal.insert(proposal.begin(), copy_from.proposal.begin(), copy_from.proposal.end());
    proposal_sd.clear();
    proposal_sd.insert(proposal_sd.begin(), copy_from.proposal_sd.begin(), copy_from.proposal_sd.end());
    proposal_lower.clear();
    proposal_lower.insert(proposal_lower.begin(), copy_from.proposal_lower.begin(), copy_from.proposal_lower.end());
    proposal_upper.clear();
    proposal_upper.insert(proposal_upper.begin(), copy_from.proposal_upper.begin(), copy_from.proposal_upper.end());
    accepted.clear();
    accepted.insert(accepted.begin(), copy_from.accepted.begin(), copy_from.accepted.end());
    rejected.clear();
    rejected.insert(rejected.begin(), copy_from.rejected.begin(), copy_from.rejected.end());
    acceptance_rate.clear();
    acceptance_rate.insert(acceptance_rate.begin(), copy_from.acceptance_rate.begin(), copy_from.acceptance_rate.end());
    params_to_estim.clear();
    params_to_estim.insert(params_to_estim.begin(), copy_from.params_to_estim.begin(), copy_from.params_to_estim.end());
    total_params = copy_from.total_params;
    curr_param_to_estimate = copy_from.curr_param_to_estimate;
    old_param_value = copy_from.old_param_value;
    optimal_acceptance = copy_from.optimal_acceptance;
    lower_acceptance = copy_from.lower_acceptance;
    upper_acceptance = copy_from.upper_acceptance;
    adapt_every = copy_from.adapt_every;
    max_adapt_times = copy_from.max_adapt_times;
    stop_adapting = copy_from.stop_adapting;
}

Parameter::Parameter(std::string filename) {
    curr_param_to_estimate = 0;
    stop_adapting = false;
    std::ifstream file(filename);
    std::string line;
    file >> total_params >> optimal_acceptance >> lower_acceptance >> upper_acceptance;
    file >> adapt_every;
    std::getline(file, line);
    max_adapt_times = std::stod(line);
    parameter_values.resize(total_params);
    parameter_names.resize(total_params);
    estimate.resize(total_params);
    transform.resize(total_params);
    prior.resize(total_params);
    prior_par_1.resize(total_params);
    prior_par_2.resize(total_params);
    prior_par_3.resize(total_params);
    proposal.resize(total_params);
    proposal_sd.resize(total_params);
    proposal_lower.resize(total_params);
    proposal_upper.resize(total_params);
    accepted.resize(total_params);
    rejected.resize(total_params);
    acceptance_rate.resize(total_params);
    int i = 0;
    while (i<total_params) {
        getline(file, line);
        std::istringstream input (line);
        for (int item=0; item!=3; ++item) {
            std::getline(input, line, ' ');
            if (item==0) parameter_values[i] = std::stod(line);
            if (item==1) parameter_names[i] = line;
            if (item==2) {
                std::string str1 ("True");
                bool estimate_or_not = (line.compare(str1)==0);
                if (!estimate_or_not) {
                    str1 = "TRUE";
                    estimate_or_not = (line.compare(str1)==0);
                    if (!estimate_or_not) {
                        str1 = "true";
                        estimate_or_not = (line.compare(str1)==0);
                    }
                }
                estimate[i] = estimate_or_not;
                if (estimate[i]) params_to_estim.push_back(i);
            }
        }
        if (estimate[i]) {
            for (int item=3; item!=12; ++item) {
                getline(input, line, ' ');
                if (item==3) transform[i] = line;
                if (item==4) prior[i] = line;
                if (item==5) prior_par_1[i] = std::stod(line);
                if (item==6) prior_par_2[i] = std::stod(line);
                if (item==7) prior_par_3[i] = std::stod(line);
                if (item==8) proposal[i] = line;
                if (item==9) proposal_sd[i] = std::stod(line);
                if (item==10) proposal_lower[i] = std::stod(line);
                if (item==11) proposal_upper[i] = std::stod(line);
            }
        }
        ++i;
    }
    if (params_to_estim.size() > 0) curr_param_to_estimate = params_to_estim.front();
    file.close();
}

bool Parameter::is_estim(int index) const {
    return(estimate[index]);
}

int Parameter::get_curr_estim() const {
    return (curr_param_to_estimate);
}

void Parameter::set_next_param() {
    bool next_found = false;
    while (!next_found) {
        ++curr_param_to_estimate;
        if (curr_param_to_estimate == estimate.size()) curr_param_to_estimate = 0;
        next_found = estimate[curr_param_to_estimate];
    }
}

int Parameter::get_total_params() const {
    return (total_params);
}

int Parameter::get_total_params_to_estim() const {
    return ((int)params_to_estim.size());
}

int Parameter::get_estim_index(int index) const {
    return (params_to_estim[index]);
}

double Parameter::get(std::string param_name) const {
    int itr=0;
    while (parameter_names[itr]!=param_name) ++itr;
    double x = parameter_values[itr];
    return(x);
}

double Parameter::get(int index) const {
    return(parameter_values[index]);
}

double Parameter::get_lower(int index) const {
    return (proposal_lower[index]);
}

double Parameter::get_upper(int index) const {
    return (proposal_upper[index]);
}

void Parameter::set(int index, double value) {
    parameter_values[index] = value;
}

std::vector <double> Parameter::get_values_vector() const {
    return (parameter_values);
}

std::vector <std::string> Parameter::get_names_vector() const {
    return (parameter_names);
}


std::string Parameter::getname(int index) const {
    return(parameter_names[index]);
}

void Parameter::reset() {
    std::fill(accepted.begin(), accepted.end(), 0.0);
    std::fill(rejected.begin(), rejected.end(), 0.0);
    std::fill(acceptance_rate.begin(), acceptance_rate.end(), 0.0);
}

void Parameter::accept() {
    ++accepted[curr_param_to_estimate];
    double total = accepted[curr_param_to_estimate] + rejected[curr_param_to_estimate];
    acceptance_rate[curr_param_to_estimate] = accepted[curr_param_to_estimate] / total;
    set_next_param();
}

void Parameter::reject() {
    ++rejected[curr_param_to_estimate];
    double total = accepted[curr_param_to_estimate] + rejected[curr_param_to_estimate];
    acceptance_rate[curr_param_to_estimate] = accepted[curr_param_to_estimate] / total;
    parameter_values[curr_param_to_estimate] = old_param_value;
    set_next_param();
}

double Parameter::get_acceptance() const {
    double average_accept = 0.0;
    for (int i=0; i!=acceptance_rate.size(); ++i) {
        average_accept += acceptance_rate[i];
    }
    average_accept /= (double)acceptance_rate.size();
    return (average_accept);
}

void Parameter::stop_adapt() {
    stop_adapting = true;
}

void Parameter::start_adapt() {
    stop_adapting = false;
    reset();
}

void Parameter::adapt() {
    if (!stop_adapting) {
        int itr=total_params;
        while (!estimate[itr]) --itr;
        if ((accepted[itr] + rejected[itr]) >= adapt_every) {
            for (int i=0; i!=estimate.size(); ++i) {
                if (estimate[i]) {
                    double curr_accept = acceptance_rate[i];
                    if ((curr_accept < lower_acceptance) | (curr_accept > upper_acceptance)) {
                        curr_accept = std::max(0.0001, curr_accept);
                        curr_accept = std::min(0.99, curr_accept);
//                        proposal_sd[i] *= qnorm(optimal_acceptance/2.0, 0.0, 1.0, 1, 0) / qnorm(curr_accept/2.0, 0.0, 1.0, 1, 0);
                        proposal_sd[i] *= exp(0.999/2.0*(curr_accept-optimal_acceptance));
                        proposal_sd[i] = std::max(proposal_sd[i], 0.000001);
                    }
                    accepted[i] = 0.0;
                    rejected[i] = 0.0;
                }
            }
            --max_adapt_times;
        }
    }
    if (max_adapt_times == 0) stop_adapt();
}

double Parameter::get_transform(double value, std::string transformation, bool reverse) const {
    double newvalue = value;
    if (transformation=="inverse") {
        newvalue = 1.0/value;
    }
    else if (transformation=="log") {
        if (!reverse) newvalue = log(value);
        else newvalue = exp(value);
    }
    else if (transformation=="logit") {
        if (!reverse) newvalue = log(value/(1-value));
        else newvalue = exp(value) / (1.0 + exp(value));
    }
    return (newvalue);
}

double Parameter::get_lognormal_sd (double MEAN, double SD) const {
    double sigma = pow(log(1.0+SD*SD/MEAN/MEAN), 0.5);
    return (sigma);
}

double Parameter::get_lognormal_mean (double MEAN, double SD) const {
    double zeta = log(MEAN/pow(1.0+(SD*SD/MEAN/MEAN), 0.5));
    return (zeta);
}

double Parameter::propose(gsl_rng * rng) {
    if (!stop_adapting) adapt();
    std::string trans = transform[curr_param_to_estimate];
    old_param_value = parameter_values[curr_param_to_estimate];
    std::string proposal_distribution = proposal[curr_param_to_estimate];
    double old = get_transform(old_param_value, trans, false);
    double curr;
    double lower = proposal_lower[curr_param_to_estimate];
    double upper = proposal_upper[curr_param_to_estimate];
    double sd = proposal_sd[curr_param_to_estimate];
    double correction_factor, tempA, tempB, tempC, tempD;
    if (proposal_distribution == "lognormal") {
        double zeta = get_lognormal_mean(old, sd);
        double sigma = get_lognormal_mean(old, sd);
        curr = exp(zeta+sigma*gsl_ran_gaussian(rng, 1.0));
        while ((curr < lower) | (curr > upper)) {
            curr = exp(zeta+sigma*gsl_ran_gaussian(rng, 1.0));
        }
        tempA = gsl_cdf_ugaussian_P((log(upper)-log(old))/sd);
        tempB = gsl_cdf_ugaussian_P((log(upper)-log(curr))/sd);
        correction_factor = (curr*tempA) / (old*tempB);
    } else {
        curr = old+gsl_ran_gaussian(rng, sd);
        while ((curr < lower) | (curr > upper)) {
            curr = old+gsl_ran_gaussian(rng, sd);
        }
        tempA = gsl_cdf_ugaussian_P((upper-old)/sd);
        tempB = gsl_cdf_ugaussian_P((lower-old)/sd);
        tempC = gsl_cdf_ugaussian_P((upper-curr)/sd);
        tempD = gsl_cdf_ugaussian_P((lower-curr)/sd);
        correction_factor = (tempA - tempB) / (tempC - tempD);
    }
    curr = get_transform(curr, trans, true);
    parameter_values[curr_param_to_estimate] = curr;
    return (log(correction_factor));
}

double Parameter::get_prior(double original_value, int index) const {
    double value = Parameter::get_transform(original_value, transform[index], false);
    if (prior[index]=="unif") {
        if ((value >= prior_par_1[index]) & (value <= prior_par_2[index])) {
            return (1.0/(prior_par_2[index]-prior_par_1[index]));
        }
        return (0.0);
    }
    else if (prior[index]=="gamma") {
        return (gsl_ran_gamma_pdf(value, prior_par_1[index], prior_par_2[index]));
    }
    else if (prior[index]=="beta") {
        return(gsl_ran_beta_pdf(value, prior_par_1[index], prior_par_2[index]));
    }
    else if (prior[index]=="lognormal") {
        double m = get_lognormal_mean(prior_par_1[index], prior_par_2[index]);
        double s = get_lognormal_sd(prior_par_1[index], prior_par_2[index]);
        return (gsl_ran_lognormal_pdf(value, m, s));
    }
    return (0.0);
}


double Parameter::get_prior_ratio() const {
    int i = curr_param_to_estimate;
    return(get_prior(parameter_values[i], i) - get_prior(old_param_value, i));
}

double Parameter::get_prior_all() const {
    double prior_value = 0.0;
    for (int i=0; i!=total_params; ++i) {
        if (estimate[i]) {
            prior_value += get_prior(parameter_values[i], i);
        }
    }
    if (!std::isfinite(prior_value)) prior_value =  -0.1*std::numeric_limits<double>::max();
    return(prior_value);
}


bool Parameter::param_exists(std::string param_name) const {
    bool found = false;
    int itr = total_params;
    while (!found & (itr > 0)) {
        --itr;
        found = parameter_names[itr] == param_name;
    }
    return(found);
}

void Parameter::transform_param(int index, bool reverse) {
    parameter_values[index] = get_transform(parameter_values[index], transform[index], reverse);
}




MCMCoptions::MCMCoptions(){}

MCMCoptions::MCMCoptions(const MCMCoptions& old_copy) {
    particles = old_copy.particles;
    iterations = old_copy.iterations;
    log_every = old_copy.log_every;
    pfilter_every = old_copy.pfilter_every;
    pfilter_threshold = old_copy.pfilter_threshold;
    which_likelihood = old_copy.which_likelihood;
    num_trees = old_copy.num_trees;
    log_filename = old_copy.log_filename;
    traj_filename = old_copy.traj_filename;
    model = old_copy.model;
    total_dt = old_copy.total_dt;
    sim_dt = old_copy.sim_dt;
    num_groups = old_copy.num_groups;
    seed = old_copy.seed;
    verbose = old_copy.verbose;
    save_traj = old_copy.save_traj;
    use_lhs = old_copy.use_lhs;
    lhs_divides = old_copy.lhs_divides;
    lhs_iterations = old_copy.lhs_iterations;
    num_threads = old_copy.num_threads;
    rng = old_copy.rng;
}

MCMCoptions::MCMCoptions(std::string filename) {
    std::ifstream file(filename);
    std::string line;
    verbose = true;
    save_traj = true;
    use_lhs = false;
    num_threads = omp_get_num_threads();
    num_groups = 1;
    heat_factor = 1.0;
    heat_length = 0;
    cool_rate = 1.0;
    while (file >> line) {
        if (line == "particles") file >> particles;
        else if (line == "iterations") file >> iterations;
        else if (line == "log_every") file >> log_every;
        else if (line == "pfilter_every") file >> pfilter_every;
        else if (line == "pfilter_threshold") file >> pfilter_threshold;
        else if (line == "which_likelihood") file >> which_likelihood;
        else if (line == "num_trees") file >> num_trees;
        else if (line == "log_filename") file >> log_filename;
        else if (line == "traj_filename") file >> traj_filename;
        else if (line == "model") file >> model;
        else if (line == "verbose") {
            int num;
            file >> num;
            if (num == 0) verbose = false;
        }
        else if (line == "save_traj") {
            int num;
            file >> num;
            if (num == 0) save_traj = false;
        }
        else if (line == "use_lhs") {
            int num;
            file >> num;
            if (num > 0) use_lhs = true;
        }
        else if (line == "lhs_divides") file >> lhs_divides;
        else if (line == "lhs_iterations") file >> lhs_iterations;
        else if (line == "num_threads") {
            file >> num_threads;
        }
        else if (line=="heat_factor") file >> heat_factor;
        else if (line=="heat_length") file >> heat_length;
        else if (line=="cool_rate") file >> cool_rate;
        else file >> line;
    }
    file.close();
}

