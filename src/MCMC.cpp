//
//  MCMC.cpp
//  EpiGenMCMC
//
//  Created by Lucy Li on 06/05/2016.
//  Copyright (c) 2016 Lucy Li, Imperial College London. All rights reserved.
//

#define MATHLIB_STANDALONE

#include <iostream>
#include "MCMC.h"



namespace EpiGenMCMC_MCMC {
    void initialise_logfile(std::string filename, Parameter model_params, MCMCoptions options, unsigned int seed1, unsigned int seed2) {
        std::chrono::time_point<std::chrono::system_clock> curr_time = std::chrono::system_clock::now();
        std::time_t start_time = std::chrono::system_clock::to_time_t(curr_time);
        std::cout << filename << std::endl;
        std::ofstream file(filename);
        file << "# Seeds: " << seed1 << ", " << seed2 << std::endl;
        file << "# Start time: " << std::ctime(&start_time) << std::endl;
        file << "# Starting parameters:" << std::endl;
        for (int i=0; i!=model_params.get_total_params(); ++i) {
            file << "#   " << model_params.getname(i) << " ";
            file << model_params.get(i) << std::endl;
        }
        file << "# MCMC details:" << std::endl;
        file << "#   " << options.iterations << " iterations" << std::endl;
        file << "#   " << options.particles << " particles" << std::endl;
        file << "#   Log parameters every " << options.log_every << " iterations" << std::endl;
        file << "#   Filter at most every " << options.pfilter_every << " time steps" << std::endl;
        file << "#   Filter if ESS < " << options.pfilter_threshold*100.0 << "% of particles" << std::endl;
        file << "#   Likelihood based on ";
        if (options.which_likelihood == 0) {
            file << "both epidemiologic and genealogical data";
            if (options.num_trees > 1) file << "(" << options.num_trees << " trees)";
        }
        else if (options.which_likelihood == 1) file << "epidemiologic data";
        else file << "genealogical data";
        file << std::endl;
        file << "state\tposterior\tlikelihood\tprior";
        for (int i=0; i!=model_params.get_total_params(); ++i) {
            if (model_params.is_estim(i)) file << "\t" << model_params.getname(i);
        }
        //        file << std::endl << "0";
        //        for (int i=0; i!=model_params.get_total_params(); ++i) {
        //            if (model_params.is_estim(i)) file << "\t" << model_params.get(i);
        //        }
        file << std::endl;
        file.close();
    }
    void mcmc_log(std::string filename, Parameter model_params, int iteration, double loglik, double logprior) {
        std::ofstream file(filename, std::ios::app);
        if (!std::isfinite(logprior)) logprior = -0.1*std::numeric_limits<double>::max();
        file << iteration << "\t" << loglik+logprior << "\t" << loglik << "\t" << logprior;
        for (int i=0; i!=model_params.get_total_params(); ++i) {
            if (model_params.is_estim(i)) file << "\t" << model_params.get(i);
        }
        file << std::endl;
        file.close();
    }
    void lhs(Parameter model_params, MCMCoptions options, std::vector <double> & lower_bounds, std::vector <double> & upper_bounds, int divides, int index_ML, bool verbose) {
        // Latin Hypercube Sampling
        double min_val;
        double max_val;
        double interval_size;
        double value;
        if (lower_bounds.size() < 1) { // First iteration of sampling
            for (int i=0; i!=model_params.get_total_params_to_estim(); ++i) {
                int j = model_params.get_estim_index(i);
                min_val = model_params.get_lower(j);
                max_val = model_params.get_upper(j);
                interval_size = (max_val - min_val) / (double) divides;
                for (int interval_i=0; interval_i!=divides; ++interval_i) {
                    value = min_val+(double)interval_i*interval_size;
                    lower_bounds.push_back(value);
                    upper_bounds.push_back(value+interval_size);
                }
            }
            if (verbose) std::cout << "Drawing a Latin Hypercube Sample for " << lower_bounds.size()/divides << " parameters ";
            if (verbose) std::cout << "each divided into " << divides << " partitions." << std::endl;
        }
        else {
            for (int i=0; i!=model_params.get_total_params_to_estim(); ++i) {
                int j = model_params.get_estim_index(i);
                min_val = lower_bounds[i*divides+index_ML];
                max_val = upper_bounds[i*divides+index_ML];
                interval_size = (max_val - min_val) / (double) divides;
                for (int interval_i=0; interval_i!=divides; ++interval_i) {
                    value = min_val+(double)interval_i*interval_size;
                    lower_bounds[i*divides+interval_i] = std::max(model_params.get_lower(j), value);
                    upper_bounds[i*divides+interval_i] = std::min(model_params.get_upper(j), value+interval_size);
                }
            }
        }
        std::vector <int> indices (divides);
        for (int i=0; i!=divides; ++i) indices[i] = i;
        int param_index;
        std::vector <double> lowers_matrix (lower_bounds.size());
        std::vector <double> uppers_matrix = lowers_matrix;
        // Construct grid of parameter value combinations to calculate the likelihood
        for (int i=0; i!=model_params.get_total_params_to_estim(); ++i) {
            gsl_ran_shuffle(options.rng[0], &indices[0], indices.size(), sizeof(int));
            for (int j=0; j!=indices.size(); ++j) {
                param_index = i*divides+indices[j];
                lowers_matrix[i*divides+j] = lower_bounds[param_index];
                uppers_matrix[i*divides+j] = upper_bounds[param_index];
            }
        }
        lowers_matrix.swap(lower_bounds);
        uppers_matrix.swap(upper_bounds);
    }
    int mcmc (Parameter & model_params, MCMCoptions & options, Trajectory & traj, TimeSeriesData & epi_data, TreeData & tree_data, MultiTreeData & multitree_data) {
        if (options.verbose) std::cout << "Begin MCMC" << std::endl;
        if (options.verbose) std::cout << "Using " << options.num_threads << " cores." << std::endl;
        unsigned int seed1 = (int)time(0);
        unsigned int seed2 = (int)time(0)+10;
        options.seed = seed1;
        const gsl_rng_type * rng_t = gsl_rng_rand;
        gsl_rng **rng;
        rng = (gsl_rng **) malloc(options.num_threads * sizeof(gsl_rng *));
        for (int i = 0; i != options.num_threads; ++i) {
            rng[i] = gsl_rng_alloc(rng_t);
            gsl_rng_set(rng[i], seed1+i*options.particles);
        }
        options.rng = rng;
        initialise_logfile(options.log_filename, model_params, options, seed1, seed2);
        if (options.verbose) std::cout << "MCMC log file initialised" << std::endl;
        double loglik;
        double newloglik;
        double prior_ratio;
        double correction_factor;
        double accept_with_prob;
        double ran_num;
        bool accept;
        Trajectory temptraj = traj;
        Particle particles(options.particles, traj);
        Model sim_model;
        if (options.use_lhs) {
            bool save_traj_in_mcmc_chain = options.save_traj;
            if (options.save_traj) options.save_traj = false;
            std::vector <double> lower_bounds;
            std::vector <double> upper_bounds;
            int index_max = 0;
            for (int i=0; i!=options.lhs_iterations; ++i) {
                if (options.verbose) std::cout << "Iteration 1" << std::endl;
                lhs(model_params, options, lower_bounds, upper_bounds, options.lhs_divides, index_max, options.verbose);
                double loglik = -std::numeric_limits<double>::max();
                double max_posterior = loglik;
                for (int j=0; j!=options.lhs_divides; ++j) {
                    if (options.verbose) std::cout << j << ") ";
                    // Set parameter values to the next combination
                    for (int param_i=0; param_i!=model_params.get_total_params_to_estim(); ++param_i) {
                        int curr_param_index = model_params.get_estim_index(param_i);
                        int grid_index = param_i*options.lhs_divides+j;
                        double curr_param_value = (lower_bounds[grid_index]+upper_bounds[grid_index])/2.0;
                        model_params.set(curr_param_index, curr_param_value);
                        model_params.transform_param(curr_param_index, true);
                        std::cout << model_params.get(curr_param_index) << " ";
                    }
                    // Likelihood of current parameter combination
//                    loglik = -std::numeric_limits<double>::max();
                    loglik = EpiGenPfilter::pfilter(sim_model, model_params, options, particles, temptraj, epi_data, tree_data, multitree_data);
                    if (options.verbose) std::cout << "loglik: " << loglik << std::endl;
                    if (j==0) {
                        max_posterior = loglik + model_params.get_prior_all();
                        index_max = 0;
                    }
                    else {
                        if (loglik + model_params.get_prior_all() > max_posterior) {
                            // If a new maximum posterior density is found, update index_ML
                            // to the index of the current parameter combination
                            index_max = j;
                            max_posterior = loglik + model_params.get_prior_all();
                        }
                    }
                }
            }
            // Set initial parameter values to those that produced the maximum posterior density
            for (int param_i=0; param_i!=model_params.get_total_params_to_estim(); ++param_i) {
                int curr_param_index = model_params.get_estim_index(param_i);
                int grid_index = param_i*options.lhs_divides+index_max;
                double curr_param_value = (lower_bounds[grid_index]+upper_bounds[grid_index])/2.0;
                model_params.set(curr_param_index, curr_param_value);
                model_params.transform_param(curr_param_index, true);
            }
            if (save_traj_in_mcmc_chain) options.save_traj = true;
        }
        std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now(), time2;
        std::chrono::duration<double> elapsed_seconds;
        if (options.verbose) time1 = std::chrono::system_clock::now();
        loglik = EpiGenPfilter::pfilter(sim_model, model_params, options, particles, traj, epi_data, tree_data, multitree_data);
        if (options.verbose) {
            time2 = std::chrono::system_clock::now();
            elapsed_seconds = time2-time1;
            time1 = time2;
            std::cout << "State TimeInterval Acceptance" << std::endl;
            std::cout << "0 " << elapsed_seconds.count() << " 0.0"<< std::endl;
        }
        if (options.save_traj) traj.print_to_file(0, options.traj_filename, options.pfilter_every, true);
        mcmc_log(options.log_filename, model_params, 0, loglik, model_params.get_prior_all());
        if (options.heat_factor > 1.0) model_params.stop_adapt();
        for (int iter = 1; iter != options.iterations; ++iter) {
            correction_factor = model_params.propose(options.rng[0]);
            particles.reset_weights();
            newloglik = EpiGenPfilter::pfilter(sim_model, model_params, options, particles, temptraj, epi_data, tree_data, multitree_data);
            prior_ratio = model_params.get_prior_ratio();
            accept_with_prob = (newloglik - loglik) + prior_ratio + correction_factor;
            if (options.heat_factor > 1.0) {
                if (iter <= (options.heat_length*model_params.get_total_params_to_estim())) accept_with_prob += log(options.heat_factor);
                else {
                    if ((iter % model_params.get_total_params_to_estim())==0) {
                        options.heat_factor -= std::min(options.heat_factor-1.0, options.cool_rate);
                        accept_with_prob += log(options.heat_factor);
                        if (options.heat_factor <= 1.0) {
                            model_params.start_adapt();
                        }
                    }
                }
            }
            accept = accept_with_prob >= 0.0;
            if (!accept) {
                ran_num = gsl_rng_uniform(options.rng[0]);
                if (log(ran_num) < accept_with_prob) accept = true;
            }
            if (accept) {
                loglik = newloglik;
                traj = temptraj;
                model_params.accept();
            }
            else {
                model_params.reject();
            }
            if ((iter)%options.log_every == 0) {
                if (options.save_traj) traj.print_to_file(iter, options.traj_filename, options.pfilter_every, true);
                mcmc_log(options.log_filename, model_params, iter, loglik, model_params.get_prior_all());
                if (options.verbose) {
                    time2 = std::chrono::system_clock::now();
                    elapsed_seconds = time2-time1;
                    time1 = time2;
                    std::cout << iter << " " << elapsed_seconds.count() << " ";
                    std::cout << model_params.get_acceptance() << std::endl;
                }
            }
        }
        particles.clear();
        return (1);
    }
}
