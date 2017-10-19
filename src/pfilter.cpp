//
//  pfilter.cpp
//  EpiGenMCMC
//
//  Created by Lucy Li on 05/05/2016.
//  Copyright (c) 2016 Lucy Li, Imperial College London. All rights reserved.
//

#include "pfilter.h"
#include <iostream>



namespace EpiGenPfilter {
    double pfilter(Model & sim_model, Parameter & model_params, MCMCoptions & options, Particle &particles, Trajectory & output_traj, TimeSeriesData &epi_data, TreeData &tree_data, MultiTreeData &multitree_data) {
        double loglik = 0.0;
        int num_groups = options.num_groups;
        int num_particles = options.particles;
        int init_seed = options.seed;
        int total_dt = options.total_dt;
        double sim_dt = options.sim_dt;
        int total_steps = ceil((double)total_dt/(double)options.pfilter_every);
        int add_dt = 0;
        double ESS_threshold = options.pfilter_threshold*(double)num_particles;
        Likelihood likelihood_calc;
//        std::vector <Parameter> values;// (options.num_threads, model_params);
//        for (int i=0; i!=options.num_threads; ++i) values.push_back(model_params);
//        for (int i=0; i!=model_params.get_total_params(); ++i) values.push_back(model_params.get(i));
        std::vector <std::vector<double> > values(options.num_threads, std::vector<double>(model_params.get_total_params(), 0.0));
        for (int i=0; i!=options.num_threads; ++i) {
            for (int j=0; j!=model_params.get_total_params(); ++j) {
                values[i][j] = model_params.get(j);
            }
        }
//        printf("Size of values = %d\n",values.size());
        double reporting_rate = 1.0;
        if (model_params.param_exists("reporting")) {
            reporting_rate = model_params.get("reporting");
        }
        std::vector <std::string> param_names = model_params.get_names_vector();
        std::vector <std::vector<std::string> > param_names_threads (options.num_threads);
        if (model_params.param_exists("time_before_data")) {
            add_dt = model_params.get("time_before_data");
        }
        if (options.save_traj) {
            if (add_dt > 0) {
                particles.start_particle_tracing(add_dt+total_dt, num_groups);
            }
            else if (add_dt < 0) {
                particles.start_particle_tracing(add_dt+total_dt, num_groups);
                total_steps = ceil((double)(total_dt+add_dt)/(double)options.pfilter_every);
            }
            else {
                particles.start_particle_tracing(total_dt, num_groups);
            }
        }
        std::vector <Model> models;
        for (int i=0; i<options.num_threads; ++i) {
            models.push_back(sim_model);
        }
        std::vector <int> add_dt_threads (options.num_threads, add_dt);
        std::vector <int> start_dt_threads (options.num_threads, 0);
        std::vector <int> end_dt_threads (options.num_threads, add_dt);
        std::vector <double> dt_threads (options.num_threads, sim_dt);
        std::vector <int> total_dt_threads(options.num_threads, total_dt);
        std::vector <double> reporting_rate_threads(options.num_threads, reporting_rate);
        std::vector <int> num_groups_threads(options.num_threads, num_groups);
        // Simulate model and calculate likelihood assuming no observed data
        if (model_params.param_exists("time_before_data")) {
            if (add_dt > 0) {
                omp_set_num_threads(options.num_threads);
//                std::vector <Trajectory *> curr_trajs;
//                for (int i=0; i!=num_particles; ++i) {
//                    curr_trajs.push_back(particles.get_traj(i));
//                }
#pragma omp parallel for shared(particles, values)
                for (int i=0; i<num_particles; i++) {
                    int tn = omp_get_thread_num();
                    gsl_rng* r = gsl_rng_alloc( gsl_rng_mt19937 );
                    gsl_rng_set( r, omp_get_thread_num() + i );
                    // Adjust length of trajectory
                    particles.get_traj(i)->resize(add_dt, num_groups);
                    models[tn].simulate(values[tn], param_names_threads[tn], particles.get_traj(i), 0, add_dt_threads[tn], dt_threads[tn], total_dt_threads[tn], r);
                    if (options.which_likelihood<2) {
                        double w = likelihood_calc.binomial_lik(reporting_rate_threads[tn], particles.get_traj(i)->get_total_traj(), add_dt_threads[tn]+total_dt_threads[tn], 0, add_dt_threads[tn], num_groups_threads[tn], false);
                        particles.set_weight(w, i, false);
                    }
                    if (options.save_traj) {
                        particles.save_traj_to_matrix(i, 0, add_dt);
                        particles.save_ancestry(i, 0, add_dt);
                    }
                }
            }
        }
        init_seed += num_particles;
        int t=0;
        int start_dt;
        int end_dt;
        for (t=0; t!=total_steps; ++t) {
//            std::vector<double> we(options.particles, 0.0), wg(options.particles, 0.0);
            start_dt = t*options.pfilter_every;
            end_dt = std::min(total_dt, (t+1)*options.pfilter_every);
            std::fill(start_dt_threads.begin(), start_dt_threads.end(), start_dt);
            std::fill(end_dt_threads.begin(), end_dt_threads.end(), end_dt);
            omp_set_num_threads(options.num_threads);
#pragma omp parallel for shared (particles, values)
            for (int i=0; i<num_particles; i++) {
                int tn = omp_get_thread_num();
                gsl_rng* r = gsl_rng_alloc( gsl_rng_mt19937 );
                gsl_rng_set( r, omp_get_thread_num() + i );
                // Adjust length of trajectory
                //                    if (tn==0) std::cout << i << ' ' << std::endl;
                particles.get_traj(i)->resize(end_dt-start_dt, options.num_groups);
                models[tn].simulate(values[tn], param_names_threads[tn], particles.get_traj(i), start_dt_threads[tn], end_dt_threads[tn], dt_threads[tn], total_dt_threads[tn], r);
                double w = 1.0;
                double temp = 0.0;
                if (options.which_likelihood<2) {
                    double A = particles.get_traj(i)->get_total_traj();
                    temp = likelihood_calc.binomial_lik(reporting_rate_threads[tn], A, epi_data.get_data_ptr(0), add_dt_threads[tn]+total_dt_threads[tn], start_dt_threads[tn], end_dt_threads[tn], add_dt_threads[tn], num_groups_threads[tn], false);
                    w *= temp;
                    //                        we[i] = log(temp);
                }
                if (options.which_likelihood != 1) {
                    temp = likelihood_calc.coalescent_lik(particles.get_traj(i)->get_traj_ptr(0, 0), particles.get_traj(i)->get_traj_ptr(1, 0),
                                                          tree_data.get_binomial_ptr(0), tree_data.get_interval_ptr(0), tree_data.get_ends_ptr(0),
                                                          start_dt_threads[tn], end_dt_threads[tn], add_dt_threads[tn], false);
                    w *= temp;
                    //                        wg[i] = log(temp);
                }
                particles.set_weight(w, i, true);
                if (options.save_traj) {
                    particles.save_traj_to_matrix(i, start_dt_threads[tn]+add_dt_threads[tn], end_dt_threads[tn]+add_dt_threads[tn]);
                    particles.save_ancestry(i, start_dt_threads[tn]+add_dt_threads[tn], end_dt_threads[tn]+add_dt_threads[tn]);
                }
            }
//            std::cout << "Epi Weight: " << std::accumulate(we.begin(), we.end(), 0.0) << " Gen Weight: " << std::accumulate(wg.begin(), wg.end(), 0.0) << " Total: " << particles.get_total_weight() << std::endl;
            double curr_ESS = particles.get_ESS();
            if (curr_ESS < ESS_threshold) {
                double total_weight = particles.get_total_weight();
                if (total_weight == 0.0) {
                    loglik += -0.1*std::numeric_limits<double>::max();
//                    std::cout << std::accumulate(epi_data.get_data_ptr(0)+start_dt, epi_data.get_data_ptr(0)+end_dt, 0.0) << " : " << particles.get_traj(0)->get_traj(0) << std::endl;
                    std::cout << "stop time: " << end_dt << std::endl;
                    break;
                } else {
                    loglik += log(total_weight) - log(num_particles);
                }
                particles.resample(options.rng[0]);
            }
            else {
                particles.reset_parents();
            }
        }
        if (options.save_traj) {
            output_traj.resize((total_dt+add_dt), num_groups);
            //if (loglik > -0.1*std::numeric_limits<double>::max()) {
                particles.retrace_traj(output_traj, options.rng[0]);
            //}
        }
        for (int i=0; i!=num_particles; ++i) {
            particles.get_traj(i)->reset();
        }
        std::vector < std::vector<double> >().swap(values);
        return (loglik);
    }
}
