//
//  particle.cpp
//  EpiGenMCMC
//
//  Created by Lucy Li on 13/04/2016.
//  Copyright (c) 2016 Lucy Li, Imperial College London. All rights reserved.
//

#include "particle.h"
#include <numeric> // std::accumulate
#include <chrono>


Particle::Particle() {}

Particle::Particle(int input_num_particles, Trajectory init_traj) {
    num_particles = input_num_particles;
    for (int i=0; i!=num_particles; ++i) {
        trajectories.push_back(new Trajectory(init_traj));
    }
    next_particles.resize(num_particles);
    weights.resize(num_particles);
    normalised_weights.resize(num_particles);
    non_zero_particles.resize(num_particles);
    cumulative_weights.resize(num_particles);
    std::fill(weights.begin(), weights.end(), 1.0);
    parents.resize(num_particles);
    for (int i=0; i!=num_particles; ++i) {
        parents[i] = i;
    }
    num_groups = 1;
}

void Particle::start_particle_tracing(int input_num_time_steps, int input_num_groups) {
    num_groups = input_num_groups;
    num_time_steps = input_num_time_steps;
    overall_traj.clear();
    particle_ancestry.clear();
    overall_traj.resize(num_time_steps*num_groups*num_particles, 0.0);
    particle_ancestry.resize(num_time_steps*num_groups*num_particles, 0);
}

void Particle::save_traj_to_matrix(int start_dt, int end_dt) {
    for (int particlei=0; particlei!=num_particles; ++particlei) {
        for (int timei=start_dt; timei!=end_dt; ++timei) {
            for (int group_id=0; group_id!=num_groups; ++group_id) {
                int index = num_time_steps*num_groups*particlei+num_time_steps*group_id+timei;
                double value = trajectories[particlei]->get_traj(0, timei-start_dt, group_id);
                overall_traj[index] = value;
            }
        }
    }
}

void Particle::save_traj_to_matrix(int particle_id, int start_dt, int end_dt) {
    for (int timei=start_dt; timei!=end_dt; ++timei) {
        for (int group_id=0; group_id!=num_groups; ++group_id) {
            int index = num_time_steps*num_groups*particle_id+num_time_steps*group_id+timei;
            double value = trajectories[particle_id]->get_traj(0, timei-start_dt, group_id);
            overall_traj[index] = value;
        }
    }
}

void Particle::save_ancestry(int start_dt, int end_dt) {
    for (int particlei=0; particlei!=num_particles; ++particlei) {
        for (int timei=start_dt; timei!=end_dt; ++timei) {
            for (int group_id=0; group_id!=num_groups; ++group_id) {
                int ancestry_i = num_time_steps*num_groups*particlei+num_time_steps*group_id+timei;
                particle_ancestry[ancestry_i] = parents[particlei];
            }
        }
    }
}

void Particle::save_ancestry(int particle_i, int start_dt, int end_dt) {
    for (int timei=start_dt; timei!=end_dt; ++timei) {
        for (int group_id=0; group_id!=num_groups; ++group_id) {
            int ancestry_i = num_time_steps*num_groups*particle_i+num_time_steps*group_id+timei;
            particle_ancestry[ancestry_i] = parents[particle_i];
        }
    }
}


void Particle::reset_parents() {
    for (int i=0; i!=parents.size(); ++i) {
        parents[i] = i;
    }
}

int Particle::get_traj_random(gsl_rng * rng) {
    int select = gsl_rng_uniform_int(rng, num_particles);
    return (select);
}

void Particle::retrace_traj(Trajectory& output_traj, gsl_rng * rng) {
    int total_steps = num_time_steps*num_groups;
    int a = get_traj_random(rng); // Randomly select a particle to trace back its ancestry
    int timei = num_time_steps-1; // Current time step (i.e. last time step of the simulation)
    int groupi; // Current group number
    int index = 0;
    double traj_size;
    for (groupi=num_groups-1; groupi>=0; --groupi) {
        index = total_steps*a+num_time_steps*groupi+timei; // Position on the overall_traj
        traj_size = overall_traj[index]; // Epidemic size for particle a at time timei
        output_traj.set_traj(0, traj_size, timei, groupi);
    }
    a = particle_ancestry[index];
    for (timei=num_time_steps-2; timei>=0; --timei) {
        for (groupi=num_groups-1; groupi>=0; --groupi) {
            index = total_steps*a+num_time_steps*groupi+timei; // Position on the overall_traj
            traj_size = overall_traj[index];
            output_traj.set_traj(0, traj_size, timei, groupi);
        }
        a = particle_ancestry[index];
    }
}

Trajectory* Particle::get_traj(int index) const {
    return (trajectories[index]);
}



void Particle::set_weight(double w, int particle_num, bool multiply) {
    if (multiply) weights[particle_num] *= w;
    else weights[particle_num] = w;
}

double Particle::get_weight(int particle_num) {
    return(weights[particle_num]);
}
double Particle::get_total_weight() {
    double total = std::accumulate(weights.begin(), weights.end(), 0.0);
    return(total);
}

void Particle::normalise() {
    double total=get_total_weight();
    for (int i=0; i!=weights.size(); ++i) {
        normalised_weights[i] = weights[i] / total;
    }
}

double Particle::get_ESS() {
    if (get_total_weight() == 0.0) return (0.0);
    normalise();
    double total = 0.0;
    for (int i=0; i!=normalised_weights.size(); ++i) {
        total += normalised_weights[i] * normalised_weights[i];
    }
    return (1.0/total);
}

void Particle::replace(int old_part, int new_part) {
    trajectories[old_part]->replace(trajectories[new_part]);
//    Trajectory * to_delete = trajectories[old_part];
//    trajectories[old_part] = new Trajectory (*trajectories[new_part]);
//    delete(to_delete);
}

void rmultinomial (int num_samples, int N, std::vector<double>&p, std::vector<int>&n) {
    n.resize(num_samples, 0);
    double total_p = 0.0;
    for (int i=0; i!=N; ++i) {
        total_p += p[i];
        p[i] = total_p;
    }
//    std::uniform_real_distribution<double> unif (0.0, total_p);
    double interval = total_p/(double)num_samples;
    std::vector <double> z;
    for (int i=0; i!=num_samples; ++i) {
        z.push_back((double)(i+1)*interval);
    }
    bool less_than = false;
    int indicator = 0;
    for (int i=0; i!=num_samples; ++i) {
        less_than = (z[i] < p[indicator]);
        while (!less_than & (indicator < (N-1))) {
            ++indicator;
            less_than = (z[i] < p[indicator]);
        }
        ++n[indicator];
    }
}

void Particle::reset_weights() {
    std::fill(weights.begin(), weights.end(), 1.0);
}

void Particle::resample(gsl_rng * rng) {
    normalise();
    std::fill(next_particles.begin(), next_particles.end(), 0);
    gsl_ran_multinomial(rng, num_particles, num_particles, &normalised_weights[0], &next_particles[0]);
    // Create a vector 'empty' that stores the indices of particles NOT surviving to the next generation
    for (int i=0; i!=next_particles.size(); ++i) {
        if (next_particles[i]==0) empty.push_back(i);
    }
    // Replace the particles dying in the current generation with particles surviving to the next generation
    for (int i=0; i!=next_particles.size(); ++i) {
        if (next_particles[i] > 1) {
            parents[i] = i;
            for (int j=1; j!=next_particles[i]; ++j) {
                replace(empty.back(), i);
                parents[empty.back()] = i;
                empty.pop_back();
            }
        }
    }
    reset_weights();
}

void Particle::clear() {
    while(!trajectories.empty()) delete trajectories.back(), trajectories.pop_back();
//    std::vector<Trajectory*>::iterator it;
//    for (it=trajectories.begin(); it!=trajectories.end(); ++it) {
//        delete (*it);
//    }
    trajectories.clear();
}