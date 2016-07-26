//
//  particle.h
//  EpiGenMCMC
//
//  Created by Lucy Li on 13/04/2016.
//  Copyright (c) 2016 Lucy Li, Imperial College London. All rights reserved.
//

#ifndef __EpiGenMCMC__particle__
#define __EpiGenMCMC__particle__

#include <vector>
#include <gsl/gsl_randist.h>
#include "trajectory.h"

class Particle {
    int num_particles;
    int num_groups;
    int num_time_steps;
    std::vector <Trajectory*> trajectories;
    std::vector <double> weights;
    std::vector <double> normalised_weights;
    std::vector <unsigned int> next_particles;
    std::vector <int> parents;
    std::vector <double> overall_traj; // P1[Group1[T1, T2...], Group2[T1, T2...]...], P2[] ...
    std::vector <int> particle_ancestry;
    std::vector <int> non_zero_particles;
    std::vector <double> cumulative_weights;
    std::vector <int> empty;
    void replace(int, int);
    void normalise();
public:
    Particle();
    Particle(int, Trajectory);
    void start_particle_tracing(int, int);
    void save_traj_to_matrix(int, int);
    void save_traj_to_matrix(int, int, int);
    void save_ancestry(int, int);
    void save_ancestry(int, int, int);
    void reset_parents();
    int get_traj_random(gsl_rng *);
    void retrace_traj(Trajectory&, gsl_rng *);
    Trajectory* get_traj(int) const;
    void set_weight(double, int, bool);
    double get_weight(int);
    double get_total_weight();
    double get_ESS();
    void reset_weights();
    void resample(gsl_rng *);
    void clear();
};


#endif /* defined(__EpiGenMCMC__particle__) */
