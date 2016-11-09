//
//  trajectory.h
//  EpiGenMCMC
//
//  Created by Lucy Li on 08/04/2016.
//  Copyright (c) 2016 Lucy Li, Imperial College London. All rights reserved.
//

#ifndef EpiGenMCMC_trajectory_h
#define EpiGenMCMC_trajectory_h

#include <fstream>
#include <string>
#include <vector>
#include <numeric>
#include <iostream>
#include <cmath>

class Trajectory{
    std::vector <double> trajectory;
    std::vector <double> trajectory2;
    std::vector <double> trajectory3;
    std::vector <double> initial_states;
    std::vector <double> states;
    std::vector <int> recoveries;
    int num_time_steps;
    int num_groups;
public:
    Trajectory();
    Trajectory(int, int);
    void initialise_file(std::string, int);
    void resize(int);
    void resize(int, int);
    void initialise_states(std::string);
    void initialise_states(std::vector<double>);
    void print_to_file (int, std::string, int, bool);
    void print_to_file (std::string, int, bool);
    void change_num_groups(int);
    void set_traj(int, double, int);
    void set_traj(int, double, int, int);
//    void set_traj2(double, int);
    void set_state(double, int);
    double get_traj(int, int) const;
    double get_traj(int, int, int) const;
//    double get_traj2(int) const;
    double get_total_traj() const;
    void resize_recoveries(int);
    void add_recovery_time(int);
    void add_recovery_time(int, int);
    int num_recover_at(int);
    int num_recover_between(int, int);
    int num_recover_after(int);
    void delete_recoveries_before(int);
    std::vector<double> get_traj_range(int, int, int) const;
    std::vector<double> get_traj_range(int, int, int, int) const;
//    std::vector<double> get_traj2_range(int, int) const;
    std::vector<double>::iterator get_traj_ptr(int, int);
//    std::vector<double>::iterator get_traj2_ptr(int);
    double * get_curr_states();
    double get_init_state(int) const;
    double get_state(int) const;
    int get_curr_states_size() const;
    void reset();
    void replace(Trajectory*);
};



#endif
