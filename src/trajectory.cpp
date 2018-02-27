//
//  trajectory.cpp
//  EpiGenMCMC
//
//  Created by Lucy Li on 13/04/2016.
//  Copyright (c) 2016 Lucy Li, Imperial College London. All rights reserved.
//

#include "trajectory.h"


void Trajectory::initialise_file (std::string filename, int every) {
    std::ofstream file(filename);
    file << "state";
    if (num_groups > 1) {
        for (int i=0; i<num_groups; ++i) {
            for (int t=0; t<num_time_steps; t+=every) {
                file << "\tT" << t << "Group" << i;
            }
        }
    }
    else {
        for (int t=0; t<num_time_steps; t+=every) {
            file << "\tIncT" << t;
        }
        for (int t=0; t<num_time_steps; t+=every) {
            file << "\tPrevT" << t;
        }
    }
    file << std::endl;
    file.close();
}

Trajectory::Trajectory() {
    trajectory.resize(1);
    trajectory2.resize(1);
    trajectory3.resize(1);
    initial_states.resize(1);
    recoveries.resize(1);
    states.resize(1);
}

Trajectory::Trajectory(int numT, int numGroups) {
    num_time_steps = numT;
    num_groups = numGroups;
    trajectory.resize(numT*numGroups);
    trajectory2.resize(numT);
    trajectory3.resize(numT);
    recoveries.resize(numT);
//    initialise_file(filename);
}


void Trajectory::resize(int length) {
    if (trajectory.size() != length) {
        trajectory.resize(length);
        trajectory2.resize(length);
        trajectory3.resize(length);
    }
    fill(trajectory.begin(), trajectory.end(), 0.0);
    fill(trajectory2.begin(), trajectory2.end(), 0.0);
    fill(trajectory3.begin(), trajectory3.end(), 0.0);
//    num_time_steps = length;
}

void Trajectory::resize(int length, int groups) {
    if (trajectory.size() != length*groups) {
        trajectory.resize(length*groups);
        trajectory2.resize(length);
        trajectory3.resize(length);
    }
    fill(trajectory.begin(), trajectory.end(), 0.0);
    fill(trajectory2.begin(), trajectory2.end(), 0.0);
    fill(trajectory3.begin(), trajectory3.end(), 0.0);
    //    num_time_steps = length;
}

void Trajectory::print_to_file(int iteration, std::string filename, int every, bool sum_across) {
    std::ofstream file(filename, std::ios::app);
    file << iteration;
    int total_steps = (int) trajectory.size();
    if (num_groups > 1) {
        for (int i=0; i<num_groups; ++i) {
            for (int t=0; t<total_steps; t+=every) {
                double count = 0.0;
                if (sum_across) {
                    for (int j=0; j<std::min(every, iteration-t); ++j) {
                        count += trajectory[i*num_time_steps+t+j];
                    }
                }
                else {
                    count = trajectory[i*num_time_steps+t];
                }
                file << "\t" << count;
            }
        }
    }
    else {
        for (int t=0; t<total_steps; t+=every) {
            double count = 0.0;
            if (sum_across) {
                for (int j=0; j<std::min(every, total_steps-t); ++j) {
                    count += trajectory[t+j];
                }
            }
            else {
                count = trajectory[t];
            }
            file << "\t" << count;
        }
        for (int t=0; t<total_steps; t+=every) {
            double count = 0.0;
            if (sum_across) {
                for (int j=0; j<std::min(every, num_time_steps-t); ++j) {
                    count += trajectory2[t+j];
                }
                count /= (double) std::min(every, num_time_steps-t);
            }
            else {
                count = trajectory2[t];
            }
            file << "\t" << std::round(count);
        }
    }
    file << std::endl;
    file.close();
}

void Trajectory::print_to_file(std::string filename, int every, bool sum_across) {
    std::ofstream file(filename, std::ios::app);
    int total_steps = (int) trajectory.size();
    if (num_groups > 1) {
        for (int i=0; i<num_groups; ++i) {
            for (int t=0; t<num_time_steps; t+=every) {
                double count = 0.0;
                if (sum_across) {
                    for (int j=0; j<std::min(every, total_steps-t); ++j) {
                        count += trajectory[i*num_time_steps+t+j];
                    }
                }
                else {
                    count = trajectory[i*num_time_steps+t];
                }
                file << "\t" << count;
            }
        }
    }
    else {
        for (int t=0; t<total_steps; t+=every) {
            double count = 0.0;
            if (sum_across) {
                for (int j=0; j<std::min(every, num_time_steps-t); ++j) {
                    count += trajectory[t+j];
                }
            }
            else {
                count = trajectory[t];
            }
            file << "\t" << count;
        }
        for (int t=0; t<total_steps; t+=every) {
            double count = 0.0;
            if (sum_across) {
                for (int j=0; j<std::min(every, num_time_steps-t); ++j) {
                    count += trajectory2[t+j];
                }
                count /= (double) std::min(every, num_time_steps-t);
            }
            else {
                count = trajectory2[t];
            }
            file << "\t" << std::round(count);
        }
    }
    file << std::endl;
    file.close();
}

void Trajectory::initialise_states(std::string filename) {
    std::ifstream file(filename);
    double num;
    while (file >> num) {
        initial_states.push_back(num);
        states.push_back(num);
    }
    file.close();
}

void Trajectory::initialise_states(std::vector<double>input_vector) {
    for (int i=0; i!=input_vector.size(); ++i) {
        states[i] = input_vector[i];
        initial_states[i] = input_vector[i];
    }
}


void Trajectory::set_traj(int i, double value, int time) {
    if (i==0) trajectory[time] = value;
    else if (i==1) trajectory2[time] = value;
    else trajectory3[time] = value;
}

void Trajectory::set_traj(int i, double value, int time, int group) {
    trajectory[group*num_time_steps+time] = value;
}

//void Trajectory::set_traj2(double value, int time) {
//    trajectory2[time] = value;
//}

void Trajectory::set_state(double value, int time) {
    states[time] = value;
}

double Trajectory::get_traj(int i, int time_index) const {
    if (i==0) return(trajectory[time_index]);
    else if (i==1) return (trajectory2[time_index]);
    else return (trajectory3[time_index]);
}

double Trajectory::get_traj(int i, int time_index, int group_id) const {
    int index = num_time_steps*group_id+time_index;
    return (get_traj(i, index));
}

//
//double Trajectory::get_traj2(int time_index) const {
//    return(trajectory2[time_index]);
//}

void Trajectory::resize_recoveries(int total_times) {
    recoveries.resize(total_times);
}

void Trajectory::add_recovery_time(int new_recovery_time) {
    ++recoveries[new_recovery_time];
}

void Trajectory::add_recovery_time(int new_recovery_time, int number) {
    recoveries[new_recovery_time] += number;
}

int Trajectory::num_recover_at(int time_to_find) {
    if (recoveries.size() < 1) return (0);
    return (recoveries[time_to_find]);
}

int Trajectory::num_recover_between(int time1, int time2) {
    if (recoveries.size() < 1) return (0);
    int total = 0;
    for (int i=time1; i!=time2; ++i) {
        total += recoveries[i];
    }
    return (total);
}

int Trajectory::num_recover_after(int time_to_find) {
    if (recoveries.size() < 1) return (0);
    int total = 0;
    for (int i=time_to_find; i!=recoveries.size(); ++i) {
        total += recoveries[i];
    }
    return (total);
}

void Trajectory::delete_recoveries_before(int time_to_find) {
   recoveries.erase(recoveries.begin(), recoveries.begin()+time_to_find);
}


std::vector <double> Trajectory::get_traj_range(int id, int start, int end) const {
    std::vector <double> output;
    if (num_groups > 1) {
        for (int group=0; group!=num_groups; ++group) {
            for (int i=start; i!=end; ++i) {
                if (id==0) output.push_back(trajectory[group*num_time_steps+i]);
                else if (id==1) output.push_back(trajectory2[group*num_time_steps+i]);
                else output.push_back(trajectory3[group*num_time_steps+i]);
            }
        }
    }
    else {
        for (int i=start; i!=end; ++i) {
            output.push_back(trajectory[i]);
        }
    }
    return(output);
}

std::vector <double> Trajectory::get_traj_range(int id, int start, int end, int group_id) const {
    std::vector <double> output;
    for (int i=start; i!=end; ++i) {
        if (id==0) output.push_back(trajectory[group_id*num_time_steps+i]);
        else if (id==1) output.push_back(trajectory2[group_id*num_time_steps+i]);
        else output.push_back(trajectory3[group_id*num_time_steps+i]);
    }
    return(output);
}

//std::vector <double> Trajectory::get_traj2_range(int start, int end) const {
//    std::vector <double> output;
//    for (int i=start; i!=end; ++i) {
//        output.push_back(trajectory2[i]);
//    }
//    return(output);
//}

std::vector<double>::iterator Trajectory::get_traj_ptr(int i, int index) {
    if (i==0) return(trajectory.begin()+index);
    else if (i==1) return(trajectory2.begin()+index);
    else return(trajectory3.begin()+index);
}

//std::vector<double>::iterator Trajectory::get_traj2_ptr(int index) {
//    return(trajectory2.begin()+index);
//}

double Trajectory::get_total_traj() const {
    double total_data = std::accumulate(trajectory.begin(), trajectory.end(), 0.0);
    return (total_data);
}

double* Trajectory::get_curr_states() {
    return (&states[0]);
}

double Trajectory::get_init_state(int index) const {
    return (states[index]);
}

double Trajectory::get_state(int index) const {
    return(states[index]);
}

int Trajectory::get_curr_states_size() const {
    return ((int)states.size());
}

void Trajectory::reset() {
//    fill(trajectory.begin(), trajectory.end(), 0.0);
//    fill(trajectory2.begin(), trajectory2.end(), 0.0);
    states = initial_states;
    fill(recoveries.begin(), recoveries.end(), 0.0);
}

void Trajectory::replace(Trajectory * new_traj) {
//    for (int i=0; i!=new_traj->trajectory.size(); ++i) {
//        trajectory[i] = new_traj->trajectory[i];
//        trajectory2[i] = new_traj->trajectory2[i];
//    }
    states.clear();
    states.insert(states.begin(), new_traj->states.begin(), new_traj->states.end());
    initial_states.clear();
    initial_states.insert(initial_states.begin(), new_traj->initial_states.begin(), new_traj->initial_states.end());
    for (int i=0; i!=trajectory.size(); ++i) {
        trajectory[i] = new_traj->get_traj(0, i);
    }
    for (int i=0; i!=trajectory2.size(); ++i) {
        trajectory2[i] = new_traj->get_traj(1, i);
    }
    for (int i=0; i!=trajectory3.size(); ++i) {
        trajectory3[i] = new_traj->get_traj(2, i);
    }
    for (int i=0; i!=recoveries.size(); ++i) {
        recoveries[i] = new_traj->recoveries[i];
    }
    num_time_steps = new_traj->num_time_steps;
}
