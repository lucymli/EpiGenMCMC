//
//  data.h
//  EpiGenMCMC
//
//  Created by Lucy Li on 07/04/2016.
//  Copyright (c) 2016 Lucy Li, Imperial College London. All rights reserved.
//

#ifndef EpiGenMCMC_data_h
#define EpiGenMCMC_data_h

#include <vector>
#include <string>
#include <fstream>
#include <sstream>


//
// TimeSeriesData class
//

class TimeSeriesData {
    std::vector <double> data;
    int first_data_time;
    int num_time_steps;
    double time_step_size;
    int num_groups;
public:
    TimeSeriesData();
    TimeSeriesData(std::string);
    int get_T() const;
    double get_dt() const;
    int get_num_group() const;
    int get_first() const;
    double get(int) const;
    double get(int, int) const;
    std::vector<double>::iterator get_data_ptr(int);
};






//
// TreeData Class
//


class TreeData {
    std::vector <double> binomials;
    std::vector <double> intervals;
    std::vector <double> starts;
    std::vector <double> ends;
    int first_data_time;
    int num_time_steps;
    double time_step_size;
public:
    TreeData();
    TreeData(std::string);
    int get_T () const;
    double get_dt () const;
    int get_first() const;
    int get_num_events(int) const;
    void get_binomial (int, std::vector <double>&) const;
    void get_binomial (int, int, std::vector <double>&) const;
    double get_binomial (int, int) const;
    std::vector<double>::iterator get_binomial_ptr(int);
    void get_interval (int, std::vector <double>&) const;
    void get_interval (int, int, std::vector <double>&) const;
    double get_interval (int, int) const;
    std::vector<double>::iterator get_interval_ptr(int);
    std::vector<double>::iterator get_starts_ptr(int);
    std::vector<double>::iterator get_ends_ptr(int);
};









//
// MultiTreeData Class
//


class MultiTreeData {
    std::vector <TreeData* > tree_ptrs;
    int first_data_time;
    int num_trees;
    int num_time_steps;
    double time_step_size;
public:
    MultiTreeData();
    MultiTreeData(std::string);
    int get_num_tree() const;
    int get_T () const;
    double get_dt () const;
    int get_first() const;
};







#endif
