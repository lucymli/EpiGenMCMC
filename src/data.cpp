//
//  data.cpp
//  EpiGenMCMC
//
//  Created by Lucy Li on 13/04/2016.
//  Copyright (c) 2016 Lucy Li, Imperial College London. All rights reserved.
//

#include "data.h"


//
// TimeSeriesData
//

TimeSeriesData::TimeSeriesData(){};

TimeSeriesData::TimeSeriesData(std::string filename) {
    std::ifstream file (filename);
    double num;
    file >> num_time_steps >> time_step_size >> num_groups;
    bool first_found = false;
    if (num_groups == 1) {
        while (file >> num) {
            data.push_back(num);
            if (!first_found & (num > 0)) {
                first_data_time = data.size()-1.0;
                first_found = true;
            }
        }
    }
    else {
        int pos=0;
        int data_pos;
        data.resize(num_time_steps*num_groups);
        while (file >> num) {
            data_pos = (pos%num_groups)*num_time_steps+pos/num_groups;
            data[data_pos] = num;
            ++pos;
            if (!first_found & (num > 0)) {
                first_data_time = pos/num_groups;
                first_found = true;
            }
        }
    }
    file.close();
}

int TimeSeriesData::get_T() const {
    return (num_time_steps);
}

double TimeSeriesData::get_dt() const {
    return (time_step_size);
}

int TimeSeriesData::get_num_group() const {
    return (num_groups);
}

int TimeSeriesData::get_first() const {
    return(first_data_time);
}

double TimeSeriesData::get(int index) const {
    return (data[index]);
}
double TimeSeriesData::get(int time_index, int group_index) const {
    return (data[group_index*num_time_steps+time_index]);
}

std::vector<double>::iterator TimeSeriesData::get_data_ptr(int index) {
    return (data.begin());
}

//
// TreeData
//

TreeData::TreeData(){}

TreeData::TreeData(std::string filename) {
    std::ifstream file (filename);
    std::string line;
    file >> num_time_steps >> time_step_size;
    bool first_found = false;
    // Read in binomial data
    for (int t=0; t!=num_time_steps; ++t) {
        std::getline(file, line);
        if (line.empty()) std::getline(file, line);
        std::istringstream input(line);
        if (t==0) starts.push_back(0);
        else {
            starts.push_back(binomials.size());
        }
        while (std::getline(input, line, ' ')) {
            binomials.push_back (std::stod(line));
            if (!first_found & (std::stod(line) > 1.0)) {
                first_data_time = t;
                first_found = true;
            }
        }
        ends.push_back(binomials.size()-1.0);
    }
    // Read in intervals data
    for (int t=0; t!=num_time_steps; ++t) {
        std::getline(file, line);
        std::istringstream input(line);
        while (std::getline(input, line, ' ')) {
            intervals.push_back (std::stod(line));
        }
    }
    file.close();
}

int TreeData::get_T() const {
    return (num_time_steps);
}

double TreeData::get_dt() const {
    return (time_step_size);
}

int TreeData::get_first() const {
    return(first_data_time);
}

int TreeData::get_num_events(int time_index) const {
    return (ends[time_index]-starts[time_index]+1);
}

void TreeData::get_binomial (int time_index, std::vector <double>& output) const {
    for (int i=starts[time_index]; i<=ends[time_index]; ++i) {
        output[i-starts[time_index]] = binomials[i];
    }
}

void TreeData::get_binomial (int start_time, int end_time, std::vector <double>& output) const {
    for (int i=starts[start_time]; i<=ends[end_time-1]; ++i) {
        output[i-starts[start_time]] = binomials[i];
    }
}

double TreeData::get_binomial (int time_index, int event_index) const {
    return (binomials[starts[time_index]+event_index]);
}

std::vector<double>::iterator TreeData::get_binomial_ptr(int index) {
    return (binomials.begin()+starts[index]);
}

void TreeData::get_interval (int time_index, std::vector <double>& output) const {
    for (int i=starts[time_index]; i<=ends[time_index]; ++i) {
        output[i-starts[time_index]] = intervals[i];
    }
}

double TreeData::get_interval (int time_index, int event_index) const {
    return (intervals[starts[time_index]+event_index]);
}

void TreeData::get_interval (int start_time, int end_time, std::vector <double>& output) const {
    for (int i=starts[start_time]; i<=ends[end_time-1]; ++i) {
        output[i-starts[start_time]] = intervals[i];
    }
}

std::vector<double>::iterator TreeData::get_interval_ptr(int index) {
    return (intervals.begin()+starts[index]);
}

std::vector<double>::iterator TreeData::get_starts_ptr(int index) {
    return(starts.begin()+index);
}

std::vector<double>::iterator TreeData::get_ends_ptr(int index) {
    return(ends.begin()+index);
}


//
// MultiTreeData
//

MultiTreeData::MultiTreeData(){}

MultiTreeData::MultiTreeData(std::string filename) {
    std::ifstream file (filename);
    std::string tree_filename;
    file >> num_trees;
    for (int i=0; i!=num_trees; ++i) {
        file >> tree_filename;
        TreeData *ptr = new TreeData(tree_filename);
        tree_ptrs.push_back(ptr);
    }
    file.close();
    num_time_steps = tree_ptrs[0]->get_T();
    time_step_size = tree_ptrs[0]->get_dt();
    for (int i=0; i!=num_trees; ++i) {
        first_data_time = std::min(first_data_time, tree_ptrs[i]->get_first());
    }
};

int MultiTreeData::get_num_tree() const {
    return (num_trees);
}

int MultiTreeData::get_T() const {
    return (num_time_steps);
}

double MultiTreeData::get_dt() const {
    return (time_step_size);
}

int MultiTreeData::get_first() const {
    return (first_data_time);
}
