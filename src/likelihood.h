//
//  likelihood.h
//  EpiGenMCMC
//
//  Created by Lucy Li on 13/04/2016.
//  Copyright (c) 2016 Lucy Li, Imperial College London. All rights reserved.
//

#ifndef __EpiGenMCMC__likelihood__
#define __EpiGenMCMC__likelihood__

#include <vector>
#include <stdio.h>
#include <limits>
#include <gsl/gsl_cdf.h>
#include "parameter.h"

class Likelihood {
public:
    Likelihood();
    double binomial_lik(double, double, int, int, int, int, bool);
    double binomial_lik(double, double, std::vector<double>::iterator, int, int, int, int, int, bool);
    double coalescent_lik(std::vector<double>::iterator, std::vector<double>::iterator, std::vector<double>::iterator, std::vector<double>::iterator, std::vector<double>::iterator, int, int, int, bool);
};

#endif /* defined(__EpiGenMCMC__likelihood__) */
