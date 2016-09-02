//
//  model.cpp
//  EpiGenMCMC
//
//  Created by Lucy Li on 09/05/2016.
//  Copyright (c) 2016 Lucy Li, Imperial College London. All rights reserved.
//

#include "model.h"



Model::Model() {
    use_deterministic = false;
    custom_prob.resize(0);
    custom_prob_alpha = 1.0;
    custom_prob_scale = 1.0;
}

void Model::set_deterministic(bool logical) {
    use_deterministic = logical;
}
