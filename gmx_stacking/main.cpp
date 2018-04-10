//
//  main.cpp
//  gmx_stacking
//
//  Created by Yiming Tang on 09/04/2018.
//  Copyright © 2018 Yiming Tang. All rights reserved.
//

#include <iostream>
#include "stacking.hpp"

int main(int argc, char * argv[]) {
    // insert code here...
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<stacking>(argc, argv);

}
