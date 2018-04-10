//
//  stacking.hpp
//  gmx_stacking
//
//  Created by Yiming Tang on 09/04/2018.
//  Copyright © 2018 Yiming Tang. All rights reserved.
//

#ifndef stacking_hpp
#define stacking_hpp

#include <stdio.h>
#include <gromacs/trajectoryanalysis.h>
#include <string>
#include <vector>

using namespace std;
using namespace gmx;

struct coordinate
{
    double x;
    double y;
    double z;
};


class stacking: public TrajectoryAnalysisModule
{
public:
    stacking();
    
    virtual void initOptions(IOptionsContainer          *options,
                             TrajectoryAnalysisSettings *settings);
    
    virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                              const TopologyInformation        &top);
    
    virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                              TrajectoryAnalysisModuleData *pdata);
    
    virtual void finishAnalysis(int nframes);
    
    virtual void writeOutput();
    
    double calDistance(ConstArrayRef<rvec> Group1, ConstArrayRef<rvec> Group2);
    double calAngle(ConstArrayRef<rvec> Group1, ConstArrayRef<rvec> Group2);
    
private:
    
    AnalysisData    data_probability_;
    
    std::string     fn_density_;
    
    SelectionList   ref_;
    SelectionList   sel_;
    
    double          cutoff_;
    
    AnalysisDataAverageModulePointer avem_probability_;
    
    t_topology     *top_;
    t_atoms         atoms_;
    
    
    std::string     fnEnergySurface_;
    std::string     fnEnergySurfaceRaw_;
    std::string     fnPropability_;
    
    double          maxDistance_;
    double          stepDistance_;
    double          maxAngle_;
    double          stepAngle_;
    
    double          temperature_;
    
};

#endif /* stacking_hpp */
