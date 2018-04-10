//
//  stacking.cpp
//  gmx_stacking
//
//  Created by Yiming Tang on 09/04/2018.
//  Copyright © 2018 Yiming Tang. All rights reserved.
//

#include "stacking.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <gromacs/trajectoryanalysis.h>
#include <gromacs/fileio/matio.h>
#include <gromacs/fileio/gmxfio.h>


double stacking::calDistance(ConstArrayRef<rvec> Group1, ConstArrayRef<rvec> Group2)
{
    ConstArrayRef<rvec>::iterator iter_group1;
    ConstArrayRef<rvec>::iterator iter_group2;
    
    //cout << Group1.size() << '\t';
    //cout << Group2.size() << endl;;
    
    double x_1 = 0.0, y_1 = 0.0, z_1 = 0.0;
    double x_2 = 0.0, y_2 = 0.0, z_2 = 0.0;
    
    for(iter_group1 = Group1.begin(); iter_group1 != Group1.end(); ++iter_group1)
    {
        x_1 += *iter_group1[0];
        y_1 += *iter_group1[1];
        z_1 += *iter_group1[2];
    }
    
    x_1 = x_1 / Group1.size();
    y_1 = y_1 / Group1.size();
    z_1 = z_1 / Group1.size();
    
    for(iter_group2 = Group2.begin(); iter_group2 != Group2.end(); ++iter_group2)
    {
        x_2 += *iter_group2[0];
        y_2 += *iter_group2[1];
        z_2 += *iter_group2[2];
    }
    
    x_2 = x_2 / Group2.size();
    y_2 = y_2 / Group2.size();
    z_2 = z_2 / Group2.size();
    
    
    //cout << x_1 << '\t' << y_1 << '\t' << z_1 << endl;
    //cout << x_2 << '\t' << y_2 << '\t' << z_2 << endl;
    //cout << sqrt( pow(x_1 - x_2, 2.0) + pow(y_1 - y_2, 2.0) + pow(z_1 - z_2, 2.0)) << endl << endl;
    return sqrt( pow(x_1 - x_2, 2.0) + pow(y_1 - y_2, 2.0) + pow(z_1 - z_2, 2.0));
    
}

double stacking::calAngle(ConstArrayRef<rvec> Group1, ConstArrayRef<rvec> Group2)
{
    ConstArrayRef<rvec>::iterator iter_group1;
    ConstArrayRef<rvec>::iterator iter_group2;
    
    if (Group1.size() < 3 || Group2.size() < 3)
    {
        cerr << "You have groups that contains less than three atoms!" << endl;
        exit(1);
    }
    
    coordinate p1, p2, p3;
    
    iter_group1 = Group1.begin();
    p1.x = *iter_group1[0];
    p1.y = *iter_group1[1];
    p1.z = *iter_group1[2];
    
    iter_group1++;
    p2.x = *iter_group1[0];
    p2.y = *iter_group1[1];
    p2.z = *iter_group1[2];
    
    iter_group1++;
    p3.x = *iter_group1[0];
    p3.y = *iter_group1[1];
    p3.z = *iter_group1[2];
    
    cout << p1.x << '\t' << p1.y << '\t' << p1.z << endl;
    cout << p2.x << '\t' << p2.y << '\t' << p2.z << endl;
    cout << p3.x << '\t' << p3.y << '\t' << p3.z << endl;
    
    double a1 = ( (p2.y-p1.y)*(p3.z-p1.z)-(p2.z-p1.z)*(p3.y-p1.y) );
    double b1 = ( (p2.z-p1.z)*(p3.x-p1.x)-(p2.x-p1.x)*(p3.z-p1.z) );
    double c1 = ( (p2.x-p1.x)*(p3.y-p1.y)-(p2.y-p1.y)*(p3.x-p1.x) );
    
    cout << a1 << '\t' << b1 << '\t' << c1 << endl;
    
    iter_group2 = Group2.begin();
    p1.x = *iter_group2[0];
    p1.y = *iter_group2[1];
    p1.z = *iter_group2[2];
    
    iter_group2++;
    p2.x = *iter_group2[0];
    p2.y = *iter_group2[1];
    p2.z = *iter_group2[2];
    
    iter_group2++;
    p3.x = *iter_group2[0];
    p3.y = *iter_group2[1];
    p3.z = *iter_group2[2];
    
    cout << p1.x << '\t' << p1.y << '\t' << p1.z << endl;
    cout << p2.x << '\t' << p2.y << '\t' << p2.z << endl;
    cout << p3.x << '\t' << p3.y << '\t' << p3.z << endl;
    
    double a2 = ( (p2.y-p1.y)*(p3.z-p1.z)-(p2.z-p1.z)*(p3.y-p1.y) );
    double b2 = ( (p2.z-p1.z)*(p3.x-p1.x)-(p2.x-p1.x)*(p3.z-p1.z) );
    double c2 = ( (p2.x-p1.x)*(p3.y-p1.y)-(p2.y-p1.y)*(p3.x-p1.x) );
    
    cout << a2 << '\t' << b2 << '\t' << c2 << endl;
    
    double angle = acos((a1*a2+b1*b2+c1*c2) / sqrt(a1*a1+b1*b1+c1*c1) / sqrt(a2*a2+b2*b2+c2*c2)) / (2 * M_PI) * 180;
    
    cout << angle << endl << endl;
    
    return min(abs(angle-0), abs(180 - angle));
}

stacking::stacking(): cutoff_(1.5), maxDistance_(1.2), maxAngle_(90.0), stepDistance_(0.01), stepAngle_(1.0)
{
    registerAnalysisDataset(&data_probability_, "probability");
}

void stacking::initOptions(gmx::IOptionsContainer *options, gmx::TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] =
    {
        "Stacking Free Energy Surface (SFES) Calculator by Yiming Tang @ Fudan. \n",
        "This Tool calculates the Joint Probabilit Density (JPD) of two molecules",
        "as a function of angle (0-90) - distance (0-1.2). The tool ouputs the free",
        "energy surface defined by E = - R * T * log(H), where H is the probability.",
        "This tool takes two selection groups named ref and sel, each should contains",
        "groups each of which contains one benzene. Angles and Distances will be calculated",
        "on each group pair between ref and sel."
    };
    
    //settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
    
    settings->setHelpText(desc);
    
    options->addOption(FileNameOption("energy").filetype(eftUnknown).legacyType(efXPM).outputFile()
                       .store(&fnEnergySurface_).defaultBasename("energy")
                       .description("Energy Surface xpm map").required());
    
    options->addOption(FileNameOption("raw-energy").filetype(eftGenericData).outputFile()
                       .store(&fnEnergySurfaceRaw_)
                       .description("Ras Energy Surface Data for further analysis"));
    
    options->addOption(FileNameOption("raw-probability").filetype(eftGenericData).outputFile()
                       .store(&fnPropability_)
                       .description("Raw Probability Data for further analysis"));
    
    options->addOption(DoubleOption("cutoff").store(&cutoff_)
                       .description("Benzenes whose minimum distance are beyond this cutoff is not considered contacted"));
    
    options->addOption(DoubleOption("maxD").store(&maxDistance_).defaultValue(1.2)
                       .description("Maxmimum Centroid Distance to output."));
    
    options->addOption(DoubleOption("stepD").store(&stepDistance_).defaultValue(0.01)
                       .description("Step of Centroid Distance to output."));
    
    options->addOption(DoubleOption("maxA").store(&maxAngle_).defaultValue(90)
                       .description("Maxmimum Centroid Angle to output."));
    
    options->addOption(DoubleOption("stepA").store(&stepAngle_).defaultValue(1)
                       .description("Step of Centroid Angle to output."));
    
    options->addOption(SelectionOption("ref")
                       .storeVector(&ref_).required().multiValue()
                       .description("Reference Groups of benzenes to calculate angle and distance"));
    
    options->addOption(SelectionOption("sel")
                       .storeVector(&sel_).required().multiValue()
                       .description("Selection Groups of benzenes to calculate angle and distance"));
    
    options->addOption(DoubleOption("temperature").store(&temperature_).defaultValue(298)
                       .description("Temperature in K for energy calculation"));
}

void stacking::initAnalysis(const gmx::TrajectoryAnalysisSettings &settings, const gmx::TopologyInformation &top)
{
    //atoms_ = top.topology()->atoms;
    //top_   = top.topology();
    
    // initial data set columns and rows
    
    data_probability_.setDataSetCount((int)(maxDistance_ / stepDistance_));
    
    for (int i = 0; i < (int)(maxDistance_ / stepDistance_); i++)
    {
        data_probability_.setColumnCount(i, (int)(maxAngle_ / stepAngle_));
    }
    
    avem_probability_.reset(new AnalysisDataAverageModule());
    data_probability_.addModule(avem_probability_);
    
}

void stacking::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc, gmx::TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle dhProbability = pdata->dataHandle(data_probability_);
    dhProbability.startFrame(frnr, fr.time);
    
    // Constructed a new temp all-zero matrix
    
    int **tempProbability = new int*[(int)(maxDistance_ / stepDistance_)];
    for(int i = 0; i < (int)(maxDistance_ / stepDistance_); i++)
    {
        tempProbability[i] = new int[(int)(maxAngle_ / stepAngle_)];
    }
    for(int i = 0 ; i < (int)(maxDistance_ / stepDistance_); i++)
    {
        for(int j = 0; j < (int)(maxAngle_ / stepAngle_); j++)
        {
            tempProbability[i][j] = 0;
        }
    }
    
    // Begin analysis
    
    for (size_t gRef = 0 ; gRef < ref_.size(); ++gRef)
    {
        for (size_t gSel = 0; gSel < sel_.size(); ++gSel)
        {
            //cout << "ref_size = " << ref_.size() << endl;
            //cout << "sel_size = " << sel_.size() << endl;
            
            
            const Selection &ref = pdata->parallelSelection(ref_[gRef]);
            const Selection &sel = pdata->parallelSelection(sel_[gSel]);
            
            if (ref.atomIndices()[0] == sel.atomIndices()[0])
            {
                continue;
            }
            
            //cout << ref.coordinates().size() << '\t' << sel.coordinates().size() << endl;
            /*
            double myDistance = calDistance(ref.coordinates(), sel.coordinates());
            double myAngle    = calAngle(ref.coordinates(), sel.coordinates());
            */
            
            ////// Calculating distance //////////////
            
            //ConstArrayRef<rvec>::iterator iter_group1;
            //ConstArrayRef<rvec>::iterator iter_group2;
            
            //cout << Group1.size() << '\t';
            //cout << Group2.size() << endl;;
            
            double x_1 = 0.0, y_1 = 0.0, z_1 = 0.0;
             double x_2 = 0.0, y_2 = 0.0, z_2 = 0.0;
            
            for (int i = 0; i < ref.coordinates().size(); i++)
            {
                x_1 += ref.coordinates()[i][0];
                y_1 += ref.coordinates()[i][1];
                z_1 += ref.coordinates()[i][2];
                //cout << ref.coordinates()[i][0] << '\t' << ref.coordinates()[i][1] << '\t' << ref.coordinates()[i][2] << endl;
            }
            /*
            for(iter_group1 = ref.coordinates().begin(); iter_group1 != ref.coordinates().end(); ++iter_group1)
            {
                x_1 += *iter_group1[0];
                y_1 += *iter_group1[1];
                z_1 += *iter_group1[2];
            }
            */
            x_1 = x_1 / ref.coordinates().size();
            y_1 = y_1 / ref.coordinates().size();
            z_1 = z_1 / ref.coordinates().size();
            //cout << "CENTROID: " << x_1 << '\t' << y_1 << '\t' << z_1 << endl;
            
            
            for (int i = 0; i < sel.coordinates().size(); i++)
            {
                x_2 += sel.coordinates()[i][0];
                y_2 += sel.coordinates()[i][1];
                z_2 += sel.coordinates()[i][2];
                //cout << sel.coordinates()[i][0] << '\t' << sel.coordinates()[i][1] << '\t' << sel.coordinates()[i][2] << endl;
            }
            
            /*
            for(iter_group2 = sel.coordinates().begin(); iter_group2 != sel.coordinates().end(); ++iter_group2)
            {
                x_2 += *iter_group2[0];
                y_2 += *iter_group2[1];
                z_2 += *iter_group2[2];
            }
            */
            
            x_2 = x_2 / sel.coordinates().size();
            y_2 = y_2 / sel.coordinates().size();
            z_2 = z_2 / sel.coordinates().size();
            //cout << "CENTROID: " << x_2 << '\t' << y_2 << '\t' << z_2 << endl;
            
            //cout << x_1 << '\t' << y_1 << '\t' << z_1 << endl;
            //cout << x_2 << '\t' << y_2 << '\t' << z_2 << endl;
            //cout << sqrt( pow(x_1 - x_2, 2.0) + pow(y_1 - y_2, 2.0) + pow(z_1 - z_2, 2.0)) << endl << endl;
            double myDistance = sqrt( pow(x_1 - x_2, 2.0) + pow(y_1 - y_2, 2.0) + pow(z_1 - z_2, 2.0));
            
            //cout << myDistance << endl << endl ;
            
            //////////////////////////////////////////
            
            /// Calculating Angle ////////////////////
            
            if (ref.coordinates().size() < 3 || sel.coordinates().size() < 3)
            {
                cerr << "You have groups that contains less than three atoms!" << endl;
                exit(1);
            }
            
            coordinate p1, p2, p3;
            
            p1.x = ref.coordinates()[0][0];
            p1.y = ref.coordinates()[0][1];
            p1.z = ref.coordinates()[0][2];
            
            p2.x = ref.coordinates()[1][0];
            p2.y = ref.coordinates()[1][1];
            p2.z = ref.coordinates()[1][2];
            
            p3.x = ref.coordinates()[2][0];
            p3.y = ref.coordinates()[2][1];
            p3.z = ref.coordinates()[2][2];
            
            //cout << p1.x << '\t' << p1.y << '\t' << p1.z << endl;
            //cout << p2.x << '\t' << p2.y << '\t' << p2.z << endl;
            //cout << p3.x << '\t' << p3.y << '\t' << p3.z << endl;
            
            double a1 = ( (p2.y-p1.y)*(p3.z-p1.z)-(p2.z-p1.z)*(p3.y-p1.y) );
            double b1 = ( (p2.z-p1.z)*(p3.x-p1.x)-(p2.x-p1.x)*(p3.z-p1.z) );
            double c1 = ( (p2.x-p1.x)*(p3.y-p1.y)-(p2.y-p1.y)*(p3.x-p1.x) );
            
           // cout << ":::: \t" << a1 << '\t' << b1 << '\t' << c1 << endl;
            
            p1.x = sel.coordinates()[0][0];
            p1.y = sel.coordinates()[0][1];
            p1.z = sel.coordinates()[0][2];
            
            p2.x = sel.coordinates()[1][0];
            p2.y = sel.coordinates()[1][1];
            p2.z = sel.coordinates()[1][2];
            
            p3.x = sel.coordinates()[2][0];
            p3.y = sel.coordinates()[2][1];
            p3.z = sel.coordinates()[2][2];
            
            //cout << p1.x << '\t' << p1.y << '\t' << p1.z << endl;
            //cout << p2.x << '\t' << p2.y << '\t' << p2.z << endl;
            //cout << p3.x << '\t' << p3.y << '\t' << p3.z << endl;
            
            double a2 = ( (p2.y-p1.y)*(p3.z-p1.z)-(p2.z-p1.z)*(p3.y-p1.y) );
            double b2 = ( (p2.z-p1.z)*(p3.x-p1.x)-(p2.x-p1.x)*(p3.z-p1.z) );
            double c2 = ( (p2.x-p1.x)*(p3.y-p1.y)-(p2.y-p1.y)*(p3.x-p1.x) );
            
            //cout << ":::: \t" << a2 << '\t' << b2 << '\t' << c2 << endl;
            
            double angle = acos((a1*a2+b1*b2+c1*c2) / sqrt(a1*a1+b1*b1+c1*c1) / sqrt(a2*a2+b2*b2+c2*c2)) / (M_PI) * 180;
            
            //cout << angle << endl << endl;
            
            double myAngle =  min(abs(angle-0), abs(180 - angle));
            
            //////////////////////////////////////////
            
            
            
            
            if (myDistance <= maxDistance_ && myAngle <= maxAngle_ )
            {
                //cout << myDistance << '\t' << myAngle << endl;
                //cout <<(int)floor(myDistance / stepDistance_) << '\t' << (int)floor(myAngle/stepAngle_) << endl << endl;
                tempProbability[(int)floor(myDistance / stepDistance_)][(int)floor(myAngle/stepAngle_)] += 1;
            }
        }
    }
    
    for(int i = 0; i < maxDistance_ / stepDistance_; i++)
    {
        for(int j = 0; j < maxAngle_ / stepAngle_; j++)
        {
            dhProbability.selectDataSet(i);
            dhProbability.setPoint(j, tempProbability[i][j]);
            if (tempProbability[i][j] != 0)
            {
            //cout << tempProbability[i][j];
            }
        }
    }
    
    dhProbability.finishFrame();
    
    
}

void stacking::finishAnalysis(int /*nframes*/)
{
    
}

void stacking::writeOutput()
{
    
    real DistanceVector[(int)(maxDistance_ / stepDistance_)];
    for(int i = 0; i < (int)(maxDistance_ / stepDistance_) ; i++)
    {
        DistanceVector[i] = i * stepDistance_;
    }
    
    real AngleVector[(int)(maxAngle_ / stepAngle_)];
    
    for(int j = 0; j < (int)(maxAngle_ / stepAngle_); j++)
    {
        AngleVector[j] = j * stepAngle_;
    }
    
    
    if       (fnEnergySurface_.empty())         {fnEnergySurface_ = "energy.xpm";}
    else if  (fnEnergySurface_.compare(".xpm")) {}
    else                                        {fnEnergySurface_ += ".xpm";}
    
    // Construct matrix probability
    real **matProbability = new real*[(int)(maxDistance_ / stepDistance_)];
    for(int i = 0; i < (int)(maxDistance_ / stepDistance_); i++)
    {
        matProbability[i] = new real[(int)(maxAngle_ / stepAngle_)];
    }
    
    real **matEnergy = new real*[(int)(maxDistance_ / stepDistance_)];
    for(int i = 0; i < (int)(maxDistance_ / stepDistance_); i++)
    {
        matEnergy[i] = new real[(int)(maxAngle_ / stepAngle_)];
    }
    
    int scale = (int)(sel_.size()) * (int)(ref_.size());
    
    for(int i = 0 ; i < (int)(maxDistance_ / stepDistance_); i++)
    {
        for(int j = 0; j < (int)(maxAngle_ / stepAngle_); j++)
        {
            matProbability[i][j] = avem_probability_->average(i, j) ;
            //cout << i << '\t' << j << '\t' << avem_probability_->average(i, j) << '\t' << avem_probability_->average(i, j) / (float)scale << '\t' << matProbability[i][j] << endl;
            matEnergy[i][j]      = - 8.31447 * temperature_ * log(matProbability[i][j] / (float)scale) * 0.0002389;
        }
    }
    
    t_rgb rlo, rhi;
    rlo.r = 0.0; rlo.g = 0.0; rlo.b = 1.0;
    rhi.r = 1.0; rhi.g = 0.0; rhi.b = 0.0;
    int nlevels  = 400;
    
    FILE *fpEnergySurface;
    fpEnergySurface = fopen(fnEnergySurface_.c_str(), "w");
    
    
    
    write_xpm(fpEnergySurface, 0, "Free Energy Surface"
              , "Contact Probability", "Distance", "angle"
              , (int)(maxDistance_ / stepDistance_), (int)(maxAngle_ / stepAngle_), DistanceVector, AngleVector
              , matEnergy, 0, 20, rlo, rhi, &nlevels);
    
    fclose(fpEnergySurface);
    
    if (!fnEnergySurfaceRaw_.empty())
    {
        if (fnEnergySurfaceRaw_.compare(".dat"))    {}
        else                                        {fnEnergySurfaceRaw_ += ".dat";}
        
        FILE *fpEnergySurfaceRaw;
        fpEnergySurfaceRaw = fopen(fnEnergySurfaceRaw_.c_str(), "w" );
        
        for(int i = 0 ; i < (int)(maxDistance_ / stepDistance_); i++)
        {
            for(int j = 0; j < (int)(maxAngle_ / stepAngle_); j++)
            {
                fprintf(fpEnergySurfaceRaw, "%f ", matEnergy[i][j]);
            }
            fprintf(fpEnergySurfaceRaw, "\n");
        }
        fclose(fpEnergySurfaceRaw);
        
        
    }
    
    if (!fnPropability_.empty())
    {
        if (fnPropability_.compare(".dat"))    {}
        else                                        {fnPropability_ += ".dat";}
        
        FILE *fpProbability;
        fpProbability = fopen(fnPropability_.c_str(), "w" );
        
        for(int i = 0 ; i < (int)(maxDistance_ / stepDistance_); i++)
        {
            for(int j = 0; j < (int)(maxAngle_ / stepAngle_); j++)
            {
                fprintf(fpProbability, "%f ", matProbability[i][j]);
            }
            fprintf(fpProbability, "\n");
        }
        fclose(fpProbability);
        
    }
    
}


















