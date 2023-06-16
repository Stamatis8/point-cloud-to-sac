#ifndef ON_SURFACE_POINT_CLOUD_TO_SAC_HPP
#define ON_SURFACE_POINT_CLOUD_TO_SAC_HPP

#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>

#include "on-surface-point-cloud-to-sac-utils.hpp"

/* prototypes */

template<typename integer = int>
std::vector<std::vector<double>> SectionalAreaXwiseYsymmetrical(std::vector<std::vector<double>> points, integer N, integer M);

template<typename integer = int>
std::vector<std::vector<double>> SectionalAreaXwiseYsymmetrical(std::vector<std::vector<double>> points, std::vector<std::vector<double>> sections);

template<typename integer = int>
std::vector<std::vector<double>> SectionalAreaXwiseYsymmetrical(std::vector<std::vector<double>> points, integer N, double lambda);

/* functions */

template<typename integer = int>
std::vector<std::vector<double>> SectionalAreaXwiseYsymmetrical(std::vector<std::vector<double>> points, integer N, double lambda){
    
    /* sort points along X coordinate and calculate Xmin, Xmax */

    std::sort(points.begin(), points.end(), [](std::vector<double> a, std::vector<double> b){
        if(a.at(0) < b.at(0)) return true;
        return false;
    });

    double Xmin = points.at(0).at(0);
    double Xmax = points.back().at(0);

    /* compose functions to make result */

    return SectionalAreaXwiseYsymmetrical(points, pc2sac::DiscreteBucketsRelative(Xmin,Xmax,N,lambda));
}

template<typename integer = int>
std::vector<std::vector<double>> SectionalAreaXwiseYsymmetrical(std::vector<std::vector<double>> points, std::vector<std::vector<double>> sections){

    /*
        points assumed to be sorted X-wise
    */

    /* create N buckets of points */

    integer N = sections.size();//number of sections

    //initialize X buckets
    std::vector<std::vector<std::vector<double>>> Xbuckets = std::vector<std::vector<std::vector<double>>>(N, std::vector<std::vector<double>>());

    //fill buckets
    integer count = 0;//count of points
    integer Np = points.size();//points size
    for(integer i = 0; i < N; i++){//Xbuckets.at(i) contains points from x.at(i) to x.at(i+1)
        
        if(points.at(count).at(0) < sections.at(i).at(0)){//skip points not belonging to any section
            i--;
            count++;
            continue;
        }

        while(count < Np && points.at(count).at(0) < sections.at(i).at(1) && points.at(count).at(0) >= sections.at(i).at(0)){//add point to bucket
            Xbuckets.at(i).push_back(points.at(count));
            count++;
        }

        // check if any bucket is empty and throw error 
        if(Xbuckets.at(i).size() == 0){
            std::cout << "Insufficient number of points: could not fill all buckets along X" << std::endl;
            throw(0.0);
        }
    }

    /* project points in the same bucket to YZ plane at their midpoint and sort Z-wise */
    for(integer i = 0; i < N; i++){

        //project
        for(integer j = 0; j < Xbuckets.at(i).size(); j++){//project
            Xbuckets.at(i).at(j).at(0) = (sections.at(i).at(0) + sections.at(i).at(1))/2.0;
        }

        //sort
        std::sort(Xbuckets.at(i).begin(), Xbuckets.at(i).end(), [](std::vector<double> a, std::vector<double> b){
            if (a.at(2) < b.at(2)) return true;
            return false;
        });

    }
    
    auto area = [](std::vector<double> a, std::vector<double> b){
        //calculate area created by line between a, b
        //  if projected on Z plane along (0,-1,0) or (0,1,0)
        return std::sqrt(std::pow(a.at(0)-b.at(0),2.0)+std::pow(a.at(2)-b.at(2),2.0))*(std::abs(b.at(1))+(a.at(1)-b.at(1))/2.0);

    };

    /* using 'Xbuckets' calculate areas */
    std::vector<std::vector<double>> areas = std::vector<std::vector<double>>(2, std::vector<double>(N, 0.0));
    for(integer i = 0; i < N; i++){
        areas.at(0).at(i) = Xbuckets.at(i).at(0).at(0);//longitudinal position of section
        for(integer j = 0; j < (Xbuckets.at(i).size()-1); j++){
            areas.at(1).at(i) += area(Xbuckets.at(i).at(j), Xbuckets.at(i).at(j+1)); 
        }
        areas.at(1).at(i) *= 2.0;
    }

    return areas;
}

template<typename integer = int>
std::vector<std::vector<double>> SectionalAreaXwiseYsymmetrical(std::vector<std::vector<double>> points, integer N, integer M){

    /* sort points according to X-coordinate */

    std::sort(points.begin(), points.end(), [](std::vector<double> a, std::vector<double> b){
        if(a.at(0) < b.at(0)) return true;
        return false;
    });

    /* determine Xmin, Xmax and generate xi, x0 = Xmin, xN = Xmax */

    double Xmin = points.at(0).at(0);
    double Xmax = points.back().at(0);
    std::vector<double> x = std::vector<double>(N+1,0.0);
    for(integer i = 0; i < N+1; i++){
        x.at(i) = Xmin + (Xmax-Xmin)*double(i)/double(N);
    }

    /* create 'points_in_buckets' (see above) */

    //initialize points_in_buckets (see above)
    std::vector<std::vector<std::vector<double>>> points_in_buckets = pc2sac::PointCloudToBucketsXwise<integer>(points, N, M);
    

    auto area = [](std::vector<double> a, std::vector<double> b){
        //calculate area created by line between a, b
        //  if projected on Z plane along (0,-1,0) or (0,1,0)
        return std::sqrt(std::pow(a.at(0)-b.at(0),2.0)+std::pow(a.at(2)-b.at(2),2.0))*(std::abs(b.at(1))+(a.at(1)-b.at(1))/2.0);

    };

    /* using 'points_in_buckets' calculate areas */
    std::vector<std::vector<double>> areas = std::vector<std::vector<double>>(2, std::vector<double>(N, 0.0));
    for(integer i = 0; i < N; i++){
        areas.at(0).at(i) = x.at(i)/2.0 + x.at(i+1)/2.0;
        for(integer j = 0; j < (M-1); j++){
            areas.at(1).at(i) += area(points_in_buckets.at(i).at(j), points_in_buckets.at(i).at(j+1)); 
        }
    }

    return areas;
}

#endif //ON_SURFACE_POINT_CLOUD_TO_SAC_HPP