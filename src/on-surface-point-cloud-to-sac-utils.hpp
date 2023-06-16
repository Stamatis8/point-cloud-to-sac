#ifndef ON_SURFACE_POINT_CLOUD_TO_SAC_UTILS_HPP
#define ON_SURFACE_POINT_CLOUD_TO_SAC_UTILS_HPP

#include <vector>
#include <iostream>
#include <algorithm>

namespace pc2sac{

/* prototypes */

template<typename integer = int>
std::vector<std::vector<std::vector<double>>> PointCloudToBucketsXwise(std::vector<std::vector<double>> points, integer N, integer M);

template<typename integer = int>
std::vector<std::vector<double>> DiscreteBucketsRelative(double a, double b, integer N, double lambda);


/* functions */

template<typename integer = int>
std::vector<std::vector<std::vector<double>>> PointCloudToBucketsXwise(std::vector<std::vector<double>> points, integer N, integer M){

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

    /* create N buckets of points */

    //initialize X buckets
    std::vector<std::vector<std::vector<double>>> Xbuckets = std::vector<std::vector<std::vector<double>>>(N, std::vector<std::vector<double>>());

    //fill buckets
    integer count = 0;//count of points
    for(integer i = 0; i < N; i++){//Xbuckets.at(i) contains points from x.at(i) to x.at(i+1)
        
        while(points.at(count).at(0) < x.at(i+1)){
            Xbuckets.at(i).push_back(points.at(count));
            count++;
        }

        // check if any bucket is empty and throw error 
        if(Xbuckets.at(i).size() == 0){
            std::cout << "Insufficient number of points: could not fill all buckets along X" << std::endl;
            throw(0.0);
        }
    }

    /* loop over buckets and create 'points_in_buckets' (see above) */

    //initialize points_in_buckets (see above)
    std::vector<std::vector<std::vector<double>>> points_in_buckets =
        std::vector<std::vector<std::vector<double>>>(N, std::vector<std::vector<double>>(M, std::vector<double>()));
    
    //initialize Z buckets
    std::vector<std::vector<std::vector<double>>> Zbuckets = std::vector<std::vector<std::vector<double>>>(M, std::vector<std::vector<double>>());

    //generate z with M+1 elements
    std::vector<double> z = std::vector<double>(M+1,0.0);

    for(integer i = 0; i < N; i++){

        /* sort Xbucket.at(i) according to Z-coordinate */

        std::sort(Xbuckets.at(i).begin(), Xbuckets.at(i).end(), [](std::vector<double> a, std::vector<double> b){
            if(a.at(2) < b.at(2)) return true;
            return false;
        });

        /* determine Zmin, Zmax and generate zi, z0 = Zmin, zM = Zmax */

        double Zmin = Xbuckets.at(i).at(0).at(2);
        double Zmax = Xbuckets.at(i).back().at(2);
        for(integer j = 0; j < M+1; j++){
            z.at(j) = Zmin + (Zmax-Zmin)*double(j)/double(M);
        }

        /* generate M buckets of points */

        count = 0;//count of points
        for(integer j = 0; j < M; j++){//Zbuckets.at(i) contains points from x.at(i) to t.at(i+1)
            
            Zbuckets.at(j).clear();
            
            while(Xbuckets.at(i).at(count).at(2) < z.at(j+1)){
                Zbuckets.at(j).push_back(Xbuckets.at(i).at(count));
                count++;
            }
        
            //throw error if any Zbucket is empty
            if(Zbuckets.at(j).size()==0){
                std::cout << "Insufficient number of points: could not fill all buckets along Z" << std::endl;
                throw(0.0);
            }
        }

        /* for each Zbucket save only element which is closest to mid */
        //  (i.e. closest to (x.at(i)/2.0 + x.at(i+1)/2.0), anything, z.at(j)/2.0 + z.at(j+1)/2.0)
        for(integer j = 0; j < M; j++){

            auto mid_metric = [i,j,x,z](std::vector<double> a, std::vector<double> b){
                if(
                    (std::abs(a.at(0) - x.at(i)/2.0 - x.at(i+1)/2.0) + std::abs(a.at(2) - z.at(j)/2.0 - z.at(j+1)/2.0)) < 
                    (std::abs(b.at(0) - x.at(i)/2.0 - x.at(i+1)/2.0) + std::abs(b.at(2) - z.at(j)/2.0 - z.at(j+1)/2.0))
                ) return true;
                return false;
            };

            points_in_buckets.at(i).at(j) = *std::min_element(Zbuckets.at(j).begin(),Zbuckets.at(j).end(), mid_metric);
        }
    }

    return points_in_buckets;
}

template<typename integer = int>
std::vector<std::vector<std::vector<double>>> PointCloudToBucketsXwise(std::vector<std::vector<double>> points, std::vector<std::vector<double>> sections){
    
    /* sort points according to X-coordinate */

    std::sort(points.begin(), points.end(), [](std::vector<double> a, std::vector<double> b){
        if(a.at(0) < b.at(0)) return true;
        return false;
    });

    /* create N buckets of points */

    integer N = sections.size();//number of sections

    //initialize X buckets
    std::vector<std::vector<std::vector<double>>> Xbuckets = std::vector<std::vector<std::vector<double>>>(N, std::vector<std::vector<double>>());

    //fill buckets
    integer count = 0;//count of points
    integer Np = points.size();//size of points
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

    return Xbuckets;
}

template<typename integer = int>
std::vector<std::vector<double>> DiscreteBucketsRelative(double a, double b, integer N, double lambda){
    /*
        Description: discretizes [a,b] to N uniform subintervals, each subinterval of length lambda*(a-b)/(N-1) and centered at
            t_i = a + (b-a)*i/(N-1) for i = 0,..., N-1. The boundary subintevals are centered about t_0 = a and t_{N-1} = b
            so their length will be lambda*(a-b)/N / 2.0
        Output:
            std::vector<std::vector<double>> discretization
                - [discretization.at(i).at(0),discretization.at(i).at(1)] is the ith subinterval
        Example:
            DiscreteBucketsRelative(0.0, 1.0, 5, 0.8) = {
                {0.0,0.1},
                {0.15,0.35},
                {0.4,0.6},
                {0.65,0.85}
                {0.9,1.0}
            }
        Notes:
            N > 1
            lambda < 1.0 will ensure that the subintervals are mutually distinct
    */
    
    //subinterval length
    double l = lambda*(b-a)/double(N-1);

    //initialize result
    std::vector<std::vector<double>> t = std::vector<std::vector<double>>(N, std::vector<double>(2,0.0));

    for(int i = 0; i < N; i++){
        t.at(i).at(0) = (a+(b-a)*double(i)/double(N-1)) - l/2.0;
        t.at(i).at(1) = t.at(i).at(0) + l;
    }

    //amend endpoints
    t.at(0).at(0) = a;
    t.back().at(1) = b;

    return t;
}


}//pc2sac namespace
#endif //ON_SURFACE_POINT_CLOUD_TO_SAC_UTILS_HPP