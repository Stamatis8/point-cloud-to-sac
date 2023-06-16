#include <vector>
#include <iostream>
#include <string>
#include <algorithm>

#include "../src/point-cloud-to-sac.hpp"

#include "../modelers/WigleyModeler.hpp"

#include "../src/utils/SimpleMonteCarlo.hpp"
#include "../src/utils/WriteToFile.hpp"
#include "../src/utils/TransposeVector.hpp"

int main(){
    /*
        Verify validity of `SectionalAreaXwiseYsymmetrical()` by comparing it with the 
        analytic sectional area curve of the Wigley hull
    */

    /* general parameters */
    
    int N = 10;//number of point to evaluate SAC at

    double lambda = 0.5;//subinterval length percentage of standard uniform subinterval

    int Ncloud = 3000;//number of points in point cloud

    std::string triangulation_filename = "test-1-wigley.stl";// analytic triangulation
    std::string cloud_filename = "test-1-point-cloud.txt";// point-cloud
    std::string cloud_reduced_filename = "test-1-point-cloud-reduced.txt";// point-cloud reduced (put into buckets)
    std::string analytic_areas_filename = "test-1-analytic.dat";// analytic sectional area curve
    std::string approx_areas_filename = "test-1-approx.dat";// approximate sectional area curve

    /* initialize modeller */

    double L = 2.0;
    double B = 0.5;
    double T = 0.2;
    double c1 = 4.3;
    double c2 = 1.3;
    double c3 = 2.3;
    WigleyModeler wigley {{0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}};
    wigley.set_design({L, B, T, c1, c2, c3});
    //triangulate
    wigley.triangulate(triangulation_filename,3000);

    /* generate point cloud */
    std::vector<std::vector<double>> points = SimpleMonteCarlo({{-1.0,1.0},{0.0,1.0}},Ncloud);
    //map to surface
    for(int i = 0; i < Ncloud; i++){
        points.at(i) = wigley.evaluate(points.at(i));
    }
    WriteToFile(points, cloud_filename);

    /* reduce point cloud to buckets and save in 'cloud_reduced_filename' */
    std::vector<std::vector<double>> sections = pc2sac::DiscreteBucketsRelative(-1.0,1.0,N,lambda);
    std::vector<std::vector<std::vector<double>>> points_in_buckets = pc2sac::PointCloudToBucketsXwise(points,sections);
    std::vector<std::vector<double>> points_in_buckets_compact;
    for(int i = 0; i < N; i++){
        points_in_buckets_compact.insert(points_in_buckets_compact.end(),points_in_buckets.at(i).begin(),points_in_buckets.at(i).end());
    }
    WriteToFile(points_in_buckets_compact, cloud_reduced_filename);

    /* evaluate approx sectional area curve and save to file */
    std::vector<std::vector<double>> approx_sac = SectionalAreaXwiseYsymmetrical(points, N, lambda);
    WriteToFile(TransposeVector(approx_sac), approx_areas_filename);

    /* evaluate analytic sectional area curve and save to file */
    std::vector<std::vector<double>> analytic_sac = std::vector<std::vector<double>>(2, std::vector<double>(1000,0.0));
    for(int i = 0; i < 1000; i++){
        analytic_sac.at(0).at(i) = -1.0 + 2.0*double(i)/999.0;//[-1,1]
        analytic_sac.at(1).at(i) = wigley.SectionalArea((analytic_sac.at(0).at(i)+1.0)/2.0);
    }
    WriteToFile(TransposeVector(analytic_sac), analytic_areas_filename);
    analytic_sac = std::vector<std::vector<double>>(2, std::vector<double>(N,0.0));
    for(int i = 0; i < N; i++){
        analytic_sac.at(0).at(i) = approx_sac.at(0).at(i);//[-1,1]
        analytic_sac.at(1).at(i) = wigley.SectionalArea((analytic_sac.at(0).at(i)+1.0)/2.0);
    }

    /* evaluate error and print */
    std::vector<double> error = std::vector<double>(N, 0.0);
    for(int i = 0; i < N; i++){
        error.at(i) = std::abs(analytic_sac.at(1).at(i) - approx_sac.at(1).at(i));
    }
    std::cout << "A point cloud of " << Ncloud << " points was used to approximate the sectional area curve of the Wigley hull";
    std::cout << " with an error of " << *std::max_element(error.begin(),error.end()) << std::endl;

    system("gnuplot test-1-gnuplot.sh");

    return 0;
}