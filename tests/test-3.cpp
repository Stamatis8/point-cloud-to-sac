#include <vector>
#include <iostream>
#include <string>
#include <algorithm>

#include "../src/point-cloud-to-sac.hpp"

#include "../modelers/KCSsimModeler.hpp"

#include "../src/smpl-triangulation/make_stl.hpp"

#include "../src/utils/SimpleMonteCarlo.hpp"
#include "../src/utils/WriteToFile.hpp"
#include "../src/utils/ReadFile.hpp"
#include "../src/utils/TransposeVector.hpp"
#include "../src/utils/simpson.hpp"
#include "../src/utils/grid.hpp"

int main(){
    /*
        Evaluate cross sectional area curve for KCSsim hull (see `KCSsimModeler.hpp`) using 
        `SectionalAreaXwiseYsymmetrical()`. Test for validity by integating and comparing with 
        model volume
    */

    /* general parameters */
    
    int N = 51;//number of point to evaluate SAC at (make it odd for simpson)

    int Ncloud = 30000;//number of points in point cloud

    double lambda = 1.0;//subinterval length percentage of standard uniform subinterval

    int Ntriangles = 50000;//number of triangles

    std::string triangulation_filename = "test-3-KCSsim.stl";// analytic triangulation
    std::string cloud_filename = "test-3-point-cloud.txt";// point-cloud
    std::string cloud_reduced_filename = "test-3-point-cloud-reduced.txt";// point-cloud reduced (put into buckets)
    std::string approx_areas_filename = "test-3-approx.dat";// approximate sectional area curve

    /* initialize modeller */

    std::string filename = "KCSsim.stl";
    std::string matlab_executable_directory = "";
    std::vector<std::vector<double>> design {
        {0.0,1.0},
        {0.8},
        {0.8},
        {0.265},
        {0.543},
        {0.0},
        {0.3},
        {0.65}, 
        {0.7},
        {0.5},
        {0.5},          
        {0.5},    
        {0.5},      
        {0.5},
        {0.5},
        {0.05},  
        {0.5},   
        {0.7},      
        {0.8},      
        {1.0},       
        {1.0},   
        {0.0},
        {0.5},
        {0.5},
        {0.5},
        {0.5},
        {0.5},
        {0.5},
        {0.5}
    };
    KCSsimModeler hull {design, double(Ntriangles)};
    hull.set_design({1.0},matlab_executable_directory);
    //save triangulation
    smpl_triangulation::make_stl(hull.get_triangulation(), triangulation_filename);

    /* generate point cloud */
    //std::vector<std::vector<double>> points = ReadFile("surface_points.txt");
    //std::vector<std::vector<double>> points = SimpleMonteCarlo({{0,1},{0,1}}, Ncloud);
    std::vector<std::vector<double>> points = grid({{0,1},{0,1}}, std::ceil(std::sqrt(double(Ncloud))), std::ceil(std::sqrt(double(Ncloud))));
    points = hull.evaluate(&points, "");
    WriteToFile(points, cloud_filename);

    /* reduce point cloud to buckets and save in 'cloud_reduced_filename' */
    std::vector<std::vector<double>> sections = pc2sac::DiscreteBucketsRelative(-1.0,0.0,N,lambda);
    std::vector<std::vector<std::vector<double>>> points_in_buckets = pc2sac::PointCloudToBucketsXwise(points,sections);
    std::vector<std::vector<double>> points_in_buckets_compact;
    for(int i = 0; i < N; i++){
        points_in_buckets_compact.insert(points_in_buckets_compact.end(),points_in_buckets.at(i).begin(),points_in_buckets.at(i).end());
    }
    WriteToFile(points_in_buckets_compact, cloud_reduced_filename);

    /* evaluate approx sectional area curve and save to file */
    std::vector<std::vector<double>> approx_sac = SectionalAreaXwiseYsymmetrical(points, N, lambda);
    WriteToFile(TransposeVector(approx_sac), approx_areas_filename);

    /* evaluate approx. and analytical volume and print */
    double approx_volume = simpson(approx_sac.at(1),approx_sac.at(0));
    double accurate_volume = hull.moment(0,0,0,false,false);
    std::cout << "Volume of mesh with " << Ntriangles << " triangles: " << accurate_volume << std::endl;
    std::cout << "Volume as integral of SAC: " << approx_volume << std::endl;
    std::cout << "Relative difference: " << std::abs(accurate_volume-approx_volume)/std::abs(accurate_volume) << std::endl;

    system("gnuplot test-3-gnuplot.sh");

    return 0;
}