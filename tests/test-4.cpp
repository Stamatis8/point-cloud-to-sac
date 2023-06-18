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
        Evaluate the derivatives of the cross-sectional area of the KCSsim hull (see `KCSsimModeler.hpp`) at the aft and at the bow using `SectionalAreaXwiseYsymmetrical()` and 
        `numerical-differentiation.hpp`.
    */

    /* general parameters sac calculation */
    
    int N = 101;//number of point to evaluate SAC at (make it odd for simpson)

    int Ncloud = 100000;//number of points in point cloud

    double lambda = 1.0;//subinterval length percentage of standard uniform subinterval

    int Ntriangles = 100000;//number of triangles

    std::string triangulation_filename = "test-4-KCSsim.stl";// analytic triangulation
    std::string cloud_filename = "test-4-point-cloud.txt";// point-cloud
    std::string cloud_reduced_filename = "test-4-point-cloud-reduced.txt";// point-cloud reduced (put into buckets)
    std::string approx_areas_filename = "test-4-approx.dat";// approximate sectional area curve
    std::string cloud_filename_aft = "test-4-point-cloud-aft.txt";// point-cloud
    std::string cloud_filename_fwd = "test-4-point-cloud-fwd.txt";// point-cloud
    std::string approx_deriv_1_filename_aft = "test-4-approx-deriv-1-aft.dat";// 1 point approx derivative at aft
    std::string approx_deriv_1_filename_fwd = "test-4-approx-deriv-1-fwd.dat";// 1 point approx derivative at fwd

    /* initialize modeller */

    std::string matlab_executable_directory = "";
    std::vector<std::vector<double>> design {{1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0},{0.0, 1.0}};
    KCSsimModeler hull {design, double(Ntriangles)};
    srand(110);
    std::cout<<SimpleMonteCarlo(std::vector<std::vector<double>>(design.begin()+1,design.end()),1).at(0).at(0)<<std::endl;
    return 0;
    hull.set_design(SimpleMonteCarlo(std::vector<std::vector<double>>(design.begin()+1,design.end()),1).at(0),matlab_executable_directory);
    //save triangulation
    smpl_triangulation::make_stl(hull.get_triangulation(), triangulation_filename);

    // std::string matlab_executable_directory = "";
    // std::vector<std::vector<double>> design {
    //     {0.0,1.0},//c1
    //     {0.8},//c2
    //     {0.8},//c3
    //     {0.265},//t1 c4
    //     {0.543},//t2 c5
    //     {0.0},//t3 c6
    //     {0.3},//t4 c7
    //     {0.65},//t5 c8
    //     {0.7},//t6 c9
    //     {0.5},//t7 c10
    //     {0.5},//t8 c11         
    //     {0.5},//t9 c12  
    //     {0.5},//t10 c13      
    //     {0.5},//t11 c14
    //     {0.5},//t12 c15
    //     {0.05},//t13 c16  
    //     {0.5},//t14 c17  
    //     {0.7},//t15 c18     
    //     {0.8},//t16 c19     
    //     {1.0},//t17 c20      
    //     {1.0},//t18 c21
    //     {0.0},//t19 c22
    //     {0.5},//t20 c23
    //     {0.5},//t21 c24
    //     {0.5},//t22 c25
    //     {0.5},//t23 c26
    //     {0.5},//t24 c27
    //     {0.5},//t25 c28
    //     {0.5}//t26 c29
    // };
    // KCSsimModeler hull {design, double(Ntriangles)};
    // hull.set_design({1.0},matlab_executable_directory);
    // //save triangulation
    // smpl_triangulation::make_stl(hull.get_triangulation(), triangulation_filename);

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

    /* evaluate aft point cloud */

    int NcloudEnds = 5000;
    std::vector<std::vector<double>> points_fwd = grid({{0.0,1.0},{0.0,0.1}}, std::ceil(std::sqrt((7.0)*double(NcloudEnds))), std::ceil(std::sqrt((1.0/7.0)*double(NcloudEnds))));
    points_fwd = hull.evaluate(&points_fwd, "");
    WriteToFile(points_fwd ,cloud_filename_fwd);
    std::vector<std::vector<double>> points_aft = grid({{0.0,1.0},{0.9,1.0}}, std::ceil(std::sqrt((7.0)*double(NcloudEnds))), std::ceil(std::sqrt((1.0/7.0)*double(NcloudEnds))));
    points_aft = hull.evaluate(&points_aft, "");
    WriteToFile(points_aft ,cloud_filename_aft);

    /* evaluate deriv tangents with 1 point and save to file */
    std::sort(points_fwd.begin(),points_fwd.end(),[](std::vector<double> a, std::vector<double> b){return a.at(0)<b.at(0);});
    double bulbous_end_location = points_fwd.back().at(0);
    double deriv_domain = 0.05;//evaluate tangent from point of tangency to +- deriv_domain
    double dSA_1pnt_aft = hull.dSectionalArea(0.0);
    double dSA_1pnt_fwd = hull.dSectionalArea(1.0);
    std::vector<std::vector<double>> approx_dsac_aft = std::vector<std::vector<double>>(1000, std::vector<double>(2,0.0));
    std::vector<std::vector<double>> approx_dsac_fwd = std::vector<std::vector<double>>(1000, std::vector<double>(2,0.0));
    double t;
    for(int i = 0; i < 1000; i++){
        t = -1.0 + deriv_domain*double(i)/999.0;
        approx_dsac_aft.at(i) = {t, hull.SectionalArea(0.0) + dSA_1pnt_aft*(t+1.0)};
        t = bulbous_end_location - deriv_domain*double(i)/999.0;
        approx_dsac_fwd.at(i) = {t, hull.SectionalArea(1.0) + dSA_1pnt_fwd*(t-bulbous_end_location)};
    }
    WriteToFile(approx_dsac_aft, approx_deriv_1_filename_aft);
    WriteToFile(approx_dsac_fwd, approx_deriv_1_filename_fwd);

    system("gnuplot test-4-gnuplot.sh");

    return 0;
}