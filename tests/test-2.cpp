#include <vector>
#include <iostream>
#include <string>
#include <algorithm>

#include "../src/point-cloud-to-sac.hpp"

#include "../src/numerical-differentiation.hpp"

#include "../modelers/WigleyModeler.hpp"

#include "../src/utils/SimpleMonteCarlo.hpp"
#include "../src/utils/WriteToFile.hpp"
#include "../src/utils/TransposeVector.hpp"

int main(){
    /*
        Test first derivative of sectional are curve generated through `SectionalAreaXwiseYsymmetrical()` 
        at boundary points, by comparing it with the analytic version of the Wigley hull
    */

    /* general parameters */
    
    int Ncloud = 1000;//number of points at each end in point cloud

    double delta = 0.00000000001;//length of aft and bow section to generate point clouds

    double lambda = 1.0;//lambda to be used to generate sections

    std::vector<std::vector<double>> sections_aft;//sections to evaluate sectional-area at in the form of
        //sub-intervals (see 'SectionalAreaXwiseYsymmetrical()'), for aft end of geometry
    sections_aft = pc2sac::DiscreteBucketsRelative(-1.0,-1.0 + delta, 6, lambda);

    std::vector<std::vector<double>> sections_fwd;//sections to evaluate sectional-area at in the form of
        //sub-intervals (see 'SectionalAreaXwiseYsymmetrical()'), for aft end of geometry
    sections_fwd = pc2sac::DiscreteBucketsRelative(1.0-delta, 1.0, 6, lambda);

    std::string triangulation_filename = "test-2-wigley.stl";// analytic triangulation
    std::string cloud_filename_aft = "test-2-point-cloud-aft.txt";// point-cloud at aft
    std::string cloud_filename_fwd = "test-2-point-cloud-fwd.txt";// point-cloud at fwd
    std::string cloud_reduced_filename_aft = "test-2-point-cloud-reduced-aft.txt";// point-cloud reduced (put into buckets) at aft
    std::string cloud_reduced_filename_fwd = "test-2-point-cloud-reduced-fwd.txt";// point-cloud reduced (put into buckets) at fwd
    std::string analytic_areas_filename = "test-2-analytic.dat";// analytic sectional area curve
    std::string analytic_deriv_filename_aft = "test-2-analytic-deriv-aft.dat";// analytic section area curve tangent at aft
    std::string analytic_deriv_filename_fwd = "test-2-analytic-deriv-fwd.dat";// analytic section area curve tangent at fwd
    std::string approx_deriv_1_filename_aft = "test-2-approx-deriv-1-aft.dat";// 1 point approx derivative at aft
    std::string approx_deriv_5_filename_aft = "test-2-approx-deriv-5-aft.dat";// 5 point approx derivative at aft
    std::string approx_deriv_1_filename_fwd = "test-2-approx-deriv-1-fwd.dat";// 1 point approx derivative at fwd
    std::string approx_deriv_5_filename_fwd = "test-2-approx-deriv-5-fwd.dat";// 5 point approx derivative at fwd
    
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

    /* generate point clouds */

    std::vector<std::vector<double>> points_aft = SimpleMonteCarlo({{-1.0,-1.0+delta},{0.0,1.0}},Ncloud);
    std::vector<std::vector<double>> points_fwd = SimpleMonteCarlo({{1.0-delta,1.0},{0.0,1.0}},Ncloud);
    //map to surface
    for(int i = 0; i < Ncloud; i++){
        points_aft.at(i) = wigley.evaluate(points_aft.at(i));
        points_fwd.at(i) = wigley.evaluate(points_fwd.at(i));
    }
    WriteToFile(points_aft, cloud_filename_aft);
    WriteToFile(points_fwd, cloud_filename_fwd);
    //sort
    std::sort(points_aft.begin(),points_aft.end(),[](std::vector<double> a, std::vector<double> b){return a.at(0)<b.at(0);});
    std::sort(points_fwd.begin(),points_fwd.end(),[](std::vector<double> a, std::vector<double> b){return a.at(0)<b.at(0);});

    /* reduce point cloud to buckets and save in 'cloud_reduced_filename' */
    std::vector<std::vector<std::vector<double>>> points_in_buckets_aft = pc2sac::PointCloudToBucketsXwise(points_aft,sections_aft);
    std::vector<std::vector<std::vector<double>>> points_in_buckets_fwd = pc2sac::PointCloudToBucketsXwise(points_fwd,sections_fwd);
    std::vector<std::vector<double>> points_in_buckets_compact_aft;
    std::vector<std::vector<double>> points_in_buckets_compact_fwd;
    for(int i = 0; i < 6; i++){
        points_in_buckets_compact_aft.insert(points_in_buckets_compact_aft.end(),points_in_buckets_aft.at(i).begin(),points_in_buckets_aft.at(i).end());
        points_in_buckets_compact_fwd.insert(points_in_buckets_compact_fwd.end(),points_in_buckets_fwd.at(i).begin(),points_in_buckets_fwd.at(i).end());

    }
    WriteToFile(points_in_buckets_compact_aft, cloud_reduced_filename_aft);
    WriteToFile(points_in_buckets_compact_fwd, cloud_reduced_filename_fwd);

    /* evaluate analytic sectional area curve and save to file */
    std::vector<std::vector<double>> analytic_sac = std::vector<std::vector<double>>(2, std::vector<double>(1000,0.0));
    for(int i = 0; i < 1000; i++){
        analytic_sac.at(0).at(i) = -1.0 + 2.0*double(i)/999.0;//[-1,1]
        analytic_sac.at(1).at(i) = wigley.SectionalArea((analytic_sac.at(0).at(i)+1.0)/2.0);
    }
    WriteToFile(TransposeVector(analytic_sac), analytic_areas_filename);

    /* evaluate analytic tangents and save to file */

    double deriv_domain = 0.3;//evaluate tangent from point of tangency to +- deriv_domain
    double dSA_analytic_aft = wigley.dSectionalArea(0.0);
    double dSA_analytic_fwd = wigley.dSectionalArea(1.0);
    double SA_analytic_aft = wigley.SectionalArea(0.0);
    double SA_analytic_fwd = wigley.SectionalArea(1.0);
    double t = 0.0;
    std::vector<std::vector<double>> analytic_dsac_aft = std::vector<std::vector<double>>(1000, std::vector<double>(2,0.0));
    std::vector<std::vector<double>> analytic_dsac_fwd = std::vector<std::vector<double>>(1000, std::vector<double>(2,0.0));
    for(int i = 0; i < 1000; i++){
        t = -1.0 + deriv_domain*double(i)/999.0;
        analytic_dsac_aft.at(i) = {t, SA_analytic_aft + dSA_analytic_aft*(t+1.0)};
        t = 1.0 - deriv_domain*double(i)/999.0;
        analytic_dsac_fwd.at(i) = {t, SA_analytic_fwd + dSA_analytic_fwd*(t-1.0)};
    }
    WriteToFile(analytic_dsac_aft, analytic_deriv_filename_aft);
    WriteToFile(analytic_dsac_fwd, analytic_deriv_filename_fwd);

    /* evaluate sectional areas */

    std::vector<std::vector<double>> areas_aft = SectionalAreaXwiseYsymmetrical(points_aft, sections_aft);
    std::vector<std::vector<double>> areas_fwd = SectionalAreaXwiseYsymmetrical(points_fwd, sections_fwd);

    areas_aft = TransposeVector(areas_aft);
    areas_fwd = TransposeVector(areas_fwd);

    /* evaluate deriv tangents with 1 point and save to file */

    double dSA_1pnt_aft = Fwd1pointNumDiff(areas_aft.at(0).at(1), areas_aft.at(1).at(1), areas_aft.at(1).at(0) - areas_aft.at(0).at(0));
    double dSA_1pnt_fwd = Bck1pointNumDiff(areas_fwd.at(5).at(1), areas_fwd.at(4).at(1), areas_fwd.at(5).at(0) - areas_fwd.at(4).at(0));
    std::vector<std::vector<double>> approx_dsac_aft = std::vector<std::vector<double>>(1000, std::vector<double>(2,0.0));
    std::vector<std::vector<double>> approx_dsac_fwd = std::vector<std::vector<double>>(1000, std::vector<double>(2,0.0));
    for(int i = 0; i < 1000; i++){
        t = -1.0 + deriv_domain*double(i)/999.0;
        approx_dsac_aft.at(i) = {t, SA_analytic_aft + dSA_1pnt_aft*(t+1.0)};
        t = 1.0 - deriv_domain*double(i)/999.0;
        approx_dsac_fwd.at(i) = {t, SA_analytic_fwd + dSA_1pnt_fwd*(t-1.0)};
    }
    WriteToFile(approx_dsac_aft, approx_deriv_1_filename_aft);
    WriteToFile(approx_dsac_fwd, approx_deriv_1_filename_fwd);

    /* evaluate deriv tangents with 2 points */

    double dSA_2pnt_aft = Fwd2pointNumDiff(
        areas_aft.at(0).at(1),
        areas_aft.at(1).at(1),
        areas_aft.at(2).at(1),
        areas_aft.at(1).at(0) - areas_aft.at(0).at(0));
    double dSA_2pnt_fwd = Bck2pointNumDiff(
        areas_fwd.at(5).at(1),
        areas_fwd.at(4).at(1),
        areas_fwd.at(3).at(1),
        areas_fwd.at(5).at(0) - areas_fwd.at(4).at(0));

    /* evaluate deriv tangents with 3 points */

    double dSA_3pnt_aft = Fwd3pointNumDiff(
        areas_aft.at(0).at(1),
        areas_aft.at(1).at(1),
        areas_aft.at(2).at(1),
        areas_aft.at(3).at(1),
        areas_aft.at(1).at(0) - areas_aft.at(0).at(0));
    double dSA_3pnt_fwd = Bck3pointNumDiff(
        areas_fwd.at(5).at(1),
        areas_fwd.at(4).at(1),
        areas_fwd.at(3).at(1),
        areas_fwd.at(2).at(1),
        areas_fwd.at(5).at(0) - areas_fwd.at(4).at(0));

    /* evaluate deriv tangents with 4 points */

    double dSA_4pnt_aft = Fwd4pointNumDiff(
        areas_aft.at(0).at(1),
        areas_aft.at(1).at(1),
        areas_aft.at(2).at(1),
        areas_aft.at(3).at(1),
        areas_aft.at(4).at(1),
        areas_aft.at(1).at(0) - areas_aft.at(0).at(0));
    double dSA_4pnt_fwd = Bck4pointNumDiff(
        areas_fwd.at(5).at(1),
        areas_fwd.at(4).at(1),
        areas_fwd.at(3).at(1),
        areas_fwd.at(2).at(1),
        areas_fwd.at(1).at(1),
        areas_fwd.at(5).at(0) - areas_fwd.at(4).at(0));

    /* evaluate deriv tangents with 5 points and save to file */

    double dSA_5pnt_aft = Fwd5pointNumDiff(
        areas_aft.at(0).at(1),
        areas_aft.at(1).at(1),
        areas_aft.at(2).at(1),
        areas_aft.at(3).at(1),
        areas_aft.at(4).at(1),
        areas_aft.at(5).at(1), 
        areas_aft.at(1).at(0) - areas_aft.at(0).at(0));
    double dSA_5pnt_fwd = Bck5pointNumDiff(
        areas_fwd.at(5).at(1),
        areas_fwd.at(4).at(1),
        areas_fwd.at(3).at(1),
        areas_fwd.at(2).at(1),
        areas_fwd.at(1).at(1),
        areas_fwd.at(0).at(1),
        areas_fwd.at(5).at(0) - areas_fwd.at(4).at(0));
    approx_dsac_aft = std::vector<std::vector<double>>(1000, std::vector<double>(2,0.0));
    approx_dsac_fwd = std::vector<std::vector<double>>(1000, std::vector<double>(2,0.0));
    for(int i = 0; i < 1000; i++){
        t = -1.0 + deriv_domain*double(i)/999.0;
        approx_dsac_aft.at(i) = {t, SA_analytic_aft + dSA_5pnt_aft*(t+1.0)};
        t = 1.0 - deriv_domain*double(i)/999.0;
        approx_dsac_fwd.at(i) = {t, SA_analytic_fwd + dSA_5pnt_fwd*(t-1.0)};
    }
    WriteToFile(approx_dsac_aft, approx_deriv_5_filename_aft);
    WriteToFile(approx_dsac_fwd, approx_deriv_5_filename_fwd);

    /* print out error */
    std::cout << "Analytic Aft Derivative - 1 point Approx Aft Derivative: " << std::abs(dSA_analytic_aft - dSA_1pnt_aft)/std::abs(dSA_analytic_aft) << std::endl;
    std::cout << "Analytic Aft Derivative - 2 point Approx Aft Derivative: " << std::abs(dSA_analytic_aft - dSA_2pnt_aft)/std::abs(dSA_analytic_aft) << std::endl;
    std::cout << "Analytic Aft Derivative - 3 point Approx Aft Derivative: " << std::abs(dSA_analytic_aft - dSA_3pnt_aft)/std::abs(dSA_analytic_aft) << std::endl;
    std::cout << "Analytic Aft Derivative - 4 point Approx Aft Derivative: " << std::abs(dSA_analytic_aft - dSA_4pnt_aft)/std::abs(dSA_analytic_aft) << std::endl;
    std::cout << "Analytic Aft Derivative - 5 point Approx Aft Derivative: " << std::abs(dSA_analytic_aft - dSA_5pnt_aft)/std::abs(dSA_analytic_aft) << std::endl;
    std::cout << "Analytic Fwd Derivative - 1 point Approx Fwd Derivative: " << std::abs(dSA_analytic_fwd - dSA_1pnt_fwd)/std::abs(dSA_analytic_fwd) << std::endl;
    std::cout << "Analytic Fwd Derivative - 2 point Approx Fwd Derivative: " << std::abs(dSA_analytic_fwd - dSA_2pnt_fwd)/std::abs(dSA_analytic_fwd) << std::endl;
    std::cout << "Analytic Fwd Derivative - 3 point Approx Fwd Derivative: " << std::abs(dSA_analytic_fwd - dSA_3pnt_fwd)/std::abs(dSA_analytic_fwd) << std::endl;
    std::cout << "Analytic Fwd Derivative - 4 point Approx Fwd Derivative: " << std::abs(dSA_analytic_fwd - dSA_4pnt_fwd)/std::abs(dSA_analytic_fwd) << std::endl;
    std::cout << "Analytic Fwd Derivative - 5 point Approx Fwd Derivative: " << std::abs(dSA_analytic_fwd - dSA_5pnt_fwd)/std::abs(dSA_analytic_fwd) << std::endl;
    system("gnuplot test-2-gnuplot.sh");

    return 0;
}