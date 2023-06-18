#ifndef KCSSIMMODELER_HPP
#define KCSSIMMODELER_HPP

#include <vector>

//triangulations utility
#include "../src/smpl-triangulation/smpl_triangulation.hpp"

//triangulation geometric moment calculator
#include "../src/geom_moments/geom_moments.hpp"

//FileIO
#include "../src/utils/WriteToFile.hpp"
#include "../src/utils/ReadFile.hpp"

class KCSsimModeler{
    /*
        This is a 29-parameter parametric modeler which is similar to the KCS hull [1]. For
            more details see [2]

        [1] http://www.simman2008.dk/KCS/container.html
        [2] "Geometric Moment-Dependent Global Sensitivity Analysis without Simulation Data: Application to Ship Hull Form Optimisation"
					Shahroz Khan, Panagiotis Kaklis, Andrea Serani, Matteo Diez, 2022

        
    */
public:

    /* constructor */

    KCSsimModeler(
        std::vector<std::vector<double>> D,
        double N
    ){
        /*
            This constructor can be used to create instances of the parametric modeler in less than
            29 parameters (i.e. specify a fixed value for some of the parameters and a range 
            for the rest)
            
            D is the input design space and fixed parameters. 

            N is the number of vertices (approximately) that the geometry triangulation 
            should have

            D must still be of length 29. However, for each fixed ith parameter the size of
            D.at(i) MUST be equal to 1 and have said fixed value. For any other jth non-fixed
            parameter its range will be [D.at(j).at(0),D.at(j).at(1)] and the size of D.at(j) 
            MUST be 2

            So, for example, if the user wants to fix the first and 14th parameters at 0.5, but have 
            all others between 0 and 1,

                D.at(0) = { 0.5 }
                D.at(13) = { 0.5 }
                D.at(j) = { 0.0, 1.0 } for j = 1...12 14...29
        */

       /* create design space and set fixed parameters */

        this->fixed_parameters_mask = std::vector<bool>(29,false);//assume initially that all parameters are non-fixed
        this->design = std::vector<double>(29,0);

        for(int i = 0; i < D.size(); i++){
            if(D.at(i).size() == 1)
            {//ith parameter is fixed
                this->fixed_parameters_mask.at(i) = true;
                this->design.at(i) = D.at(i).at(0);
            }
            else if (D.at(i).size() == 2)
            {//ith parameter is not fixed
                this->design_space_attribute.push_back(D.at(i));
            }
        }

        /* create triangulation preimage */

        //surface domain is set to be [-1,1]\times[0,1] with [0,1]\times[0,1] referring to the
        //  port side of the hull and [-1,0]\times[0,1] referring to the starboard side of the 
        //  hull
    
        this->triangulation_preimage = smpl_triangulation::PlanarTriangulationStitchLeftRight({-1.0,1.0},{0.0,1.0}, N);
        
        //constructing 'triangulation_preimage_nodes_in_surface_domain'
        this->triangulation_preimage_nodes_in_surface_domain = this->triangulation_preimage.nodes;
        for(int i = 0; i < this->triangulation_preimage.nodes.size(); i++){
            if(this->triangulation_preimage_nodes_in_surface_domain.at(i).at(0) < 0.0){
                this->triangulation_preimage_nodes_in_surface_domain.at(i).at(0) *= -1.0;
            }
        }
    };

    /* Set Design */

    void set_design(std::vector<double> design_in, const std::string matlab_executable_directory = ""){

        //matlab_executable_directory must be in backslashes. For example
        //   matlab_executable_directory = "folder1\\folder1_1\\folder_1_1_1\\"

        /* update this->design */

        //design_in.size() must be equal to this->design_space_attribute.size()
        int j = 0;//design_in index
        for(int i = 0; i < 29; i++){
            if(!this->fixed_parameters_mask.at(i)){//ith parameter is not fixed
                this->design.at(i) = design_in.at(j);
                j += 1;
            }
        }

        /* evaluate triangulation */

        //generate dummy triangulation as copy of planar triangulation
        smpl_triangulation::Triangulation<std::vector<double>> dummy_triangulation = this->triangulation_preimage;
        
        //map onto surface
        dummy_triangulation.nodes = this->evaluate(&this->triangulation_preimage_nodes_in_surface_domain, matlab_executable_directory);

        //translate vertices appropriately to cover entire hull
        for(int i = 0; i < dummy_triangulation.nodes.size(); i++){
            if(this->triangulation_preimage.nodes.at(i).at(0) < 0.0){
                //mirror on-surface vertex to starboard side
                dummy_triangulation.nodes.at(i).at(1) *= -1.0;
            }
        }

        //write to this->triangulation
        this->triangulation = dummy_triangulation;

        /* SAaft, SAbow, dSAaft, dSAbow  */
        
    }

    /* Get Design Space */
	
	std::vector<std::vector<double>> design_space() const{
		return this->design_space_attribute;
	}

    /* Get Particulars */

	std::vector<double> get_particulars() const{
		return this->design;
	}

    /* Moment Evaluator */

    double moment(int p, int q, int r, bool is_translation_invariant = false, bool is_scaling_invariant = false) const{
		/*
			Description:
				  - Calculates the s = p + q + r order geometric moment of the current design in modeler.
				  - is_translation_invariant  == true calculates the translation invariant of said moment
				  - is_scaling_invariant == true calculates scaling invariant of said moment
		*/
        smpl_triangulation::Triangulation<std::vector<double>> T = this->triangulation;//passing by copy
		return MomentSthOrder(T,p,q,r,p+q+r,is_translation_invariant,is_scaling_invariant);
	}

    double SectionalArea(double x) const{
        /*
            Description: returns cross sectional area from x = 0.0 (aft) to x = 1.0 (bow)
        */

        return 0.0;
    }

    double dSectionalArea(double x) const{
        /*
            Description: returns derivative of cross sectional area curve from x = 0.0 (aft) to x = 1.0 (bow)
        */
        
        return 0.0;
    }

    smpl_triangulation::Triangulation<std::vector<double>> get_triangulation(){
        return this->triangulation;
    }

    std::vector<std::vector<double>> evaluate(
        std::vector<std::vector<double>> * parametric_points, 
        std::string matlab_executable_directory) const{
        /*
            Description: pass a pointer to a number of list of parametric points, and map them onto the surface. 

            Input:
                - std::vector<std::vector<double>> * parametric_points
                    - pointer to parametric_points
                    - all parametric_points must be in [0,1]^2
                    - v parametric-direction is longitudinal-direction, with (0,0) -> aft and (0,1) -> bow
                - std::string matlab_executable_directory
                    - directory were the relevant matlab executable is located at
                    - all auxiliary files (design.txt, planar_points.txt, surface_points.txt) are created and read
                        from this directory     
                    - directory must by in backslahes, example 
                        matlab_executable_directory = "folder1\\folder1_1\\folder_1_1_1\\"
            Output:
                - std::vector<std::vector<double>> surface_points
                    - parametric_points mapped onto surface
                    - they are all mapped to one side of the hull
        */

        //write design to a file
        WriteToFile(this->design,std::string(matlab_executable_directory+"design.txt").c_str(),false);

        //write all preimage vertices to a file
        WriteToFile(*parametric_points,std::string(matlab_executable_directory+"planar_points.txt").c_str());

        //create new file with preimage points mapped to on-surface points through matlab script
        //system("matlab -batch \"run(\'src/modelers/KCSsimModelerEvaluate.m\')\"");
        system(std::string(matlab_executable_directory+"KCSsimModelerEvaluateLocal.exe").c_str());
    
        return ReadFile(std::string(matlab_executable_directory+"surface_points.txt").c_str());

    }

protected:

    std::vector<std::vector<double>> design_space_attribute;// ith parameter \in [design_space.at(i).at(0), design_space.at(i).at(1)]
	
    std::vector<bool> fixed_parameters_mask;// 29-item vector. If mask.at(i) = true then ith parameter is fixed

	std::vector<double> design;// current design. Automatically loaded with fixed parameters from initialization

    std::vector<std::vector<double>> aft_point_cloud;// point cloud used to evaluate sectional area derivative at aft

    std::vector<std::vector<double>> bow_point_cloud;// point cloud used to evaluate sectional area derivative at bow


    /* everything below here is specific to the current design */

    double SAaft;// sectional area at aft

    double SAbow;// sectional area at bow

    double dSAaft;// derivative of sectional area curve at aft

    double dSAbow;// derivative of sectional area curve at bow

    /* everything below here is related to the geometry triangulation */

    //std::vector<std::vector<std::vector<double>>> triangulation;// the triangulation (see smpl_triangulation)
    smpl_triangulation::Triangulation<std::vector<double>> triangulation;// the triangulation (see smpl_triangulation)

    smpl_triangulation::Triangulation<std::vector<double>> triangulation_preimage;// triangulation pre-image at surface parametric plane (this does not change between designs)

    std::vector<std::vector<double>> triangulation_preimage_nodes_in_surface_domain;// not all nodes in 'triangulation_preimage' are in
                                                                                    // the surface domain [0,1]^2. Instead, they are mapped
                                                                                    // to [0,1]^2 and saved at 
                                                                                    // 'triangulation_preimage_nodes_in_surface_domain'
                                                                                    // so that for every design, when the actual on-surface
                                                                                    // triangulation is generated, the correct on-surface
                                                                                    // vertices are mirrored to the other side of the vessel

}; //KCSsimModeler

#endif //KCSSIMMODELER_HPP