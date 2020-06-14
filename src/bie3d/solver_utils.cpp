#include "solver_utils.hpp"
#include "fmm3d/pvfmm_bis_interface.hpp"
#include "markgrid.hpp"
#include "bdry3d/p4est_interface.hpp"

// MJM TODO remove these
BEGIN_EBI_NAMESPACE

//aux functions
int dim() { return 3; } 
//-----------------------------------------------------------------------------
// Petsc Options shortcuts
//-----------------------------------------------------------------------------

// TODO replace annoying calls to PetscOptionsGetXXX(...) with simple accessors
// by accesing options locally, it's making the code messy and obnoxious when
// there's no need

//-----------------------------------------------------------------------------
// Petsc Vector shortcuts
//-----------------------------------------------------------------------------

// Returns number of local points on a single processor
int  num_local_points(Vec pos) { 
    int64_t tmp;
    VecGetLocalSize(pos, &tmp);
    return tmp/dim(); 
}

// Returns global number of points across all processors
int  num_global_points(Vec pos) {
    int64_t tmp;
    VecGetSize(pos, &tmp);
    return tmp/dim();
}

// Store in beg and end the indices owned by the current processor
void local_index_range(Vec pos, int& beg, int& end) {
    int64_t petsc_beg, petsc_end;
    VecGetOwnershipRange(pos, &petsc_beg, &petsc_end);
    beg = petsc_beg/dim();
    end = petsc_end/dim();
}


int denscale(int sd, Vec constant, Vec den, Vec scaled_density)
{
  ebiFunctionBegin;
  double* constant_ptr;
  double* scaled_den_ptr;
  double* density_ptr;

  iC( VecGetArray(constant, &constant_ptr) );
  iC( VecGetArray(den, &density_ptr) );
  iC( VecGetArray(scaled_density, &scaled_den_ptr) );
  int64_t num_local_pts;
  iC( VecGetLocalSize(constant, &num_local_pts) );  
  for(int k=0; k<num_local_pts; k++)
	 for(int d=0; d<sd; d++){
		scaled_den_ptr[k*sd+d] = density_ptr[k*sd+d] * constant_ptr[k];
     }

  iC( VecRestoreArray(constant, &constant_ptr) );
  iC( VecRestoreArray(den, &density_ptr) );
  iC( VecRestoreArray(scaled_density, &scaled_den_ptr) );

  ebiFunctionReturn(0);
}

int denscale(int sd, vector<Vec>& constant_vecs, Vec den, Vec scaled_density)
{
  ebiFunctionBegin;
  vector<double*> constant_ptrs;
  constant_ptrs.resize(constant_vecs.size());

  for(int i=0; i<constant_vecs.size(); i++) 
	 iC( VecGetArray(constant_vecs[i], &(constant_ptrs[i])) );

  double* density_ptr;
  double* scaled_den_ptr;
  iC( VecGetArray(den, &density_ptr) );
  iC( VecGetArray(scaled_density, &scaled_den_ptr) );

  int64_t num_local_pts;
  iC( VecGetLocalSize(constant_vecs[0], &num_local_pts) );

  for(int k=0; k<num_local_pts; k++) {
	 double prod = 1;	 
     for(int i=0; i<constant_vecs.size(); i++)
         prod *= constant_ptrs[i][k];

	 for(int d=0; d<sd; d++)
		scaled_den_ptr[k*sd+d] = density_ptr[k*sd+d] * prod;
  }
  for(int i=0; i<constant_vecs.size(); i++) {
	 iC( VecRestoreArray(constant_vecs[i], &(constant_ptrs[i])) );
  }

  iC( VecRestoreArray(den, &density_ptr) );
  iC( VecRestoreArray(scaled_density, &scaled_den_ptr) );
  ebiFunctionReturn(0);
}



/*
void pvfmm_evaluation(Vec sources, Vec normals, Vec targets, Kernel3d kernel,
        Vec density, Vec potential){
    PvFMM* fmm = new PvFMM(sources, normals, targets, kernel);
    fmm->evaluate(density, potential);
    delete fmm;
}
*/

void find_closest_on_surface_points(Vec near_targets, PatchSamples* patch_samples, 
        Vec& closest_points, Vec& closest_face_points){

    int num_local_targets = num_local_points(near_targets);
    int num_local_samples = patch_samples->local_num_sample_points();

    double* near_targets_ptr;
    double* closest_points_ptr;
    double* closest_face_points_ptr;
    VecGetArray(near_targets, &near_targets_ptr);
    VecGetArray(closest_points, &closest_points_ptr);
    VecGetArray(closest_face_points, &closest_face_points_ptr);
    
    DblNumMat near_targets_local(DIM, num_local_targets, false, near_targets_ptr);
    DblNumMat closest_points_local(DIM, num_local_targets, false, closest_points_ptr);
    DblNumMat closest_face_points_local(patch_samples->face_point_size_in_doubles(),
                    num_local_targets, false, closest_face_points_ptr);

    for(int ti = 0; ti < num_local_targets; ti++){
        Point3 current_target(near_targets_local.clmdata(ti));
        Point3 closest_sample;
        double closest_sample_distance = DBL_MAX;
        int closest_sample_index = 0;
        int patch_containing_closest_sample = 0;

        // Find the closest collocation point to the ti-th target
        // MJM TODO replace this loop with a call to PvFMM tree if possible?
        // MJM TODO move to a function on PatchSamples
        for(int pi = 0; pi < patch_samples->bdry()->patches().size(); pi++){

            DblNumMat sample_point_3d_position         = patch_samples->sample_point_3d_position(pi);
            DblNumMat sample_point_parametric_preimage = patch_samples->sample_point_parametric_preimage(pi);

            for(int si = 0; si < num_local_samples; si++){
                Point3 current_sample(sample_point_3d_position.clmdata(si));
                double current_sample_distance = (current_sample - current_target).length();

                if(current_sample_distance < closest_sample_distance){
                    closest_sample          = current_sample;
                    closest_sample_distance = current_sample_distance;
                    closest_sample_index    = si;
                    patch_containing_closest_sample = pi;
                }
            }
        }

        // Start searching for the closest on-surface point using the nearest
        // collocation point as a starting point
        
        // MJM NOTE this code is aggressively adapted from the relavent section
        // of markgrid(); may or may not be optimal...
        
        Point3 closest_sample_point = patch_samples->
            sample_point_3d_position(patch_containing_closest_sample).
            clmdata(closest_sample_index);

        Point2 closest_sample_point_parametric_preimage = patch_samples->
            sample_point_parametric_preimage(patch_containing_closest_sample).
            clmdata(closest_sample_index);
        Patch* current_patch = patch_samples->bdry()->patches()[patch_containing_closest_sample];

        Point3 position_and_derivs[3];
        current_patch->xy_to_patch_coords(
                closest_sample_point_parametric_preimage,
                PatchSamples::EVAL_VL|PatchSamples::EVAL_FD,
                (double*) position_and_derivs);

        Point3 target_to_closest = current_target - position_and_derivs[0];
        Point3 deriv_x = position_and_derivs[1];
        Point3 deriv_y = position_and_derivs[2];
        Point3 normal = cross(deriv_x, deriv_y);
        normal /= normal.length();

        double p_dot_n = dot(target_to_closest, normal);
        bool is_valid = true;

        while(abs(p_dot_n) < 0.975*target_to_closest.length() &&
                target_to_closest.length() > 1e-12){
            Point3 pp = target_to_closest - p_dot_n*normal;
            Point3 n_cross_deriv_x = cross(normal, deriv_x);
            Point3 n_cross_deriv_y = cross(normal, deriv_y);
            Point2 parametric_update(
                    dot(n_cross_deriv_y, pp)/dot(n_cross_deriv_y, deriv_x),
                    dot(n_cross_deriv_x, pp)/dot(n_cross_deriv_x, deriv_y));

            // And let's not forget...
            // LEXING, 0.25 VERY IMPORTANT
            closest_sample_point_parametric_preimage += (1.0/64.0)*parametric_update;
            current_patch->is_xy_valid(closest_sample_point_parametric_preimage.array(), is_valid); 
            if(!is_valid){
                cerr << "sample not on patch: " << closest_sample_point_parametric_preimage <<" " << parametric_update << endl;
                assert(is_valid);
            } else {
                current_patch->xy_to_patch_coords(
                        closest_sample_point_parametric_preimage,
                        PatchSamples::EVAL_VL|PatchSamples::EVAL_FD,
                        (double*) position_and_derivs);

                Point3 target_to_closest = current_target - position_and_derivs[0];
                Point3 deriv_x = position_and_derivs[1];
                Point3 deriv_y = position_and_derivs[2];
                Point3 normal = cross(deriv_x, deriv_y);
                normal /= normal.length();
            }
        }
        Point3 closest_point = position_and_derivs[0];
        for(int d =0; d < DIM; d++){
            closest_points_local(d,ti) = closest_point(d);
        }

        FacePointOverlapping* face_point_ptr = 
            (FacePointOverlapping*) (closest_face_points_local.clmdata(ti));
        current_patch->xy_to_face_point(closest_sample_point_parametric_preimage, face_point_ptr);
    }
    VecRestoreArray(near_targets, &near_targets_ptr);
    VecRestoreArray(closest_points, &closest_points_ptr);
    VecRestoreArray(closest_face_points, &closest_face_points_ptr);
}

// Write sample point values to a tab-separated text file for parsing in Python
void write_to_text_file(
        string filename,
        Vec sample_points, // size DIM x num_local_samples
        int num_local_samples,
        Vec density, // size source_dof() x num_local_samples
        int source_dof,
        Vec potential, // size target_dof() x num_local_samples
        int target_dof,
        Vec val1, //size val1_stride x num_local_samples
        int val1_stride,
        Vec val2, //size val2_stride x num_local_samples
        int val2_stride,
        Vec val3, //size val3_stride x num_local_samples
        int val3_stride){
    
    DblNumMat sample_points_local   = get_local_vector(DIM, num_local_samples, sample_points);
    DblNumMat density_local         = get_local_vector(source_dof, num_local_samples, density);
    DblNumMat potential_local       = get_local_vector(target_dof, num_local_samples, potential);
    DblNumMat val1_local            = get_local_vector(val1_stride, num_local_samples, val1);
    DblNumMat val2_local            = get_local_vector(val2_stride, num_local_samples, val2);
    DblNumMat val3_local            = get_local_vector(val3_stride, num_local_samples, val3);
    
    //Dump these arrays to a file
    //Note this doesn't work in parallel
    ofstream out;
    out.open(filename.c_str());
    out.precision(16);

    // Write the column headers for each vector
    out << "pos_1" << "\t" 
        << "pos_2" << "\t" 
        << "pos_3" << "\t";
    for(int i = 0; i < source_dof; i++){
        out << "density_" << i << "\t";
    }
    
    for(int i = 0; i < target_dof; i++){
        out << "potential_" << i << "\t";
    }

    for(int i = 0; i < val1_stride; i++){
        out << "val1_" << i << "\t";
    }
    
    for(int i = 0; i < val2_stride; i++){
        out << "val2_" << i << "\t";
    }

    for(int i = 0; i < val3_stride; i++){
        out << "val3_" << i << "\t";
    }
    out << "\n";

    // Write the data for each sample point
    for(int si = 0; si < num_local_samples; si++){
        for(int i = 0; i < DIM; i++)
            out << sample_points_local(i,si) << "\t";

        for(int i = 0; i < source_dof; i++)
            out << density_local(i,si) << "\t";
        
        for(int i = 0; i < target_dof; i++)
            out << potential_local(i,si) << "\t";

        for(int i = 0; i < val1_stride; i++)
            out << val1_local(i,si) << "\t";
        
        for(int i = 0; i < val2_stride; i++)
            out << val2_local(i,si) << "\t";

        for(int i = 0; i < val3_stride; i++)
            out << val3_local(i,si) << "\t";

        out << "\n";
    }
    out.close();

    /*
    sample_points_local.restore_local_vector(); // MJM BUG causes seg fault, FIXME
    density_local.restore_local_vector();
    potential_local.restore_local_vector();
    if(val1_stride != 0){
        val1_local.restore_local_vector();
    }
    if(val2_stride != 0){
        val2_local.restore_local_vector();
    }
    if(val3_stride != 0){
        val3_local.restore_local_vector();
    }
    */

}


void write_to_file(
        string filename,
        DblNumMat m,
        DblNumVec v){

    ofstream out;
    out.open(filename.c_str());
    out.precision(16);
    assert(m.n() == v.m());
    if(m.n() != 0 && m.m() != 0){
        for(int d = 0; d < m.m(); d++){
            out << "m_" << d << "\t";
        }
    }
    if(v.m() != 0){
        out << "v" << endl;
    }
    for(int i = 0; i < m.n(); i++){
        for(int d = 0; d < m.m(); d++){
            out << m(d,i) << "\t";
        }
            out << v(i) << endl;
    }
    out.close();
}


void evaluate_singularity_solution(
        Kernel3d kernel,
        Vec singularity_positions,      // DIM x num_singularities
        Vec singularity_densities,      // source_dof x num_singularities
        Vec target_points,              // DIM x num_targets
        Vec target_potential){          // targe_dof x num_targets
    

    // MJM TODO add asserts to verify correct size vectors based on kernel DOF's
   /* 

    cout << "before fmm initialize" <<  endl;
    PvFMM* fmm = new PvFMM(
            singularity_positions,
            singularity_normals,
            target_points,
            kernel);
    Vec dummy_normals;
    VecCreateMPI(MPI_COMM_WORLD, num_local_samples*3, PETSC_DETERMINE, &dummy_normals);
    VecSet(dummy_normals, 0.);
    cout << "after fmm initialize" <<  endl;
    
    cout << "before fmm evaluate" <<  endl;
    fmm->evaluate(singularity_densities, target_potential);
    cout << "after fmm evaluate" <<  endl;
    delete fmm;
*/
}

void axis_aligned_singularities(
        MPI_Comm comm,
        Kernel3d problem_kernel,
        Vec& singularity_positions,
        Vec& singularity_normals,
        Vec& singularity_densities){

    // Place singularities at +/-R in the x-, y-, and z- dimensions, with
    // constant density 4pi and unit "normal" direction pointing toward the
    // origin
    cout << "creating singularity solution" << endl;
    int source_dof = problem_kernel.get_sdof();
    int target_dof = problem_kernel.get_tdof();
    //Vec singularity_positions, singularity_normals, singularity_densities;
    int num_singularities = 6;
    VecCreateMPI(comm, num_singularities*3, PETSC_DETERMINE, &singularity_positions);
    VecCreateMPI(comm, num_singularities*3, PETSC_DETERMINE, &singularity_normals);
    VecCreateMPI(comm, num_singularities*source_dof, PETSC_DETERMINE, &singularity_densities);

    double zero = 0.;
    VecSet(singularity_positions, zero);
    VecSet(singularity_normals, zero);
    VecSet(singularity_densities, zero);

    double* singularity_position_ptr;
    double* singularity_normal_ptr;
    double* singularity_density_ptr;
    VecGetArray(singularity_positions, &singularity_position_ptr);
    VecGetArray(singularity_normals, &singularity_normal_ptr);
    VecGetArray(singularity_densities, &singularity_density_ptr);
   
    DblNumMat positions(DIM, num_singularities, false, singularity_position_ptr);
    DblNumMat   normals(DIM, num_singularities, false, singularity_normal_ptr);
    DblNumMat densities(source_dof, num_singularities, false, singularity_density_ptr);

    double R = 5.;
    double fourpi = 4*M_PI;
    for(int i=0; i < num_singularities/2; i++){
        // positive unit vector in each dimension
        positions(i, i) = 1.;
        // negative unit vector in each dimension
        positions(i, i + num_singularities/2) = -1.;

        // normals pointing toward the center
        normals(i, i) = -1.;
        normals(i, i + num_singularities/2) = 1.;

        densities(0, i) = -1; // 
        //densities(2, i) = 1; // 
    }
    for(int i=0; i < num_singularities; i++)
        densities(0, i) = -1;
    
    VecRestoreArray(singularity_positions, &singularity_position_ptr);
    VecRestoreArray(singularity_normals, &singularity_normal_ptr);
    VecRestoreArray(singularity_densities, &singularity_density_ptr);


    // Scale points and density appropriately
    VecScale(singularity_positions, R); 
    VecScale(singularity_densities, fourpi);
    cout << "created singularity solution" << endl;

}

Vec evaluate_singularities_along_basis(
        MPI_Comm comm,
        Kernel3d problem_kernel, 
        Vec target_points){

    // Place singularities at +/-R in the x-, y-, and z- dimensions, with
    // constant density 4pi and unit "normal" direction pointing toward the
    // origin
    cout << "creating singularity solution" << endl;
    int source_dof = problem_kernel.get_sdof();
    int target_dof = problem_kernel.get_tdof();
    Vec singularity_positions, singularity_normals, singularity_densities;
    int num_singularities = 6;
    VecCreateMPI(comm, num_singularities*3, PETSC_DETERMINE, &singularity_positions);
    VecCreateMPI(comm, num_singularities*3, PETSC_DETERMINE, &singularity_normals);
    VecCreateMPI(comm, num_singularities*source_dof, PETSC_DETERMINE, &singularity_densities);

    double zero = 0.;
    VecSet(singularity_positions, zero);
    VecSet(singularity_normals, zero);
    VecSet(singularity_densities, zero);

    double* singularity_position_ptr;
    double* singularity_normal_ptr;
    double* singularity_density_ptr;
    VecGetArray(singularity_positions, &singularity_position_ptr);
    VecGetArray(singularity_normals, &singularity_normal_ptr);
    VecGetArray(singularity_densities, &singularity_density_ptr);
   
    DblNumMat positions(DIM, num_singularities, false, singularity_position_ptr);
    DblNumMat   normals(DIM, num_singularities, false, singularity_normal_ptr);
    DblNumMat densities(source_dof, num_singularities, false, singularity_density_ptr);

    double R = 5.;
    double fourpi = 4*M_PI;
    for(int i=0; i < num_singularities/2; i++){
        // positive unit vector in each dimension
        positions(i, i) = 1.;
        // negative unit vector in each dimension
        positions(i, i + num_singularities/2) = -1.;

        // normals pointing toward the center
        normals(i, i) = -1.;
        normals(i, i + num_singularities/2) = 1.;

        //densities(0, i) = -1; // 
        //densities(2, i) = 1; // 
    }
    for(int i=0; i < num_singularities; i++)
        densities(0, i) = -1;
    
    VecRestoreArray(singularity_positions, &singularity_position_ptr);
    VecRestoreArray(singularity_normals, &singularity_normal_ptr);
    VecRestoreArray(singularity_densities, &singularity_density_ptr);


    // Scale points and density appropriately
    VecScale(singularity_positions, R); 
    VecScale(singularity_densities, fourpi);
    cout << "created singularity solution" << endl;


    // Create target potential vector
    Vec target_potential;
    int num_local_targets = num_local_points(target_points);
    VecCreateMPI(comm, num_local_targets*target_dof, PETSC_DETERMINE, &target_potential);
    
    // Evaluate potential induced by singularities
    evaluate_singularity_solution(
            problem_kernel,
            singularity_positions,      
            singularity_normals,        
            singularity_densities,      
            target_points,              
            target_potential);
    cout << "inside solver_utils" << endl;
    return target_potential;
}
void evaluate_singularity_solution(
        Kernel3d kernel,
        Vec singularity_positions,      // DIM x num_singularities
        Vec singularity_normals,        // DIM x num_singularities
        Vec singularity_densities,      // source_dof x num_singularities
        Vec target_points,              // DIM x num_targets
        Vec target_potential){          // targe_dof x num_targets
    

    // MJM TODO add asserts to verify correct size vectors based on kernel DOF's
    
    cout << "before fmm initialize" <<  endl;
    PvFMM* fmm = new PvFMM(
            singularity_positions,
            singularity_normals,
            target_points,
            kernel);
    cout << "after fmm initialize" <<  endl;
    
    cout << "before fmm evaluate" <<  endl;
    fmm->evaluate(singularity_densities, target_potential);
    cout << "after fmm evaluate" <<  endl;
    delete fmm;

}


Vec evaluate_solution_x(MPI_Comm comm, Kernel3d problem_kernel, Vec target_points){
    int num_local_targets = num_local_points(target_points);
    int target_dof = problem_kernel.get_tdof();
    Vec target_potential;
    VecCreateMPI(comm,
            num_local_targets*target_dof, 
            PETSC_DETERMINE,
            &target_potential);
    
    DblNumMat target_potential_local(target_dof, num_local_targets, target_potential);
    DblNumMat target_points_local(DIM, num_local_targets, target_points);
    for(int i = 0; i < num_local_targets; i++){
        for(int d = 0; d < target_dof; d++){
            target_potential_local(d,i) = d == 1 ? target_points_local(0,i) : 0.;
        }
    }
    target_potential_local.restore_local_vector();
    target_points_local.restore_local_vector();
    return target_potential;
}

Vec evaluate_solution_zxy(MPI_Comm comm, Kernel3d problem_kernel, Vec target_points){
    int num_local_targets = num_local_points(target_points);
    int target_dof = problem_kernel.get_tdof();
    Vec target_potential;
    VecCreateMPI(comm,
            num_local_targets*target_dof, 
            PETSC_DETERMINE,
            &target_potential);
    
    DblNumMat target_potential_local(target_dof, num_local_targets, target_potential);
    DblNumMat target_points_local(DIM, num_local_targets, target_points);
    for(int i = 0; i < num_local_targets; i++){
        target_potential_local(0,i) =  target_points_local(2,i);
        target_potential_local(1,i) =  target_points_local(0,i);
        target_potential_local(2,i) =  target_points_local(1,i);
    }
    target_potential_local.restore_local_vector();
    target_points_local.restore_local_vector();
    return target_potential;
}

void copy_values_to_subvec( vector<int64_t> vec_indicies_to_copy, int64_t stride, 
        CopyDirection copy_dir,
        Vec& subvec, Vec& vec){
    int subvec_size = Petsc::get_vec_local_size(subvec)/stride;
    int vec_size = Petsc::get_vec_local_size(vec)/stride;
    assert(vec_size >= subvec_size);
    assert(subvec_size == vec_indicies_to_copy.size());
    DblNumMat subvec_local = get_local_vector(stride, subvec_size, subvec);
    DblNumMat vec_local = get_local_vector(stride, vec_size, vec);

    for(int i = 0; i < vec_indicies_to_copy.size(); i++){
        int index = vec_indicies_to_copy[i];
        for(int d = 0; d < stride; d++){
            if(copy_dir == COPY_TO){
                subvec_local(d,i) = vec_local(d,index);
            } else if (copy_dir == COPY_FROM){
                vec_local(d,index) = subvec_local(d,i);
            }
        }
    }

    subvec_local.restore_local_vector();
    vec_local.restore_local_vector();

}




END_EBI_NAMESPACE
