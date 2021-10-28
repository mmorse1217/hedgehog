#ifndef __SOLVER_UTILS__
#define __SOLVER_UTILS__

#include "common/ebi.hpp"
#include "common/kernel3d.hpp"
#include "bdry3d/patch_samples.hpp"
//#include <p4est.h>
#include "bdry3d/p4est_interface.hpp"
BEGIN_EBI_NAMESPACE
//#define DIM 3
#include <tuple>


//aux functions
int dim();
int mpiSize();

// Returns number of local points on a single processor
int  num_local_points(Vec pos);

// Returns global number of points across all processors
int  num_global_points(Vec pos);

void local_index_range(Vec pos, int& beg, int& end);

//--------------------------
int denscale(int dof, Vec cst, Vec den, Vec fnl); //scale density by alf or wcb
int denscale(int dof, vector<Vec>& csts, Vec den, Vec fnl);

//Vec pvfmm_evaluation(Vec sources, Vec normals, Vec targets, Kernel3d kernel, Vec density);

// Write sample point values to a tab-separated text file for parsing in Python
void write_to_text_file(
        string filename,
        Vec sample_points, // size DIM x num_local_samples
        int num_local_samples,
        Vec density, // size source_dof() x num_local_samples
        int source_dof,
        Vec potential, // size target_dof() x num_local_samples
        int target_dof,
        Vec val1=NULL, //size val1_stride x num_local_samples
        int val1_stride=0,
        Vec val2=NULL, //size val2_stride x num_local_samples
        int val2_stride=0,
        Vec val3=NULL, //size val3_stride x num_local_samples
        int val3_stride=0);


void write_to_file(
        string filename,
        DblNumMat m,
        DblNumVec v);

void axis_aligned_singularities(
        MPI_Comm comm,
        Kernel3d problem_kernel,
        Vec& singularity_positions,
        Vec& singularity_normals,
        Vec& singularity_densities); 

void evaluate_singularity_solution(
        Kernel3d kernel,
        Vec singularity_positions,      // DIM x num_singularities
        Vec singularity_densities,      // source_dof x num_singularities
        Vec target_points,              // DIM x num_targets
        Vec target_potential);          // targe_dof x num_targets

void evaluate_singularity_solution(
        Kernel3d kernel,
        Vec singularity_positions,      // DIM x num_singularities
        Vec singularity_normals,        // DIM x num_singularities
        Vec singularity_densities,      // source_dof x num_singularities
        Vec target_points,              // DIM x num_targets
        Vec target_potential);          // targe_dof x num_targets

Vec evaluate_singularities_along_basis(
        MPI_Comm comm,
        Kernel3d problem_kernel, 
        Vec target_points);

Vec evaluate_solution_x(MPI_Comm comm, Kernel3d problem_kernel, Vec target_points);
Vec evaluate_solution_zxy(MPI_Comm comm, Kernel3d problem_kernel, Vec target_points);

enum CopyDirection{
   COPY_TO= 1,
   COPY_FROM = 2
};

void copy_values_to_subvec( vector<int64_t> subvec_indicies, int64_t stride, 
        CopyDirection copy_dir,
        Vec& subvec, Vec& vec);
END_EBI_NAMESPACE

#endif
