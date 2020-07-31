#ifndef _UTILS_HPP_
#define _UTILS_HPP_
#include "ebi.hpp"
#include "kernel3d.hpp"
#include "vec3t.hpp"
using hedgehog::DblNumMat;
//class PatchSamples;
namespace ArgParse {

    int parse_kernel(char*  arg);
    int parse_layer(char*  arg);
    void read_command_line_args(int argc, char** argv, string options_file);
};


namespace Test {
    void compute_relative_error(Vec true_potential, Vec computed_potential, 
            MPI_Comm comm, int num_local_targets, int target_dof,
            string eval_type);
    Vec create_scaled_surface_target_points(vector<double> scale_factors, 
            //PatchSamples* patch_samples);
        Vec sample_points);

    Vec generate_random_vector(int size, double upper=0, double lower=1);
    Vec compute_dirichlet_boundary_data(
            MPI_Comm comm,
            hedgehog::Kernel3d kernel,
            Vec singularity_positions,
            Vec singularity_densities,
            Vec target_positions,
            Vec target_normals);
    Vec compute_neumann_boundary_data(
            MPI_Comm comm,
            hedgehog::Kernel3d kernel,
            Vec singularity_positions,
            Vec singularity_densities,
            Vec target_positions,
            Vec target_normals);

    Point3 cartesian_to_spherical(Point3 p);
double associated_legendre_function(int m, int n, double x);
double spherical_harmonic(int m, int n, Point3 x);
            Vec compute_spherical_harmonic_bc(MPI_Comm comm, 
                    Vec target_positions);
            Vec compute_spherical_harmonic_density(MPI_Comm comm, 
                    Vec target_positions);

    string get_domain();
    string get_kernel();
};
// TODO make sub-namespace of Petsc
namespace Options {
    double get_double_from_petsc_opts(string opt_name);
    bool get_bool_from_petsc_opts(string opt_name);
    int get_int_from_petsc_opts(string opt_name);
    string get_string_from_petsc_opts(string opt_name);
    void set_value_petsc_opts(string opt_name, string value);
};


namespace Petsc{
    int get_vec_size(Vec v);
    int get_vec_local_size(Vec v);
    void create_mpi_vec(MPI_Comm comm, int64_t size, Vec& v);
    void destroy_vec(Vec& v);
    // Concatenates two vectors A and B into vector C = [A B]
    Vec concatenate(Vec A, Vec B);
};

namespace Debug{

    void print_mat(DblNumMat m);
    void save_mat(DblNumMat m, string file);
    void load_mat(DblNumMat& m, string file);
    string load_file_to_string(string file);
};
#endif
