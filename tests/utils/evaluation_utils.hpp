#ifndef __EVALUATION_UTILS_HPP__
#define __EVALUATION_UTILS_HPP__
#include "bie3d/solver_gmres_double_layer.hpp"
#include "bdry3d/patch_surf_blended.hpp"
#include "common/kernel3d.hpp"
BEGIN_EBI_NAMESPACE

// Type of boundary data
enum class BoundaryDataType {
    CONSTANT_DENSITY = 0,
    HARMONIC = 1,
    RANDOM = 2,
    SPHERICAL_HARMONIC = 3,
    ANALYTIC = 4
};
// Type of singulariy distribution, if boundary data is harmonic
enum class SingularityType {
    // singularities on +/-c*e_i, i=1,2,3 where e_i is the standard basis
    AXIS_ALIGNED = 0,
    // singularities on the sphere of radius r, 
    // uniformly sampled on [0,2*\pi]x[0,\pi] 
    SPHERE = 1,
    // a single point singularity
    SINGLE = 2
};
enum class AnalyticSolutionType {
    ONE_PLUS_XYZ = 0
};
// which evaluation scheme we want to test
enum class EvaluationScheme{
    // Check value of the density
    NONE = -1,
    SMOOTH_QUAD = 0,
    NEAR_SINGULAR = 1,
    NEAR_QBKIX = 2,
    ON_SINGULAR = 3,
    ON_QBKIX = 4,
    ON_QBKIX_AVERAGE = 5,
    // Evaluate solution via Green's Identity 
    GREENS_IDENTITY =6,
    AUTOEVAL_QBKIX=7,
    CHECK_DENSITY = 8
};

// Target point distribution
enum class TargetType {
    GRID = 0,
    COLLOCATION_POINTS = 1,
    SPHERE_FAR = 2,
    SINGLE = 3,
    LINE = 4,
    RING = 5,
    ON_SURFACE_NON_COLLOCATION = 6,
    PLANE= 7
};
enum class SolverMatvec{
    ON_SINGULAR = 0,
    ON_QBKIX = 1
};

// how to solve (1/2*I+D)\phi = f
enum class SolutionScheme {
    // Explicit density is provided, simply evaluate at targets
    EXPLICIT_DENSITY = 0,
    // Solve for the density from given boundary data via GMRES
    GMRES_SOLVE = 1,
    NONE = 2
};

class TestConfig{
    public:
        string error_filename;

        // which type of matvec to use inside GMRES in the solver
        EvaluationType solver_matvec_type;
        // how we are computing the solution to the pde
        SolutionScheme solution_scheme;
        // which evaluation scheme(s) we want to test once the density is
        // computed
        EvaluationScheme evaluation_scheme;

        TargetType target_type;
        BoundaryDataType bc_type;
        AnalyticSolutionType analytic_solution_type;
        SingularityType singularity_type;
        // location of singularity for SingularityType::SINGLE
        Point3 single_singularity_location;
        // radius of sphere of singularities for SingularityType::SPHERE
        double sphere_radius_bc;
        // radius of sphere of target points for TargetType::SPHERE_FAR
        double sphere_radius_targets;
        //3d point to use as target
        Point3 single_target_point;
        
        Point3 target_plane_point;
        Point3 target_plane_vec1;
        Point3 target_plane_vec2;
        int num_targets=25;
        // how far along +/-e_i to place singularities in
        // SingularityType::AXIS_ALIGNED
        double c;

        // use FMM for far marking in point marking
        bool compute_far_marking;
        bool time_coarse_fmm;
        bool dump_values;
        TestConfig():TestConfig(""){;}
        TestConfig(string s):error_filename(s){;}
};

void check_error(DblNumMat true_values, DblNumMat computed_values, 
        double eps_abs, double eps_rel);

Vec compute_error_petsc_vec(Vec true_value, Vec computed_potential);

Vec generate_targets_in_cube(int num_samples_1d, Cube target_domain);

void tear_down(Vec targets, Vec true_potential, Vec computed_potential, 
        TestConfig test,   int target_dof,
        double eps_abs, double eps_rel, bool dump_values);

void setup_singularities(TestConfig test, unique_ptr<SolverGMRESDoubleLayer>& solver,
        Vec& positions, Vec& strengths);

void setup_target_data(unique_ptr<SolverGMRESDoubleLayer>& solver, 
        TestConfig test_type, Vec& targets, 
        NumVec<OnSurfacePoint>& closest_points ,
        Vec& computed_potential, Vec& true_potential);

void setup_problem_data(unique_ptr<SolverGMRESDoubleLayer>& solver, 
        TestConfig test_type, Vec& boundary_data, Vec& solved_density);

void solve_and_evaluate(unique_ptr<SolverGMRESDoubleLayer>& solver,
        Vec boundary_data, Vec targets, NumVec<OnSurfacePoint> closest_points,
        Vec& solved_density, Vec& computed_potential,
        TestConfig test_type);

void run_test(PatchSurf* surface,
        TestConfig test_type, bool dump_values=true);

void test_far_eval(unique_ptr<SolverGMRESDoubleLayer>& solver,
        TestConfig test_type, bool dump_values);

/*
 * TODO delete
 */
void test_far_eval_constant_density(unique_ptr<SolverGMRESDoubleLayer>& solver,
        bool dump_values);
void test_singular_quadrature_constant_density(unique_ptr<SolverGMRESDoubleLayer>& solver,
        bool dump_values);
void test_solver_constant_density(unique_ptr<SolverGMRESDoubleLayer>& solver,
        bool dump_values);
void test_qbkix_constant_density(unique_ptr<SolverGMRESDoubleLayer>& solver,
        bool dump_values);

void test_constant_boundary_data(PatchSurf* surface, bool dump_values=false);

void test_qbkix_extrapolation_no_quad(PatchSurf* surface, 
        double (*true_function)(Point3),
        bool dump_values=false);
/*
 * TODO delete
 */

Vec compute_analytic_boundary_data(MPI_Comm comm,
        AnalyticSolutionType analytic_solution_type,
        Vec target_positions);

unique_ptr<SolverGMRESDoubleLayer> setup_solver(PatchSurf* surface, 
        EvaluationType eval_type, bool compute_refined_surface=true);

Vec compute_3d_position_of_on_surface_points(
        NumVec<OnSurfacePoint> on_surface_points,
        PatchSurfFaceMap* face_map);

void check_potential(unique_ptr<PatchSurfBlended>& surface, 
        EvaluationType solver_matvec_type,  string eval_type, 
        string output_folder, int i=0);

double eval_x(Point3 p);
Point2 grad_exp_sin_cos(Point2 s);
Point3 grad_exp_sin_sin_sin(Point3 s);
double eval_exp_sin_sin(Point3 p);
unique_ptr<PatchSurfFaceMap> setup_face_map(PatchSurfFaceMap::SurfaceType surface_type, 
        int num_uniform_refinement_lvls, string ref_type=string("uniform"));
void setup_and_run_face_map_convergence_test(
        int patch_order,
        int patch_refinement_factor,
        int kernel_enum,
        PatchSurfFaceMap::SurfaceType surface_type,
        string domain,
        string output_folder,
        void (*convergence_test)(PatchSurfFaceMap*, string, int),
        void (*options_init)(),
        string polynomial_patch_filename = string(),
        int num_iterations=5,
        string ref_type=string("uniform"));


void test_gmres_solve_far_field_eval(
        PatchSurfFaceMap* surface,
        string output_folder,
        int i);
void test_gmres_solve_near_eval(
        PatchSurfFaceMap* surface,
        string output_folder,
        int i);

void gmres_test_base_options();
void solver_test_base_options();

void laplace_singluarity_cube(Vec samples, int dof,Vec& potential);
void laplace_singluarity_flat_patch(Vec samples, int dof,Vec& potential);
void laplace_singluarity_propeller(Vec samples, int dof,Vec& potential);

END_EBI_NAMESPACE
#endif
