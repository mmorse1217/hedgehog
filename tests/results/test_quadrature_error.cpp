#include "../catch.hpp"
#include "bie3d/solver_gmres_double_layer.hpp"
#include "bie3d/evaluator_far.hpp"
#include "bie3d/markgrid.hpp"
#include "common/stats.hpp"
#include "common/utils.hpp"
#include "bie3d/error_estimate.hpp"
#include "common/vtk_writer.hpp"
#include "../utils/evaluation_utils.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "bdry3d/p4est_refinement.hpp"
#include <math.h>
using namespace Ebi;
void laplace_dl_potential(Vec targets, int dof,Vec& potential){
    Vec source;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, 1*DIM, source);
    PetscInt idx[3] = {0,1,2};
    double v[3] = {0., 0.,.1};
    VecSetValues(source, 3, idx, v, INSERT_VALUES);
    VecAssemblyBegin(source);
    VecAssemblyEnd(source);
    Vec src_normal;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, 1*DIM, src_normal);
    v[2] =1.;
    VecSetValues(src_normal, 3, idx, v, INSERT_VALUES);
    VecAssemblyBegin(src_normal);
    VecAssemblyEnd(src_normal);
    
    Kernel3d kernel(121, vector<double>(2,1.));

    unique_ptr<PvFMM> fmm(new PvFMM());
    fmm->initialize_fmm(source, src_normal, targets, kernel);

    VecSet(potential, 0.);
    Vec density;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, 1, density);
    VecSet(density, 1.);

    // Evaluate
    fmm->evaluate_direct(density, potential);
}

//static double dist_to_patch = .01;
//double compute_integral(PatchSurfFaceMap* face_map, double dist_to_patch){
double compute_integral(PatchSurfFaceMap* face_map, Point3 target_point){
    //setup sampling and quadrature evaluator with current spacing
    unique_ptr<PatchSamples> samples (new PatchSamples("", ""));
    int num_patches = face_map->num_patches();
    vector<int> patch_partition(num_patches, 0);
    samples->bdry() = face_map;
    samples->patch_partition() = patch_partition;
    samples->setup();
    
    Vec target;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, 1*DIM, target);
    PetscInt idx[3] = {0,1,2};
    //double v[3] = {0., 0.,dist_to_patch};
    VecSetValues(target, 3, idx, target_point.array(), INSERT_VALUES);
    VecAssemblyBegin(target);
    VecAssemblyEnd(target);
    
    Kernel3d kernel(121, vector<double>(2,1.));
    unique_ptr<EvaluatorFar> evaluator(new EvaluatorFar("", ""));
    evaluator->knl()                  = kernel;
    evaluator->target_3d_position()   = target;
    evaluator->set_surface_discretization(samples.get());
    evaluator->setFromOptions();
    evaluator->setup();
    
    Vec density;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, samples->local_num_sample_points(), density);
    VecSet(density, 1.);
    Vec potential;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, 1, potential);
    VecSet(potential, 0.);
    // compute the integral
    evaluator->eval(density,  potential);
    double integral = 0;
    {
        DblNumMat p_local(1, potential);
        integral = p_local(0,0);
    }
    VecDestroy(&target);
    VecDestroy(&density);
    VecDestroy(&potential);
    return integral;
}

double error_estimate(int q, double delta, double phi, double max_jacobian){
    //double phi = 1.;
    //double max_jacobian = 2.;
    //double r= dist_to_patch; // TODO fix
    //double delta = max_jacobian*r; // TODO fix
    //double delta = dist_to_patch; // TODO fix
    double r = 1./(max_jacobian)*delta;
    double temp = 1./log(q);
    if(max_jacobian > 1.)
        temp *= pow(max_jacobian,3);
    else
        temp *= 1./pow(max_jacobian,3);
    //temp *= 3;
    //return .5*phi*max_jacobian*sqrt(M_PI)/tgamma(4)*sqrt(2*q*2*q)*delta*exp(-2*q*r);
    return phi*max_jacobian*sqrt(M_PI)/tgamma(4)*sqrt(2*q*2*q)*delta*exp(-4*q*r)*temp; // works for flat patch and both sides of curved patch
    //return phi*max_jacobian*sqrt(M_PI)/tgamma(4)*sqrt(2*q*2*q)*delta*exp(-4*q*r)*temp*1./(pow(max_jacobian,3)); 
}

void perturb_flat_patch(unique_ptr<PatchSurfFaceMap>& face_map, double tx, double ty){
    assert(face_map->num_patches() == 1);
    auto coefficients= face_map->face_map().control_points_write(0);
    // 4x4 grid of coefficients
    Point3 dx(0.,0.,tx);
    Point3 dy(0.,0.,ty);
    for (int i = 0; i < 4; i+=3) {
        for (int j = 1; j < 3; j++) {
            (*coefficients)(i*4 + j) += dx; 
        }
    }
    for (int i = 1; i < 3; i++) {
        for (int j = 0; j < 4; j+=3) {
            (*coefficients)(i*4 + j) += dy; 
        }
    }
    for (int i = 0; i < 4; i+=3) {
        for (int j = 0; j < 4; j+=3) {
            (*coefficients)(i*4 + j) += (dy+dx)*.5; 
            
        }
    }
        
    
}


void run_quad_estimate_test(string domain, Point3 target_point, int num_steps){
   
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "3");
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/"+domain+".wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/"+domain+".wrl");
    Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_files/poly/"+domain+".poly");
    Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "12");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".15");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".15");
    Options::set_value_petsc_opts("-bis3d_pts_max", "1000000000");

    // TODO need higher multipole order>16to evaluate beyond 9 digits below
    //Options::set_value_petsc_opts("-bis3d_np", "16");
    int qbkix_order = Options::get_int_from_petsc_opts("-near_interpolation_num_samples");

    // setup and adaptively refine true face-map to serve as a proxy for the
    // true integral value
    unique_ptr<PatchSurfFaceMap> refined_face_map(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    refined_face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
    refined_face_map->_coarse = true;

    refined_face_map->setFromOptions();
    int dx =+2.;
    int dy =-.5;
    refined_face_map->setup();
    perturb_flat_patch(refined_face_map, dx, dy);

    //refined_face_map->refine_test();
    refined_face_map->refine_uniform(3);
    auto f = refined_face_map.get();
    double true_integral = compute_integral(refined_face_map.get(),target_point);

    // setup single patch surface
    unique_ptr<PatchSurfFaceMap> face_map(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
    face_map->_coarse = true;
    face_map->setFromOptions();
    face_map->setup();
    perturb_flat_patch(face_map, dx, dy);
    face_map->refine_test();
    
    Markgrid::NearFieldMap closest_points = 
        Markgrid::compute_closest_points_on_patches(target_point, 
                0, 
                face_map.get(),
                vector<uint>(1,0));
    auto closest_point = closest_points[0][0];
    cout << "distance from patch to point:" << closest_point.distance_from_target << endl;


    auto patch = face_map->subpatch(0);

    vector<double> integrals(num_steps);
    
    for (int qi = 0; qi < num_steps; qi++) {
        double step = 1./(1.+double(qi));
        Options::set_value_petsc_opts("-bis3d_spacing", to_string(step));
        double approx_integral = compute_integral(face_map.get(), target_point);
        integrals[qi] = approx_integral;
    }
    DblNumMat data(3, num_steps);


    for (int i = 0; i < num_steps; i++) {
        double integral = integrals[i];
        double error = fabs(integral - true_integral)/fabs(true_integral);
        double phi = 1.;
        Options:: set_value_petsc_opts("-bis3d_spacing", to_string(1./double(i+1)));
        int q = floor(1./Options::get_double_from_petsc_opts("-bis3d_spacing")) +1;
        
        unique_ptr<PatchSamples> samples (new PatchSamples("", ""));
        vector<int> patch_partition(face_map->num_patches(), 0);
        samples->bdry() = face_map.get();
        samples->patch_partition() = patch_partition;
        samples->setup();
        double L = patch->characteristic_length();

        DblNumMat density_values(1, samples->local_num_sample_points());
        setvalue(density_values, 1.);

        DblNumMat uv_values = samples->sample_point_parametric_preimage(0);

        double estimate = ErrorEstimate::evaluate_error_estimate(patch, target_point, 
                closest_point, q, uv_values, density_values);

        cout << "computed error:" << error << ", estimated error: " << estimate << ", ratio: " << error/estimate << endl;
        CHECK(estimate/error <= 5.);

        data(0,i) = 1./(1.+double(i));
        data(1,i) = error; 
        data(2,i) = estimate;
    }
        Debug::save_mat(data, "quadrature_error_estimate_"+domain+".txt");

}


void error_sweep_over_curvature(string domain, int num_dist_to_patch, int num_quad_steps, int num_curvature_steps, double min_t, double max_t){

    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/"+domain+".wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/"+domain+".wrl");
    Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_files/poly/"+domain+".poly");
    Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "3");

    Options::set_value_petsc_opts("-near_interpolation_num_samples", "12");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".15");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".15");
    Options::set_value_petsc_opts("-bis3d_pts_max", "1000000000");

    // TODO need higher multipole order>16to evaluate beyond 9 digits below
    //Options::set_value_petsc_opts("-bis3d_np", "16");
    int qbkix_order = Options::get_int_from_petsc_opts("-near_interpolation_num_samples");

    double min_d = 1e-2;
    double max_d = 1.;
    // setup and adaptively refine true face-map to serve as a proxy for the
    int iter =0;
    for (int ti = 0; ti < num_curvature_steps; ti++) {
        for (int tj = 0; tj < num_curvature_steps; tj++) {
            for (int ci = 0; ci < num_dist_to_patch; ci++) {

                unique_ptr<PatchSurfFaceMap> refined_face_map(new PatchSurfFaceMap("BD3D_", "bd3d_"));
                refined_face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
                refined_face_map->_coarse = true;

                refined_face_map->setFromOptions();
                refined_face_map->setup();
                double tx = (max_t-min_t)/double(num_curvature_steps-1)*ti + min_t;
                double ty = (max_t-min_t)/double(num_curvature_steps-1)*tj + min_t;
                cout << "tx: " << tx << ", ty: " << ty << endl;
                double dist_to_patch  = (max_d-min_d)/double(num_dist_to_patch-1)*ci + min_d;

                perturb_flat_patch(refined_face_map, tx, ty);
                //refined_face_map->refine_test();
                refined_face_map->refine_uniform(3);
                auto f = refined_face_map.get();
                //resolve_function(refined_face_map->_p4est, f, &laplace_dl_potential, 1, 1e-11);


                // setup single patch surface
                unique_ptr<PatchSurfFaceMap> face_map(new PatchSurfFaceMap("BD3D_", "bd3d_"));
                face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
                face_map->_coarse = true;

                face_map->setFromOptions();
                face_map->setup();
                perturb_flat_patch(face_map, tx, ty);
                face_map->refine_test();

                auto patch = face_map->subpatch(0);
                vector<Point3> sample_point(3, Point3());
                Point2 uv(.5, .5);
                patch->xy_to_patch_coords(uv.array(), PatchSamples::EVAL_VL|PatchSamples::EVAL_FD, (double*)sample_point.data());
                Point3 normal = cross(sample_point[1], sample_point[2]);
                normal /= normal.length();
                Point3 target_point = sample_point[0] - dist_to_patch*normal;
                /*
                   Markgrid::NearFieldMap closest_points = 
                   Markgrid::compute_closest_points_on_patches(target_point, 
                   0, 
                   face_map.get(),
                   vector<uint>(1,0));
                   auto closest_point = closest_points[0][0];
                   dist_to_patch = closest_point.distance_from_target;
                   cout << "dist_to_patch: " << dist_to_patch << endl;
                   */

                double true_integral = compute_integral(refined_face_map.get(),target_point);

                //int num_steps = 25;
                double max_jacobian;
                patch->estimate_jacobian(&max_jacobian);
                double L = patch->characteristic_length();
                double mean_curvature = patch->mean_curvature(uv);
                double gaussian_curvature = patch->gaussian_curvature(uv);
                cout << "length:" << L << endl;
                cout <<  "max_jacobian: " << max_jacobian << endl;
                cout << "mean_curvature: " << mean_curvature << endl;
                cout  << "gaussian_curvature: " << gaussian_curvature << endl;

                vector<double> integrals(num_quad_steps);
                for (int qi = 0; qi < num_quad_steps; qi++) {
                    double step = 1./(1.+double(qi));
                    Options::set_value_petsc_opts("-bis3d_spacing", to_string(step));
                    double approx_integral = compute_integral(face_map.get(), target_point);
                    integrals[qi] = approx_integral;
                }

                DblNumMat data(7+6, num_quad_steps);
                for (int i = 0; i < num_quad_steps; i++) {
                    double integral = integrals[i];
                    double error = fabs(integral - true_integral)/fabs(true_integral);

                    data(0,i) = i+1; // quad order
                    data(1,i) = error; 
                    data(2,i) = mean_curvature;
                    data(3,i) = gaussian_curvature;
                    data(4,i) = max_jacobian;
                    data(5,i) = dist_to_patch;
                    data(6,i) = L;
                    for (int di = 0; di < DIM; di++) {
                        data(7+di, i) = sample_point[0](di);
                        data(7+DIM+di, i) = target_point(di);
                    }
                }
                Debug::save_mat(data, "output/quad_estimate_data/"+domain+"quad_error_curve_param_"+to_string(iter++)+".txt");
            }
        }
    }
}

void test_bie_quad(string domain, Point3 target_point, int num_steps){

    /*
       Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/flat_patch.wrl");
       Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/flat_patch.wrl");
       Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_files/poly/flat_patch.poly");
       Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/parabaloid.wrl");
       Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/parabaloid.wrl");
       Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_files/poly/parabaloid.poly");
       */
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "10");
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/"+domain+".wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/"+domain+".wrl");
    Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "8");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".15");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".15");
    Options::set_value_petsc_opts("-bis3d_pts_max", "1000000000");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");




    // setup boundary geometry
    unique_ptr<PatchSurfFaceMap> face_map(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    face_map->_surface_type = PatchSurfFaceMap::BLENDED;
    face_map->_coarse = true;
    face_map->setFromOptions();
    face_map->setup();
    face_map->refine_test();
    write_face_map_patches_to_vtk(DblNumMat(0,0), 
            vector<int>(face_map->num_patches(),0),
            face_map.get(), 0, "output/");
    write_face_map_patch_bounding_boxes_to_vtk(
            vector<int>(face_map->num_patches(),0),
            face_map.get(), 0, "output/",false);
    
    auto f = face_map.get();
    // adaptively upsample geometry
    /*refine_patches_for_fixed_qbkix_points(face_map->_p4est, f);
    face_map->patches() = p4est_to_face_map_subpatches(face_map->_p4est, f);*/

    NumVec<OnSurfacePoint> closest_points = 
        Markgrid::mark_target_points(DblNumMat(3,1,false, target_point.array()), face_map.get(), false);

    /*Markgrid::NearFieldMap closest_points = 
        Markgrid::compute_closest_points_on_patches(target_point, 
                0, 
                face_map.get(),
                vector<uint>(1,0));*/
    auto closest_point = closest_points(0);
    double dist_to_patch = closest_point.distance_from_target;
    auto patch = face_map->subpatch(closest_point.parent_patch);
    

    Options::set_value_petsc_opts("-bis3d_spacing", ".0333");
    double true_integral = compute_integral(f,target_point);


    Point2 uv = closest_point.parametric_coordinates;
    double L = patch->characteristic_length();
    //dist_to_patch /= L;
    double mean_curvature = patch->mean_curvature(uv);
    double gaussian_curvature = patch->gaussian_curvature(uv);
    Point3 ret[3];
    patch->xy_to_patch_coords(uv.array(), PatchSamples::EVAL_VL|PatchSamples::EVAL_FD, (double*)ret);
    double J = cross(ret[1], ret[2]).length();

    double hmax = (2*mean_curvature + sqrt(4*mean_curvature*mean_curvature - 4*gaussian_curvature))/(2*gaussian_curvature);
    hmax = fabs(hmax);
    
    DblNumMat data(3, num_steps);
    double target_accuracy = 1e-6;
    vector<double> integrals(num_steps);
    for (int qi = 0; qi < num_steps; qi++) {
        double step = 1./(1.+double(qi));
        Options::set_value_petsc_opts("-bis3d_spacing", to_string(step));

        double approx_integral = compute_integral(face_map.get(), target_point);
        integrals[qi] = approx_integral;
        double integral = integrals[qi];
        double error = fabs(integral - true_integral)/fabs(true_integral);
        //cout << "h: " << 1./(1.+double(i)) << ", value: " << integral << ", error:" << error << ", estimated error: " << error_estimate(i+1) << endl;
        //double estimate = error_estimate(i+1);
        double phi = 1.;

        //double estimate = error_estimate2(i+1, dist_to_patch, phi, max_jacobian );
        //double estimate = error_estimate3(i+1, dist_to_patch, max_jacobian, pullback);
        
        cout << "params:" <<  dist_to_patch/L << ", " << L << ", " << -L*mean_curvature << ", " << J << endl;
        double estimate = error_estimate_fit(double(qi+1), dist_to_patch/L, phi, -L*mean_curvature);
        //double estimate = error_estimate2(qi+1, dist_to_patch, phi, J);
        //double actual_near_zone_size = near_zone_approx_size(i+1, phi, max_jacobian, target_accuracy);
        //double approx_near_zone_size = near_zone_approx_size(i+1, phi, max_jacobian, estimate);
        cout << "computed error:" << error << ", estimated error: " << estimate << ", ratio: " << error/estimate << endl;
        CHECK(error < estimate);
        /*cout << "actual near-zone size: " <<actual_near_zone_size << 
            ", approx near-zone size: " <<approx_near_zone_size <<  
            ", ratio: " <<approx_near_zone_size/actual_near_zone_size <<  
            ", is quad accurate:" << is_quad_accurate_at_target(i+1, dist_to_patch, phi, 
                    max_jacobian, target_accuracy) << 
            endl;*/
        data(0,qi) = 1./(1.+double(qi));
        data(1,qi) = error; 
        data(2,qi) = estimate;
    }
    if(dist_to_patch >0)
        Debug::save_mat(data, "quadrature_error_estimate_"+domain+".txt");
    else
        Debug::save_mat(data, "quadrature_error_estimate_"+domain+"_m.txt");



}

TEST_CASE("Test quadrature error estimate", "[flat-patch][quad-error][results]"){
    Point3 target(0., 0., .25);
    run_quad_estimate_test("flat_patch",target, 15);
    //run_quad_estimate_test("flat_patch", target, 25);
    //target.z() = .0275;
    //run_quad_estimate_test("parabaloid", target, 15);
    //run_quad_estimate_test("parabaloid_flip", target, 25);
}
TEST_CASE("Test quadrature error w.r.t curvature", "[results][quad-error][curvature]"){
    error_sweep_over_curvature("flat_patch", 10, 20 , 10, -1., 1.);
}

TEST_CASE("Test error estimate vs. computed error", "[quad-error][bie]"){
    Point3 target(0., 0., .75);
    test_bie_quad("cube", target, 12);
}



TEST_CASE("Test bounding box inflation", "[bbox][flat-patch]"){
    /*
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/flat_patch.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/flat_patch.wrl");
    Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_files/poly/flat_patch.poly");
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/parabaloid.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/parabaloid.wrl");
    Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_files/poly/parabaloid.poly");
    */
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/new_ppp.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/new_ppp.wrl");
    Options::set_value_petsc_opts("-bis3d_spacing", ".05");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "12");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".15");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".15");
    Options::set_value_petsc_opts("-bis3d_pts_max", "1000000000");
    Options::set_value_petsc_opts("-target_accuracy", "1e-12");

    
    // setup single patch surface
    unique_ptr<PatchSurfFaceMap> face_map(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    face_map->_surface_type = PatchSurfFaceMap::BLENDED;
    face_map->_coarse = true;
    face_map->setFromOptions();
    face_map->setup();
    face_map->refine_test();
    for (int i = 0; i < 3; i++) {
        vector<int> pids;
        face_map->refine_uniform(1);
        for (int pi = 0; pi < face_map->num_patches(); pi++) {
            pids.push_back(pi);
        }
        write_face_map_patch_bounding_boxes_to_vtk(
                pids,
                face_map.get(), i, "output/bbox_test");
        write_face_map_patches_to_vtk(DblNumMat(0,0),pids,
                face_map.get(), i, "output/");

    }
    //assert(face_map->num_patches() == 1);
    /*write_face_map_patch_bounding_boxes_to_vtk(
            pids,
            face_map.get(), 100, "output/bbox_test");
    face_map->refine_uniform(3);
    pids.clear();
    for (int pi = 0; pi < face_map->num_patches(); pi++) {
        pids.push_back(pi);
    }
write_face_map_patches_to_vtk(DblNumMat(0,0),pids,
         face_map.get(), 100, "output/");*/


}
