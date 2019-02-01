#include "markgrid.hpp"
#include "bie3d/solver_gmres_double_layer.hpp"
#include <sampling.hpp>
#include <string>
#include <ctime>
#include <set>
#include "bdry3d/geogram_interface.hpp"

// 3/27/18 MJM - geogram defines the FAR macro; when including into markgrid, it
// explodes, so we need to undefine it. hopefully i'll remove the enums
// someday...
#ifdef FAR
#undef FAR
#endif
BEGIN_EBI_NAMESPACE
using std::pair;
using std::min;
using std::max;
using std::abs;
using std::cerr;
using std::ofstream;
using namespace Markgrid;
// ---------------------------------------------------------------------- 
// Begin Generic implementation of markgrid
// ---------------------------------------------------------------------- 

int mpiSize(){

  int mpisize;
  MPI_Comm comm;
  comm = PETSC_COMM_WORLD;
  MPI_Comm_size(comm, &mpisize);
  return mpisize;
}


NumVec<OnSurfacePoint> Markgrid::mark_target_points(
        DblNumMat input_points, 
        PatchSurfFaceMap* face_map,
        bool compute_far_marking){
    // make a list of empty OnSurfacePoint's for each input point
    NumVec<OnSurfacePoint> on_surface_points(input_points.n());

    // points are UNMARKED and NOWHERE by default
    setvalue(on_surface_points, OnSurfacePoint(-1, DBL_MAX, Point2(-1,-1), UNMARKED, -1)); 
    if(compute_far_marking)
        // Mark all points via FMM as far + in/out and unmarked for further processing
        mark_far_field(input_points, face_map, on_surface_points);


    // Mark all remaining unmarked points by near marking
    mark_near_field(input_points, face_map, on_surface_points);

    check_all_points_are_marked(input_points, on_surface_points);
    return  on_surface_points;
}

void Markgrid::mark_target_points(
        DblNumMat input_points, 
        PatchSurfFaceMap* face_map, 
        NearFieldMap& on_surface_point_map){
    // make a list of empty OnSurfacePoint's for each input point
    //NumVec<OnSurfacePoint> on_surface_points(input_points.n());

    // points are UNMARKED and NOWHERE by default
    //setvalue(on_surface_points, OnSurfacePoint(-1, DBL_MAX, Point2(-1,-1), UNMARKED, -1)); 
    //if(compute_far_marking)
        // Mark all points via FMM as far + in/out and unmarked for further processing
        //mark_far_field(input_points, face_map, on_surface_points);


    // Mark all points by near marking algorithm
    mark_near_field(input_points, face_map, on_surface_point_map);

    check_all_points_are_marked(input_points, on_surface_point_map);
    //return  on_surface_points;
}



void Markgrid::mark_far_field(DblNumMat input_points, PatchSurfFaceMap* face_map, 
        NumVec<OnSurfacePoint>& on_surface_points){
    
    assert(input_points.n() == on_surface_points.m());

    Vec potential = evaluate_fmm_with_constant_density(input_points, face_map);
    DblNumMat potential_local = get_local_vector(1, input_points.n(), potential);
    int num_near_points = 0;
    cout << "marking far field" << endl;
    for(int i = 0; i < input_points.n(); i++){
        double potential_at_target = potential_local(0,i);
        if(fabs(potential_at_target - 1.)  <= 1e-4){
            on_surface_points(i).inside_domain = INSIDE;
            on_surface_points(i).region = FAR;
            
        } else if(fabs(potential_at_target) <= 1e-4){
            on_surface_points(i).inside_domain = OUTSIDE;
            on_surface_points(i).region = FAR;
        } else {
            num_near_points++;
        }
        // Near points are still UNMARKED and NOWHERE
        //on_surface_points(i) = on_surface_point;

    }
    cout << "marked far field" << endl;
    cout << "points left to mark: " << num_near_points << endl;

    potential_local.restore_local_vector();
    VecDestroy(&potential);
}

void Markgrid::mark_near_field(DblNumMat input_points, PatchSurfFaceMap* face_map, 
        NumVec<OnSurfacePoint>& on_surface_points){

    int NN = input_points.n();
    //unique_ptr<SpatialGrid> grid(new SpatialGrid(face_map, input_points));
    unique_ptr<AABBTree> aabb_tree(new AABBTree(face_map));
    
    vector<int> qbkix_indices(NN);
    for (int i = 0; i < NN; i++)
        qbkix_indices[i] = i;

    std::srand ( unsigned ( std::time(0) ) );
    std::random_shuffle ( qbkix_indices.begin(), qbkix_indices.end() );

    double total_mark_time = 0;
cout << "marking qbkix points" << endl;
#pragma omp parallel for
    for (int j = 0; j < NN; j++){
        int i = qbkix_indices[j]; 
        if(long(100*j/NN) % 20 == 0){

            cout << long(100*j/NN) <<"%"<< endl;
        }
    
        Point3 target_point(input_points.clmdata(i));
        //NearFieldMap patches_near_target_point = 
            //find_patches_near_point(target_point, i, surface, grid);
        if(on_surface_points(i).region != FAR && on_surface_points(i).inside_domain != OUTSIDE_SPATIAL_GRID){
            double start_time = omp_get_wtime();
            NearFieldMap patches_near_target_point = 
            find_closest_patch_to_point_aabb_tree(target_point, i, face_map, aabb_tree.get());
#pragma omp critical 
            total_mark_time += omp_get_wtime() - start_time;
            //find_closest_patch_to_point_via_bfs(target_point, i, face_map, grid.get());
                //find_patches_closest_to_point(target_point, i, face_map, grid);

            int num_possible_closest_points = patches_near_target_point[i].size();
            if(num_possible_closest_points == 0){
                // this means that FMM missed marking some points. revisit
                assert(0);
            } else { 
                on_surface_points(i) = Markgrid::find_closest_on_surface_point_in_list(patches_near_target_point[i]);
                /*
                for(int si= 0; si < num_possible_closest_points; si++){
                    OnSurfacePoint possible_closest_point = 
                        patches_near_target_point[i][si];

                    on_surface_points(i) = 
                        possible_closest_point.distance_from_target < on_surface_points(i).distance_from_target ?
                        possible_closest_point : on_surface_points(i);
                }*/
            }
        }
        assert(on_surface_points(i).inside_domain != NOWHERE);
    }
    cout << "marked near field" << endl;
    cout << "total time spent marking: "<< total_mark_time << endl;
}

void Markgrid::mark_near_field(DblNumMat input_points, PatchSurfFaceMap* face_map, 
        NearFieldMap& on_surface_point_map){

    int NN = input_points.n();
    //SpatialGrid* grid = new SpatialGrid(face_map, input_points);
    //unique_ptr<SpatialGrid> grid(new SpatialGrid(face_map, input_points));
    unique_ptr<AABBTree> aabb_tree(new AABBTree(face_map));
    // Need near_patches_buffer for thread safety in loop below
    vector<NearFieldMap> near_field_map_buffer(NN);

    assert(on_surface_point_map.empty()); // <=== is this even right? only use 2c refinement!
    cout << "marking near field" << endl;

    // make a list of qbkix indices to paralellize over. shuffle the processing
    // order for better average case load balancing.
    vector<int> qbkix_indices(NN);
    for (int i = 0; i < NN; i++)
        qbkix_indices[i] = i;

    std::srand ( unsigned ( std::time(0) ) );
    std::random_shuffle ( qbkix_indices.begin(), qbkix_indices.end() );
cout << "marking qbkix points" << endl;
#pragma omp parallel for
    for (int j = 0; j < NN; j++){
        int i = qbkix_indices[j]; 
        //cout << j << endl;
        if(long(100*j/NN) % 20 == 0){

            cout << long(100*j/NN) <<"%"<< endl;
        }
        
        Point3 target_point(input_points.clmdata(i));
        
        // Find the closest on-surface points to the target point on face_map 
        // This should return 1/2/4 points depending if the closest point to the
        // target is inside the patch, along an edge or a vertex.
        NearFieldMap near_field_map= 
            find_closest_patch_to_point_aabb_tree(target_point, i, face_map, aabb_tree.get());
            //find_patches_closest_to_point(target_point, i, face_map, grid);
            //find_closest_patch_to_point_via_bfs(target_point, i, face_map, grid.get());
            //find_patches_closest_to_point( target_point, i, face_map, NULL);
        
        near_field_map_buffer[i] = near_field_map;
    }
    // aggregate into on_surface_point_map
    for (int i = 0; i < NN; i++){
        NearFieldMap ith_near_field_map = near_field_map_buffer[i];
        on_surface_point_map.insert(ith_near_field_map.begin(), ith_near_field_map.end());
    }
    
    cout << "marked near field" << endl;
}

void Markgrid::check_all_points_are_marked(DblNumMat input_points, 
        NumVec<OnSurfacePoint>& on_surface_points){
#pragma omp parallel for
    for (int i = 0; i < input_points.n(); i++){
        if(on_surface_points(i).region != FAR && on_surface_points(i).inside_domain != OUTSIDE_SPATIAL_GRID){
            assert(on_surface_points(i).distance_from_target < DBL_MAX);
            assert(on_surface_points(i).target_index == i);
            assert(on_surface_points(i).parent_patch != -1);
        }
        assert(on_surface_points(i).inside_domain != NOWHERE);
        assert(on_surface_points(i).region != UNMARKED);
    }
}

void Markgrid::check_all_points_are_marked(DblNumMat input_points, 
        NearFieldMap& on_surface_point_map){
    cout << on_surface_point_map.size()  << ", " << input_points.n() << endl;
    int num_points = input_points.n();
    int num_final_on_surface_points = on_surface_point_map.size();
    assert(num_points == num_final_on_surface_points);
#pragma omp parallel for
    for (int i = 0; i < input_points.n(); i++){
        vector<OnSurfacePoint> closest_on_surface_points_to_target = on_surface_point_map[i];
        //for(int j = 0; j < closest_on_surface_points_to_target.size(); j++){
        for(auto p : closest_on_surface_points_to_target){
            //OnSurfacePoint p = closest_on_surface_points_to_target[j];
            if(p.region != FAR && p.inside_domain != OUTSIDE_SPATIAL_GRID){
                assert(p.distance_from_target < DBL_MAX);
                assert(p.target_index == i);
                assert(p.parent_patch != -1);
            }
            assert(p.inside_domain != NOWHERE);
            assert(p.region != UNMARKED);
        }
    }
}



/*NumVec<OnSurfacePoint> Markgrid::mark_targets_region(
        DblNumMat input_points, 
        PatchSamples* patch_samples
        ){

    NumVec<OnSurfacePoint> closest_on_surface_points = 
                            compute_closest_on_surface_points(
                                    input_points, 
                                    (PatchSurfFaceMap*) patch_samples->bdry());

    return  closest_on_surface_points;
}*/

Vec Markgrid::evaluate_fmm_with_constant_density(DblNumMat input_points, PatchSurfFaceMap* face_map){
    double multipole_order = Options::get_int_from_petsc_opts("-bis3d_np");
    double equation = Options::get_int_from_petsc_opts("-kt");
    double spacing= Options::get_double_from_petsc_opts("-bis3d_spacing");
    double refined_spacing= Options::get_double_from_petsc_opts("-bis3d_rfdspacing");

    // increase the multipole order to mark more points
    Options::set_value_petsc_opts("-bis3d_np", "12");
    // using laplace regardless
    Options::set_value_petsc_opts("-kt", "111");
    
    // use refined spacing
    ostringstream stream;
    stream << refined_spacing;
    Options::set_value_petsc_opts("-bis3d_spacing", stream.str());
    stream.str( std::string() );
    stream.clear();



    unique_ptr<SolverGMRESDoubleLayer> solver(new SolverGMRESDoubleLayer("BIS3D_", "bis3d_"));
    solver->bdry() = dynamic_cast<PatchSurf*>(face_map);
    vector<int> patch_partition(face_map->num_patches(), 0);
    solver->patch_partition() = patch_partition;
    solver->dom() = Options::get_int_from_petsc_opts("-dom");
    solver->eqcoefs() = vector<double>(2, 1.);
    solver->setFromOptions(); 
    solver->_compute_refined_surface = false;
    solver->setup();
    
    int source_degrees_of_freedom = 1;
    int target_degrees_of_freedom = 1;
    int num_targets = input_points.n();
    int num_samples = solver->patch_samples()->local_num_sample_points();
    
    Vec targets;
    VecCreateMPI(
            PETSC_COMM_WORLD,
            num_targets*DIM,
            PETSC_DETERMINE,
            &targets);

    vector<int64_t> to_fill(num_targets*DIM, 0);
    // TODO c++11 update
    for(int64_t i =0; i < num_targets*DIM; i++)
        to_fill[i] = i;
    

    VecSetValues(targets,           // vec to insert in 
            num_targets*DIM,        // number of points to insert
            to_fill.data(),        //indices of targets to insert ith value of y into
            input_points.data(),   // data to copy
            INSERT_VALUES);         // insert mode
    //DblNumMat targets_local = get_local_vector(DIM, num_targets, targets);
    VecAssemblyBegin(targets);
    VecAssemblyEnd(targets);
    //targets_local.restore_local_vector();

    Vec density;
    VecCreateMPI(
            PETSC_COMM_WORLD,
            num_samples*source_degrees_of_freedom,
            PETSC_DETERMINE,
            &density);

    VecSet(density, 1.); // set constant density

    Vec potential;
    VecCreateMPI(
            PETSC_COMM_WORLD,
            num_targets*target_degrees_of_freedom,
            PETSC_DETERMINE,
            &potential);

    VecSet(potential, 0.); // set constant density

    solver->fareval(targets, VAR_U, density, potential);
    cout <<"point marking final potential" << endl;
    //VecView(potential, PETSC_VIEWER_STDOUT_WORLD);
    VecDestroy(&density);
    VecDestroy(&targets);



    // restore old multipole order and equation type
    stream << multipole_order;
    Options::set_value_petsc_opts("-bis3d_np", stream.str());
    stream.str( std::string() );
    stream.clear();

    stream << equation;
    Options::set_value_petsc_opts("-kt", stream.str());
    stream.str( std::string() );
    stream.clear();
    
    stream << spacing;
    Options::set_value_petsc_opts("-bis3d_spacing", stream.str());
    stream.str( std::string() );
    stream.clear();
    return potential;
}

NumVec<OnSurfacePoint> Markgrid::compute_closest_on_surface_points(
        DblNumMat input_points, 
        PatchSurfFaceMap* surface){

    int NN = input_points.n();
    SpatialGrid* grid = new SpatialGrid(surface);
   
    // point all targets into the grid
    for (int i = 0; i < NN; i++){
        Point3 target_point(input_points.clmdata(i));
        grid->insert(i, target_point);
    }
    
    NumVec<OnSurfacePoint> closest_points(NN);
    setvalue(closest_points, OnSurfacePoint());
#pragma omp parallel for
    for (int i = 0; i < NN; i++){
        //cout << "i: " << i << endl;
        
        Point3 target_point(input_points.clmdata(i));
        //NearFieldMap patches_near_target_point = 
            //find_patches_near_point(target_point, i, surface, grid);
        NearFieldMap patches_near_target_point = 
            find_patches_closest_to_point(target_point, i, surface, grid);

        int num_possible_closest_points = patches_near_target_point[i].size();
        if(num_possible_closest_points == 0){
            //cout << target_point << endl;
            OnSurfacePoint dummy_point;
            dummy_point.region = FAR;
            closest_points(i) = dummy_point;
        } else { 
            for(int si= 0; si < num_possible_closest_points; si++){
                OnSurfacePoint possible_closest_point = 
                    patches_near_target_point[i][si];
                
                closest_points(i) = 
                    possible_closest_point.distance_from_target < closest_points(i).distance_from_target ?
                    possible_closest_point : closest_points(i);
            }
        }

    }
    return closest_points;
}

OnSurfacePoint Markgrid::closest_point_on_patch_to_target(Point3 target_point, 
        FaceMapSubPatch* patch){
    
    return closest_point_on_patch_to_target(target_point, Point2(.5, .5), patch);
}

OnSurfacePoint Markgrid::closest_point_on_patch_to_target(Point3 target_point, 
        Point2 initial_uv_guess,
        FaceMapSubPatch* patch){
    OnSurfacePoint closest_point;
    
    OnSurfacePoint initial_sample;

    //initial_sample.parent_patch = patch->V();
    initial_sample.parent_patch = patch->_id;
    initial_sample.parametric_coordinates = initial_uv_guess;

    // Find the closest point on the infinite extension of the polynomial on the 
    // face of the patch
    OnSurfacePoint closest_point_on_polynomial_face = 
        closest_on_surface_point_on_extended_patch(
                initial_sample,
                target_point,
                patch,
                PLANE);
    //min_distance = closest_point_on_polynomial_face.distance_from_target;
    bool is_valid;
    patch->is_xy_valid(
            closest_point_on_polynomial_face.parametric_coordinates.array(),
            is_valid);
    // TODO ensure that geometry implies that if closest point is outside
    // domain, then the closest point is along an edge (same for vertices).
    // if not add a distance check
    if(!is_valid) {
        // if the point closest to the patch is outside the defined domain bounds 
        // [0,1]^2 of the patch, find the closest point along each edge
        OnSurfacePoint closest_point_on_polynomial_edge = 
            closest_point_along_patch_edges(
                    target_point,
                    patch);//,closest_point_on_polynomial_face.parametric_coordinates);

        patch->is_xy_valid(
                closest_point_on_polynomial_edge.parametric_coordinates.array(), 
                is_valid);
        if(!is_valid) {
            // if the point closest to the patch is outside of the [0,1] on any of the 
            // edges, compute the distance from the four vertices to the target
            OnSurfacePoint closest_polynomial_patch_corner = 
                closest_point_patch_corners(
                        target_point,
                        patch);

            patch->is_xy_valid(
                    closest_polynomial_patch_corner.parametric_coordinates.array(), 
                    is_valid);

            ebiAssert(is_valid);
            closest_point = closest_polynomial_patch_corner;

        } else {
            // The closest point is contained in one of the patch edges,
            // we're done.
            closest_point = closest_point_on_polynomial_edge;
        }
    } else {
        // The closest point we found is inside the patch face, all done
        closest_point = closest_point_on_polynomial_face;
    }
    return closest_point;
}

OnSurfacePoint Markgrid::select_closest_point_biased_toward_target_patch(
        vector<OnSurfacePoint> closest_on_surface_points, 
        FaceMapSubPatch* target_patch, Point3 target_point,
        PatchSurfFaceMap* face_map){
    cout.precision(16);
    OnSurfacePoint closest_point = 
        find_closest_on_surface_point_in_list(closest_on_surface_points);
    int patch_id = target_patch->V();
    /*
    cout << "computing closest point on " << patch_id << " and target " << closest_point.target_index << endl;
    cout << " absolute closest point to point " << closest_point.target_index<< ": " 
        << closest_point.parent_patch << ", " << closest_point.distance_from_target << endl;
    */
    vector<OnSurfacePoint> nearby_closest_points;
    OnSurfacePoint point_on_target_patch;
    for(auto on_surface_point : closest_on_surface_points){
        //double threshold = target_patch->characteristic_length()*1e-4 + 1e-13;
        double threshold = target_patch->characteristic_length()*1e-1 + 1e-13;

        //cout << "target parent patch: " << on_surface_point.parent_patch<< endl;
        //cout << "target patch dist: " << on_surface_point.distance_from_target << endl;

        double distance_difference = 
            fabs(on_surface_point.distance_from_target 
                    - closest_point.distance_from_target);

        //cout << "target dist from closest: " << distance_difference << ", " <<threshold << endl;
        //cout << "u,v coords: " << closest_point.parametric_coordinates << endl;

        if(distance_difference < threshold){
            nearby_closest_points.push_back(on_surface_point);
        }/* else
            if(on_surface_point.parent_patch == target_patch->V()){
                cout << "TARGET PATCH " << target_patch->V() << " SKIPPED: " << endl;
                cout << "closest dist: " << closest_point.distance_from_target << endl;
                cout << "target patch dist: " << on_surface_point.distance_from_target << endl;
            }*/
    }
    /*cout << "nearby_closest_points: " << endl;
    for(const auto& osp : nearby_closest_points){
        cout << "closest point to point " << osp.target_index<< ": " 
            << osp.parent_patch << ", " << osp.distance_from_target << endl;

    }*/
    //target 30, patch 109, bad point is from patch 54
    if(nearby_closest_points.size() > 1){ // we're on an edge or corner
        OnSurfacePoint closest_point_intermed = closest_point;
        for(auto on_surface_point : nearby_closest_points){

            auto current_patch =  face_map->subpatch(on_surface_point.parent_patch);
            auto closest_patch = face_map->subpatch(closest_point_intermed.parent_patch);
            Point3 current_position;
            Point3 closest_position;

            // find the distance between the closest point and
            // current point
            //current_patch->eval_unsafe(
            current_patch->xy_to_patch_coords(
                    on_surface_point.parametric_coordinates.array(),
                    PatchSamples::EVAL_VL,
                    current_position.array());
            //closest_patch->eval_unsafe(
            closest_patch->xy_to_patch_coords(
                    closest_point_intermed.parametric_coordinates.array(),
                    PatchSamples::EVAL_VL,
                    closest_position.array());

            double dist_closest_to_current =
                (closest_position-current_position).length();
            double dist_closest_to_target =
                (closest_position-target_point).length();
            double dist_current_to_target =
                (current_position-target_point).length();
            double mean_dist_to_target = (dist_current_to_target + dist_closest_to_target)/2.;

            assert(fabs(dist_current_to_target - on_surface_point.distance_from_target) <= 1e-14);
            assert(fabs(dist_closest_to_target - closest_point_intermed.distance_from_target) <= 1e-14);
            /*cout << "dist_closest_to_current:" << dist_closest_to_current << endl;
            cout << "dist_closest_to_target: " << dist_closest_to_target  << endl;
            cout << "dist_current_to_target: " << dist_current_to_target  << endl;
            cout << "mean_dist_to_target:    " << mean_dist_to_target     << endl; 
            cout << mean_dist_to_target*1e-3 + 1e-14 << endl;
            cout << on_surface_point.parent_patch << ", " << target_patch->_id << ", " 
                        << closest_point_intermed.parent_patch  << endl;*/
            // The distance between the closest point and the
            // current is negligible compared to the distance
            // between the target and the closest point
            if(dist_closest_to_current <= mean_dist_to_target*1e-1+ 1e-14){
                //cout << "first cond" <<endl;
                // if this point is on the patch we need,
                // replace the closest point
                if(on_surface_point.parent_patch == target_patch->_id &&
                        closest_point_intermed.parent_patch != target_patch->_id){
                //cout << "second cond" <<endl;
                    closest_point_intermed = on_surface_point;
                }
            }

        }
        closest_point = closest_point_intermed;
    }
    /*cout << "final closest point: " << closest_point.parent_patch << 
             ", " << closest_point.distance_from_target << endl;*/
    return closest_point;
}



NearFieldMap Markgrid::find_points_near_patch(
        DblNumMat input_points,
        FaceMapSubPatch* patch){
    
    NearFieldMap near_points_to_patch;
    int patch_id = patch->V();
    near_points_to_patch[patch_id] = vector<OnSurfacePoint>();

    double min_distance = DBL_MAX;
    cout << input_points.n()<< endl;
#pragma omp parallel for
    for(int i = 0; i < input_points.n(); i++){
        Point3 target_point(input_points.clmdata(i));
        
        OnSurfacePoint closest_point = closest_point_on_patch_to_target(target_point, patch);
        closest_point.parent_patch = patch_id;
        closest_point.target_index = i;


        // if this distance is < .5*characteristic length of the patch, the point is
        // near to the surface. Keep it 
            //cout << "distance from patch: "<< closest_point.distance_from_target << "near-zone: " << .5*patch->characteristic_length() << endl;
        if( closest_point.distance_from_target <  .5*patch->characteristic_length()){
            closest_point.region = NEAR;
            near_points_to_patch[patch_id].push_back(closest_point);
        } 
    }
    return near_points_to_patch;
}

// TODO need patch neighbor information first
//
Markgrid::NearFieldMap Markgrid::find_patches_near_point(
        Point3 target_point,
        int target_index,
        PatchSurfFaceMap* face_map,
        SpatialGrid* grid){
    NearFieldMap near_patches_to_point;
    near_patches_to_point[target_index] = vector<OnSurfacePoint>();


    // TODO update loop over all patches to a grid query for near patches and
    // loop over the result
    vector<uint> near_patches;
    if(grid != NULL){
        near_patches = (*grid->boxes_near_points())[target_index];
        //near_patches = grid->boxes_closest_to_point(target_index);
    } else {
        //near_patches = vector<uint>(face_map->num_patches());
        for(int i = 0; i < face_map->num_patches(); i++){
            near_patches.push_back(i);
        }
    }
    for(int npi = 0; npi < near_patches.size(); npi++){
        int pi = near_patches[npi];
        FaceMapSubPatch* patch = (FaceMapSubPatch*) face_map->patches()[pi];

        OnSurfacePoint closest_point = 
            closest_point_on_patch_to_target(target_point, patch);
        closest_point.parent_patch = pi;// patch->V();
        closest_point.target_index = target_index;


        Point3 ret[3];
        patch->xy_to_patch_coords(closest_point.parametric_coordinates,
                PatchSamples::EVAL_VL|PatchSamples::EVAL_FD,
                (double*) ret);
        Point3 normal = cross(ret[1], ret[2]).dir();
        if (dot(normal, target_point-ret[0]) < 0){
            //cout << "inside" << endl;
            closest_point.inside_domain = INSIDE; 
        }  else {
            //cout << "outside" << endl;
            closest_point.inside_domain = OUTSIDE; 
        }

        if( closest_point.distance_from_target <  .5*patch->characteristic_length()){
            // We've found a near point; we need to know if it's inside or
            // outside; NOTE need more robust in/out test here

            closest_point.region = NEAR;
            near_patches_to_point[target_index].push_back(closest_point);
        } 
    }
    //cout << "number of near_patches: " << near_patches.size() << endl;

    return near_patches_to_point;
}
Markgrid::NearFieldMap Markgrid::find_closest_patch_to_point_aabb_tree(
        Point3 target_point,
        int target_index,
        PatchSurfFaceMap* face_map,
        AABBTree* aabb_tree){
    NearFieldMap near_patches_to_point;
    near_patches_to_point[target_index] = vector<OnSurfacePoint>();

    // find the closest point from the target(distance d away) in the tree
    double start = omp_get_wtime();
    NearFieldMap closest_point_map = 
        find_patches_closest_to_point(target_point, target_index, face_map, aabb_tree);
    
    OnSurfacePoint closest_point = find_closest_on_surface_point_in_list(closest_point_map[target_index]);
    
    // should have returned a single closest patch
    //assert(closest_point_map[target_index].size() == 1);

    // create a bounding box for the sphere of radius d containing the target
    // NOTE we search a little bit further in order to avoid some adversarial
    // edge cases. 
    Point3 bbox_max = target_point + Point3(closest_point.distance_from_target*1.1);
    Point3 bbox_min = target_point - Point3(closest_point.distance_from_target*1.1);
    
    vector<uint> additional_patches = aabb_tree->patches_intersecting_bbox(bbox_min, bbox_max);

    NearFieldMap additional_patches_to_check= 
        compute_closest_points_on_patches(target_point, 
                target_index, 
                face_map,
                additional_patches);
    
    assert(closest_point_map.size() == 1);
    assert(additional_patches_to_check.size() == 1);
    
    for(const auto on_surface_point : additional_patches_to_check[target_index]){
        closest_point_map[target_index].push_back(on_surface_point);
    }

    return closest_point_map;

}

Markgrid::NearFieldMap Markgrid::find_closest_patch_to_point_via_bfs(
        Point3 target_point,
        int target_index,
        PatchSurfFaceMap* face_map,
        SpatialGrid* grid){
    NearFieldMap near_patches_to_point;
    near_patches_to_point[target_index] = vector<OnSurfacePoint>();

    assert(grid != NULL); // we can't do this without the patches in a grid already!

    // Find the closest patch bounding boxes to the target via bfs, terminating at the first
    // level containing >= 1 bounding box.
    NearFieldMap closest_point_map = 
        find_patches_closest_to_point(target_point, target_index, face_map, grid);


    // of the patch bounding boxes found, determine the closest one
    OnSurfacePoint closest_point = find_closest_on_surface_point_in_list(closest_point_map[target_index]);

    /*
     * Since the closest bounding box doesn't necessarily contain the closest
     * patch (e.g. consider an axis aligned flat patch two grid levels away, and a patch 
     * of higher curvature whose bounding box intersects grid level zero or one), 
     * we need determine how many more levels of grid searching need to be
     * done. 
     *
     * This can be accelerated in the future (only need to search
     * ceil(d/box_size - level of bfs termination) ), but for now, we'll search
     * levels i=0, ..., grid_levels_to_search =
     * distance_to_closest_point/grid->box_size().
     */
    int grid_levels_to_search = int(ceil(closest_point.distance_from_target/grid->box_size()))+1;

    // determine which patches that we've already searched in 
    // find_patches_closest_to_point() above
    vector<uint> patches_processed;
    vector<OnSurfacePoint> closest_on_surface_points = closest_point_map[target_index];
    for(int i =0; i < closest_on_surface_points.size(); i++){
        OnSurfacePoint point = closest_on_surface_points[i];
        patches_processed.push_back(point.parent_patch);
    }

    for(int grid_level = 0; grid_level < grid_levels_to_search; grid_level++){
        // Get all bounding boxes on current grid_level
        vector<uint> patch_ids = 
            grid->bounding_boxes_level_l_from_target_point(target_index, grid_level);
        sort(patch_ids.begin(), patch_ids.end());
        
        
        // only compute the new patches
        vector<uint> new_patch_ids;
        std::set_difference(patch_ids.begin(), patch_ids.end(), 
                patches_processed.begin(), patches_processed.end(),
                    std::inserter(new_patch_ids, new_patch_ids.begin()));
        //checking...
        for(auto new_id : new_patch_ids){
            // new id should be in returned by the grid search
            assert(std::find(patch_ids.begin(), patch_ids.end(), new_id) != patch_ids.end());
            // new id should NOT have already been processed.
            assert(std::find(patches_processed.begin(), patches_processed.end(), new_id) == patches_processed.end());
        }
        
      
        // add the ids we've found into the set of patch ids we've already seen
        std::copy(patch_ids.begin(),patch_ids.end(),std::inserter(patches_processed, patches_processed.end()));
        for(auto new_id : new_patch_ids){
            // new id should have been added to patches_processed
            assert(std::find(patches_processed.begin(), patches_processed.end(), new_id) != patches_processed.end());
        }
        sort(patches_processed.begin(), patches_processed.end());
        
        // Compute the closest point on each patch (with a bounding box on this
        // level) to the target.
        NearFieldMap near_patches_to_point_current_level = 
            compute_closest_points_on_patches(
                    target_point, 
                    target_index, 
                    face_map,
                    new_patch_ids);

        // List of on-surface points of nearby patches to target_point
        vector<OnSurfacePoint> on_surface_points_current_level = 
            near_patches_to_point_current_level[target_index];

        // Compute the closest on surface point to target_point
        //OnSurfacePoint closest_point_current_level = 
        //    find_closest_on_surface_point_in_list(on_surface_points_current_level);
        
        for(int i = 0; i < on_surface_points_current_level.size(); i++){
            closest_point_map[target_index].push_back(on_surface_points_current_level[i]);
        }

   }

    return closest_point_map;
}



OnSurfacePoint Markgrid::find_closest_on_surface_point_in_list(
        vector<OnSurfacePoint> on_surface_points){

    OnSurfacePoint closest_point;
    closest_point.distance_from_target = DBL_MAX;
    
    for(int i = 0; i < on_surface_points.size(); i++){
        OnSurfacePoint on_surface_point = on_surface_points[i];
        if(on_surface_point.distance_from_target < closest_point.distance_from_target){
            closest_point = on_surface_point;
        }
    }
    return closest_point;
}

vector<OnSurfacePoint> Markgrid::collect_nearby_on_surface_points(
        vector<OnSurfacePoint> on_surface_points, 
        OnSurfacePoint closest_point, 
        double eps){

    vector<OnSurfacePoint> nearby_points;
    // Find the on-surface points that have an O(eps) absolute difference in
    // distance_from_target from closest_point.distance_from_target
    for(int i = 0; i < on_surface_points.size(); i++){
        OnSurfacePoint on_surface_point = on_surface_points[i];
        double distance_from_target_abs_error = 
            fabs(on_surface_point.distance_from_target - closest_point.distance_from_target);
        double distance_from_target_rel_error = distance_from_target_abs_error/fabs(closest_point.distance_from_target);
        if( distance_from_target_rel_error <= eps){
            nearby_points.push_back(on_surface_point);
        }
    }

    return nearby_points;
}

Markgrid::NearFieldMap Markgrid::find_patches_closest_to_point(
        Point3 target_point,
        int target_index,
        PatchSurfFaceMap* face_map,
        AABBTree* aabb_tree){

    // TODO update loop over all patches to a grid query for near patches and
    // loop over the result
    vector<uint> near_patches;
    if(aabb_tree != NULL){
        uint nearest_patch =aabb_tree->closest_patch_to_point(target_point);
        near_patches.push_back(nearest_patch);

    } else {
        for(int i = 0; i < face_map->num_patches(); i++){
            near_patches.push_back(i);
        }
    }
    double start = omp_get_wtime();
    NearFieldMap near_patches_to_point = 
        compute_closest_points_on_patches(target_point, 
                target_index, 
                face_map,
                near_patches);

    return near_patches_to_point;
}

Markgrid::NearFieldMap Markgrid::find_patches_closest_to_point(
        Point3 target_point,
        int target_index,
        PatchSurfFaceMap* face_map,
        SpatialGrid* grid){
    // TODO update loop over all patches to a grid query for near patches and
    // loop over the result
    vector<uint> near_patches;
    if(grid != NULL){
        near_patches = grid->bounding_boxes_closest_to_point(target_index);

    } else {
        //near_patches = vector<uint>(face_map->num_patches());
        for(int i = 0; i < face_map->num_patches(); i++){
            near_patches.push_back(i);
        }
    }
    NearFieldMap near_patches_to_point = 
        compute_closest_points_on_patches(target_point, 
                target_index, 
                face_map,
                near_patches);


    return near_patches_to_point;
}

Markgrid::NearFieldMap Markgrid::compute_closest_points_on_patches(
        Point3 target_point,
        int target_index,
        PatchSurfFaceMap* face_map,
        vector<uint> near_patches){
    
    return compute_closest_points_on_patches(
            target_point,
            target_index,
            face_map,
            vector<Point2>(near_patches.size(), Point2(.5, .5)),
            near_patches);
}

Markgrid::NearFieldMap Markgrid::compute_closest_points_on_patches(
        Point3 target_point,
        int target_index,
        PatchSurfFaceMap* face_map,
        vector<Point2> initial_uv_guesses,
        vector<uint> near_patches){
    stats.result_plus_equals("total number of points processed", 1);
    assert(initial_uv_guesses.size() == near_patches.size());

    NearFieldMap near_patches_to_point;
    near_patches_to_point[target_index] = vector<OnSurfacePoint>();

    for(int npi = 0; npi < near_patches.size(); npi++){
        int pi = near_patches[npi];
        FaceMapSubPatch* patch = (FaceMapSubPatch*) face_map->patches()[pi];

        double start = omp_get_wtime();
        OnSurfacePoint closest_point = 
            closest_point_on_patch_to_target(target_point, patch);
        stats.result_plus_equals("patch opt time",  omp_get_wtime()- start );
        
        populate_closest_point_data(
                target_index,
                target_point,
                patch,
                closest_point);

        near_patches_to_point[target_index].push_back(closest_point);
    }
    return near_patches_to_point;
}

void Markgrid::populate_closest_point_data(
        int target_index,
        Point3 target_point,
        FaceMapSubPatch* patch,
        OnSurfacePoint& closest_point){
    
    closest_point.parent_patch =  patch->_id;
    closest_point.target_index = target_index;


    Point3 ret[3];
    patch->xy_to_patch_coords(closest_point.parametric_coordinates,
            PatchSamples::EVAL_VL|PatchSamples::EVAL_FD,
            (double*) ret);
    Point3 normal = cross(ret[1], ret[2]).dir();
    double residual_proj_onto_normal = dot(normal, target_point-ret[0]);
    double threshold = dynamic_cast<PatchSurfFaceMap*>(patch->_face_map_patch->bdry())->_on_surface_threshold;
    
    // mark inside/outside
    if (fabs(residual_proj_onto_normal) < threshold){
        closest_point.inside_domain = ON_SURFACE; 
    } else if(dot(normal, target_point-ret[0]) < 0){
        closest_point.inside_domain = INSIDE; 
    } else {
        closest_point.inside_domain = OUTSIDE; 
    }

    double target_accuracy = Options::get_double_from_petsc_opts("-target_accuracy");
    double spacing = Options::get_double_from_petsc_opts("-bis3d_spacing");
    double mean_curvature = patch->mean_curvature(closest_point.parametric_coordinates);
    int n = int(floor(1./spacing))+1;
    // mark near/far
    if(!is_quad_accurate_at_target_fit(n, closest_point.distance_from_target,
            1., mean_curvature, target_accuracy)){
        // We've found a near point; we need to know if it's inside or
        // outside; NOTE need more robust in/out test here

        closest_point.region = NEAR;
    }  else {
        closest_point.region = FAR;
    }

}





OnSurfacePoint Markgrid::closest_point_patch_corners(
        Point3 target_point,
        FaceMapSubPatch* patch){
    
    OnSurfacePoint closest_point_corners;
    closest_point_corners.parent_patch = patch->_id;
    //closest_point_corners.parent_patch = patch->V();
    closest_point_corners.distance_from_target = DBL_MAX;
    
    // Load a vector of OnSurfacePoint's with values at the domain corners
    vector<OnSurfacePoint> corners(4, OnSurfacePoint());
    corners[0].parametric_coordinates = Point2(0., 0.);
    corners[1].parametric_coordinates = Point2(0., 1.);
    corners[2].parametric_coordinates = Point2(1., 0.);
    corners[3].parametric_coordinates = Point2(1., 1.);
    
    for(int i = 0; i < corners.size(); i++){
        OnSurfacePoint corner = corners[i];
        
        Point3 position[1];
        patch->xy_to_patch_coords(corner.parametric_coordinates,
                PatchSamples::EVAL_VL,
                (double*) position);
        
        double corner_to_target_distance = (position[0] - target_point).length();
        if(closest_point_corners.distance_from_target > corner_to_target_distance){
            closest_point_corners = corner;
            closest_point_corners.distance_from_target = corner_to_target_distance;
        }
    }
    return closest_point_corners;
}

OnSurfacePoint Markgrid::closest_point_along_patch_edges(Point3 target_point,
                    FaceMapSubPatch* patch,Point2 initial_guess){
    OnSurfacePoint initial_sample;
    //initial_sample.parent_patch = patch->V();
    initial_sample.parent_patch = patch->_id;

    OnSurfacePoint closest_point_on_any_edge;
    //closest_point_on_any_edge.parent_patch = patch->V();
    closest_point_on_any_edge.parent_patch = patch->_id;
    closest_point_on_any_edge.distance_from_target = DBL_MAX;
   
    vector<OnSurfacePoint> initial_samples(4, initial_sample);
    // y = 0
    initial_samples[0].parametric_coordinates = Point2(.5, 0);
    // y = 1
    initial_samples[1].parametric_coordinates = Point2(.5, 1);
    // x = 0
    initial_samples[2].parametric_coordinates = Point2(0, .5);
    // x = 1
    initial_samples[3].parametric_coordinates = Point2(1, .5);

    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            DescentType descent_type = (DescentType) (i+1); 
            int index = 2*i+j;
            
            OnSurfacePoint closest_point_on_single_edge = 
                closest_on_surface_point_on_extended_patch(
                        initial_samples[index],
                        target_point,
                        patch,
                        descent_type);
            
            if(closest_point_on_single_edge.distance_from_target < 
                    closest_point_on_any_edge.distance_from_target){
                closest_point_on_any_edge = closest_point_on_single_edge;
            }
        }
    }
    return closest_point_on_any_edge;
}


Point2 Markgrid::newton_direction(
        Point3 X_u, Point3 X_v,
        Point3 X_uu, Point3 X_uv, Point3 X_vv,
        Point3 n, Point3 p, DescentType descent_type){

    Point2 grad(-2*dot(p, X_u), -2*dot(p,X_v));
    // We have second derivatives so we can take a Newton step.
    // TODO write out derivative of step direction
    double A = -2.*( dot(p, X_uu) - dot(X_u, X_u));
    double B = -2.*( dot(p, X_uv) - dot(X_u, X_v));
    double C = -2.*( dot(p, X_vv) - dot(X_v, X_v));
    //Hessian  is [A B]
    //            [B C]
    // so the inverse is 
    //      H^-1 = 1./(AC-B^2)*[C -B]
    //                         [-B A]
    // applied to a vector v = (v1,v2), we have
    // H^-1*v = 1./(AC-B^2)*(C*v1 - B*v2, -B*v1 + A*v2)
    
    // For newton's method we need to solve
    // H \delta = -G
    // for hessian H and gradient G.
    
    grad *= -1.;
    double du = C*grad.x() - B*grad.y();
    double dv = -B*grad.x() + A*grad.y();
    
    du *= 1./(A*C-B*B);
    dv *= 1./(A*C-B*B);
    Point2 xyd;
    switch(descent_type){
        case PLANE: 
            xyd = Point2(du, dv);
            break;

        case X_DIRECTION:
            xyd = Point2(grad.x()/A, 0.);
            break;

        case Y_DIRECTION:
            xyd = Point2(0., grad.y()/C);
            break;

        default:
            assert(0);
    }
    return xyd;
}


Point2 Markgrid::gradient_direction(Point3 X_u, Point3 X_v, Point3 n,
        Point3 p, 
        DescentType descent_type){

    double pn = dot(p, n);
    Point3 pp = p - pn*n; // projection of p onto tangent plane 
    Point3 up = cross(n,X_u);
    Point3 vp = cross(n,X_v);
    Point2 xyd;
    switch(descent_type){
        case PLANE: 
            // If we're finding the closest point on a patch face
            // map point into parametric domain

            xyd = Point2(dot(vp,pp)/dot(vp,X_u), dot(up,pp)/dot(up,X_v)); 
            break;
        case X_DIRECTION:
            // If we're finding the closest point on a patch edge in the
            // (1,0) direction
            xyd = Point2(dot(vp,pp)/dot(vp,X_u), 0.); 
            break;
        case Y_DIRECTION:
            // If we're finding the closest point on a patch edge in the
            // (0,1) direction
            xyd = Point2(0., dot(up,pp)/dot(up,X_v)); 
            break;
        default:
            assert(0);
    }
    return xyd;
}

bool Markgrid::stopping_criteria(Point3 u, Point3 v, Point3 p, 
        DescentType descent_type){

    switch(descent_type){
        case PLANE: 
            return (abs(dot(u/u.length(), p)) > 1e-15|| abs(dot(v/v.length(), p)) > 1e-15);
        case X_DIRECTION:
            return abs(dot(u/u.length(), p)) > 1e-15;
        case Y_DIRECTION:
            return abs(dot(v/v.length(), p)) > 1e-15;
        default:
            assert(0);
    }
}

bool is_xy_in_neighborhood_of_domain(Point2 xy){
    return xy.x() >= -.2 && 
           xy.y() >= -.2 && 
           xy.x() <= 1.2 && 
           xy.y() <= 1.2;
}

OnSurfacePoint Markgrid::closest_on_surface_point_on_extended_patch(
        OnSurfacePoint closest_on_surface_sample,
        Point3 target_point,
        FaceMapSubPatch* patch, 
        DescentType descent_type){	 
    OnSurfacePoint closest_on_surface_point;

    // MJM BUG if the closest on surface sample is on the edge or corner of a
    // patch, then there are 2 or 4 samples at the same location; the last one 
    // that was processed may be chosen as the "closest." If this point is not
    // on the same patch as the actual closest on-surface point, then walking
    // along the surface will cause the code to break. Need to fix.
    //int pi = closest_on_surface_sample.parent_patch;
    closest_on_surface_point.parent_patch= closest_on_surface_sample.parent_patch;
    Point2 xyc = closest_on_surface_sample.parametric_coordinates;
    //Copied comment from marknode()
    //p = vector pointing from the ith grid point of regular sampling of 
    //    the FMM box to  closest sample points
    /*        
     *         |          
     *  x ---> o <--- x 
     *         |
     *    <=== |
     *         |
     *
     *  x = possible grid point locations
     *  o = nearest boundary sample point to x
     *  --> = p
     *  <== = surface normal at o
     *
     * Below: n = normal at o.
     * If: p \dot n > 0, x is outside the boundary hence outside
     *     p \dot n = 0, p is perp to n, hence x is closer to the 
     *                   surface at another point assuming geometry 
     *                   is nice. NOTE this is important to 
     *                   consider for complex geometry as an edge case;
     *                   but this will never be exactly zero in practice
     *     p \dot n < 0, x is in the interior
     */

    //MJM BUG if p might need to be negated in order for logic to
    //work geometrically?

    Point3 ret[6];
    iC( patch->xy_to_patch_coords(xyc, 
                PatchSamples::EVAL_VL|
                PatchSamples::EVAL_FD|
                PatchSamples::EVAL_2ND_DERIV, 
                (double*)ret) );

    Point3 p = target_point - ret[0];
    Point3 u = ret[1];
    Point3 v = ret[2]; 
    
    Point3 X_uu = ret[3];
    Point3 X_uv = ret[4];
    Point3 X_vv = ret[5];

    Point3 n = cross(u, v);  n /= n.l2();

    double pn = dot(p, n);
    bool is_valid = true;
    int max_iter = 250;
    
    int iter = 0;
    double eval_time =omp_get_wtime();
    while((iter < max_iter) &&  
            stopping_criteria(u, v, p, descent_type) &&  
            is_xy_in_neighborhood_of_domain(xyc)) { 
        
        Point2 xyd = newton_direction(u, v, X_uu, X_uv, X_vv, n, p, descent_type);


        double step_size = Markgrid::backtracking_line_search(patch, xyc, xyd, target_point); 
        
        // When searching patch interior, optimization can sometimes converge
        // too quickly in one component but not the other, i.e. u is resolved to
        // 1e-13, but v is O(1) wrong. The newton direction for the 2d
        // optimization in this case is correct but too small in magnitude. 
        // We check to see if the step in either direction is too small; if so,
        // so we begin to perform the 1d newton search.
        if(descent_type == PLANE){
            if( fabs(step_size*xyd[0]) <= 1e-10 || fabs(step_size*xyd[1]) <= 1e-10)
            {
                Point2 xydx = newton_direction(u, v, X_uu, X_uv, X_vv, n, p, X_DIRECTION);
                Point2 xydy = newton_direction(u, v, X_uu, X_uv, X_vv, n, p, Y_DIRECTION);
                
        
                double step_sizex = Markgrid::backtracking_line_search(patch, xyc, xydx, target_point); 
                double step_sizey = Markgrid::backtracking_line_search(patch, xyc, xydy, target_point); 

       
                if(fabs(step_sizex*xydx.x()) > 1e-10){ 
                    xyd = xydx; 
                    step_size = step_sizex;
                }
                else if(fabs(step_sizey*xydy.y()) > 1e-10){
                    xyd = xydy; 
                    step_size = step_sizey;
                }
                
            }
        } 
        xyc += step_size * xyd; 

        patch->eval_unsafe(xyc, 
                PatchSamples::EVAL_VL|
                PatchSamples::EVAL_FD|
                PatchSamples::EVAL_2ND_DERIV, 
                (double*)ret);
        p = target_point - ret[0]; 
        u = ret[1];
        v = ret[2];
        X_uu = ret[3];
        X_uv = ret[4];
        X_vv = ret[5];
        n = cross(u, v);
        n /= n.l2(); 
        pn = dot(p, n);
        iter++;
        if( (step_size*xyd).l2()<1e-13){
            break;
        }

    }
    eval_time = omp_get_wtime() - eval_time;

    closest_on_surface_point.parametric_coordinates= xyc;
    closest_on_surface_point.distance_from_target = (ret[0] - target_point).l2();

    return closest_on_surface_point;
}

// Working backtracking line-search
double Markgrid::backtracking_line_search(FaceMapSubPatch* patch, 
        Point2 xyc, 
        Point2 xyd, 
        Point3 target_point){
    double step_size = 1.;
    Point3 ret[3];
    patch->eval_unsafe(xyc, PatchSamples::EVAL_VL|PatchSamples::EVAL_FD, (double*)ret);
    Point3 X_u = ret[1];
    Point3 X_v = ret[2];
    Point3 p = ret[0]-target_point;
    Point2 grad(-2*dot(p, X_u), -2*dot(p,X_v));
    double init_distance = p.l2();
    double min_distance = init_distance;
    double armijo_approx = -dot(xyd,grad);
    double armijo_c = .2;
    
    double best_step = 1e-16;

    Point2 xy_line = xyc + step_size*xyd;
    bool is_valid = true;

    bool break_flag = false;
    int count = 0;
    while(step_size > 1e-16){
        patch->eval_unsafe(
                xy_line,
                PatchSamples::EVAL_VL,
                (double*)ret);
        //stats.result_plus_equals("eval_unsafe", eend-estart);

        // if the current step size satisfies armijo line search criteria,
        // update the best step_size
        double distance = (ret[0] - target_point).l2();
        min_distance = init_distance + armijo_c*armijo_approx*step_size;
        
        if(distance <  min_distance){
            best_step = step_size;
            min_distance = distance;
            break_flag = true;
        } 

        count++;
        if(break_flag)
            break;

        step_size /= 2.;

        xy_line = xyc + step_size* xyd;
    }
    return best_step;
    
}


END_EBI_NAMESPACE
