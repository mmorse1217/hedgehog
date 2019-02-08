#include "../catch.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "bdry3d/patch_samples.hpp"
#include "bdry3d/p4est_interface.hpp"
#include "bie3d/markgrid.hpp"
#include "bie3d/spatial_grid.hpp"
#include <algorithm>

using namespace Ebi;
using namespace Markgrid;


TEST_CASE("Spatial grid", "[geom][grid]"){
    const Point3 spatial_min(-1., -1., -1.);
    const Point3 spatial_max(1., 1., 1.);
    double box_length = .1; // =>10x10x10 boxes covering [-1,1]^3

    SpatialGrid grid(box_length, spatial_min, spatial_max);

    int num_boxes = grid.num_boxes(); 
    
     SECTION("Test spatial grid construction"){

        for(uint i = 0; i < grid.grid_boxes().size(); i++){
            Box* grid_box = grid.grid_boxes(i);

            // ensure a grid box is of the proper size
            CHECK( fabs((grid_box->max().x() - grid_box->min().x()) - box_length) <= 1e-15);
            CHECK( fabs((grid_box->max().y() - grid_box->min().y()) - box_length) <= 1e-15);
            CHECK( fabs((grid_box->max().z() - grid_box->min().z()) - box_length) <= 1e-15);

            // ensure it's actually in the domain of interest
            CHECK(grid_box->max() <= grid.max());
            CHECK(grid_box->min() >= grid.min());
        }

        // Randomly select points in space, see what grid box they land in, and
        // ensure that the grid box that is returned actually contains the point
        int num_tests = 100; 
        for(int i = 0; i < num_tests; i++){
            // Random point in [0,1]^3
            Point3 random_point(drand48(), drand48(), drand48());
            
            // Rescale to [spatial_min, spatial_max] \subset R^3
            random_point.x() *= (spatial_max.x() - spatial_min.x());
            random_point.y() *= (spatial_max.y() - spatial_min.y());
            random_point.z() *= (spatial_max.z() - spatial_min.z());

            random_point += spatial_min;
            
            // see where this point actually lands in the grid
            Box* box = grid(random_point);
            CHECK(box->max().x() >= random_point.x());
            CHECK(box->min().x() <= random_point.x());
            
            CHECK(box->max().y() >= random_point.y());
            CHECK(box->min().y() <= random_point.y());
            CHECK(box->max().z() >= random_point.z());
            CHECK(box->min().z() <= random_point.z());

        }
        cout << "test 1" << endl;
    }

    SECTION("Test immediate grid box neighbors"){
        // Check the neighbors of each grid box and ensure that they cover the
        // proper spatial region
        for(int i= 0; i < num_boxes; i++){
            for(int j= 0; j < num_boxes; j++){
                for(int k= 0; k < num_boxes; k++){
                    Box* box = grid(i,j,k);

                    // find immediate neighbors i.e. one box length away from
                    // box
                    vector<Box*> neighbor_boxes = grid.neighbors(box);
                    
                    for(uint ni =0 ; ni < neighbor_boxes.size(); ni++){
                        
                        double distance_from_grid_box = box_length;
                        Box* neighbor = neighbor_boxes[ni];

                        // the box is only a neighbor if the distance from the
                        // min point to the target box's min point is equal to
                        // the box length
                        Point3 distance = abs(box->min() - neighbor->min());
                        bool is_neighbor = 
                            fabs(distance.x() - box_length) < 1e-15 ||
                            fabs(distance.y() - box_length) < 1e-15 ||
                            fabs(distance.z() - box_length) < 1e-15;
                        
                        CHECK(is_neighbor);
                    }
                }
            }
        }
        cout << "test 2" << endl;
    }
    
    SECTION("Test further grid box neighbors"){
        // Check the neighbors of each grid box and ensure that they cover the
        // proper spatial region
        for(int l = 0; l < 4; l++){
            for(int i= 0; i < num_boxes; i++){
                for(int j= 0; j < num_boxes; j++){
                    for(int k= 0; k < num_boxes; k++){
                        Box* box = grid(i,j,k);
                        // get the neighbors to box on level l 
                        // l=0: self,
                        // l=1: immediate neighbors
                        // l>=2: neighbors distance l from box in grid
                        vector<Box*> neighbor_boxes = grid.neighbors(box, l);

                        for(uint ni =0 ; ni < neighbor_boxes.size(); ni++){

                            Box* neighbor = neighbor_boxes[ni];
                            double distance_from_grid_box = double(l)*box_length;


                            // the box is only a neighbor if the distance from its
                            // min point to box's min point is equal to
                            // the grid box length
                            Point3 distance = abs(box->min() - neighbor->min());
                            bool is_neighbor = 
                                fabs(distance.x() - distance_from_grid_box) < 1e-15 ||
                                fabs(distance.y() - distance_from_grid_box) < 1e-15 ||
                                fabs(distance.z() - distance_from_grid_box) < 1e-15;

                            CHECK(is_neighbor);
                        }
                    }
                }
            }
        }
        cout << "test 3" << endl;
    }
    
    SECTION("Test spatial grid point insertion"){
        // insert a point in each box, then ask for the box containing that
        // point and check that the point got stored in the right box
        DblNumMat target_points(3, num_boxes*num_boxes*num_boxes);
        for(int i = 0; i < num_boxes; i++){
            for(int j = 0; j < num_boxes; j++){
                for(int k = 0; k < num_boxes; k++){
                    int index = num_boxes*num_boxes*i + num_boxes*j + k;
                    double step = 2./double(num_boxes);
                    
                    target_points(0,index) = -1. + box_length/2. + i*step;
                    target_points(1,index) = -1. + box_length/2. + j*step;
                    target_points(2,index) = -1. + box_length/2. + k*step;
                    
                    grid.insert(index, Point3(target_points.clmdata(index)));
                }
            }
        }

        // Now make sure the point actually got placed in the box you think it
        // should have
        for(int i = 0; i < target_points.n(); i++){ 
            Point3 target(target_points.clmdata(i));
            Box* box = grid(target);

            CHECK(box->_points[0] == i);
            CHECK(box->_points.size() == 1);

            // make sure the point actually lives in this box
            CHECK(box->min() < target);
            CHECK(box->max() > target);
        }
        
        cout << "test 4" << endl;
    }
    
    SECTION("Test bounding box insertion"){
        // insert a bounding box slightly smaller than each grid box, then ask for the grid boxes containing that
        // bounding box and check that the bounding box got stored in the right grid boxes
        const double eps = 1e-8;
        for(uint i = 0; i < grid.grid_boxes().size(); i++){
            Box* box = grid.grid_boxes(i);

            // slightly shrink the bounding box to fit inside the grid box
            Point3 bounding_box_max(box->max().x() - eps, 
                                    box->max().y() - eps,
                                    box->max().z() - eps);
            Point3 bounding_box_min(box->min().x() + eps, 
                                    box->min().y() + eps,
                                    box->min().z() + eps);
            BoundingBox* bbox = new BoundingBox(bounding_box_min, bounding_box_max);

            grid.insert(i, bbox);
        }

        // get the ith grid box and make sure it contains the ith bounding box
        for(uint i = 0; i < grid.grid_boxes().size(); i++){
            Box* box = grid.grid_boxes(i);
            vector<uint> bounding_boxes = box->_bounding_boxes;
            CHECK(bounding_boxes[0] == i);   
            CHECK(bounding_boxes.size() == 1);   
        }
        cout << "test 5" << endl;
    }
    
    SECTION("Test neighbor index generation"){
        // Check that the indices on level l are actually that far away from
        // (i,j,k)
        for(int level = 1; level < 4; level++){
            int d = 2*level+1;
            for(int i = 0; i < num_boxes; i++){
                for(int j = 0; j < num_boxes; j++){
                    for(int k = 0; k < num_boxes; k++){
                        Index3 index(i, j, k);
                        vector<Index3> neighbor_indices = 
                            grid.generate_neighbor_indices(index, level);
                        for(size_t m = 0; m < neighbor_indices.size(); m++){
                            Index3 neighbor = neighbor_indices[m];
                            Index3 delta = abs(index  - neighbor);

                            int max = INT_MIN;
                            for(int d=0; d < 3; d++)
                                max = max < delta(d) ? delta(d) : max;
                            CHECK(max == level);
                        }

                    }
                }
            }
        }
        cout << "test 5.5" << endl;
    }
    
    SECTION("Test neighbor box generation"){
        // TODO kill this test it's redundant
        int min_index = 0;
        int max_index = num_boxes;
        for(int level = 1; level < 4; level++){
            for(int i = min_index; i < max_index; i++){
                for(int j = min_index; j < max_index; j++){
                    for(int k = min_index; k < max_index; k++){
                        Index3 index(i, j, k);
                        Box* box = grid(i, j, k);

                        vector<Box*> neighbors= 
                            grid.neighbors(box, level);
                        for(size_t m = 0; m < neighbors.size(); m++){
                            Box* neighbor= neighbors[m];
                            Index3 delta = abs(index  - neighbor->index());

                            int max = INT_MIN;
                            for(int d=0; d < 3; d++)
                                max = max < delta(d) ? delta(d) : max;
                            CHECK(max == level);
                        }

                    }
                }
            }
        }
        cout << "test 6" << endl;
    }

    SECTION("Test large bounding box insertion"){
        // Insert bounding boxes generated by applying a buffer to bbox
        // bounds of size k*box_length to an existing grid box. Then, grid boxes
        // that are neighbors on level <= k will have the same bo
        // 
        vector<int> box_ids(num_boxes);
        for(int i = 0; i < num_boxes; i++)
            box_ids[i] = i;
       
        int k = 2; 
        // get a bounding box 
        for(int it = 0; it < box_ids.size(); it++){
            int index = box_ids[it];
            
            //get the (index, index,index)th box
            Box* box = grid(index,index,index);

            Point3 bounding_box_min = box->min();
            Point3 bounding_box_max = box->max();
            
            Point3 box_center = (bounding_box_max + bounding_box_min)*.5;
           
            // shift to origin
            bounding_box_min -= box_center;
            bounding_box_max -= box_center;
            
            // scale by k
            bounding_box_min *= double(k);
            bounding_box_max *= double(k);
           
            // shift back to box center
            bounding_box_min += box_center;
            bounding_box_max += box_center;

            // make a new bounding box centered on the it-th grid box k times
            // larger than the original grid box.
            BoundingBox* bounding_box = new BoundingBox(bounding_box_min, bounding_box_max);
            
            // ... and put it into the grid
            grid.insert(it, bounding_box);
            for(int level = 0; level < k; level++){
                // now, find the neighbors of the grid box that generated the bounding box
                vector<Box*> neighbors= 
                    grid.neighbors(box, level);

                // make sure that the bounding box id ended up in those grid
                // boxes
                for(int ni = 0; ni < neighbors.size(); ni++){
                    Box* neighbor_box = neighbors[ni];
                    vector<uint> bounding_boxes = neighbor_box->_bounding_boxes;
                    CHECK(find(bounding_boxes.begin(), bounding_boxes.end(), it) != bounding_boxes.end());
                }
            }

            // check all other grid boxes that the bounding box didn't land
            // in those by mistake
            for(int level = k; level < num_boxes; level++){
                // find the neighbors of the grid box that generated the 
                // bounding box, that are strictly further away than the
                // bounding box size
                vector<Box*> neighbors= 
                    grid.neighbors(box, level);

                // make sure that the bounding box id isn't in those grid boxes
                for(int ni = 0; ni < neighbors.size(); ni++){
                    Box* neighbor_box = neighbors[ni];
                    vector<uint> bounding_boxes = neighbor_box->_bounding_boxes;
                    CHECK(find(bounding_boxes.begin(), bounding_boxes.end(), it) == bounding_boxes.end());
                }
            }

        }
        cout << "test 6" << endl;

    }

    SECTION("Test grid containing axis-aligned plane"){

        // Sample an axis aligned plane of the grid and
        // insert the points into the grid (sampling rate = grid split rate, so
        // samples will live on grid box corners). bounding box is known, so we
        // can ask for points intersecting bounding box.
      
        NumMatrix sample_points(3, pow(num_boxes,3));
        setvalue(sample_points, 0.);

        for(int d = 0; d < 3; d++){
            for(int i = 0; i < num_boxes; i++){
                for(int j = 0; j < num_boxes; j++){
                    int index = d*num_boxes*num_boxes + i*num_boxes + j;

                    Point3 sample(sample_points.clmdata(index));
                    switch(d){
                        case 0: // x = 0
                            sample.y() = i*box_length - 1. + 1e-15;
                            sample.z() = j*box_length - 1. + 1e-15;
                            break;
                        case 1: // y = 0
                            sample.x() = i*box_length - 1. + 1e-15;
                            sample.z() = j*box_length - 1. + 1e-15;
                            break;
                        case 2: // z = 0
                            sample.x() = i*box_length - 1. + 1e-15;
                            sample.y() = j*box_length - 1. + 1e-15;
                            break;
                        default:
                            assert(0);
                    }

                    grid.insert(index, sample);

                }
            }
            
            Point3 bounding_box_min(-1); 
            Point3 bounding_box_max(1);
                switch(d){
                    case 0: // x = 0
                        bounding_box_min.x() = -1e-14;
                        bounding_box_max.x() = 1e-14;
                        break;
                    case 1: // y = 0
                        bounding_box_min.y() = -1e-14;
                        bounding_box_max.y() = 1e-14;
                        break;
                    case 2: // z = 0
                        bounding_box_min.z() = -1e-14;
                        bounding_box_max.z() = 1e-14;
                        break;
                    default:
                        assert(0);
                }
            auto plane_bounding_box = new BoundingBox(bounding_box_min, bounding_box_max);
            
            grid.insert(d, plane_bounding_box);

            map<uint, vector<uint> > boxes_containing_samples = grid.points_near_boxes();
            //map<uint, vector<uint> >* samples_in_boxes = grid.boxes_near_points(true);
            map<uint, vector<uint> >* samples_in_boxes = grid.boxes_near_points(true);
            
            assert(samples_in_boxes);
           
            vector<uint> samples_in_bounding_box = boxes_containing_samples[d];
            for(int i = 0; i < num_boxes; i++){
                for(int j = 0; j < num_boxes; j++){
                    int index = d*num_boxes*num_boxes + i*num_boxes + j;
                    
                    Point3 sample(sample_points.clmdata(index));

                    // Ensure that the bounding box of the plane actually contains all the
                    // points we sampled from it, both in terms of the grid
                    // query that asks for points inside a bounding box and in
                    // physical space.
                    CHECK(find(samples_in_bounding_box.begin(), samples_in_bounding_box.end(), index) 
                            != samples_in_bounding_box.end());
                    
                    CHECK(sample < plane_bounding_box->max());
                    CHECK(sample > plane_bounding_box->min());

                    // Also check that the index-th sample point lives in the
                    // same grid box as its bounding box (the dual of the above
                    // test);
                    vector<uint> bounding_boxes_near_samples = (*samples_in_boxes)[index];
                    CHECK(find(bounding_boxes_near_samples.begin(), bounding_boxes_near_samples.end(), d) 
                            != bounding_boxes_near_samples.end());

                }
            }
        }
        cout << "test 7 " << endl;
        
    }

    SECTION("Test Point and patch bounding box insertion"){
        // Make a face-map.  Sample each patch at it's corners .
        // Insert each patch and samples into the grid. check to make sure the
        // points sampled each patch fall into the same grid box as the generating  
        // patch (i.e. the patch bounding box for the generating patch should
        // intersect the same grid boxes containing samples).

        Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl");
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
        //Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");
        Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
        Options::set_value_petsc_opts("-bd3d_facemap_fit_accuracy", "1e-3");
        Options::set_value_petsc_opts("-bis3d_spacing", "1");
        Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
        Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor", "0.");

        PatchSurfFaceMap* face_map = new PatchSurfFaceMap("BD3D_", "bd3d_");
        face_map->_surface_type = PatchSurfFaceMap::BLENDED;
        face_map->setFromOptions();
        face_map->setup();
        face_map->refine_test();


        // make a grid using face-map patches and insert all near-zone 
        // bounding boxes
        SpatialGrid grid(face_map);

        // sample the face-map
        PatchSamples* patch_samples = new PatchSamples("", "");
        patch_samples->bdry() = face_map;
        vector<int> patch_partition(face_map->patches().size(), 0);
        patch_samples->patch_partition() = patch_partition;
        patch_samples->setup();

        int num_samples = Petsc::get_vec_size(patch_samples->sample_point_3d_position())/DIM;
        DblNumMat sample_points_local = 
            get_local_vector(DIM, num_samples, patch_samples->sample_point_3d_position());

        // insert each sample
        for(int i = 0; i < num_samples; i++){
            Point3 sample(sample_points_local.clmdata(i));
            grid.insert(i, sample);
        }

        // check that the samples on a given patch get hashed to the same grid
        // box as the bounding box for the patch that generated them.
        double spacing = Options::get_double_from_petsc_opts("-bis3d_spacing");
        int num_samples_per_patch = floor(1./spacing)+1;
        num_samples_per_patch *= num_samples_per_patch;

        for(int pi = 0; pi < face_map->patches().size(); pi++){
            for(int i = 0; i < num_samples_per_patch; i++){
                int index = pi*num_samples_per_patch + i;

                Point3 sample(sample_points_local.clmdata(index));
                Box* grid_box = grid(sample);
                vector<uint> bounding_boxes_in_grid_box = grid_box->_bounding_boxes;

                CHECK(find(bounding_boxes_in_grid_box.begin(),
                            bounding_boxes_in_grid_box.end(),
                            pi) != bounding_boxes_in_grid_box.end());
            }
        }
        sample_points_local.restore_local_vector(); 




    }

    /*SECTION("Point and patch bounding box insertion"){
        // Make a face-map.  Sample each patch at it's corners .
        // Insert each patch and samples into the grid. check to make sure the
        // points sampled each patch fall into the same grid box as the generating  
        // patch (i.e. the patch bounding box for the generating patch should
        // intersect the same grid boxes containing samples).
       
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl");
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
        //Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "0");
        Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
        Options::set_value_petsc_opts("-bis3d_spacing", ".1");

        PatchSurfFaceMap* face_map = new PatchSurfFaceMap("BD3D_", "bd3d_");
        face_map->_surface_type = PatchSurfFaceMap::BLENDED;
        face_map->setFromOptions();
        face_map->setup();
        face_map->refine_test();
        

        // make a grid using face-map patches and insert all near-zone 
        // bounding boxes
        SpatialGrid grid(0.0938037, Point3(-1,-1,-1), Point3(1,1,1));
        
        // sample the face-map
        PatchSamples* patch_samples = new PatchSamples("", "");
        patch_samples->bdry() = face_map;
        vector<int> patch_partition(face_map->patches().size(), 0);
        patch_samples->patch_partition() = patch_partition;
        patch_samples->setup();
        int patch_id = 659 ;
        int num_samples = patch_samples->sample_point_3d_position(patch_id).n();
        DblNumMat sample_points_local = 
            patch_samples->sample_point_3d_position(patch_id);
        cout << sample_points_local << endl;
        
        //FaceMapSubPatch* p = (FaceMapSubPatch*) face_map->patches()[patch_id];
        FaceMapSubPatch* p = (FaceMapSubPatch*) face_map->patches()[patch_id];
        Point3 patch_min, patch_max;
        p->bounding_box(patch_min, patch_max);
        BoundingBox* b = new BoundingBox(patch_min, patch_max);
        grid.insert(patch_id, b);
        // insert each sample
        for(int i = 0; i < num_samples; i++){
            Point3 sample(sample_points_local.clmdata(i));
            grid.insert(i, sample);
        }

        
        // check that the samples on a given patch get hashed to the same grid
        // box as the bounding box for the patch that generated them.
        double spacing = Options::get_double_from_petsc_opts("-bis3d_spacing");
        int num_samples_per_patch = floor(1./spacing)+1;
        num_samples_per_patch *= num_samples_per_patch;
        for(int pi = patch_id; pi < patch_id+1; pi++){
            for(int i = 0; i < num_samples_per_patch; i++){
                //int index = pi*num_samples_per_patch + i;
                int index = i;
                
                Point3 sample(sample_points_local.clmdata(index));
                Box* grid_box = grid(sample);
                vector<uint> bounding_boxes_in_grid_box = grid_box->_bounding_boxes;
                vector<uint> _points_in_grid_box = grid_box->_points;
                cout << "sample" << endl;
                cout << sample<< endl;
                cout << "grid box" << endl;
                cout << grid_box->index()<< endl;
                cout << grid_box->min()<< endl;
                cout << grid_box->max()<< endl;
                cout << "patch bbox" << endl;
                cout << patch_min<< endl;
                cout << patch_max<< endl;

                cout<< "patch " << pi << ": ";
                for(int k = 0; k < bounding_boxes_in_grid_box.size(); k++)
                    cout << bounding_boxes_in_grid_box[k] << ", " ;
                cout << endl;
                cout<< "patch " << pi << ": ";
                for(int k = 0; k < _points_in_grid_box.size(); k++)
                    cout << _points_in_grid_box[k] << ", " ;
                cout << endl << endl; 
                if (find(bounding_boxes_in_grid_box.begin(),
                            bounding_boxes_in_grid_box.end(),
                            pi) == bounding_boxes_in_grid_box.end()){
                    cout << "False," << endl;
                } else {
                    cout << "True," << endl;

                }
                CHECK(find(bounding_boxes_in_grid_box.begin(),
                            bounding_boxes_in_grid_box.end(),
                            pi) != bounding_boxes_in_grid_box.end());
            }
        }
        sample_points_local.restore_local_vector(); 
        cout << "grid box size: " << grid.box_size()<< endl;
    }*/
}

TEST_CASE("Test grid search", "[geom][grid][search]"){

    // Make a face-map.  Sample each patch at it's corners .
    // Insert each patch and samples into the grid. check to make sure the
    // points sampled each patch fall into the same grid box as the generating  
    // patch (i.e. the patch bounding box for the generating patch should
    // intersect the same grid boxes containing samples).

    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    //Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");
    Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
    Options::set_value_petsc_opts("-bd3d_facemap_fit_accuracy", "1e-3");
    Options::set_value_petsc_opts("-bis3d_spacing", "1");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor", "0.");

    PatchSurfFaceMap* face_map = new PatchSurfFaceMap("BD3D_", "bd3d_");
    face_map->_surface_type = PatchSurfFaceMap::BLENDED;
    face_map->setFromOptions();
    face_map->setup();
    face_map->refine_test();


    // make a grid using face-map patches and insert all near-zone 
    // bounding boxes
    SpatialGrid grid(face_map);

    // sample the face-map
    PatchSamples* patch_samples = new PatchSamples("", "");
    patch_samples->bdry() = face_map;
    vector<int> patch_partition(face_map->patches().size(), 0);
    patch_samples->patch_partition() = patch_partition;
    patch_samples->setup();

    int num_samples = Petsc::get_vec_size(patch_samples->sample_point_3d_position())/DIM;
    DblNumMat sample_points_local = 
        get_local_vector(DIM, num_samples, patch_samples->sample_point_3d_position());

    // insert each sample
    for(int i = 0; i < num_samples; i++){
        Point3 sample(sample_points_local.clmdata(i));
        grid.insert(i, sample);
    }

    // check that the samples on a given patch get hashed to the same grid
    // box as the bounding box for the patch that generated them.
    double spacing = Options::get_double_from_petsc_opts("-bis3d_spacing");
    int num_samples_per_patch = floor(1./spacing)+1;
    num_samples_per_patch *= num_samples_per_patch;
    int num_patches = face_map->patches().size();
    for(int pi = 0; pi < num_patches; pi++){
        for(int i = 0; i < num_samples_per_patch; i++){
            int index = pi*num_samples_per_patch + i;

            Point3 sample(sample_points_local.clmdata(index));
            Box* grid_box = grid(sample);
            vector<uint> bounding_boxes_in_grid_box = grid_box->_bounding_boxes;

            CHECK(find(bounding_boxes_in_grid_box.begin(),
                        bounding_boxes_in_grid_box.end(),
                        pi) != bounding_boxes_in_grid_box.end());
        }
    }
    sample_points_local.restore_local_vector(); 


    // generate qbkix targets from samples
    vector<int> indices(1, Options::get_int_from_petsc_opts("-near_interpolation_num_samples") - 1);
    vector<int> patch_ids;

    for(int i =0; i < num_patches; i++)
        patch_ids.push_back(i);

    Vec qbkix_points = patch_samples->generate_qbkix_points_from_sample_points(indices, patch_ids);

    int64_t num_qbkix_samples;
    VecGetSize(qbkix_points, &num_qbkix_samples);
    num_qbkix_samples /= DIM;
    DblNumMat qbkix_points_local = get_local_vector(DIM, num_qbkix_samples, qbkix_points);
    grid.insert_points(qbkix_points_local);
    //qbkix_points_local.restore_local_vector();
    //for(int i = 0; i < qbkix_points_local.n(); i++){
    vector<Markgrid::NearFieldMap> temp(qbkix_points_local.n());
#pragma omp parallel for
    for(int pi = 0; pi < num_patches;  pi++){
        for(int i = 0; i < num_samples_per_patch;  i++){
            int index = pi*num_samples_per_patch + i;
            cout << "i: " << index << endl;

            Point3 qbkix_point(qbkix_points_local.clmdata(index));
            Markgrid::NearFieldMap closest_patches = 
                Markgrid::find_closest_patch_to_point_via_bfs(
                        qbkix_point, index, face_map, &grid);
            temp[index] = closest_patches;
        }
    }
    // aggregate
    Markgrid::NearFieldMap closest_points;
    for(const auto& map : temp){
        closest_points.insert(map.begin(), map.end());
    }

    // check that the parent patch is found in the grid search
    for(int pi = 0; pi < num_patches;  pi++){
        for(int i = 0; i < num_samples_per_patch;  i++){
            int index = pi*num_samples_per_patch + i;
            vector<OnSurfacePoint> closest_patches = closest_points[index];

            vector<int> patches_found;
            for(auto point : closest_patches){
                patches_found.push_back(point.parent_patch);
            }
            CHECK(std::find(patches_found.begin(), patches_found.end(), pi) != patches_found.end());
        }
    }

    // check that the parent patch is found in the grid search
    NumVec<OnSurfacePoint> out (num_samples_per_patch*num_patches);
    vector<int> pids(num_patches);
    for(int pi = 0; pi < num_patches;  pi++){
        pids[pi] = pi;
        for(int i = 0; i < num_samples_per_patch;  i++){
            int index = pi*num_samples_per_patch + i;

            auto patch = FaceMapSubPatch::as_subpatch(face_map->patch(pi));
            Point3 qbkix_point(qbkix_points_local.clmdata(index));
            
            vector<OnSurfacePoint> closest_patches = closest_points[index];
            OnSurfacePoint closest_point =
                Markgrid::select_closest_point_biased_toward_target_patch(
                        closest_patches, 
                        patch, qbkix_point,
                        face_map);

            cout << "near patches: " ;
            for(const auto& o: closest_patches){
                cout << "(" << o.parent_patch << ", " << o.distance_from_target << "), ";
            }
            cout << "computed closest point: " << closest_point.parent_patch << ", " << closest_point.distance_from_target << endl;
            cout << endl;
            CHECK(closest_point.parent_patch == pi);
            if(closest_point.parent_patch != pi){
                closest_point.parent_patch = -200;
            }
            out(index) = closest_point;
        }
    }
            dump_vtk_data_for_paraview(qbkix_points_local,
                    out, 1000,
                    pids, face_map, "data/");



}

