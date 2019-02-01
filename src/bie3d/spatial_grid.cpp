#include "spatial_grid.hpp"

BEGIN_EBI_NAMESPACE

using namespace Markgrid;

void SpatialGrid::discretize(double box_size){
    // discretize [-1,1]^3 into ((1- -1)/box_size)^3 number of boxes
    box_size = std::max(box_size, .005);
    _box_size = box_size;
    cout << "grid box size: " << box_size << endl;
    _num_boxes = int(ceil((1. - (-1.))/box_size));
    _num_boxes = _num_boxes % 2 == 1 ? _num_boxes + 1 : _num_boxes;
    cout << "_num_boxes: " << _num_boxes<< endl;
    //_grid_boxes.reserve(pow(_num_boxes,3));
    /*for(int i = 0; i < _num_boxes; i++){
        for(int j = 0; j < _num_boxes; j++){
            for(int k = 0; k < _num_boxes; k++){
                //int index = i*(_num_boxes*_num_boxes)+j*_num_boxes + k;
                Point3 box_min(
                        index_to_coord(i), 
                        index_to_coord(j),
                        index_to_coord(k));
                Point3 box_max(
                        index_to_coord(i+1), 
                        index_to_coord(j+1),
                        index_to_coord(k+1));
                _grid_boxes.push_back(new Box(i,j,k, box_min, box_max));
                //_grid_boxes[index] =new Box(i,j,k, box_min, box_max);
            }
        }
    }*/

}

void SpatialGrid::initialize_grid_box(int i,int j, int k, Box*& box){
    assert(box == nullptr); // only initialize an empty box pointer
    Point3 box_min(
            index_to_coord(i), 
            index_to_coord(j),
            index_to_coord(k));
    Point3 box_max(
            index_to_coord(i+1), 
            index_to_coord(j+1),
            index_to_coord(k+1));
    box = new Box(i,j,k,box_min, box_max);

}
Box* SpatialGrid::grid_boxes(uint i){
    Box* b = _grid_box_map[i];
    if(!b){
        int ai = i % _num_boxes;
        i /= _num_boxes;
        int aj = i % _num_boxes;
        i /= _num_boxes;
        int ak= i;
        initialize_grid_box(ai,aj,ak,b);
    }
    return b;
}

SpatialGrid::SpatialGrid(vector<Patch*> patches):
    SpatialGrid(patches, DblNumMat(3,0))
{

}
SpatialGrid::SpatialGrid(PatchSurfFaceMap* face_map):
    SpatialGrid(face_map->patches(), DblNumMat(3,0))
{
}

SpatialGrid::SpatialGrid(DblNumMat points):
    SpatialGrid(vector<Patch*>(), points)
{
}
SpatialGrid::SpatialGrid(PatchSurfFaceMap* face_map, DblNumMat points, int num_boxes):
    SpatialGrid(face_map->patches(), points, num_boxes)
{

}
SpatialGrid::SpatialGrid(vector<Patch*> patches, DblNumMat points, int num_boxes):
    _max(Point3(1.,1., 1.)), 
    _min(Point3(-1., -1., -1.)), 
    _bounding_boxes_near_points(NULL)
{
    if(!patches.empty()){
        // discretize with grid size proportional to patch size if we can
        insert_face_map_patches(patches);
    } else {
        // discrete with either default number of boxes per dim or
        // user-specified
        _box_size = (1. - (-1.))/double(num_boxes);
        discretize(_box_size);

    }

    insert_points(points);
}

void SpatialGrid::insert_points(DblNumMat points){
    for (int i = 0; i < points.n(); i++){
        Point3 point(points.clmdata(i));
        // make sure the target point to insert actually lands in the grid box
        // TODO make omp parallel loop
        if(point < _max && point > _min){
            insert(i, point);
        } else {
            assert(0);
        }
    }
} 

void SpatialGrid::insert_bounding_boxes(vector<BoundingBox*> bounding_boxes){
    int num_bounding_boxes = bounding_boxes.size();
    vector<vector<int> > grid_boxes_intersecting_bounding_boxes(num_bounding_boxes);
    cout << "inserting bounding boxes" << endl;
#pragma omp parallel for
    for(int pi= 0; pi < num_bounding_boxes; pi++){
        double eps = 4e-16;
        BoundingBox* box = bounding_boxes[pi];
        int key = pi;
        Point3 bbox_min = box->min();
        Point3 bbox_max = box->max();

        if(bbox_min.x() <= -1.)
            bbox_min.x() = -1. + eps;
        if(bbox_min.y() <= -1.)
            bbox_min.y() = -1. + eps;
        if(bbox_min.z() <= -1.)
            bbox_min.z() = -1. + eps;

        if(bbox_max.x() >= 1.)
            bbox_max.x() = 1. - eps;
        if(bbox_max.y() >= 1.)
            bbox_max.y() = 1. - eps;
        if(bbox_max.z() >= 1.)
            bbox_max.z() = 1. - eps;

        Index3 min_box = point_to_index(bbox_min); 
        Index3 max_box = point_to_index(bbox_max); 

        //_stored_bounding_boxes[key] = box;

        //vector<int> local_grid_boxes= grid_boxes_intersecting_bounding_boxes[pi];
        for(int i = min_box.x()-1; i <= max_box.x(); i++){
            for(int j = min_box.y()-1; j <= max_box.y(); j++){
                for(int k = min_box.z()-1; k <= max_box.z(); k++){
                    // patch pi intersects box i,j,k
                    int index = hash(i,j,k); 
                    grid_boxes_intersecting_bounding_boxes[pi].push_back(index);

                }
            }
        }

    }

    // build inverted index from grid box ids to bounding box ids that intersect
    // the ith grid box
    unordered_map<int, vector<int>, index_hash, index_equal> bounding_boxes_in_grid_box;
    //set<int> grid_box_id_set;
    cout << "num bounding_boxes: " << num_bounding_boxes << endl;
    int counter = 0;
    // dump all grid boxes intersecting bounding_boxes into the hash map
    // with an empty value
    for(int pi= 0; pi < num_bounding_boxes; pi++){
        for(const auto& grid_box_id : grid_boxes_intersecting_bounding_boxes[pi]){
            // search for the grid box in the map
            auto it = bounding_boxes_in_grid_box.find(grid_box_id);

            // if we don't find it, give it a default value of an empty bbox
            // list
            if(it == bounding_boxes_in_grid_box.end()){
                bounding_boxes_in_grid_box[grid_box_id] = vector<int>();
                counter++;
            } 
        }
        //cout << "bounding_boxes_in_grid_box size: "<< counter << endl;
    }
    cout << "num grid boxes: " << bounding_boxes_in_grid_box.size() << endl;
    counter = 0;
    for(int pi= 0; pi < num_bounding_boxes; pi++){
        for(const auto& grid_box_id : grid_boxes_intersecting_bounding_boxes[pi]){
            bounding_boxes_in_grid_box[grid_box_id].push_back(pi);
            counter++;
        }
        //cout << "bounding_boxes_in_grid_box size: "<< counter << endl;
    }
    
    cout << "copying keys" << endl;
    vector<int> grid_box_ids;
    grid_box_ids.reserve(bounding_boxes_in_grid_box.size());
    for(const auto& kv: bounding_boxes_in_grid_box){
        grid_box_ids.push_back(kv.first);
    }

    //vector<Box*> grid_boxes(grid_box_ids.size());
    cout << "parallel GridMaps" << grid_box_ids.size() << endl;
    vector<GridMap> grid_box_maps(grid_box_ids.size());
    cout << "openmp" << endl;
    /*
    // build inverted index from grid box ids to bounding box ids that intersect
    // the ith grid box
    unordered_map<int, vector<int> > bounding_boxes_in_grid_box;
    set<int> grid_box_id_set;
    cout << "num bounding_boxes: " << num_bounding_boxes << endl;
    int counter = 0;
    for(int pi= 0; pi < num_bounding_boxes; pi++){
        for(const auto& grid_box_id : grid_boxes_intersecting_bounding_boxes[pi]){
            bounding_boxes_in_grid_box[grid_box_id].push_back(pi);
            grid_box_id_set.insert(grid_box_id);
            counter++;
        }
        cout << "bounding_boxes_in_grid_box size: "<< counter << endl;
    }
    vector<int> grid_box_ids(grid_box_id_set.begin(), grid_box_id_set.end());
    vector<Box*> grid_boxes(grid_box_ids.size());
    cout << "num grid boxes: " << grid_box_ids.size() << endl;
    vector<GridMap> grid_box_maps(grid_box_ids.size());
    */
    // initialize all distinct grid boxes we've found
#pragma omp parallel for
    for(int bi=0; bi < grid_box_ids.size(); bi++){
        int grid_box_id = grid_box_ids[bi];
    //for(const auto& grid_box_id : grid_box_ids){
        Box* b=nullptr;// = grid_box_maps[bi][grid_box_id];
        //cout << b << endl;
        if(!b){
            int i = grid_box_id;
            int ai = i % _num_boxes;
            i /= _num_boxes;
            int aj = i % _num_boxes;
            i /= _num_boxes;
            int ak= i;
            initialize_grid_box(ai,aj,ak,b);
        //cout << b << endl;
        }
        for(const auto& pi : bounding_boxes_in_grid_box[grid_box_id])
            b->add_bounding_box(pi);
        grid_box_maps[bi][grid_box_id] = b;
        //bounding_boxes_in_grid_box.erase(grid_box_id);
    }
    for(int bi=0; bi < grid_box_ids.size(); bi++){
        int box_id = grid_box_ids[bi];
        GridMap grid_map = grid_box_maps[bi];
        /*cout << "box: " << bi << endl;
        for(const auto& kv : grid_map){
        cout << kv.first << ": " << kv.second << endl;
        }*/
        _grid_box_map.insert(grid_map.begin(), grid_map.end());
        //_grid_box_map[box_id] = grid_box_map[box_id];
    }

}

void SpatialGrid::insert_face_map_patches(PatchSurfFaceMap* face_map){
    insert_face_map_patches(face_map->patches());
}
void SpatialGrid::insert_face_map_patches(vector<Patch*> patches){
    vector<FaceMapSubPatch*> face_map_subpatches;
    face_map_subpatches.reserve(patches.size());
    int i=0;
    cout << "insertin patches" << endl;
    for(const auto& patch : patches){
        assert(patch);
        auto p = dynamic_cast<FaceMapSubPatch*>(patch);
        cout << i++ << ": " << p << endl;
        face_map_subpatches.push_back(p);
    }
    insert_face_map_patches(face_map_subpatches);
}
void SpatialGrid::insert_face_map_patches(vector<FaceMapSubPatch*> patches){
    
    // compute average patch size
    double average_patch_size = 0.;
    
    int num_patches = patches.size();
    for(int pi= 0; pi < num_patches; pi++){
        auto patch = patches[pi];
        average_patch_size += patch->characteristic_length();
    }
    
    average_patch_size /= double(num_patches);
    
    //_box_size = average_patch_size;
    // discretize based on patch size
    int num_boxes = int(ceil((1. - (-1.))/average_patch_size));
    if(num_boxes % 2 == 1){
        num_boxes += 1;
        //_box_size = (1. - (-1.))/double(num_boxes+1);
    }
    //cout << "constructor init:" << endl;
    //cout << num_boxes << endl;
    _box_size = (1. - (-1.))/double(num_boxes);

    discretize(_box_size);
    
    vector<BoundingBox*> bounding_boxes(num_patches);
#pragma omp parallel for
    for(int pi= 0; pi < num_patches; pi++){
        auto patch = dynamic_cast<FaceMapSubPatch*>(patches[pi]);
        Point3 min, max;
        patch->bounding_box(min, max);
        auto box = new BoundingBox(min, max);
        bounding_boxes[pi] = box;
        //_stored_bounding_boxes[pi] = box;
    }
    for(int pi= 0; pi < num_patches; pi++)
        _stored_bounding_boxes[pi] = bounding_boxes[pi];
    insert_bounding_boxes(bounding_boxes);
/*
    // insert bounding boxes for each patch
    for(int pi= 0; pi < face_map->patches().size(); pi++){
        FaceMapSubPatch* patch = (FaceMapSubPatch*) face_map->patches()[pi];
        Point3 min, max;
        patch->bounding_box(min, max);
        BoundingBox* patch_bounding_box = new BoundingBox(min, max);
        insert(pi, patch_bounding_box);
    }*/
}


SpatialGrid::SpatialGrid(double box_size, Point3 spatial_min, Point3 spatial_max):
    _box_size(box_size),  _max(spatial_max),_min(spatial_min), _bounding_boxes_near_points(nullptr)
{
    discretize(box_size);

}
SpatialGrid::~SpatialGrid(){
cout << "destroying spatial grid" << endl;
cout << "num boxes: " << _grid_boxes.size() << " , " << _grid_box_map.size() << endl;
    for(auto box : _grid_boxes){
        cout << "deleting Box [" << box << "]" << endl;
        delete box;
    }
    for(auto index_box_pair: _grid_box_map){
        delete index_box_pair.second;
    }
    if(_bounding_boxes_near_points != nullptr){
        delete _bounding_boxes_near_points;
    }
}

Box*& SpatialGrid::SpatialGrid::operator()(int i, int j, int k){
    assert(i < _num_boxes && j < _num_boxes && k < _num_boxes);
    assert(i >= 0 && j >= 0 && k >= 0);
    //cout << i << ", " << j << ", " << k << endl;
    //cout << hash(i,j,k) << endl;
    Box* box = _grid_box_map[hash(i,j,k)];
    if(!box){
        initialize_grid_box(i,j,k,box);
        _grid_box_map[hash(i,j,k)] = box;
    }
    return box;
    //return _grid_boxes[hash(i, j, k)];
}

Box*& SpatialGrid::SpatialGrid::operator()(Point3 point){
    return (*this)(
            coord_to_index(point.x()),
            coord_to_index(point.y()),
            coord_to_index(point.z())
            );
}


Index3 SpatialGrid::SpatialGrid::point_to_index(Point3 point){
    return Index3(
                coord_to_index(point.x()),
                coord_to_index(point.y()),
                coord_to_index(point.z())
            );
}

int  SpatialGrid::SpatialGrid::coord_to_index(double x){
    //return floor((x + 1.)/_box_size);
    //// MJM WARNING CHANGED FLOOR TO CEILING
    //AND REFINEMENT BECAME SYMMETRIC
    //
    //return ceil((x + 1.)/_box_size);
    return floor((x + 1.)/_box_size);
}
double  SpatialGrid::SpatialGrid::index_to_coord(int i){
    return i*_box_size -1.;
}

void SpatialGrid::insert(int key, Point3 point){
    assert(point.x() <= 1.);
    assert(point.y() <= 1.);
    assert(point.z() <= 1.);
    assert(point.x() >= -1.);
    assert(point.y() >= -1.);
    assert(point.z() >= -1.);
    /*cout << point << endl;
            cout <<coord_to_index(point.x()) << "," << 
            coord_to_index(point.y())<< "," <<
            coord_to_index(point.z()) << "," << endl;*/
    Box* box_containing_point = (*this)(point);

    box_containing_point->add_point(key);
    _stored_points[key] = point;
}

void SpatialGrid::insert(int key, BoundingBox* box){
    // trim bounding box to fit in 
    double eps = 4e-16;
   Point3 bbox_min = box->min();
   Point3 bbox_max = box->max();
    
   if(bbox_min.x() <= -1.)
         bbox_min.x() = -1. + eps;
    if(bbox_min.y() <= -1.)
         bbox_min.y() = -1. + eps;
    if(bbox_min.z() <= -1.)
         bbox_min.z() = -1. + eps;

    if(bbox_max.x() >= 1.)
         bbox_max.x() = 1. - eps;
    if(bbox_max.y() >= 1.)
         bbox_max.y() = 1. - eps;
    if(bbox_max.z() >= 1.)
         bbox_max.z() = 1. - eps;

   Index3 min_box = point_to_index(bbox_min); 
   Index3 max_box = point_to_index(bbox_max); 

   _stored_bounding_boxes[key] = box;

   for(int i = min_box.x(); i <= max_box.x(); i++){
       for(int j = min_box.y(); j <= max_box.y(); j++){
           for(int k = min_box.z(); k <= max_box.z(); k++){
               Box* grid_box_intersecting_new_bounding_box = (*this)(i,j,k);
               grid_box_intersecting_new_bounding_box->add_bounding_box(key);

           }
       }
   }
  //dump_grid(); 
}

map<uint, vector<uint> > init_empty_map(uint num_items){
    map<uint, vector<uint> > near_map;
    
    for(uint i =0; i < num_items; i++){
        near_map[i] = vector<uint>();
    }

    return  near_map;
}

vector<uint> get_points_in_box(Box* box){
    return box->_points;
}

vector<uint> get_bounding_boxes_in_box(Box* box){
    return box->_bounding_boxes;
}

// returns a map from object B ids to a list of object A ids that share grid boxes with that object B id
map<uint, vector<uint> > temp(
        vector<Box*> _grid_boxes,  // list of grid boxes
        vector<uint> (*get_object_a_from_box)(Box* box), // extracts lists of type A objects from a grid box
        vector<uint> (*get_object_b_from_box)(Box* box), // extracts lists of type B objects from a grid box 
        map<uint, vector<uint> >& box_to_near_points){   // an empty input map with (# of type B objects contained in the grid ) empty elements

    // For each grid box 
    for(uint i = 0; i < _grid_boxes.size(); i++){
        Box* box = _grid_boxes[i];

        // Get the list of type A and type B objects contained in the ith grid box
        vector<uint> object_a_in_box = get_object_a_from_box(box);
        vector<uint> object_b_in_box = get_object_b_from_box(box);

        // Each object A in ith grid box [obja_1, obja_2, ... obja_n] 
        // shares a grid box with each object B  [objb_1, objb_2, ... objb_k]
        // in the ith grid box
        // We need to add the pairs (obja_i: [objb_1, ..., objb_k]) to the map
        // for i= 1,..., n
        
        // for i=1...n
        for(uint k = 0; k < object_a_in_box.size(); k++){
            //obja_i
            uint obj_a_id= object_a_in_box[k];
            // for j=1, ... k
            for(uint j =0; j < object_b_in_box.size(); j++){
                // objb_j
                uint obj_b_id= object_b_in_box[j];
                // append objb_j to the list of near points to obja_i
                box_to_near_points[obj_a_id].push_back(obj_b_id);

            }
        }
    }
    return box_to_near_points;

}
vector<Box*> SpatialGrid::collect_grid_boxes(GridMap grid_box_map){
    vector<Box*> grid_boxes;
    grid_boxes.reserve(grid_box_map.size());
    for(const auto& grid_box : grid_box_map){
        grid_boxes.push_back(grid_box.second);

    }
    return grid_boxes;

}

map<uint, vector<uint> >* SpatialGrid::boxes_near_points(bool update){
    // for every grid box
    //if(_bounding_boxes_near_points.empty()){
    if(_bounding_boxes_near_points == nullptr){
        map<uint, vector<uint> >* box_to_near_points = new map<uint, vector<uint> >();
        temp(collect_grid_boxes(_grid_box_map), 
                &get_points_in_box,
                &get_bounding_boxes_in_box,
                *box_to_near_points);
         _bounding_boxes_near_points = box_to_near_points ;
    } else if(update){
        delete _bounding_boxes_near_points;
        _bounding_boxes_near_points  = new map<uint, vector<uint> >();
        temp(collect_grid_boxes(_grid_box_map), 
                &get_points_in_box,
                &get_bounding_boxes_in_box,
                *_bounding_boxes_near_points );
         
    }
    return _bounding_boxes_near_points;
    //return box_to_near_points;

}

map<uint, vector<uint> > SpatialGrid::points_near_boxes(){
    map<uint, vector<uint> > point_to_near_boxes;// = init_empty_map(boxes.size());
        temp(collect_grid_boxes(_grid_box_map), 
            &get_bounding_boxes_in_box,
            &get_points_in_box,
            point_to_near_boxes);
    
    return point_to_near_boxes;
}

vector<uint> SpatialGrid::bounding_boxes_level_l_from_target_point(uint target_index, int l){
    Box* box_containing_target_point = (*this)(_stored_points[target_index]);
    vector<Box*> neighbors_l_away_from_target = 
        this->neighbors(box_containing_target_point, l);

    //cout << "level: " << l << endl;
    //cout << "current level: " << neighbors_l_away_from_target.size() << endl;
    
    vector<Box*> neighbors_l_plus_one_away_from_target = 
        this->neighbors(box_containing_target_point, l+1);
    //cout << "next level: " << neighbors_l_plus_one_away_from_target.size() << endl;

    neighbors_l_away_from_target.insert(
            neighbors_l_away_from_target.end(),
            neighbors_l_plus_one_away_from_target.begin(),
            neighbors_l_plus_one_away_from_target.end());

    vector<uint> bounding_boxes_intersecting_neighbors;

    for(vector<Box*>::iterator it = neighbors_l_away_from_target.begin();
            it != neighbors_l_away_from_target.end();
            it++){

        Box* neighbor = *it;
        bounding_boxes_intersecting_neighbors.insert(
                bounding_boxes_intersecting_neighbors.end(),
                neighbor->_bounding_boxes.begin(),
                neighbor->_bounding_boxes.end());

    }
    sort(bounding_boxes_intersecting_neighbors.begin(), bounding_boxes_intersecting_neighbors.end());
    bounding_boxes_intersecting_neighbors.erase(
            unique(bounding_boxes_intersecting_neighbors.begin(), 
                   bounding_boxes_intersecting_neighbors.end()),
            bounding_boxes_intersecting_neighbors.end());
    return bounding_boxes_intersecting_neighbors;
}

vector<uint> SpatialGrid::bounding_boxes_closest_to_point(uint target_index){

    vector<Box*> neighbors; 

    // check the box that the target point lives in for intersecting bounding
    // boxes
    vector<uint> nearby_bounding_boxes = (*boxes_near_points())[target_index];

    // If there's no bounding boxes in the target point's box...
    if(nearby_bounding_boxes.empty()){
        //cout << "searching for more boxes :( " << endl;
        Box* box_containing_target_point = (*this)(_stored_points[target_index]);
        for(size_t l = 1; l < _num_boxes; l++){
            // get neighbors distance l from the box containing the target point
            
            vector<uint> bounding_boxes_intersecting_neighbors = 
                bounding_boxes_level_l_from_target_point(target_index, l);
            /*
            vector<Box*> neighbors_l_away_from_target = 
                this->neighbors(box_containing_target_point, l);
            vector<Box*> neighbors_l_plus_one_away_from_target = 
                this->neighbors(box_containing_target_point, l+1);
            neighbors_l_away_from_target.insert(
                    neighbors_l_away_from_target.end(),
                    neighbors_l_plus_one_away_from_target.begin(),
                    neighbors_l_plus_one_away_from_target.end());
            vector<Box*> neighbors_l_plus_two_away_from_target = 
                this->neighbors(box_containing_target_point, l+2);
            neighbors_l_away_from_target.insert(
                    neighbors_l_away_from_target.begin(),
                    neighbors_l_plus_two_away_from_target.begin(),
                    neighbors_l_plus_two_away_from_target.end());

            // collect all the bounding boxes intersecting the neighbors
            vector<uint> bounding_boxes_intersecting_neighbors;
            for(vector<Box*>::iterator it = neighbors_l_away_from_target.begin();
                    it != neighbors_l_away_from_target.end();
                    it++){
                Box* neighbor = *it;
                bounding_boxes_intersecting_neighbors.insert(
                        bounding_boxes_intersecting_neighbors.end(),
                        neighbor->_bounding_boxes.begin(),
                        neighbor->_bounding_boxes.end());

            }
                    */

            // if there are any bounding boxes, return them, 
            if(!bounding_boxes_intersecting_neighbors.empty()) {
                nearby_bounding_boxes = bounding_boxes_intersecting_neighbors;
                break;
            }
            // otherwise, keep searching on the next level
        }
        // remove duplicates?
        /*
        sort(nearby_bounding_boxes.begin(), nearby_bounding_boxes.end());
        nearby_bounding_boxes.erase(
                unique(nearby_bounding_boxes.begin(), 
                       nearby_bounding_boxes.end()),
               nearby_bounding_boxes.end());*/
     
    }/* else {
        Box* box_containing_target_point = (*this)(_stored_points[target_index]);
        vector<Box*> neighbors_one_away_from_target = 
            this->neighbors(box_containing_target_point, 1);
            
        vector<uint> bounding_boxes_intersecting_neighbors;
        for(vector<Box*>::iterator it = neighbors_one_away_from_target.begin();
                it != neighbors_one_away_from_target.end();
                it++){
            Box* neighbor = *it;
            bounding_boxes_intersecting_neighbors.insert(
                    bounding_boxes_intersecting_neighbors.end(),
                    neighbor->_bounding_boxes.begin(),
                    neighbor->_bounding_boxes.end());

        }
        nearby_bounding_boxes.insert(
                nearby_bounding_boxes.end(),
                bounding_boxes_intersecting_neighbors.begin(),
                bounding_boxes_intersecting_neighbors.end());

    }*/

    // sort and remove duplicates
    sort( nearby_bounding_boxes.begin(), nearby_bounding_boxes.end() );
    nearby_bounding_boxes.erase( 
            unique( 
                nearby_bounding_boxes.begin(), 
                nearby_bounding_boxes.end() ), 
            nearby_bounding_boxes.end() );

    return nearby_bounding_boxes;

}

void SpatialGrid::enumerate_indices_within_bounds(int i_min, int i_max,
                int j_min, int j_max, 
                int k_min, int k_max, 
                vector<Index3>& neighbor_indices){
    for(int i = i_min; i <= i_max; i++){
        for(int j = j_min; j <= j_max; j++){
            for(int k = k_min; k <= k_max; k++){
                neighbor_indices.push_back(Index3(i,j,k));        
            }
        }
    }

}

vector<Index3> SpatialGrid::generate_neighbor_indices(Index3 index, size_t level){

    /**
     * To enumerate neighboring indices from "index" that live on the cube of
     * indices of
     * distance "level," it's enough to just them list off. Let d = 2*level+1,
     * or the length of the edge of the cube containing the neighbors. In the 
     * case of an interior (index - level > 0) box, there are 
     * d^2 + 4*(d-1)^2  + (d-2)^2 boxes on the cube (pick an arbitrary starting
     * face and include all four edges, this face has d^2 boxes. The four faces
     * that share an edge with this face each have (d-1)^2 boxes, when one edge
     * is excluded from each (since that edge is included in another face)
     * [Equivalently, there are two faces here of size d*(d-1) boxes and two more
     * of (d-1)*(d-2) boxes]. the remaining face only contributes (d-2)^2 
     * uncounted boxes (note that all edges have been already counted) 
     * These are the boxes we need to list.
     * 
     * We will enumerate all boxes on the cube of width 2*level+1, centered at
     * index and assume index is an "interior" point, in the sense that 
     * 0 <= index \pm level < num_boxes.
     * We will do a final pass over the boxes and remove indices that don't
     * correspond to grid boxes.
     */
    vector<Index3> neighbor_indices;
    int xi = index.x(); 
    int yi = index.y(); 
    int zi = index.z(); 
    int i_fixed, j_fixed, k_fixed; 
    int d = 2*level + 1; // neighbor cube edge width

    // Choose the plane x = xi - level as the starting face
    i_fixed = xi - level; 
    enumerate_indices_within_bounds(
            i_fixed, i_fixed,       // stay in the plane x = xi - level
            yi - level, yi + level, // iterate over all y, 2*level+1 steps
            zi - level, zi + level, // iterate over all z, 2*level+1 steps
            neighbor_indices);             // total of (2*level+1)^2 new indices added

  
    assert(neighbor_indices.size() == d*d);
    // For the next four loops, we skip the boxes in the plane x = xi -level
    // enumerate plane z = zi - level
    k_fixed = zi - level; 
    enumerate_indices_within_bounds(
            xi-level+1, xi+level,       // iterate over x, 2*level steps
            yi - level+1, yi + level, // iterate over y, 2*level steps
            k_fixed, k_fixed,       // stay in the plane z = zi - level
            neighbor_indices);             // total of (2*level)^2 new indices added

    assert(neighbor_indices.size() == d*d+ (d-1)*(d-1));

    // enumerate plane y = yi + level 
    j_fixed = yi + level; 
    enumerate_indices_within_bounds(
            xi-level+1, xi+level,       // iterate over x, 2*level steps
            j_fixed, j_fixed,       // stay in the plane y = yi + level
            zi - level+1, zi + level, // iterate over  z, 2*level steps
            neighbor_indices);             // total of (2*level)^2 new indices added
    assert(neighbor_indices.size() == d*d+ 2*(d-1)*(d-1));


    // enumerate plane z = zi + level
    k_fixed = zi + level; 
    enumerate_indices_within_bounds(
            xi-level+1, xi+level,       // iterate over x, 2*level steps
            yi - level, yi + level - 1, // iterate over y, 2*level steps
            k_fixed, k_fixed,       // stay in the plane z = zi + level
            neighbor_indices);             // total of (2*level)^2 new indices added
    assert(neighbor_indices.size() == d*d+ 3*(d-1)*(d-1));

    // enumerate plane y = yi - level 
    j_fixed = yi - level; 
    enumerate_indices_within_bounds(
            xi-level+1, xi+level,       // iterate over x, 2*level steps
            j_fixed, j_fixed,       // stay in the plane y = yi - level
            zi - level, zi + level-1, // iterate over z, 2*level steps
            neighbor_indices);             // total of (2*level)^2 new indices added
    assert(neighbor_indices.size() == d*d+ 4*(d-1)*(d-1));

    // enumerate plane y = yi - level 
    i_fixed = xi + level; 
    enumerate_indices_within_bounds(
            i_fixed, i_fixed,                   // stay in the plane x = xi + level
            yi - level + 1, yi + level - 1,     // iterate over y, 2*level-1 steps
            zi - level + 1, zi + level-1,       // iterate over z, 2*level-1 steps
            neighbor_indices);                         // total of (2*level-1)^2 new indices added
    assert(neighbor_indices.size() == d*d+ 4*(d-1)*(d-1) + (d-2)*(d-2));
    return neighbor_indices;
}



vector<Box*> SpatialGrid::neighbors(Box* box, size_t level){
    //assert(level > 0);
    if(level == 0){
        vector<Box*> self;
        self.push_back(box);
        return  self;
    }else if(level == 1){
        return neighbors(box);
    } else {
        vector<Box*> neighbors; 
        Index3 index = box->index();

        vector<Index3> neighbor_indices = generate_neighbor_indices(index, level);
        for(size_t i = 0; i < neighbor_indices.size(); i++){
            Index3 neighbor_index = neighbor_indices[i];

            if(this->is_valid(neighbor_index)){
                Box* neighbor = (*this)(neighbor_index.x(), neighbor_index.y(), neighbor_index.z());
                neighbors.push_back(neighbor);
            }
        }
        return neighbors;
    }

}

vector<Box*> SpatialGrid::neighbors(Box* box){
    Index3 index = box->index();

    vector<Box*> neighbors; 
    // enumerate boxes that share face/edge/corner with box
    for(int i =-1; i <=1; i++){
        for(int j =-1; j <=1; j++){
            for(int k =-1; k <=1; k++){
                if(i == 0 && j == 0 && k == 0){
                    // this is the same index as box; skip it
                    continue;
                }

                // compute neighbor index    
                int ni = index.x() + i;
                int nj = index.y() + j;
                int nk = index.z() + k;
                // If the neighbor is inside the grid bounds, keep it.
                if(( ni >= 0 &&  ni< _num_boxes) && 
                        ( nj >= 0 &&  nj< _num_boxes) && 
                        ( nk >= 0 &&  nk< _num_boxes)){
                    
                    Box* neighbor_box = (*this)(ni, nj, nk);
                    neighbors.push_back(neighbor_box);
                }
            }
        }
    }
    int num_neighbors = neighbors.size();

    // make sure we found the right number of neighbors:
    assert( num_neighbors == 7  || // corner box
            num_neighbors == 11 || // box along an edge
            num_neighbors == 17 || // box along a face
            num_neighbors == 26); // interior box  

    return neighbors;
}

void SpatialGrid::dump_grid(){
    for(uint gi = 0; gi < _grid_boxes.size(); gi++){
        Box* grid_box = _grid_boxes[gi];
        cout << "grid box " << gi << " (min, max): (" <<  grid_box->min() <<"," << grid_box->max() << ")" <<endl;
        cout << "contained bounding boxes: " << endl;
        for(uint bi = 0; bi < grid_box->_bounding_boxes.size(); bi++){
           cout << grid_box->_bounding_boxes[bi] << ", "; 
        }
        cout << endl << "contained points: " << endl;
        for(uint pi = 0; pi < grid_box->_points.size(); pi++){
           cout << grid_box->_points[pi] << ": " << _stored_points[grid_box->_points[pi]] << endl; 
        }
        cout << endl << endl;
    }
}

END_EBI_NAMESPACE
