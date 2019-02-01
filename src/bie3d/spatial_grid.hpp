#ifndef __SPATIAL_GRID_HPP__
#define __SPATIAL_GRID_HPP__

#include "bdry3d/face_map_subpatch.hpp"
#include "common/ebiobject.hpp"
#include "bdry3d/patch_surf.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "bdry3d/patch_samples.hpp"
#include <unordered_map>
#include <tuple>
BEGIN_EBI_NAMESPACE
namespace Markgrid {
    // code yanked from 
    // https://stackoverflow.com/questions/11408934/using-a-stdtuple-as-key-for-stdunordered-map
    
    typedef tuple<int, int, int> Index;
    struct index_hash : public std::unary_function<int, std::size_t> {
        std::size_t operator()(const int& k) const {
            return k;
            //return std::get<0>(k) ^ std::get<1>(k) ^ std::get<2>(k);
            //return i*_num_boxes*_num_boxes +j*_num_boxes + k;
        }
    };

    struct index_equal : public std::binary_function<int, int, bool> {
        bool operator()(const int& v0, const int& v1) const {
            return v0 == v1;
        }
    };

    class Box{
        protected:
        int _i;
        int _j; 
        int _k;
        Point3 _min;
        Point3 _max;
        bool _init;

        public:
        vector<uint> _points;  //id of 3d points contained in this box
        vector<uint> _bounding_boxes; //id of bounding boxes contained in this box
        Box():
            _min(Point3(0.)), _max(Point3(0.)), _init(false){;}

        Box(Point3 min, Point3 max):
            _min(min), _max(max), _init(false){;}

        Box(int i, int j, int k):
            _i(i), _j(j), _k(k), _init(false) {;}

        Box(int i, int j, int k, Point3 min, Point3 max):
            _i(i), _j(j), _k(k), _min(min), _max(max), _init(true){;}
        Box(Index index, Point3 min, Point3 max):
            _i(std::get<0>(index)), 
            _j(std::get<1>(index)), 
            _k(std::get<2>(index)), _min(min), _max(max), _init(true){;}
        
        Point3 min(){
            return _min;
        }
        Point3 max(){
            return _max;
        }
        bool init() { return _init;}
        void update(int i, int j, int k, Point3 min, Point3 max){
            _i = i; 
            _j = j;
            _k = k;
            _min = min;
            _max = max;
            _init = true;
        }
        void add_point(uint key){
            _points.push_back(key);
        }
        void add_bounding_box(uint key){
            _bounding_boxes.push_back(key);
        }
        Index3 index(){
            return Index3(_i, _j, _k);
        }
        Box& operator=(const Box& b){
            _i = b._i;
            _j = b._j;
            _k = b._k;
            _min = b._min;
            _max = b._max;
            _init = b._init;
            _points = b._points;
            _bounding_boxes = b._bounding_boxes;
            return *this;
        }
        Box(Box const & b):
            _i(b._i),
            _j(b._j),
            _k(b._k),
            _min(b._min),
            _max(b._max),
            _init(b._init),
            _points(b._points),
            _bounding_boxes(b._bounding_boxes){ 
                cout << "CALLED COPY CONSTRCUTOR" << endl;
        }

    };
    
    class GridBox: public Box {
        GridBox(): Box(){;}

        GridBox(Point3 min, Point3 max): Box(min, max)
            {;}

        GridBox(int i, int j, int k): Box(i,j,k) {;}

        GridBox(int i, int j, int k, Point3 min, Point3 max):
            Box(i,j,k, min, max)
            {;}

    };

    class BoundingBox: public Box{
        public:
        BoundingBox(): Box(){;}

        BoundingBox(Point3 min, Point3 max): Box(min, max)
            {;}

        BoundingBox(int i, int j, int k): Box(i,j,k) {;}

        BoundingBox(int i, int j, int k, Point3 min, Point3 max):
            Box(i,j,k, min, max)
            {;}

    };

    class SpatialGrid {
        typedef unique_ptr<Box> BoxPtr;
        // Consider defining custom iterators over corner, exterior, and
        // interior grid boxes?
    typedef unordered_map<int, Box*, index_hash, index_equal> GridMap;
        private:
            int _num_boxes;
            double _box_size;
            map<uint, Point3> _stored_points;
            vector<Box*> _grid_boxes;
            GridMap _grid_box_map;
            Point3 _max;
            Point3 _min;
            map<uint, vector<uint> >* _bounding_boxes_near_points;
            map<uint, vector<uint> > _points_near_bounding_boxes;
        public:
            map<uint, BoundingBox*> _stored_bounding_boxes;
            SpatialGrid(PatchSurfFaceMap* face_map);
            SpatialGrid(vector<Patch*> patches);
            SpatialGrid(DblNumMat points);
        
            SpatialGrid(PatchSurfFaceMap* face_map, DblNumMat points, int num_boxes=12);
            SpatialGrid(vector<Patch*> patches, DblNumMat points, int num_boxes=12);
        SpatialGrid(double box_size, Point3 spatial_min, Point3 spatial_max);
        ~SpatialGrid();

        void initialize_grid_box(int i, int j, int k, Box*& box);
        
        void insert_points(DblNumMat points);
        void insert_face_map_patches(PatchSurfFaceMap* face_map);
        void insert_face_map_patches(vector<Patch*> patches);
        void insert_face_map_patches(vector<FaceMapSubPatch*> patches);
        void insert_bounding_boxes(vector<BoundingBox*> bounding_boxes);
        //Box*& operator()(int i, int j, int k);
        //Box*& operator()(Point3 point);
        Box*& operator()(int i, int j, int k);
        Box*& operator()(Point3 point);
        Index3 point_to_index(Point3 point);
        vector<Box*> collect_grid_boxes(GridMap grid_box_map);

        void discretize(double box_size);

        int hash(int i, int j, int k){
            return i*_num_boxes*_num_boxes +j*_num_boxes + k;
        }
        static int hash(int i, int j, int k, int _num_boxes){
            return i*_num_boxes*_num_boxes +j*_num_boxes + k;
        }
        int coord_to_index(double x);
        double index_to_coord(int i);
        void insert(int key, Point3 point);
        void insert(int key, BoundingBox* box);
        map<uint, vector<uint> >* boxes_near_points(bool update=false);
        map<uint, vector<uint> > points_near_boxes();
        // Return a list of GridBoxes in the volume grid that 
        //vector<GridBox*> intersect(QueryBox query_box);
        void dump_grid();
        const Point3 min(){
            return _min;
        }
        const Point3 max(){
            return _max;
        }
        double box_size() {
            return _box_size;
        }
        int num_boxes(){
            return _num_boxes;
        }
        const vector<Box*> grid_boxes(){
            return _grid_boxes;
        }
        Box* grid_boxes(uint i);
        vector<Box*> neighbors(Box* box);
        vector<Box*> neighbors(Box* box, size_t level);

        vector<uint> bounding_boxes_closest_to_point(uint target_index);
        vector<uint> bounding_boxes_level_l_from_target_point(uint target_index, int l);

        vector<Index3> generate_neighbor_indices(Index3 index, size_t level);
        void enumerate_indices_within_bounds(
                int i_min, int i_max,
                int j_min, int j_max, 
                int k_min, int k_max, 
                vector<Index3>& neighbor_indices);
        bool is_valid(Index3 index){
            int i = index.x();
            int j = index.y();
            int k = index.z();
            return ( i >= 0 &&  i< _num_boxes) && 
                ( j >= 0 &&  j< _num_boxes) && 
                ( k >= 0 &&  k< _num_boxes);
        
        }
        

    };
}
END_EBI_NAMESPACE
#endif 
