#include "bdry3d/geogram_interface.hpp"
#include "common/vtk_writer.hpp"
#include "bdry3d/face_map_subpatch.hpp"
#include <geogram/delaunay/delaunay_3d.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
BEGIN_EBI_NAMESPACE
using GEO::index_t;
using GEO::signed_index_t;
void AABBTree::initialize_mesh(DblNumMat vertices, IntNumMat faces){
    //_mesh = std::move(unique_ptr<GEO::Mesh>(new GEO::Mesh(3,false)));
    _mesh = new GEO::Mesh(3,false);
    int num_vertices = vertices.n();
    int num_faces = faces.n();

    assert(faces.m() == 3); // triangles only
    _mesh->vertices.create_vertices(num_vertices);
    _mesh->facets.create_triangles(num_faces);
    
    // copy vertices
    for(int i =0; i < num_vertices; i++){
        GEO::vec3& p = _mesh->vertices.point(i);
        for (int d = 0; d < DIM; d++) {
           p[d] = vertices(d,i); 
        }
    }
    // copy faces
    const int num_verts_per_tri = 3;
    for(int i =0; i < num_faces; i++){
        for (int f = 0; f < num_verts_per_tri; f++) {
           _mesh->facets.set_vertex(i, f, faces(f,i));
        }

    }
}

AABBTree::AABBTree(PatchSamples* samples){

}

void AABBTree::mesh_patches_and_build_tree(vector<Patch*> patches){
    // here, we're meshing patch bounding boxes and passing them to geogram
    //int num_patches = surface->num_patches();
    int num_patches = patches.size();
    //const int verts_per_bbox = 8;
    //const int triangles_per_bbox = 12;
    //double spacing = 0.090909;
    double spacing = .25;
    int num_samples_1d = int(round(1./spacing))+1;
    const int verts_per_bbox = num_samples_1d*num_samples_1d;
    const int triangles_per_bbox = 2*(num_samples_1d-1)*(num_samples_1d-1);
    //const int verts_per_bbox = 25;
    //const int triangles_per_bbox = 32;
    int num_total_vertices = num_patches*verts_per_bbox;
    int num_total_triangles = num_patches*triangles_per_bbox; // 

    // global surface mesh
    _surface_vertices.resize(DIM, num_total_vertices);
    _surface_triangles.resize(DIM, num_total_triangles);
    
    // mapping from patch id -> list of triangle ids associated with that patch
    vector<vector<uint> > patch_to_triangle_id_map(num_patches, vector<uint>());
#pragma omp parallel for
    for (int pi = 0; pi < num_patches; pi++) {
        auto patch = FaceMapSubPatch::as_subpatch(patches[pi]);
        
        // mesh the bounding box
        DblNumMat vertices;
        IntNumMat triangles;
        //patch->mesh_bounding_box(vertices, triangles);
        patch->mesh_patch(spacing, Rectangle(Interval(0.,1.), Interval(0.,1.)), vertices, triangles);

        assert(vertices.m() == DIM);
        assert(vertices.n() == verts_per_bbox);
        assert(triangles.m() == 3);
        assert(triangles.n() == triangles_per_bbox);

        // copy vertices and triangle into global mesh matrices
        for (int i = 0; i < verts_per_bbox; i++) {
            for (int d = 0; d < DIM; d++) {
                _surface_vertices(d, pi*verts_per_bbox+ i) = 
                    vertices(d,i);
            }
        }
        
        vector<uint> triangles_on_patch(triangles_per_bbox);
        for (int i = 0; i < triangles_per_bbox; i++) {
            for (int d = 0; d < 3; d++) {
                _surface_triangles(d, pi*triangles_per_bbox+ i) = 
                    pi*verts_per_bbox + triangles(d,i);
            }
            // append triangle ids to pi's indexing
            triangles_on_patch[i] = pi*triangles_per_bbox+i;

        }
        // save the mapping
        patch_to_triangle_id_map[pi] = triangles_on_patch;
        // sanity check
        setvalue(vertices, 0.);
        setvalue(triangles, 0);
    }
    // compute inverted indexing of triangle ids -> patch id
    _triangle_ids_to_patch_ids.resize(num_total_triangles);
    for(int pi =0; pi < num_patches; pi++){
        const auto& triangles_on_patch = patch_to_triangle_id_map[pi];
        for(const int triangle_id : triangles_on_patch){
            _triangle_ids_to_patch_ids[triangle_id] = pi;
        }
    }
    // dump triangle mesh for visualization
    //vector<int> pids;
    //pids.assign(_triangle_ids_to_patch_ids.begin(), _triangle_ids_to_patch_ids.end());
    //write_triangle_mesh_to_vtk( _surface_vertices, _surface_triangles, 0, "data/surface", pids);

    initialize_mesh(_surface_vertices, _surface_triangles);
    //_tree = std::move(unique_ptr<GEO::MeshFacetsAABB>(new GEO::MeshFacetsAABB(*(_mesh.get()))));
    _tree = new GEO::MeshFacetsAABB(*_mesh,false);

}
void AABBTree::init(){
    //GEO::CmdLine::set_arg("algo:delaunay","default");
    //GEO::CmdLine::import_arg_group("standard");
    //GEO::   CmdLine::import_arg_group("algo");

        /*GEO::CmdLine::declare_arg(
            "convex_hull", false,
            "compute just the convex hull of the points"
        );

        GEO::CmdLine::declare_arg(
            "dimension", 3, "3 for 3D, 2 for 2D"
        );
        GEO::Process::enable_multithreading(true);
*/
        GEO::Process::set_max_threads(omp_get_max_threads());
    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("pre");
    GEO::CmdLine::import_arg_group("remesh");
    GEO::CmdLine::import_arg_group("algo");
    GEO::CmdLine::import_arg_group("post");
    GEO::CmdLine::import_arg_group("opt");
    GEO::CmdLine::import_arg_group("co3ne");
    GEO::CmdLine::import_arg_group("tet");
    GEO::CmdLine::import_arg_group("poly");
    GEO::CmdLine::declare_arg(
            "nl:MKL", true,
            "Use Intel Math Kernel Library (if available in the system)"
            );
    /*GEO::CmdLine::declare_arg(
         "algo:parallel", true,
         "Use parallel standard algorithms"
     );*/
        cout << "num geogram threads: " << GEO::Process::maximum_concurrent_threads() << endl; 
}


void AABBTree::mesh_near_zone_and_build_tree(vector<Patch*> patches){
    // here we're forming a single triangle that spans the bounding box of 
    // each patch. Then, we can ask for an infintesimal bounding box of a target
    // point to compute the bounding boxes quickly
    int num_patches = patches.size();
    const int verts_per_patch= 3;
    const int triangles_per_patch = 1;
    int num_total_vertices = num_patches*verts_per_patch;
    int num_total_triangles = num_patches*triangles_per_patch; // 

    // global surface mesh
    _surface_vertices.resize(DIM, num_total_vertices);
    _surface_triangles.resize(DIM, num_total_triangles);
    
    //vector<vector<uint> > patch_to_triangle_id_map(num_patches, vector<uint>());
    _triangle_ids_to_patch_ids.resize(num_total_triangles);
#pragma omp parallel for
    for (int pi = 0; pi < num_patches; pi++) {
        auto patch = FaceMapSubPatch::as_subpatch(patches[pi]);
        
        DblNumMat vertices;
        IntNumMat triangles;
        patch->single_bounding_box_triangle(vertices, triangles);

        assert(vertices.m() == DIM);
        assert(vertices.n() == verts_per_patch);
        assert(triangles.m() == 3);
        assert(triangles.n() == triangles_per_patch);

        // copy vertices and triangle into global matrices
        for (int i = 0; i < verts_per_patch; i++) {
            for (int d = 0; d < DIM; d++) {
                _surface_vertices(d, pi*verts_per_patch+ i) = 
                    vertices(d,i);
            }
        }
        
        //vector<uint> triangles_on_patch(triangles_per_patch);
        for (int i = 0; i < triangles_per_patch; i++) {
            for (int d = 0; d < 3; d++) {
                _surface_triangles(d, pi*triangles_per_patch+ i) = 
                    pi*verts_per_patch + triangles(d,i);
            }
            //triangles_on_patch[i] = pi*triangles_per_patch+i;
        }
        //patch_to_triangle_id_map[pi] = triangles_on_patch;
        setvalue(vertices, 0.);
        setvalue(triangles, 0);
        _triangle_ids_to_patch_ids[pi] = pi;
    }
    // compute inverted indexing
    /*_triangle_ids_to_patch_ids.resize(num_total_triangles);
    for(int pi =0; pi < num_patches; pi++){
        const auto& triangles_on_patch = patch_to_triangle_id_map[pi];
        for(const int triangle_id : triangles_on_patch){
            _triangle_ids_to_patch_ids[triangle_id] = pi;
        }
    }*/
    //vector<int> pids;
    //pids.assign(_triangle_ids_to_patch_ids.begin(), _triangle_ids_to_patch_ids.end());
    //write_triangle_mesh_to_vtk( _surface_vertices, _surface_triangles, 0, "data/", pids);

    initialize_mesh(_surface_vertices, _surface_triangles);
    //_tree = std::move(unique_ptr<GEO::MeshFacetsAABB>(new GEO::MeshFacetsAABB(*(_mesh.get()))));
    _tree = new GEO::MeshFacetsAABB(*_mesh,false);

}

AABBTree::AABBTree(vector<Patch*> patches, bool near_zone){

    init();
    if(near_zone){
        cout << "making near-zone triangles." << endl;
        mesh_near_zone_and_build_tree(patches);
    } else {
        mesh_patches_and_build_tree(patches);
    }
}

AABBTree::AABBTree(PatchSurfFaceMap* surface, bool near_zone){
    init();
    if(near_zone){
        mesh_near_zone_and_build_tree(surface->patches());
    } else {
        mesh_patches_and_build_tree(surface->patches());
    }
}
AABBTree::AABBTree(PatchSurf* surface){
    int num_patches = surface->num_patches();
    double spacing = .2;
    auto p = surface->patch(0);
    int num_total_vertices = num_patches*p->mesh_num_vertices(spacing);
    int num_total_triangles = num_patches*p->mesh_num_triangles(spacing);

    // global surface mesh
    _surface_vertices.resize(DIM, num_total_vertices);
    _surface_triangles.resize(DIM, num_total_triangles);
    
    vector<vector<uint> > patch_to_triangle_id_map(surface->num_patches(), vector<uint>());
#pragma omp parallel for
    for (int pi = 0; pi < surface->num_patches(); pi++) {
        const auto& patch = surface->patch(pi);
        
        DblNumMat vertices;
        IntNumMat triangles;
        
        patch->mesh_patch(
                spacing, Rectangle(Interval(0.,1.), Interval(0.,1.)),
                vertices, triangles);
        int num_vertices_per_patch_1d = patch->mesh_num_vertices(spacing);
        int num_vertices_per_patch = patch->mesh_num_vertices(spacing);
        int num_triangles_per_patch = patch->mesh_num_triangles(spacing);

        assert(vertices.m() == DIM);
        assert(vertices.n() == num_vertices_per_patch);
        assert(triangles.m() == 3);
        assert(triangles.n() == num_triangles_per_patch);

        // copy vertices and triangle into global matrices
        for (int i = 0; i < num_vertices_per_patch; i++) {
            for (int d = 0; d < DIM; d++) {
                _surface_vertices(d, pi*num_vertices_per_patch + i) = 
                    vertices(d,i);
            }
        }
        
        vector<uint> triangles_on_patch(num_triangles_per_patch);
        for (int i = 0; i < num_triangles_per_patch; i++) {
            for (int d = 0; d < 3; d++) {
                _surface_triangles(d, pi*num_triangles_per_patch+ i) = 
                    pi*num_vertices_per_patch+triangles(d,i);
            }
            triangles_on_patch[i] = pi*num_triangles_per_patch +i;
        }
        patch_to_triangle_id_map[pi] = triangles_on_patch;
        setvalue(vertices, 0.);
        setvalue(triangles, 0);
    }
    // compute inverted indexing
    _triangle_ids_to_patch_ids.resize(num_total_triangles);
    for(int pi =0; pi < surface->num_patches(); pi++){
        const auto& triangles_on_patch = patch_to_triangle_id_map[pi];
        for(const int triangle_id : triangles_on_patch){
            _triangle_ids_to_patch_ids[triangle_id] = pi;
        }
    }
    //write_triangle_mesh_to_vtk( _surface_vertices, _surface_triangles, 0, "data/");

    initialize_mesh(_surface_vertices, _surface_triangles);
    //_tree = std::move(unique_ptr<GEO::MeshFacetsAABB>(new GEO::MeshFacetsAABB(*(_mesh.get()))));
    _tree = new GEO::MeshFacetsAABB(*_mesh);
}


AABBTree::AABBTree(DblNumMat vertices, IntNumMat faces){
    initialize_mesh(vertices, faces);
    _tree = new GEO::MeshFacetsAABB(*_mesh);
    //_tree = std::move(unique_ptr<GEO::MeshFacetsAABB>(new GEO::MeshFacetsAABB(*(_mesh.get()))));
}

uint AABBTree::closest_patch_to_point(Point3 query_point){
    GEO::vec3 geogram_query_point(query_point.x(), query_point.y(), query_point.z());
    double distance;
    GEO::vec3 closest_point;
    GEO::index_t closest_triangle_id = _tree->nearest_facet(geogram_query_point, closest_point, distance);

    return _triangle_ids_to_patch_ids[closest_triangle_id];


}
struct TrianglesNearBBox{
    public:
    vector<GEO::index_t> triangle_ids;
    TrianglesNearBBox() {}

    void operator()(GEO::index_t triangle_id){
        triangle_ids.push_back(triangle_id);
    }
};
vector<uint> AABBTree::patches_intersecting_bbox(Point3 query_box_min,Point3 query_box_max){
    GEO::Box bounding_box;
    memcpy(bounding_box.xyz_min, query_box_min.array(), sizeof(double)*3);
    memcpy(bounding_box.xyz_max, query_box_max.array(), sizeof(double)*3);
    
    /*DblNumMat vertices;
    IntNumMat triangles;
    vertices.resize(3,8);
    for (int i = 0; i < 2; i++) { // 
        for (int j = 0; j < 2; j++) { // y
            for (int k = 0; k < 2; k++) { // z
                int index = 4*i+2*j+k;
                vertices(0,index) = i ? query_box_min(0) : query_box_max(0);
                vertices(1,index) = j ? query_box_min(1) : query_box_max(1);
                vertices(2,index) = k ? query_box_min(2) : query_box_max(2);
            }
        }
    }

    // sorry future me
    triangles.resize(3,12);

    // explicitly list indices that make up quad faces of bbox as enumerated
    // above
    vector<vector<int> > bbox_quads(6);
    bbox_quads[0]= {0,2,3,1};
    bbox_quads[1]= {1,3,7,5};
    bbox_quads[2]= {4,0,1,5};
    bbox_quads[3]= {4,6,2,0};
    bbox_quads[4]= {5,7,6,4};
    bbox_quads[5]= {2,6,7,3};
    for (int i = 0; i < 6; i++) {
        const auto& quad = bbox_quads[i];

        // index set for the two triangles that make up each quad
        vector<int> t1 = {0,1,2};
        vector<int> t2 = {0,2,3};
        for (int d = 0; d < 3; d++) {
            // copy them into the output array
            triangles(d,2*i) = quad[t1[d]];
            triangles(d,2*i+1) = quad[t2[d]];
        }
    }
    cout << vertices << endl;
    cout << triangles<< endl;
    write_triangle_mesh_to_vtk(
            vertices, triangles, it, "data/query_box");
    */

    TrianglesNearBBox triangle_aggregator;
    _tree->compute_bbox_facet_bbox_intersections(bounding_box, triangle_aggregator);
    set<uint> patch_ids_set;
    for(const auto& index : triangle_aggregator.triangle_ids){
        patch_ids_set.insert(_triangle_ids_to_patch_ids[index]);
    }
    vector<uint> unique_patch_ids(patch_ids_set.begin(), patch_ids_set.end());
    /*
    cout << "patches near bbox " << it << ": ";
    for(const auto pid: unique_patch_ids){
        cout << pid << ", ";
        
    }
    cout << endl;
    cout << "triangles near bbox " << it << ": ";
    for(const auto tid: triangle_aggregator.triangle_ids){
        cout << tid << ", ";
        
    }
    cout << endl;
    it++;*/
    return unique_patch_ids;

}

vector<uint> AABBTree::patches_near_point(Point3 query_point){
    GEO::Box bounding_box;
    Point3 bbox_min = query_point - Point3(1e-12);
    Point3 bbox_max = query_point + Point3(1e-12);
    memcpy(bounding_box.xyz_min, bbox_min.array(), sizeof(double)*3);
    memcpy(bounding_box.xyz_max, bbox_max.array(), sizeof(double)*3);

    TrianglesNearBBox triangle_aggregator;

    _tree->compute_bbox_facet_bbox_intersections(bounding_box, triangle_aggregator);

    set<uint> patch_ids_set;
    for(const auto& index : triangle_aggregator.triangle_ids){
        patch_ids_set.insert(_triangle_ids_to_patch_ids[index]);
    }
    /*{
    DblNumMat vertices;
    IntNumMat triangles;
    vertices.resize(3,8);
    for (int i = 0; i < 2; i++) { // 
        for (int j = 0; j < 2; j++) { // y
            for (int k = 0; k < 2; k++) { // z
                int index = 4*i+2*j+k;
                vertices(0,index) = i ? bbox_min(0) : bbox_max(0);
                vertices(1,index) = j ? bbox_min(1) : bbox_max(1);
                vertices(2,index) = k ? bbox_min(2) : bbox_max(2);
            }
        }
    }

    // sorry future me
    triangles.resize(3,12);

    // explicitly list indices that make up quad faces of bbox as enumerated
    // above
    vector<vector<int> > bbox_quads(6);
    bbox_quads[0]= {0,2,3,1};
    bbox_quads[1]= {1,3,7,5};
    bbox_quads[2]= {4,0,1,5};
    bbox_quads[3]= {4,6,2,0};
    bbox_quads[4]= {5,7,6,4};
    bbox_quads[5]= {2,6,7,3};
    for (int i = 0; i < 6; i++) {
        const auto& quad = bbox_quads[i];

        // index set for the two triangles that make up each quad
        vector<int> t1 = {0,1,2};
        vector<int> t2 = {0,2,3};
        for (int d = 0; d < 3; d++) {
            // copy them into the output array
            triangles(d,2*i) = quad[t1[d]];
            triangles(d,2*i+1) = quad[t2[d]];
        }
    }
    write_triangle_mesh_to_vtk(
            vertices, triangles, 50, "data/query_box");

    }*/

    vector<uint> unique_patch_ids(patch_ids_set.begin(), patch_ids_set.end());
    return unique_patch_ids;
}



/*AABBTree::~AABBTree(){
}*/

//void convex_hull(FaceMapSubPatch* patch, GEO::Mesh& mesh){
 void convex_hull(FaceMapSubPatch* patch, DblNumMat& vertices, IntNumMat& faces){
     GEO::Mesh mesh(3);
    DblNumMat control_points = patch->control_points();
    //unique_ptr<GEO::Delaunay3d> delaunay(new GEO::Delaunay3d(3));
    GEO::Delaunay_var delaunay = GEO::Delaunay::create(3);
	    delaunay->set_keeps_infinite(true);
    delaunay->set_vertices(control_points.n(), control_points.data());
   /**
    * \brief Saves a Delaunay triangulation to a file.
    * \param[in] delaunay a pointer to the Delaunay triangulation.
    * \param[in] filename the name of the file to be saved.
    *  -If the example was compiled with the Geogram library, then any
    *  mesh file handled by Geogram can be used.
    *  if the example was compiled with Delaunay_psm (single file), then
    *  the points and vertices of the triangulation are output in ASCII.
    * \param[in] convex_hull_only if true, then only the triangles on the
    *  convex hull are output.
    void save_Delaunay(
	Delaunay* delaunay, const std::string& filename,
    ) 
    */
	bool convex_hull_only = true;
    GEO::vector<index_t> tri2v;
	
	if(convex_hull_only) {
	    
	    // The convex hull can be efficiently traversed only if infinite
	    // tetrahedra are kept.
	    geo_assert(delaunay->keeps_infinite());
	    
	    // The convex hull can be retrieved as the finite facets
	    // of the infinite cells (note: it would be also possible to
	    // throw away the infinite cells and get the convex hull as
	    // the facets adjacent to no cell). Here we use the infinite
	    // cells to show an example with them.
	    
	    
	    // This block is just a sanity check
	    {
		for(index_t t=0; t < delaunay->nb_finite_cells(); ++t) {
		    geo_debug_assert(delaunay->cell_is_finite(t));
		}
		
		for(index_t t=delaunay->nb_finite_cells();
		    t < delaunay->nb_cells(); ++t) {
		    geo_debug_assert(delaunay->cell_is_infinite(t));
		}
	    }
	    
	    // This iterates on the infinite cells
	    for(
		index_t t = delaunay->nb_finite_cells();
		t < delaunay->nb_cells(); ++t
	     ) {
		for(index_t lv=0; lv<4; ++lv) {
		    signed_index_t v = delaunay->cell_vertex(t,lv);
		    if(v != -1) {
			tri2v.push_back(index_t(v));
		    }
		}
	    }
	}
	
	// Using Geogram mesh I/O: copy Delaunay into a Geogram
	// mesh and save it to disk.
	
    GEO::vector<double> pts(delaunay->nb_vertices() * 3);
	for(index_t v = 0; v < delaunay->nb_vertices(); ++v) {
	    pts[3 * v] = delaunay->vertex_ptr(v)[0];
	    pts[3 * v + 1] = delaunay->vertex_ptr(v)[1];
	    pts[3 * v + 2] =
		(delaunay->dimension() >= 3) ? delaunay->vertex_ptr(v)[2] : 0.0;
	}
	
	if(convex_hull_only) {
	    mesh.facets.assign_triangle_mesh(3, pts, tri2v, true);
	} else if(delaunay->dimension() == 3) {
        GEO::vector<index_t> tet2v(delaunay->nb_cells() * 4);
	    for(index_t t = 0; t < delaunay->nb_cells(); ++t) {
	    tet2v[4 * t] = index_t(delaunay->cell_vertex(t, 0));
	    tet2v[4 * t + 1] = index_t(delaunay->cell_vertex(t, 1));
	    tet2v[4 * t + 2] = index_t(delaunay->cell_vertex(t, 2));
	    tet2v[4 * t + 3] = index_t(delaunay->cell_vertex(t, 3));
	    }
	    mesh.cells.assign_tet_mesh(3, pts, tet2v, true);
	} else if(delaunay->dimension() == 2) {
	    tri2v.resize(delaunay->nb_cells() * 3);
	    for(index_t t = 0; t < delaunay->nb_cells(); ++t) {
		tri2v[3 * t] = index_t(delaunay->cell_vertex(t, 0));
		tri2v[3 * t + 1] = index_t(delaunay->cell_vertex(t, 1));
		tri2v[3 * t + 2] = index_t(delaunay->cell_vertex(t, 2));
	    }
	    mesh.facets.assign_triangle_mesh(3, pts, tri2v, true);
	}
	mesh.show_stats();

    geogram_mesh_to_arrays(mesh, vertices, faces);

}

void geogram_mesh_to_arrays(GEO::Mesh& mesh, DblNumMat& vertices, IntNumMat& faces){
    int num_vertices = mesh.vertices.nb();
    int num_faces = mesh.facets.nb();
    vertices.resize(DIM, num_vertices);
    for (int i = 0; i < num_vertices; ++i) {
        GEO::vec3 p = mesh.vertices.point(i);
        for (int d = 0; d < DIM; d++) {
           vertices(d,i) = p[d]; 
        }
    }
    assert(mesh.facets.are_simplices());
    faces.resize(3, num_faces);
    for (int c = 0; c < num_faces; ++c) {
        for (int lv = 0; lv < 3; ++lv) {
            faces(lv,c) = mesh.facets.vertex(c, lv);
        }
    }
}


END_EBI_NAMESPACE
