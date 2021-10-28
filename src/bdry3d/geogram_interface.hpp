#ifndef __GEOGRAM_INTERFACE_HPP__
#define __GEOGRAM_INTERFACE_HPP__

#include "common/ebiobject.hpp"
#include "bdry3d/patch_surf.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "bdry3d/patch_samples.hpp"
#include <unordered_map>
#include <tuple>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>


BEGIN_EBI_NAMESPACE
class FaceMapSubPatch;
 void convex_hull(FaceMapSubPatch* patch, DblNumMat& vertices, IntNumMat& faces);
 void geogram_mesh_to_arrays(GEO::Mesh& mesh, DblNumMat& vertices, IntNumMat& faces);
class AABBTree {
    private:
        GEO::MeshFacetsAABB* _tree;
        GEO::Mesh* _mesh;

        void initialize_mesh(DblNumMat vertices, IntNumMat faces);
        int it=0;

        void mesh_patches_and_build_tree(vector<Patch*> patches);
        void mesh_near_zone_and_build_tree(vector<Patch*> patches);
        void init();

    public:
        int _num_vertices;
        int _num_faces;
        DblNumMat _surface_vertices;
        IntNumMat _surface_triangles;
        vector<uint> _triangle_ids_to_patch_ids;

        AABBTree(PatchSamples* samples);
        AABBTree(PatchSurf* surface);
        AABBTree(PatchSurfFaceMap* surface, bool near_zone=false);
        AABBTree(vector<Patch*> patches, bool near_zone=false);
        AABBTree(DblNumMat vertices, IntNumMat faces);
        ~AABBTree(){
            delete _tree;
            delete _mesh;
        }
        uint closest_patch_to_point(Point3 query_point);

        vector<uint> patches_intersecting_bbox(Point3 query_box_min,Point3 query_box_max);
        vector<uint> patches_near_point(Point3 query_point);



};


END_EBI_NAMESPACE
#endif
