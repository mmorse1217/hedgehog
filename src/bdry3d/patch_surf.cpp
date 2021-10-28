#include "patch_surf.hpp"
#include "patch_surf_face_map.hpp"
#include <bdsurf.hpp>
BEGIN_EBI_NAMESPACE

PatchSurf::~PatchSurf()
{

  for(uint i=0; i<_patches.size(); i++)
      delete _patches[i];

  _patches.resize(0);
}

// TODO cache this
NumVec<Point2> Patch::sample_patch(SamplingPattern sampling_pattern){
    double spacing = Options::get_double_from_petsc_opts("-bis3d_spacing");
    int num_samples = floor(1./spacing)+1;

    NumVec<Point2> on_surface_patch_samples(num_samples*num_samples);

    for(int si = 0; si < num_samples; si++){
        for(int sj = 0; sj < num_samples; sj++){
            int index = si*num_samples + sj;
            
            if(sampling_pattern == CHEBYSHEV){
                on_surface_patch_samples(index) = 
                    Point2(
                            (cos(si*M_PI/double(num_samples-1))+1.)/2., 
                            (cos(sj*M_PI/double(num_samples-1))+1.)/2.
                          );
            } else if(sampling_pattern == EQUISPACED){
                on_surface_patch_samples(index) = Point2(si*spacing, sj*spacing);
            } else {
                assert(0);
            }
        }
    }
    return on_surface_patch_samples;

}
DblNumMat Patch::generate_qbkix_samples(NumVec<Point2> on_surface_patch_samples){
    
    double boundary_distance_ratio = Options::get_double_from_petsc_opts("-boundary_distance_ratio");
    double interpolation_spacing_ratio = Options::get_double_from_petsc_opts("-interpolation_spacing_ratio");
    double qbkix_order = Options::get_double_from_petsc_opts("-near_interpolation_num_samples");
    double spacing = Options::get_double_from_petsc_opts("-bis3d_spacing");
    int num_total_samples = on_surface_patch_samples.m();

    DblNumMat qbkix_points(DIM, num_total_samples*qbkix_order);
    for(int si = 0; si < num_total_samples; si++){
        Point2 uv = on_surface_patch_samples(si);
        Point3 position_and_derivs[3];
        // TODO could introduce a bug here by calling virtual class methods
        // (blendsurf vs polynomials)...
        xy_to_patch_coords(uv.array(), EVAL_VALUE|EVAL_1ST_DERIV, (double*)position_and_derivs);
        Point3 position = position_and_derivs[0];
        Point3 normal = cross(position_and_derivs[1], position_and_derivs[2]).dir();

        for(int qi = 0; qi < qbkix_order; qi++){
            
            Point3 qbkix_point = position 
                - boundary_distance_ratio*spacing*characteristic_length()*normal 
                - normal * double(qi)*interpolation_spacing_ratio*spacing*characteristic_length();

            for(int d = 0; d < DIM; d++)
                qbkix_points(d, si*qbkix_order + qi) = qbkix_point(d);

        }
    }
    return qbkix_points;
}


double Patch::characteristic_length(DblNumMat normal, DblNumVec quad_weight){

    double surface_area = compute_surface_area(normal, quad_weight);
    _characteristic_length = sqrt(surface_area);
    return _characteristic_length;
}

double Patch::compute_surface_area(DblNumMat normals, DblNumVec quadrature_weight){

    // make sure we have n x n quadrature weights and n x n normal vectors at
    // each point on the patch
    assert(quadrature_weight.m() == normals.n());
    int num_samples = int(floor(sqrt(normals.n())));
    double surface_integral = 0.;
    for(int i =0; i < num_samples; i++){
        for(int j =0; j < num_samples; j++){
            int index = i*num_samples + j;
            Point3 n(normals.clmdata(index));
            double w_ij = quadrature_weight(index);
            surface_integral += n.length()*w_ij;
        }
    }
    return surface_integral;
}


void Patch::mesh_patch(double spacing,  // sampling spacing
        Rectangle sample_domain,
        DblNumMat& vertices, //3 x num_vertices matrix of positions
        IntNumMat& faces //3 x num_faces matrix of positions
        ){
    int n = mesh_num_vertices_1d(spacing);
    int num_vertices = mesh_num_vertices(spacing);
    // sample target domain
    DblNumMat sample_uv_coords(2, num_vertices, true,
            Sampling::sample_2d<Sampling::equispaced>(n, sample_domain).data()
            );
    
    // resize output array
    vertices.resize(DIM, num_vertices);

    // evaluate vertex positions
    for (int i = 0; i < num_vertices; i++) {
        Point3 vertex_position;
        Point2 vertex_uv(sample_uv_coords.clmdata(i));
        xy_to_patch_coords(vertex_uv.array(), EVAL_VALUE, vertex_position.array());
        for (int d = 0; d < DIM; d++) {
            vertices(d,i) = vertex_position(d);    
        }

    }
    mesh_patch_no_sampling(spacing, faces);
}

void Patch::mesh_patch_no_sampling(  
        double spacing, //num vertices in 1d along a regular sampling
        IntNumMat& faces //3 x num_faces matrix of positions
        ){
    const int num_verts_per_tri = 3;
    int num_triangles = mesh_num_triangles(spacing);
    faces.resize(num_verts_per_tri, num_triangles);
    int n = mesh_num_vertices_1d(spacing);
    // compute triangles. assume they are of the shape:
    // 6---7---8
    // | / | / |
    // 3---4---5
    // | / | / |
    // 0---1---2
    int findex = 0;
    for(int i = 0; i < n-1; i++){
        for (int j = 0; j < n-1; j++) {
            int vindex= (n)*i + j;
            
            // for triangle
            // 2---3
            // | / |
            // 0---1
            // lower triangle = (0,1,3)
            // upper triangle = (0,3,2)
            // for general triangle
            //
            // k+n---k+n+1
            // |   /    |
            // k  ---  k+1
            // lower triangle = (k,k+1,k+n+1)
            // upper triangle = (k,k+n+1,k+n)
            
            // lower triangle
            faces(0,findex)   = vindex;
            faces(1,findex)   = vindex+1;
            faces(2,findex)   = vindex+n+1;
            findex++;

            // upper triangle
            faces(0,findex) = vindex;
            faces(1,findex) = vindex+n+1;
            faces(2,findex) = vindex+n;
            findex++;
            
        }
    }

}

END_EBI_NAMESPACE
