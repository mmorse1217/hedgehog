#include "vtk_writer.hpp"
#include <sampling.hpp>

#include "vec3t.hpp"
#include <lean_vtk.hpp>
#include "bdry3d/patch_samples.hpp"
#include "bdry3d/face_map_subpatch.hpp"


/**
 * Use the accepted answer below to actually find the libs to link in... vtk is
 * silly
 * https://stackoverflow.com/questions/29448204/which-libraries-should-be-linked-to-this-vtk-example
 *
 */
BEGIN_EBI_NAMESPACE
enum Filetype {
    POINTS = 0,
    PATCHES = 1,
    LINES =2,
    BOXES =3,
    TRIANGLES =4
};
string build_filename(int iteration, Filetype filetype, string prefix){
    std::ostringstream ss;
    ss << iteration;

    //string file_name = Options::get_string_from_petsc_opts("-bd3d_meshfile");
    //file_name.erase(file_name.begin(), file_name.begin()+10); // chop off "wrl_files/" from file name
    //file_name.erase(file_name.end()-4, file_name.end()); // chop off ".wrl" from file name
    string file_name = Test::get_domain();
    cout << file_name << endl;
    if(prefix.empty()){
    file_name = "data/" + file_name;
    } else {
    file_name = prefix +"_"+ file_name;
    }
    if(filetype == POINTS){
        file_name +=  "_qbkix_points_";
    } else if(filetype == PATCHES){
        file_name +=  "_ref_patches_";
    } else if(filetype == LINES){
        file_name +=  "_qbkix_dirs_";
    } else if(filetype == BOXES){
        file_name +=  "_qbkix_bboxes_";
    } else if(filetype == TRIANGLES){
        file_name +=  "_mesh_";
    }
    file_name += "it";
    file_name +=  ss.str();
    file_name += filetype == POINTS ? ".vtp" : ".vtu";
    return file_name;
}

void mesh_and_save_points(DblNumMat points, DblNumMat values, string name){
    int root_n = int(sqrt(points.n()));


    IntNumMat faces;
    const int num_verts_per_tri = 3;
    int n = root_n;
    int num_triangles = 2*(n-1)*(n-1);
    faces.resize(num_verts_per_tri, num_triangles);
    cout << "num_triangles "<<  num_triangles << endl;
    cout << "n "<<  n << endl;
    
    // compute triangles. assume they are of the shape:
    // 6---7---8
    // | / | / |
    // 3---4---5
    // | / | / |
    // 0---1---2
    int findex = 0;
    for(int i = 0; i < n-1; i++){
        for (int j = 0; j < n-1; j++) {
            //int stride = 2*(n-1);
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

    write_triangle_mesh_to_vtk(points, faces, 0, name, vector<int>(), values);
}


void write_structured_data_to_vtk(vector<Vec> data_vecs, 
        int n, 
        vector<int> degrees_of_freedom, 
        vector<string> vec_names,
        string filename){
    int num_vecs = data_vecs.size();
    assert(num_vecs == degrees_of_freedom.size());
    assert(num_vecs == vec_names.size());
    for(auto vec: data_vecs){
        
    }
}

void write_general_points_to_vtk(Vec point_positions, int degrees_of_freedom,
                                 string filename, Vec point_values,
                                 string file_prefix) {

  leanvtk::VTUWriter writer;
  vector<double> points;
  vector<double> values;

  // int num_values = Petsc::get_vec_size(values)/degrees_of_freedom;
  // int num_points= Petsc::get_vec_size(points)/DIM;
  int64_t num_values;
  int64_t num_points;
  VecGetLocalSize(point_values, &num_values);
  VecGetLocalSize(point_positions, &num_points);
  num_values /= degrees_of_freedom;
  num_points /= DIM;
  assert(num_points == num_values);
  DblNumMat points_local = get_local_vector(DIM, num_points, point_positions);
  DblNumMat values_local =
      get_local_vector(degrees_of_freedom, num_values, point_values);

  for (int i = 0; i < num_points; i++) {
    for (int d = 0; d < DIM; d++) {
      points.push_back(points_local(d, i));
    }
    for (int d = 0; d < degrees_of_freedom; d++) {
      values.push_back(values_local(d, i));
    }
  }

  string file_name = build_filename(0, POINTS, filename + "lean_vtk");
  if (degrees_of_freedom == 1) {
    writer.add_scalar_field("Values", values);
  } else {
    writer.add_vector_field("Values", values, degrees_of_freedom);
  }
  writer.write_point_cloud(file_name, DIM, points);
}

void write_qbkix_points_to_vtk(
    DblNumMat qbkix_points,
    NumVec<OnSurfacePoint> final_closest_on_surface_points, int iteration,
    string file_prefix) {

  leanvtk::VTUWriter writer;
  vector<double> points;
  vector<double> qbkix_point_patch_ids;
  vector<double> relative_patch_ids;

  for (int i = 0; i < qbkix_points.n(); i++) {
    for (int d = 0; d < DIM; d++) {
      points.push_back(qbkix_points(d, i));
    }
    qbkix_point_patch_ids.push_back(
        final_closest_on_surface_points(i).parent_patch);
    relative_patch_ids.push_back(i / 4);
  }

  string file_name = build_filename(iteration, POINTS, file_prefix);
  writer.add_scalar_field("Patch Id", qbkix_point_patch_ids);
  writer.add_scalar_field("Relative Patch Id", relative_patch_ids);
  writer.write_point_cloud(file_name, DIM, points);
}

void write_triangle_mesh_to_vtk(DblNumMat vertices, IntNumMat faces,
                                int iteration, string file_prefix,
                                vector<int> corresponding_patches,
                                DblNumMat point_values) {

  leanvtk::VTUWriter writer;
  vector<double> points;
  vector<double> patch_ids;
  vector<double> corr_patch_ids;
  vector<int> triangle_ids;
  vector<double> values;

  int num_vertices = vertices.n();
  int num_faces = faces.n();
  int values_stride = point_values.m();
  assert(corresponding_patches.size() == num_faces);
  assert(num_vertices == point_values.n());

  const int num_verts_per_tri = 3;

  points.reserve(DIM * num_vertices);
  patch_ids.reserve(num_faces);
  triangle_ids.reserve(num_verts_per_tri* num_faces);
  values.reserve(values_stride * num_faces);
  corr_patch_ids.reserve(num_faces);
  assert(vertices.m() == DIM);
  assert(faces.m() == num_verts_per_tri);

  // 3 x (4*num_patches) double array to contain patch corner locations

  for (int i = 0; i < num_vertices; i++) {
    for (int d = 0; d < DIM; d++) {
      points.push_back(vertices(d, i));
    }
  }
  if (point_values.n() == num_vertices) {
    for (int i = 0; i < num_vertices; i++) {
      for (int d = 0; d < values_stride; d++) {
        values.push_back(point_values(d, i));
      }
    }
  }
  for (int f = 0; f < num_faces; f++) {
    vector<long long> face_ids(num_verts_per_tri);

    for (int d = 0; d < num_verts_per_tri; d++) {
      face_ids[d] = faces(d, f);
      triangle_ids.push_back(faces(d, f));
    }
    patch_ids.push_back(f);
    corr_patch_ids.push_back(corresponding_patches[f]);
  }

  string file_name = build_filename(iteration, TRIANGLES, file_prefix);

  writer.add_cell_scalar_field("Triangle Id", patch_ids);
  writer.add_cell_scalar_field("Patch Id", corr_patch_ids);
  writer.add_field("Values", values, values_stride);
  writer.write_surface_mesh(file_name, DIM, num_verts_per_tri, points, triangle_ids);
}

void write_face_map_patches_to_vtk(DblNumMat qbkix_points,
                                   vector<int> patches_refined_relative_ids,
                                   PatchSurfFaceMap *face_map, int iteration,
                                   string file_prefix) {
  leanvtk::VTUWriter writer;
  vector<double> points;
  vector<int> quad_ids;
  vector<double> patch_ids;
  vector<double> patch_relative_ids;

  // uv coordinates to evaluate at the patch corners... (counter clockwise
  // order);
  vector<Point2> uv_coords;
  uv_coords.push_back(Point2(0., 0.));
  uv_coords.push_back(Point2(1., 0.));
  uv_coords.push_back(Point2(1., 1.));
  uv_coords.push_back(Point2(0., 1.));

  int num_patches = face_map->patches().size();
  const int num_verts_per_quad = 4;

  for (int pi = 0; pi < num_patches; pi++) {

    FaceMapSubPatch *patch = (FaceMapSubPatch *)face_map->patches()[pi];

    // make a list of point ids defining the patch
    vector<long long> point_ids(num_verts_per_quad);
    // evaluate the patch at its corners
    for (int i = 0; i < uv_coords.size(); i++) {
      int index = num_verts_per_quad * pi + i;
      point_ids[i] = index;
      Point3 position(0.);
      patch->xy_to_patch_coords(uv_coords[i].array(), PatchSamples::EVAL_VL,
                                position.array());

      // save the position
      for (int d = 0; d < DIM; d++) {
        points.push_back(position(d));
      }
    }
    // save the quad id list
    for (int d = 0; d < num_verts_per_quad; d++) {
      quad_ids.push_back(point_ids[d]);
    }
    patch_ids.push_back(pi);
    patch_relative_ids.push_back(patches_refined_relative_ids.at(pi));
  }

  string file_name = build_filename(iteration, PATCHES, file_prefix);
  writer.add_cell_scalar_field("Patch Id", patch_ids);
  writer.add_cell_scalar_field("Relative Patch Id", patch_ids);
  writer.write_surface_mesh(file_name, DIM, num_verts_per_quad, points, quad_ids);
}

void write_face_map_mesh_to_vtk(PatchSurfFaceMap *face_map, int iteration,
                                string file_prefix, int num_samples_1d) {

  leanvtk::VTUWriter writer;
  vector<double> points;
  vector<double> patch_ids;
  vector<int> triangle_ids;

  // int num_samples_1d = 20;
  double step = 1. / (num_samples_1d - 1);
  const int num_samples = num_samples_1d * num_samples_1d;
  const int triangles_per_bbox =
      2 * (num_samples_1d - 1) * (num_samples_1d - 1);
  const int num_verts_per_tri = 3;
  DblNumMat vertices(DIM, num_samples);
  Rectangle domain = Sampling::base_domain;
  IntNumMat faces(num_verts_per_tri,
                  triangles_per_bbox); // 3 x num_faces matrix of positions
  // list of vtkQuads to view

  // uv coordinates to evaluate at the patch corners... (counter clockwise
  // order);

  int num_patches = face_map->patches().size();

  // 3 x (4*num_patches) double array to contain patch corner locations

  int num_faces = 0;
  for (int pi = 0; pi < num_patches; pi++) {

    FaceMapSubPatch *patch = face_map->subpatch(pi);

    // make a list of point ids defining the patch

    patch->mesh_patch(step, domain, vertices, faces);
    // evaluate the patch at its corners
    for (int i = 0; i < num_samples; i++) {
      long long index = i + pi * num_samples;
      Point3 position(vertices.clmdata(i));
      // save the position
      for (int d = 0; d < DIM; d++) {
        points.push_back(position(d));
      }
    }

    // for each triangle in on a patch ...
    for (int i = 0; i < triangles_per_bbox; i++) {
      // for each vertex of  a triangle...
      for (int d = 0; d < num_verts_per_tri; d++) {
        int vertex_id = faces(d, i) + pi * num_samples;
        triangle_ids.push_back(vertex_id);
      }
      // save the quad id list
      patch_ids.push_back(pi);
    }
  }

  string file_name = build_filename(iteration, TRIANGLES, file_prefix);

  writer.add_cell_scalar_field("Patch Id", patch_ids);
  writer.write_surface_mesh(file_name, DIM, num_verts_per_tri, points,
                            triangle_ids);
}

void write_face_map_patch_bounding_boxes_to_vtk(
        vector<int> patches_refined_relative_ids,
        PatchSurfFaceMap* face_map, int iteration, string file_prefix,
        bool inflate){

  leanvtk::VTUWriter writer;
  vector<double> points;
  vector<int> quad_ids;
  vector<double> patch_ids;
  vector<double> relative_patch_ids;

  int num_patches = face_map->patches().size();
  const int num_vertices_per_bbox = 8;
  const int num_faces_per_bbox = 6;
  const int num_vertices_per_face = 4;
  points.resize(DIM * num_patches * num_vertices_per_bbox);
  quad_ids.reserve(num_vertices_per_face * num_faces_per_bbox * num_patches);
  patch_ids.reserve(num_vertices_per_face * num_faces_per_bbox * num_patches);
  relative_patch_ids.reserve(num_vertices_per_face * num_faces_per_bbox * num_patches);

  // 3 x (4*num_patches) double array to contain patch corner locations
  cout << "bbox ids: " << endl;
  for (int pi = 0; pi < num_patches; pi++) {

    auto patch = face_map->subpatch(pi);

    Point3 bbox_min;
    Point3 bbox_max;
    if (inflate) {
      patch->inflated_bounding_box(bbox_min, bbox_max);
    } else {
      patch->bounding_box(bbox_min, bbox_max);
    }

    int index = DIM * num_vertices_per_bbox * pi;
    vector<Point3> bounding_box_corners = {
        Point3(bbox_min),
        Point3(bbox_min(0), bbox_min(1), bbox_max(2)),
        Point3(bbox_min(0), bbox_max(1), bbox_min(2)),
        Point3(bbox_max(0), bbox_min(1), bbox_min(2)),
        Point3(bbox_min(0), bbox_max(1), bbox_max(2)),
        Point3(bbox_max(0), bbox_min(1), bbox_max(2)),
        Point3(bbox_max(0), bbox_max(1), bbox_min(2)),
        Point3(bbox_max)};

    for (int corner = 0; corner < 8; corner++) {
      Point3 bounding_box_corner = bounding_box_corners[corner];
      for (int d = 0; d < DIM; d++) {
        points[index + DIM * corner + d] = bounding_box_corner[d];
      }
    }

    vector<long long> ids(4, 0);
    ids[0] = 8 * pi;
    ids[1] = 8 * pi + 1;
    ids[2] = 8 * pi + 5;
    ids[3] = 8 * pi + 3;
    quad_ids.insert(quad_ids.end(), ids.begin(), ids.end());
    patch_ids.push_back(pi);
    relative_patch_ids.push_back(patches_refined_relative_ids.at(pi));

    ids[0] = 8 * pi;
    ids[1] = 8 * pi + 2;
    ids[2] = 8 * pi + 6;
    ids[3] = 8 * pi + 3;
    quad_ids.insert(quad_ids.end(), ids.begin(), ids.end());
    patch_ids.push_back(pi);
    relative_patch_ids.push_back(patches_refined_relative_ids.at(pi));

    ids[0] = 8 * pi;
    ids[1] = 8 * pi + 1;
    ids[2] = 8 * pi + 4;
    ids[3] = 8 * pi + 2;
    quad_ids.insert(quad_ids.end(), ids.begin(), ids.end());
    patch_ids.push_back(pi);
    relative_patch_ids.push_back(patches_refined_relative_ids.at(pi));

    ids[0] = 8 * pi + 3;
    ids[1] = 8 * pi + 6;
    ids[2] = 8 * pi + 7;
    ids[3] = 8 * pi + 5;
    quad_ids.insert(quad_ids.end(), ids.begin(), ids.end());
    patch_ids.push_back(pi);
    relative_patch_ids.push_back(patches_refined_relative_ids.at(pi));

    ids[0] = 8 * pi + 1;
    ids[1] = 8 * pi + 5;
    ids[2] = 8 * pi + 7;
    ids[3] = 8 * pi + 4;
    quad_ids.insert(quad_ids.end(), ids.begin(), ids.end());
    patch_ids.push_back(pi);
    relative_patch_ids.push_back(patches_refined_relative_ids.at(pi));

    ids[0] = 8 * pi + 2;
    ids[1] = 8 * pi + 4;
    ids[2] = 8 * pi + 7;
    ids[3] = 8 * pi + 6;
    quad_ids.insert(quad_ids.end(), ids.begin(), ids.end());
    patch_ids.push_back(pi);
    relative_patch_ids.push_back(patches_refined_relative_ids.at(pi));
    }

    writer.add_cell_scalar_field("Patch Id", patch_ids);
    writer.add_cell_scalar_field("Relative Patch Id", relative_patch_ids);
    string file_name = build_filename(iteration, BOXES, file_prefix);
    writer.write_surface_mesh(file_name, DIM, 4, points, quad_ids);
}

void write_lines_from_qbkix_to_closest_point(
    DblNumMat qbkix_points,
    NumVec<OnSurfacePoint> final_closest_on_surface_points,
    PatchSurfFaceMap *face_map, int iteration, string file_prefix) {

  leanvtk::VTUWriter writer;
  vector<double> points;
  vector<double> point_directions;

  int num_qbkix_points = qbkix_points.n();
  assert(qbkix_points.n() == final_closest_on_surface_points.m());
  cout << num_qbkix_points << endl;

  points.reserve(DIM * num_qbkix_points);
  point_directions.reserve(DIM * num_qbkix_points);

  int num_patches = face_map->patches().size();

  for (int qi = 0; qi < num_qbkix_points; qi++) {
    Point3 qbkix_point(qbkix_points.clmdata(qi));

    // get the patch containing the closest on-surface point to the qbkix
    // point
    OnSurfacePoint on_surface_point = final_closest_on_surface_points(qi);

    int patch_conataining_closest_point = on_surface_point.parent_patch;
    Point3 dir;
    if (patch_conataining_closest_point > -1) {
      FaceMapSubPatch *patch =
          (FaceMapSubPatch *)face_map->patches()[on_surface_point.parent_patch];

      // compute the position of closest on_surface_point
      Point3 on_surface_point_position(0.);
      patch->xy_to_patch_coords(on_surface_point.parametric_coordinates.array(),
                                PatchSamples::EVAL_VL,
                                on_surface_point_position.array());

      // compute direction to closest on_surface_point_position
      dir = on_surface_point_position - qbkix_point;
    } else {
      dir = Point3(0, 0, .25);
    }

    for (int d = 0; d < DIM; d++) {
      points.push_back(qbkix_point(d));
      point_directions.push_back(dir(d));
    }
  }

  string file_name = build_filename(iteration, LINES, file_prefix);

  writer.add_vector_field("qbkix direction", point_directions, DIM);
  writer.write_point_cloud(file_name, DIM, points);
}

END_EBI_NAMESPACE
