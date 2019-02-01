#include "vtk_writer.hpp"
#include <sampling.hpp>

#include <vtkVersion.h>
#include <vtkVersion.h>
#include <vtkObjectBase.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkAppendFilter.h>
#include <vtkPolyData.h>
#include <vtkQuad.h>
#include <vtkAppendFilter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include "vec3t.hpp"
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

void write_general_points_to_vtk(Vec points, int degrees_of_freedom, 
        string filename, Vec values, string file_prefix){

    //int num_values = Petsc::get_vec_size(values)/degrees_of_freedom;
    //int num_points= Petsc::get_vec_size(points)/DIM;
    int64_t num_values;
    int64_t num_points;
    VecGetLocalSize(values, &num_values);
    VecGetLocalSize(points, &num_points);
    num_values  /=degrees_of_freedom;
    num_points/=DIM;
    cout << num_points << ", " << num_values << endl;
    assert(num_points == num_values);
    cout << "getting local vector 1" << endl;
    DblNumMat points_local = get_local_vector(DIM, num_points, points);
    cout << "getting local vector 1" << endl;
    DblNumMat values_local = get_local_vector(degrees_of_freedom, num_values, values);


    vtkSmartPointer<vtkPoints> vtk_points = vtkSmartPointer<vtkPoints>::New();
    
    vtkSmartPointer<vtkDoubleArray> vtk_values= 
        vtkSmartPointer<vtkDoubleArray>::New();

    vtk_values->SetNumberOfComponents(degrees_of_freedom);
    vtk_values->SetName("Values");
    


    for(int i =0; i < num_points; i++){
        Point3 point(points_local.clmdata(i));
        vtk_points->InsertNextPoint(point.x(), point.y(), point.z());

        vtk_values->InsertNextTuple(values_local.clmdata(i));
        /*DblNumVec value(degrees_of_freedom, false, values_loca.clmdata(i));
        for(int d =0; d < degrees_of_freedom; d++){
            vtk_values->InsertNextValue(value(d));
        }*/
    }
    vtkSmartPointer<vtkPolyData> polydata = 
        vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(vtk_points);

    polydata->GetPointData()->AddArray(vtk_values);
    polydata->GetPointData()->SetActiveScalars("Values");
 
  // Write the file
  vtkSmartPointer<vtkXMLPolyDataWriter> writer =  
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  
  string f = file_prefix + filename;
  writer->SetFileName(f.c_str());
  writer->SetInputData(polydata);
 
  // Optional - set the mode. The default is binary.
  //writer->SetDataModeToBinary();
  writer->SetDataModeToAscii();
 
  writer->Write();

}


void write_qbkix_points_to_vtk(DblNumMat qbkix_points, 
        NumVec<OnSurfacePoint> final_closest_on_surface_points, int iteration, string file_prefix){


    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    
    
    vtkSmartPointer<vtkIntArray> qbkix_point_patch_ids= 
        vtkSmartPointer<vtkIntArray>::New();
    qbkix_point_patch_ids->SetNumberOfComponents(1);
    qbkix_point_patch_ids->SetName("Patch Id");
    
    vtkSmartPointer<vtkIntArray> patch_relative_ids= 
        vtkSmartPointer<vtkIntArray>::New();
    patch_relative_ids->SetNumberOfComponents(1);
    patch_relative_ids->SetName("Refined relative Patch Id");


    for(int i =0; i < qbkix_points.n(); i++){
        Point3 qbkix_point(qbkix_points.clmdata(i));
        points->InsertNextPoint(qbkix_point.x(), qbkix_point.y(), qbkix_point.z());
        qbkix_point_patch_ids->InsertNextValue(final_closest_on_surface_points(i).parent_patch);
        patch_relative_ids->InsertNextValue(i/4);
    }
 vtkSmartPointer<vtkPolyData> polydata = 
    vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);

    polydata->GetPointData()->AddArray(qbkix_point_patch_ids);
    polydata->GetPointData()->AddArray(patch_relative_ids);
    polydata->GetPointData()->SetActiveScalars("Patch Id");
 
  // Write the file
  vtkSmartPointer<vtkXMLPolyDataWriter> writer =  
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  
  string file_name = build_filename(iteration, POINTS, file_prefix);
  
  writer->SetFileName(file_name.c_str());
  writer->SetInputData(polydata);
 
  // Optional - set the mode. The default is binary.
  //writer->SetDataModeToBinary();
  writer->SetDataModeToAscii();
 
  writer->Write();

}

void write_triangle_mesh_to_vtk(
        DblNumMat vertices, IntNumMat faces, int iteration, string file_prefix,
        vector<int> corresponding_patches){


    // list of vtkQuads to view 
    vtkSmartPointer<vtkCellArray> cellArray =
        vtkSmartPointer<vtkCellArray>::New();

    // list of vtkPoints to defining triangles
    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();
    int num_vertices = vertices.n();
    int num_faces = faces.n();
    assert(vertices.m() == DIM);
    assert(faces.m() == 3);

    // 3 x (4*num_patches) double array to contain patch corner locations
    vtkSmartPointer<vtkDoubleArray> point_position = 
        vtkSmartPointer<vtkDoubleArray>::New();
    point_position->SetNumberOfComponents(3);
    point_position->SetNumberOfTuples(num_faces*3);

    vtkSmartPointer<vtkCellArray> triangle_ids = 
        vtkSmartPointer<vtkCellArray>::New();
    
    vtkSmartPointer<vtkIntArray> patch_ids= 
        vtkSmartPointer<vtkIntArray>::New();
    patch_ids->SetNumberOfComponents(1);
    patch_ids->SetName("Triangle Id");
    vtkSmartPointer<vtkIntArray> corr_patch_ids= 
        vtkSmartPointer<vtkIntArray>::New();
    corr_patch_ids->SetNumberOfComponents(1);
    corr_patch_ids->SetName("Patch Id");

    for (int i = 0; i < num_vertices; i++) {
        Point3 vertex_position(vertices.clmdata(i));
        point_position->SetTuple(i, vertex_position.array());
    }
    for (int f = 0; f < num_faces; f++) {
        vector<long long> face_ids(3);
        
        for (int d = 0; d < 3; d++) {
            face_ids[d] = faces(d,f);
        }

        triangle_ids->InsertNextCell(3, face_ids.data()); 
        patch_ids->InsertNextValue(f);
        if(corresponding_patches.size() == num_faces){
            corr_patch_ids->InsertNextValue(corresponding_patches[f]);
        }
        
    }
    
    points->SetData(point_position);
    //cout << "number of cells... " << cellArray->GetNumberOfCells() << endl;
    
    // actual format to save to a file
    vtkSmartPointer<vtkPolyData> polydata =
        vtkSmartPointer<vtkPolyData>::New();
    // insert points + patches
    polydata->SetPoints(points);
    polydata->SetPolys(triangle_ids);
    polydata->GetCellData()->AddArray(patch_ids);
    if(corresponding_patches.size() == num_faces){
        polydata->GetCellData()->AddArray(corr_patch_ids);
    polydata->GetCellData()->SetActiveScalars("Patch Id");
    }
    polydata->GetCellData()->SetActiveScalars("Triangle Id");
    
    
    // save
    vtkSmartPointer<vtkXMLPolyDataWriter> writer =  
        vtkSmartPointer<vtkXMLPolyDataWriter>::New();

    string file_name = build_filename(iteration, TRIANGLES, file_prefix);

    writer->SetFileName(file_name.c_str());
    writer->SetInputData(polydata);


    writer->SetDataModeToAscii();
    writer->Write();

}



void write_face_map_patches_to_vtk(DblNumMat qbkix_points, 
        vector<int> patches_refined_relative_ids,
        PatchSurfFaceMap* face_map, int iteration, string file_prefix){


    // list of vtkQuads to view 
    vtkSmartPointer<vtkCellArray> cellArray =
        vtkSmartPointer<vtkCellArray>::New();

    // list of vtkPoints to defining quads 
    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();

    // uv coordinates to evaluate at the patch corners... (counter clockwise
    // order);
    vector<Point2> uv_coords;
    uv_coords.push_back(Point2(0.,0.));
    uv_coords.push_back(Point2(1.,0.));
    uv_coords.push_back(Point2(1.,1.));
    uv_coords.push_back(Point2(0.,1.));

    int num_patches = face_map->patches().size();

    // 3 x (4*num_patches) double array to contain patch corner locations
    vtkSmartPointer<vtkDoubleArray> point_position = 
        vtkSmartPointer<vtkDoubleArray>::New();
    point_position->SetNumberOfComponents(3);
    point_position->SetNumberOfTuples(num_patches*4);

    vtkSmartPointer<vtkCellArray> quad_ids = 
        vtkSmartPointer<vtkCellArray>::New();
    
    vtkSmartPointer<vtkIntArray> patch_ids= 
        vtkSmartPointer<vtkIntArray>::New();
    patch_ids->SetNumberOfComponents(1);
    patch_ids->SetName("Patch Id");


    
    vtkSmartPointer<vtkIntArray> patch_relative_ids= 
        vtkSmartPointer<vtkIntArray>::New();
    patch_relative_ids->SetNumberOfComponents(1);
    patch_relative_ids->SetName("Refined relative Patch Id");

    //point_position->SetNumberOfTuples(num_patches);
    
    for(int pi = 0; pi <num_patches; pi++){

        FaceMapSubPatch* patch = (FaceMapSubPatch*) face_map->patches()[pi];

        // make a list of point ids defining the patch 
        vector<long long > point_ids(4);
        // evaluate the patch at its corners
        for(int i = 0; i < uv_coords.size(); i++){
            int index = 4*pi+i;
            point_ids[i] = index;
            Point3 position(0.);
            patch->xy_to_patch_coords(uv_coords[i].array(), 
                    PatchSamples::EVAL_VL,
                    position.array());
            
            // save the position
            point_position->SetTuple(index, position.array());
        }
        // save the quad id list
        quad_ids->InsertNextCell(4, point_ids.data());
        patch_ids->InsertNextValue(pi);
        patch_relative_ids->InsertNextValue(patches_refined_relative_ids.at(pi));
    }
    points->SetData(point_position);
    //cout << "number of cells... " << cellArray->GetNumberOfCells() << endl;
    
    // actual format to save to a file
    vtkSmartPointer<vtkPolyData> polydata =
        vtkSmartPointer<vtkPolyData>::New();
    // insert points + patches
    polydata->SetPoints(points);
    polydata->SetPolys(quad_ids);
    //polydata->GetCellData()->SetScalars(patch_ids);
    polydata->GetCellData()->AddArray(patch_relative_ids);
    polydata->GetCellData()->AddArray(patch_ids);
    polydata->GetCellData()->SetActiveScalars("Patch Id");
    
    
    // save
    vtkSmartPointer<vtkXMLPolyDataWriter> writer =  
        vtkSmartPointer<vtkXMLPolyDataWriter>::New();

    string file_name = build_filename(iteration, PATCHES, file_prefix);

    writer->SetFileName(file_name.c_str());
    writer->SetInputData(polydata);


    writer->SetDataModeToAscii();
    writer->Write();

}

void write_face_map_mesh_to_vtk( 
        PatchSurfFaceMap* face_map, int iteration,
        string file_prefix, int num_samples_1d){


    //int num_samples_1d = 20;
    double step = 1./(num_samples_1d-1);
    const int num_samples= num_samples_1d*num_samples_1d;
    const int triangles_per_bbox = 2*(num_samples_1d-1)*(num_samples_1d-1);
    DblNumMat vertices(3, num_samples);
    Rectangle domain = Sampling::base_domain;
    IntNumMat faces(3, triangles_per_bbox); //3 x num_faces matrix of positions
    // list of vtkQuads to view 
    vtkSmartPointer<vtkCellArray> cellArray =
        vtkSmartPointer<vtkCellArray>::New();

    // list of vtkPoints to defining quads 
    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();

    // uv coordinates to evaluate at the patch corners... (counter clockwise
    // order);

    int num_patches = face_map->patches().size();

    // 3 x (4*num_patches) double array to contain patch corner locations
    vtkSmartPointer<vtkDoubleArray> point_position = 
        vtkSmartPointer<vtkDoubleArray>::New();
    point_position->SetNumberOfComponents(3);
    point_position->SetNumberOfTuples(num_patches*num_samples);

    vtkSmartPointer<vtkCellArray> triangles= 
        vtkSmartPointer<vtkCellArray>::New();
    
    vtkSmartPointer<vtkIntArray> patch_ids= 
        vtkSmartPointer<vtkIntArray>::New();
    patch_ids->SetNumberOfComponents(1);
    patch_ids->SetName("Patch Id");


    /*
    vtkSmartPointer<vtkIntArray> patch_relative_ids= 
        vtkSmartPointer<vtkIntArray>::New();
    patch_relative_ids->SetNumberOfComponents(1);
    patch_relative_ids->SetName("Refined relative Patch Id");
    */
    //point_position->SetNumberOfTuples(num_patches);
    
    int num_faces=0;
    for(int pi = 0; pi <num_patches; pi++){

        FaceMapSubPatch* patch =face_map->subpatch(pi);

        // make a list of point ids defining the patch 

        patch->mesh_patch(step, domain, vertices, faces);
        // evaluate the patch at its corners
        for(int i = 0; i < num_samples; i++){
            long long index = i + pi*num_samples;
            Point3 position(vertices.clmdata(i));
            // save the position
            point_position->SetTuple(index, position.array());
        }
        for(int i = 0; i < triangles_per_bbox; i++){
            vector<long long > ids(3,0);
            for (int d = 0; d < 3; d++) {
                ids[d] = faces(d,i) + pi*num_samples;
            }
            triangles->InsertNextCell(3, ids.data());
        // save the quad id list
        patch_ids->InsertNextValue(pi);
        }
    }
    points->SetData(point_position);
    //cout << "number of cells... " << cellArray->GetNumberOfCells() << endl;
    
    // actual format to save to a file
    vtkSmartPointer<vtkPolyData> polydata =
        vtkSmartPointer<vtkPolyData>::New();
    // insert points + patches
    polydata->SetPoints(points);
    polydata->SetPolys(triangles);
    //polydata->GetCellData()->SetScalars(patch_ids);
    polydata->GetCellData()->AddArray(patch_ids);
    polydata->GetCellData()->SetActiveScalars("Patch Id");
    
    
    // save
    vtkSmartPointer<vtkXMLPolyDataWriter> writer =  
        vtkSmartPointer<vtkXMLPolyDataWriter>::New();

    string file_name = build_filename(iteration, TRIANGLES, file_prefix);

    writer->SetFileName(file_name.c_str());
    writer->SetInputData(polydata);


    writer->SetDataModeToAscii();
    writer->Write();

}



void write_face_map_patch_bounding_boxes_to_vtk(
        vector<int> patches_refined_relative_ids,
        PatchSurfFaceMap* face_map, int iteration, string file_prefix,
        bool inflate){


    int num_patches = face_map->patches().size();
    // list of vtkQuads to view 
    vtkSmartPointer<vtkCellArray> cellArray =
        vtkSmartPointer<vtkCellArray>::New();
    cellArray->SetNumberOfCells(num_patches*6);

    // list of vtkPoints to defining quads 
    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();

    // uv coordinates to evaluate at the patch corners... (counter clockwise
    // order);
    vector<Point2> uv_coords;
    uv_coords.push_back(Point2(0.,0.));
    uv_coords.push_back(Point2(1.,0.));
    uv_coords.push_back(Point2(1.,1.));
    uv_coords.push_back(Point2(0.,1.));


    // 3 x (4*num_patches) double array to contain patch corner locations
    vtkSmartPointer<vtkDoubleArray> point_position = 
        vtkSmartPointer<vtkDoubleArray>::New();
    point_position->SetNumberOfComponents(3);
    point_position->SetNumberOfTuples(num_patches*8);

    vtkSmartPointer<vtkCellArray> quad_ids = 
        vtkSmartPointer<vtkCellArray>::New();
    //quad_ids->SetNumberOfCells(num_patches*6);
    
    vtkSmartPointer<vtkIntArray> patch_ids= 
        vtkSmartPointer<vtkIntArray>::New();
    patch_ids->SetNumberOfComponents(1);
    patch_ids->SetName("Patch Id");


    
    vtkSmartPointer<vtkIntArray> patch_relative_ids= 
        vtkSmartPointer<vtkIntArray>::New();
    patch_relative_ids->SetNumberOfComponents(1);
    patch_relative_ids->SetName("Refined relative Patch Id");

    //point_position->SetNumberOfTuples(num_patches);
    cout << "bbox ids: "<< endl;
    for(int pi = 0; pi <num_patches; pi++){

        auto patch = face_map->subpatch(pi);

        Point3 bbox_min;
        Point3 bbox_max;
        if(inflate){
            patch->inflated_bounding_box(bbox_min, bbox_max);
        } else {
            patch->bounding_box(bbox_min, bbox_max);

        }

        /*double inflation_factor = 
            Options::get_double_from_petsc_opts("-adaptive_upsampling_bbox_inflation_factor");
        // Increase by the size of the near zone
        //cout << "bbox_min: " << bbox_min << endl;
        //cout << "bbox_max: " << bbox_max << endl;
        bbox_max += Point3(inflation_factor*patch->characteristic_length());
        bbox_min -= Point3(inflation_factor*patch->characteristic_length());*/

        // I hate me.
        point_position->SetTuple(8*pi, Point3(bbox_min).array());
        point_position->SetTuple(8*pi+1, Point3(bbox_min(0), bbox_min(1), bbox_max(2)).array());
        point_position->SetTuple(8*pi+2, Point3(bbox_min(0), bbox_max(1), bbox_min(2)).array());
        point_position->SetTuple(8*pi+3, Point3(bbox_max(0), bbox_min(1), bbox_min(2)).array());
        point_position->SetTuple(8*pi+4, Point3(bbox_min(0), bbox_max(1), bbox_max(2)).array());
        point_position->SetTuple(8*pi+5, Point3(bbox_max(0), bbox_min(1), bbox_max(2)).array());
        point_position->SetTuple(8*pi+6, Point3(bbox_max(0), bbox_max(1), bbox_min(2)).array());
        point_position->SetTuple(8*pi+7, Point3(bbox_max).array());

        vector<long long> ids(4,0);// = {8*pi, 8*pi+1, 8*pi+5, 8*pi+3};
        ids[0] = 8*pi; ids[1] = 8*pi+1; ids[2] =  8*pi+5; ids[3] =  8*pi+3;
        quad_ids->InsertNextCell(4, ids.data());
        patch_ids->InsertNextValue(pi);
        patch_relative_ids->InsertNextValue(patches_refined_relative_ids.at(pi));
         
        ids[0] = 8*pi; ids[1] = 8*pi+2; ids[2] =  8*pi+6; ids[3] =  8*pi+3;
        /*for(const auto& i : ids){
            cout << i << ", ";
        }
        cout << endl;*/
        quad_ids->InsertNextCell(4, ids.data());
        patch_ids->InsertNextValue(pi);
        patch_relative_ids->InsertNextValue(patches_refined_relative_ids.at(pi));
         
        ids[0] = 8*pi; ids[1] = 8*pi+1; ids[2] =  8*pi+4; ids[3] =  8*pi+2;
        quad_ids->InsertNextCell(4, ids.data());
        patch_ids->InsertNextValue(pi);
        patch_relative_ids->InsertNextValue(patches_refined_relative_ids.at(pi));
        
        ids[0] = 8*pi+3; ids[1] = 8*pi+6; ids[2] =  8*pi+7; ids[3] =  8*pi+5;
        quad_ids->InsertNextCell(4, ids.data());
        patch_ids->InsertNextValue(pi);
        patch_relative_ids->InsertNextValue(patches_refined_relative_ids.at(pi));
         
        ids[0] = 8*pi+1; ids[1] = 8*pi+5; ids[2] =  8*pi+7; ids[3] =  8*pi+4;
        quad_ids->InsertNextCell(4, ids.data());
        patch_ids->InsertNextValue(pi);
        patch_relative_ids->InsertNextValue(patches_refined_relative_ids.at(pi));
         
        ids[0] = 8*pi+2; ids[1] = 8*pi+4; ids[2] =  8*pi+7; ids[3] =  8*pi+6;
        quad_ids->InsertNextCell(4, ids.data());
        patch_ids->InsertNextValue(pi);
        patch_relative_ids->InsertNextValue(patches_refined_relative_ids.at(pi));

    }
    points->SetData(point_position);
    //cout << "number of cells... " << cellArray->GetNumberOfCells() << endl;
    
    // actual format to save to a file
    vtkSmartPointer<vtkPolyData> polydata =
        vtkSmartPointer<vtkPolyData>::New();
    // insert points + patches
    polydata->SetPoints(points);
    polydata->SetPolys(quad_ids);
    //polydata->GetCellData()->SetScalars(patch_ids);
    polydata->GetCellData()->AddArray(patch_relative_ids);
    polydata->GetCellData()->AddArray(patch_ids);
    polydata->GetCellData()->SetActiveScalars("Patch Id");
    
    
    // save
    vtkSmartPointer<vtkXMLPolyDataWriter> writer =  
        vtkSmartPointer<vtkXMLPolyDataWriter>::New();

    string file_name = build_filename(iteration, BOXES, file_prefix);

    writer->SetFileName(file_name.c_str());
    writer->SetInputData(polydata);


    writer->SetDataModeToAscii();
    writer->Write();

    /*
       vtkSmartPointer<vtkAppendFilter> appendFilter =
       vtkSmartPointer<vtkAppendFilter>::New();
       appendFilter->AddInputData(polydata);
       appendFilter->Update();
       */ 
    //vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
    //vtkSmartPointer<vtkUnstructuredGrid>::New();
    //writer->SetInput(idf->GetOutput());
    //writer->SetFileName(fname.c_str());
    //vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
    //vtkSmartPointer<vtkUnstructuredGrid>::New();
    //unstructuredGrid->ShallowCopy(appendFilter->GetOutput());
    // Write file
    //vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
    //vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
}



void write_lines_from_qbkix_to_closest_point(DblNumMat qbkix_points, 
        NumVec<OnSurfacePoint> final_closest_on_surface_points, 
        PatchSurfFaceMap* face_map, int iteration, string file_prefix){

    int num_qbkix_points = qbkix_points.n();
    assert(qbkix_points.n() == final_closest_on_surface_points.m());
    cout << num_qbkix_points << endl;
    // list of lines from qbkix points to closest on-surface points to view 
    vtkSmartPointer<vtkCellArray> lines =
        vtkSmartPointer<vtkCellArray>::New();

    // list of qbkix + on-surface points
    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();

    int num_patches = face_map->patches().size();

    // 3 x (2*num_qbkix_points) double array to contain qbkix points and
    // on-surface points
    vtkSmartPointer<vtkDoubleArray> point_position = 
        vtkSmartPointer<vtkDoubleArray>::New();
    point_position->SetNumberOfComponents(3);
    point_position->SetNumberOfTuples(num_qbkix_points);

    vtkSmartPointer<vtkDoubleArray> point_vector = 
        vtkSmartPointer<vtkDoubleArray>::New();
    point_vector->SetNumberOfComponents(3);
    point_vector->SetNumberOfTuples(num_qbkix_points);
    point_vector->SetName("qbkix direction");
    for(int qi = 0; qi <num_qbkix_points; qi++){
        // insert qbkix point location
        Point3 qbkix_point(qbkix_points.clmdata(qi)); 
        point_position->SetTuple(qi, qbkix_point.array());

        // get the patch containing the closest on-surface point to the qbkix
        // point
        OnSurfacePoint on_surface_point = final_closest_on_surface_points(qi);
        //cout << "patch conataining closest point: " << on_surface_point.parent_patch << endl;
        int patch_conataining_closest_point = on_surface_point.parent_patch;
        Point3 dir;
        if(patch_conataining_closest_point > -1){
            FaceMapSubPatch* patch = (FaceMapSubPatch*) face_map->patches()[on_surface_point.parent_patch];

            // compute the position of closest on_surface_point
            Point3 on_surface_point_position(0.);
            patch->xy_to_patch_coords(on_surface_point.parametric_coordinates.array(), 
                    PatchSamples::EVAL_VL,
                    on_surface_point_position.array());
            //cout << "qi: " << qi << ", " << on_surface_point_position << endl;
            // insert closest on_surface_point_position
            dir = on_surface_point_position - qbkix_point;
        } else {
            dir = Point3(0,0,.25);
        }
        point_vector->SetTuple(qi, dir.array());
        //vector<long long> point_ids(2,0);
        //point_ids[0] = 2*qi;
        //point_ids[1] = 2*qi+1;
        // save point ids defining the line
        //lines->InsertNextCell(2, point_ids.data());
    }
    points->SetData(point_position);
    
    // actual format to save to a file
    vtkSmartPointer<vtkPolyData> polydata =
        vtkSmartPointer<vtkPolyData>::New();
    // insert points + patches
    polydata->SetPoints(points);
    polydata->GetPointData()->SetVectors(point_vector);
    //polydata->SetPolys(lines);

 
  // Write the file
  vtkSmartPointer<vtkXMLPolyDataWriter> writer =  
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  
  string file_name = build_filename(iteration, LINES, file_prefix);
  
  writer->SetFileName(file_name.c_str());
  writer->SetInputData(polydata);
 
  // Optional - set the mode. The default is binary.
  //writer->SetDataModeToBinary();
  writer->SetDataModeToAscii();
 
  writer->Write();

}


END_EBI_NAMESPACE
