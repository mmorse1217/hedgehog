#include "fmm_interface.hpp"
#include "common/ebi_petsc.hpp"
#include "bie3d/solver_utils.hpp"
BEGIN_EBI_NAMESPACE

// ----------------------------------------------------------------------
FMM::FMM(){
}
FMM::~FMM(){
}
void FMM::initialize_fmm(Vec source_positions, Vec source_normals,
        Vec target_positions, Kernel3d kernel, bool call_setup){



   _max_level  =         Options::get_int_from_petsc_opts("-bis3d_maxlevel"); 
   _multipole_order =    Options::get_int_from_petsc_opts("-bis3d_np");
   _max_points_per_box = Options::get_int_from_petsc_opts("-bis3d_ptsmax");



    // Gather global source points/normals and target points for all processes

   // MJM NOTE probably breaks with mpi
  
  DblNumMat srcPos = get_local_vector(DIM, num_local_points(source_positions), source_positions);
  DblNumMat srcNor = get_local_vector(DIM, num_local_points(source_normals), source_normals);
  DblNumMat trgPos = get_local_vector(DIM, num_local_points(target_positions), target_positions);
  

  DblNumMat* srcPos_ptr = &srcPos;
  DblNumMat* srcNor_ptr = &srcNor;
  DblNumMat* trgPos_ptr = &trgPos;
  
  set_kernel(kernel);
  set_src_positions(srcPos_ptr);
  set_src_normals(srcNor_ptr);
  set_trg_positions(trgPos_ptr);
  setFromOptions(); 

  if(call_setup){
      setup();
  }

  srcPos.restore_local_vector();
  srcNor.restore_local_vector();
  trgPos.restore_local_vector();
}


END_EBI_NAMESPACE
