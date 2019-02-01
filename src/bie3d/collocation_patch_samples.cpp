#include "collocation_patch_samples.hpp"

BEGIN_EBI_NAMESPACE

// ----------------------------------------------------------------------

int CollocationPatchSamples::distribute_collocation_points(
        Vec _target_3d_position, Vec _target_as_face_point, 
        PatchSamples* patch_samples, int sdof, int tdof){
        //CollocationPointData& _collocation_data){
  ebiFunctionBegin;
  
  //int sdof = this->source_dof();  
  //int tdof = this->target_dof();
  
  //1. build local collocation list
  build_local_collocation_lists(_target_as_face_point,
                                 this->colloc_point_global_indices(), 
                                 this->num_colloc_point_in_patch(), 
                                 this->colloc_point_starting_index(),
                                 patch_samples);

  //2. create _colloc_point_as_face_point, _colloc_point_3d_position
  int face_point_size_in_doubles = patch_samples->face_point_size_in_doubles();
  int local_num_collocation_points = this->colloc_point_global_indices().size();

  VecCreateMPI(mpiComm(),
          local_num_collocation_points * face_point_size_in_doubles, 
          PETSC_DETERMINE,
          &(this->colloc_point_as_face_point()));
  VecCreateMPI(mpiComm(), 
          local_num_collocation_points * dim(),   
          PETSC_DETERMINE, 
          &(this->colloc_point_3d_position()));

  // MJM: this is terribly sketchy need to convert ints to int64_t's
  // TODO remove
  vector<int64_t> tmp_vec; //= this->colloc_point_global_indices();
  for(int i = 0; i < int(this->colloc_point_global_indices().size()); i++){
      tmp_vec.push_back(this->colloc_point_global_indices()[i]);
  }
  vector<int64_t>& tmpvec = tmp_vec;

  IS colloc_index_set;
  IS target_index_set;

  Vec _colloc_point_as_face_point = this->colloc_point_as_face_point();
  Vec _colloc_point_3d_position = this->colloc_point_3d_position();

  //index set 0..length(colloc_point_global_indices)*size of facepoint in doubles  
  //@BUG DZ perhaps this should not start from zero, but rather
  // from the global index offset of the part of target_as_face_point in the current partition?
  ISCreateStride(mpiComm(), 
          tmpvec.size() * face_point_size_in_doubles,
          0, 
          1,
          &colloc_index_set);

  // block size = facepoint in doubles
  // # of blocks = #of colloc points in the current partition
  // index set =  indices of colloc. points in the global vector of target
  // points*facepoint in doubles
  ISCreateBlock( mpiComm(), 
          face_point_size_in_doubles, 
          tmpvec.size(),
          &(tmpvec[0]), 
          PETSC_COPY_VALUES, 
          &target_index_set);

  VecScatter target_to_colloc_face_point;
  // scatter from target facepoints to collocation facepoints
  // with target_index_set and colloc_index_set index sets
  // Possible problem: target_index_set seems to be local and start  from zero,
  // colloc_index_set are global
  // the idea is that this will copy target_as_face_point to locations in colloc_point_as_face_point
  // defined by indices from c2tmap

  VecScatterCreate(_target_as_face_point, 
          target_index_set, 
          _colloc_point_as_face_point, 
          colloc_index_set, 
          &target_to_colloc_face_point);

  VecScatterBegin( target_to_colloc_face_point, 
          _target_as_face_point,  
          _colloc_point_as_face_point,  
          INSERT_VALUES, 
          SCATTER_FORWARD);

  VecScatterEnd( target_to_colloc_face_point,   
          _target_as_face_point,  
          _colloc_point_as_face_point,  
          INSERT_VALUES,  
          SCATTER_FORWARD);

  iC( VecScatterDestroy(&target_to_colloc_face_point) );
  iC( ISDestroy(&colloc_index_set) );
  iC( ISDestroy(&target_index_set) );

  // do the same for target_3d_position -> colloc_point_3d_position

  ISCreateStride(mpiComm(), tmpvec.size() * dim(), 0, 1, &colloc_index_set);
  
  ISCreateBlock( mpiComm(), 
          dim(), 
          tmpvec.size(), 
          &(tmpvec[0]), 
          PETSC_COPY_VALUES, 
          &target_index_set);
  
  VecScatter target_to_colloc_pos;
  VecScatterCreate(_target_3d_position, 
          target_index_set, 
          _colloc_point_3d_position, 
          colloc_index_set,
          &target_to_colloc_pos);

  VecScatterBegin( target_to_colloc_pos, 
          _target_3d_position,  
          _colloc_point_3d_position,  
          INSERT_VALUES,  
          SCATTER_FORWARD);
  
  VecScatterEnd( target_to_colloc_pos,   
          _target_3d_position,  
          _colloc_point_3d_position,  
          INSERT_VALUES,  
          SCATTER_FORWARD);

  VecScatterDestroy(&target_to_colloc_pos);
  ISDestroy(&colloc_index_set);
  ISDestroy(&target_index_set);
  
  //3. create scatters to be used in the evaluation 
  // from collocation point to targets

  
  int num_local_colloc_points = num_local_points(_colloc_point_3d_position);
  int num_local_targets = num_local_points(_target_3d_position);

  //sdof
  // dummy vectors just to create the Scatter objects for future use
  // they just have to have correct sizes 
  Vec coltmp, trgtmp;

  iC( VecCreateMPI(mpiComm(), num_local_colloc_points*sdof, PETSC_DETERMINE, &coltmp) );  
  iC( VecCreateMPI(mpiComm(), num_local_targets*sdof, PETSC_DETERMINE, &trgtmp) );
  iC( ISCreateStride(mpiComm(), tmpvec.size()*sdof, 0, 1, &colloc_index_set) );

  // use c2map to create the index set to go from colloc points to target points
 
  iC( ISCreateBlock( mpiComm(), sdof, tmpvec.size(), &(tmpvec[0]), PETSC_COPY_VALUES, &target_index_set) );
  
  //CREATION OF THE SCATTER
  iC( VecScatterCreate(coltmp, colloc_index_set, trgtmp, target_index_set, &(this->colloc_to_target_density_scatter())) ); 

  iC( ISDestroy(&colloc_index_set) );
  iC( ISDestroy(&target_index_set) );
  iC( VecDestroy(&coltmp) );
  iC( VecDestroy(&trgtmp) );

  // now create a scatter for values, dimension may be different at target
  // using tdof instead of sdof

  iC( VecCreateMPI(mpiComm(), num_local_colloc_points*tdof, PETSC_DETERMINE, &coltmp) );
  iC( VecCreateMPI(mpiComm(), num_local_targets*tdof, PETSC_DETERMINE, &trgtmp) );

  iC( ISCreateStride(mpiComm(), tmpvec.size()*tdof, 0, 1, &colloc_index_set) );
  
  iC( ISCreateBlock( mpiComm(), tdof, tmpvec.size(), &(tmpvec[0]), PETSC_COPY_VALUES, &target_index_set) );
  
  //CREATION OF THE SCATTER
  iC( VecScatterCreate(coltmp, colloc_index_set, trgtmp, target_index_set, &(this->colloc_to_target_value_scatter())) ); 

  iC( ISDestroy(&colloc_index_set) );
  iC( ISDestroy(&target_index_set) );
  iC( VecDestroy(&coltmp) );
  iC( VecDestroy(&trgtmp) );
  
  ebiFunctionReturn(0);
}


int CollocationPatchSamples::build_local_collocation_lists(
        Vec target_as_face_point, 
        vector<int>& colloc_point_global_indices,
        vector<int>& num_colloc_point_in_patch, 
        vector<int>& colloc_point_starting_index,
        PatchSamples* patch_samples) {
  ebiFunctionBegin;
  
  vector<Patch*>& patches = patch_samples->bdry()->patches(); //get patches
  int num_patches = patches.size();

  int face_point_size_in_doubles = patch_samples->face_point_size_in_doubles();
  
  //1. target_to_colloc_map
  vector< vector<int> > target_to_colloc_mapvec; 
  target_to_colloc_mapvec.resize( mpiSize() );
  
  // get the number of target points in the partition
  int64_t lcltrgnum;
  iC(VecGetLocalSize(target_as_face_point,&lcltrgnum)); 
  lcltrgnum = lcltrgnum/face_point_size_in_doubles;
  
  // start/end global indices
  int64_t lcltrgbeg, lcltrgend;
  iC(VecGetOwnershipRange(target_as_face_point, &lcltrgbeg, &lcltrgend)); 
  lcltrgbeg = lcltrgbeg/face_point_size_in_doubles; 
  lcltrgend = lcltrgend/face_point_size_in_doubles;

  // build tcmapvec: tcmapvec[partition#] = list of pairs (patch_id, global sample point index)
  // a pair is created for each patch containing the point, i.e. 4 patches per point
  // and assigned to the list corresponding to the partition owning the patch
  // this is done for all points owned by the current partition
  double* target_as_face_point_arr;
  iC( VecGetArray(target_as_face_point, &target_as_face_point_arr) );

  for(int ti=0; ti<lcltrgnum; ti++) {
	 FacePointOverlapping* cur_face_point =
         (FacePointOverlapping*)(target_as_face_point_arr +
                 ti*face_point_size_in_doubles);

	 int curidx = lcltrgbeg + ti;
	 vector<int> piv; //patch id
	 patch_samples->bdry()->patches_containing_face_point(cur_face_point, piv);

	 for(int ii=0; ii<piv.size(); ii++) {
		target_to_colloc_mapvec[ patch_samples->patch_partition()[piv[ii]] ].push_back( piv[ii] );
		target_to_colloc_mapvec[ patch_samples->patch_partition()[piv[ii]] ].push_back( curidx );
	 }
  }
  iC( VecRestoreArray(target_as_face_point, &target_as_face_point_arr) );

  // sending sizes of lists per partition we have figured out  
  vector<int> sendszvec(mpiSize());
  for(int pi=0; pi<mpiSize(); pi++) {
	 sendszvec[pi] = target_to_colloc_mapvec[pi].size();
  }
  // receiving the sizes of the lists for the current partition, computed 
  // by other processes 
  vector<int> recvszvec(mpiSize());

  // on completion, recvszvec contains the list of sizes of lists from all
  // processes for current process
  iC( MPI_Alltoall( &(sendszvec[0]), 1, MPI_INT, &(recvszvec[0]), 1,
              MPI_INT, mpiComm() ) );

  // compute offsets of the individual lists in the concatenation of all 
  // lists in target_to_colloc_ampvec needed for communication
  vector<int> sendofvec(mpiSize());
  int sendofveccnt = 0;
  for(int pi=0; pi<mpiSize(); pi++) {
	 sendofvec[pi] = sendofveccnt;
	 sendofveccnt += sendszvec[pi];
  }
  int sendcnt = sendofveccnt; //ebiAssert((sendcnt%2)==0);

  // based on the size information we have received, compute offsets in the 
  // receive vector to be used for all2all output

  vector<int> recvofvec(mpiSize());
  int recvofveccnt = 0;
  for(int pi=0; pi<mpiSize(); pi++) {
	 recvofvec[pi] = recvofveccnt;
	 recvofveccnt += recvszvec[pi]; //MJM how is this possible? we only recieved
                                    //  one element from the alltoall?
  }
  int recvcnt = recvofveccnt; //ebiAssert((recvcnt%2)==0);

  // concatenate all lists in target_to_colloc_mapvec into the send buffer
  vector<int> sendbuf;
  for(int pi=0; pi<mpiSize(); pi++) {
	 sendbuf.insert( sendbuf.end(), target_to_colloc_mapvec[pi].begin(), target_to_colloc_mapvec[pi].end() );
  }
  assert( sendbuf.size() == sendcnt );
  vector<int> recvbuf( recvcnt );


  // if this were executed, on completion, recvbuf would contain the concatendate lists of points compiled for the current process
  // by other processes 
  //iC( MPI_Alltoallv( &(sendbuf[0]), &(sendszvec[0]), &(sendofvec[0]), MPI_INT,
  //					&(recvbuf[0]), &(recvszvec[0]), &(recvofvec[0]), MPI_INT, mpiComm() ) );

  // it seems this can only work for a single process, it only makes sense there is no distribution
  for (int ti=0; ti< recvcnt; ti++)
	 {
		recvbuf[ti] = sendbuf[ti];
	 }


  // construct a map patch id ->  list of global indices of points associated with this patch in the combined list received
  // in the recvbuf in all2all
  // for each patch, this consists of sample points of the patch itself as well as all patches overlapping this patch
  
  vector< vector<int> > c2tmapvec;
  c2tmapvec.resize( num_patches );
  for(int ti=0; ti<recvcnt; ti=ti+2) {
	 int pi = recvbuf[ti];
	 int ui = recvbuf[ti+1];
	 c2tmapvec[pi].push_back(ui);
  }
  // number of collocation points for the current partition 
  int local_num_collocation_points = recvcnt/2;
  int lclcolstt;
  // sums up numbers of collocation points for partitions with mpiRank <= current into lclcolstt
  iC( MPI_Scan( &local_num_collocation_points, &lclcolstt, 1, MPI_INT, MPI_SUM, mpiComm() ) );
  // in the global vector of collocation points, this is the start of the subvector for collocation points for this partition
  lclcolstt -= local_num_collocation_points;
  
  // Finally, build num_colloc_point_in_patch, colloc_point_starting_index and colloc_point_global_indices
  // num_colloc_point_in_patch: total number of collocation points for a patch (i.e., its own sample points, and points of overlapping patches that
  // are on this patch's faces
  // colloc_point_starting_index:  offsets of the lists of collocation points for a patch in the global vector of collocation points
  num_colloc_point_in_patch.resize(num_patches, 0); 
  colloc_point_starting_index.resize(num_patches, 0);

  // compute num_colloc_point_in_patch and colloc_point_starting_index for patches in the current partition 
  for(int pi=0; pi<num_patches; pi++) {
	 if(patch_samples->patch_partition()[pi]!=mpiRank()) {
		ebiAssert( c2tmapvec[pi].size()==0 );
	 } else {
		num_colloc_point_in_patch[pi] = c2tmapvec[pi].size();
		colloc_point_starting_index[pi] = lclcolstt;
		lclcolstt += num_colloc_point_in_patch[pi];
	 }
  }

  vector<int> tmp_num_colloc_point_in_patch(num_colloc_point_in_patch);
  vector<int> tmp_colloc_point_starting_index(colloc_point_starting_index);
  
  // get num_colloc_point_in_patch and colloc_point_starting_index from all processes; for each patch only a single process (the one owning it) has nonzero entries
  // in these arrays, so the sums have only single nonzer terms
  iC( MPI_Allreduce( &(tmp_num_colloc_point_in_patch[0]),
              &(num_colloc_point_in_patch[0]),
              num_patches,
              MPI_INT,
              MPI_SUM,
              mpiComm() ) );

  iC( MPI_Allreduce( &(tmp_colloc_point_starting_index[0]),
              &(colloc_point_starting_index[0]),
              num_patches, 
              MPI_INT, 
              MPI_SUM, 
              mpiComm() ) );
  //
  colloc_point_global_indices.clear();

  // for all patches in the current partition concatenate the lists from c2tmapvec,
  // colloc_point_global_indices is the list of collocation points for the current partition, all patches together 
  for(int pi=0; pi<num_patches; pi++) {
	 if( patch_samples->patch_partition()[pi] == mpiRank()) {
		vector<int>& tmp = c2tmapvec[pi];
		colloc_point_global_indices.insert( colloc_point_global_indices.end(), tmp.begin(), tmp.end() );
	 }
  }
  ebiAssert(colloc_point_global_indices.size()==local_num_collocation_points);
  
  ebiFunctionReturn(0);
}




// ---------------------------------------------------------------------------
int CollocationPatchSamples::num_collocation_points(int pi) {
  return this->num_colloc_point_in_patch()[pi];
}

int CollocationPatchSamples::collocation_points_starting_index(int pi) {
  return this->colloc_point_starting_index()[pi];
}

int CollocationPatchSamples::get_patch_offset(int pi){
  int a, b;
  local_index_range(colloc_point_3d_position(), a, b);

  int s = collocation_points_starting_index(pi);
  assert(s>=a && s<=b);

  return s-a;
}

double* get_pointer_to_local_vec(Vec v){
  double* arr;     
  VecGetArray(v, &arr);

  double* buf=arr;
  VecRestoreArray(v, &arr);

  return buf;
}

DblNumMat CollocationPatchSamples::colloc_point_as_face_point(int pi, int face_point_size_in_doubles) {
  int dif = get_patch_offset(pi);
  double* buf = get_pointer_to_local_vec(colloc_point_as_face_point());

  //int face_point_size_in_doubles = _patch_samples->face_point_size_in_doubles();
  return DblNumMat(face_point_size_in_doubles, num_collocation_points(pi),
          false, buf + dif * face_point_size_in_doubles);
}

DblNumMat CollocationPatchSamples::colloc_point_3d_position(int pi) {
  int dif = get_patch_offset(pi);
  double* buf = get_pointer_to_local_vec(colloc_point_3d_position());

  return DblNumMat(dim(), num_collocation_points(pi), false, buf + dif*dim());
}

DblNumMat CollocationPatchSamples::coldat(int pi, Vec dat, int dof) {
  int dif = get_patch_offset(pi);
  double* buf = get_pointer_to_local_vec(dat);

  return DblNumMat(dof, num_collocation_points(pi), false, buf + dif*dof);
}

// ---------------------------------------------------------------------- 
#undef __FUNCT__
END_EBI_NAMESPACE
