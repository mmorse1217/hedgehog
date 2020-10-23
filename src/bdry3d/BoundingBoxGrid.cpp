#include "bdry3d/BoundingBoxGrid.hpp"
#include "common/ebi_petsc.hpp"

BEGIN_EBI_NAMESPACE

template<typename Real_t>
BoundingBoxGrid<Real_t>::BoundingBoxGrid(Real_t box_size, MPI_Comm c, Real_t cell_size){
    // init periodic box size
    box_size_=box_size;
    // init communication
    comm=c;

    // number of processes
    MPI_Comm_size(comm, &np_);
    // process rank
    MPI_Comm_rank(comm, &rank_);
    // number of openmp threads
    omp_p_ = omp_get_max_threads();

    r_near_ = cell_size;
}

template<typename Real_t>
void BoundingBoxGrid<Real_t>::SetBoundingBoxGrid(const Real_t *pos_bd, const int num_points_per_patch_1d, const int num_patches)
{
    N_bbox_ = num_patches;
    
    // calculate the bounding boxes for vesicles with start, end position and min_sep
    BB_min_.ReInit(N_bbox_*DIM);
    BB_max_.ReInit(N_bbox_*DIM);

    // boundary patches bounding box
    #pragma omp parallel for
    for(size_t tid=0; tid<omp_p_; tid++){
        size_t a = ((tid+0)*num_patches)/omp_p_;
        size_t b = ((tid+1)*num_patches)/omp_p_;
        for(size_t i=a; i<b; i++){
            Real_t *mini = &BB_min_[i*DIM];
            Real_t *maxi = &BB_max_[i*DIM];

            Real_t min_tmp[3];
            Real_t max_tmp[3];
            const Real_t *first_point = &pos_bd[i*(num_points_per_patch_1d*num_points_per_patch_1d)];
            for(size_t k=0; k<DIM; k++)
            {
                min_tmp[k] = first_point[k];
                max_tmp[k] = first_point[k];
            }
            for(int i1d=0; i1d<num_points_per_patch_1d; i1d++)
            {
                for(int j1d=0; j1d<num_points_per_patch_1d; j1d++)
                {
                    int point_id = num_points_per_patch_1d*i1d+j1d;
                    for(size_t k=0; k<DIM; k++)
                    {
                        min_tmp[k] = std::min(first_point[point_id*DIM+k], min_tmp[k]);
                        max_tmp[k] = std::max(first_point[point_id*DIM+k], max_tmp[k]);
                    }
                }
            }
            for(size_t k=0; k<DIM; k++)
            {
                mini[k] = min_tmp[k]  - 1e-8;
                maxi[k] = max_tmp[k]  + 1e-8;
            }
        }
    }

}

template<typename Real_t>
void BoundingBoxGrid<Real_t>::SetBoundingBoxGrid(const PVFMMVec_t& BB_min, const PVFMMVec_t& BB_max)
{
    ASSERT(BB_min.Dim() == BB_max.Dim(), "bounding box min, max dim doesn't match");
    BB_min_ = BB_min;
    BB_max_ = BB_max;
    N_bbox_ = BB_min_.Dim()/DIM;
}

template<typename Real_t>
void BoundingBoxGrid<Real_t>::SetTrgCoord(Real_t* trg_coord, size_t N)
{
    T_.ReInit(N*DIM, trg_coord);
}

template<typename Real_t>
void BoundingBoxGrid<Real_t>::ConstructBoundingBoxGrid()
{
    ASSERT(BB_min_.Dim()>0, "empty bounding boxes min"); ASSERT(BB_max_.Dim()>0, "empty bounding boxes max");
    ASSERT(N_bbox_>0, "empty bounding boxes N_bbox_");

#ifdef BOUNDING_BOX_GRID_DEBUG
    std::cout<<"rank: "<<rank_<<" of "<<np_<<" processes has "<<omp_p_<<" threads\n";
#endif

    pvfmm::Profile::Tic("GetContactBBPair", &comm, true);
    //bool prof_state=pvfmm::Profile::Enable(false);

    // construct BB_let
    pvfmm::Profile::Tic("BBLET", &comm, true);
    SetTreeParams();

#ifdef BOUNDING_BOX_GRID_DEBUG
    std::cout<<"bbox_: "<<bbox_[0]<<","<<bbox_[1]<<","<<bbox_[2]<<","<<bbox_[3]<<".\n";
    std::cout<<"r_near_: "<<r_near_<<"\n";
    std::cout<<"tree_depth_: "<<tree_depth_<<"\n";
#endif
        
    GenerateBBPoints();
    ConstructLocalTree(BB_let_);

#ifdef BOUNDING_BOX_GRID_DEBUG
    std::cout<<"size of mid: "<<BB_let_.mid.Dim()<<"\n";
    std::cout<<BB_let_.mid<<"\n";
    std::cout<<"size of pt_cnt: "<<BB_let_.pt_cnt.Dim()<<"\n";
    std::cout<<BB_let_.pt_cnt<<"\n";
    std::cout<<"size of pt_dsp: "<<BB_let_.pt_dsp.Dim()<<"\n";
    std::cout<<BB_let_.pt_dsp<<"\n";
    std::cout<<"size of mins: "<<BB_let_.mins.Dim()<<"\n";
    std::cout<<BB_let_.mins<<"\n";
    std::cout<<"size of pt_id: "<<BB_let_.pt_id.Dim()<<"\n";
    std::cout<<BB_let_.pt_id<<"\n";
#endif
    
    pvfmm::Profile::Toc();
    // end of construct BB_let

    //pvfmm::Profile::Enable(prof_state);
    pvfmm::Profile::Toc();
}

template<typename Real_t>
void BoundingBoxGrid<Real_t>::SetTreeParams()
{
    pvfmm::Profile::Tic("TreeParams", &comm, true);

    Real_t* bbox = bbox_;
    Real_t& r_near = r_near_;
    size_t& tree_depth = tree_depth_;

    tree_depth = 0;
    assert(N_bbox_ > 0);
        
    // determine r_bbox
    std::vector<Real_t> r2_max_mp(omp_p_);
    std::vector<Real_t> r2_min_mp(omp_p_);
    #pragma omp parallel for
    for(size_t tid=0; tid<omp_p_; tid++){
        size_t a = ((tid+0)*N_bbox_)/omp_p_;
        size_t b = ((tid+1)*N_bbox_)/omp_p_;

        Real_t r2_max = 0;
        Real_t r2_min = std::numeric_limits<Real_t>::max();
        for(size_t i=a; i<b; i++){ // compute r2_bbox
            const Real_t* max_i = &BB_max_[i*DIM];
            const Real_t* min_i = &BB_min_[i*DIM];

            Real_t dx = max_i[0] - min_i[0];
            Real_t dy = max_i[1] - min_i[1];
            Real_t dz = max_i[2] - min_i[2];
            Real_t r2_box = dx*dx + dy*dy + dz*dz;

            r2_max = std::max(r2_max, r2_box);
            r2_min = std::min(r2_min, r2_box);
        }
        r2_max_mp[tid] = r2_max;
        r2_min_mp[tid] = r2_min;
    }

    // determine r_near (global max)
    double r_max_loc = 0; double r_max_glb = 0;
    double r_min_loc = std::numeric_limits<double>::max(); double r_min_glb = std::numeric_limits<double>::max();
    for(size_t tid=0; tid<omp_p_; tid++){
        r_max_loc = std::max(r2_max_mp[tid], r_max_loc);
        r_min_loc = std::min(r2_min_mp[tid], r_min_loc);
    }
    r_max_loc = std::sqrt(r_max_loc);
    r_min_loc = std::sqrt(r_min_loc);

    MPI_Allreduce(&r_max_loc, &r_max_glb, 1, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(&r_min_loc, &r_min_glb, 1, MPI_DOUBLE, MPI_MIN, comm);

    if(r_near < 0 || r_near > (r_min_glb/2))
        r_near = r_min_glb/2;

    if(box_size_>0 && r_max_glb > box_size_){
        std::cout<<"Domain too small for bounding box size\n";
        assert(false);
        exit(0);
    }

    if(box_size_<=0){ // determine the global bbox, tree_depth
        Real_t scale_x, shift_x[DIM];
        Real_t scale_tmp;
   
        // determine global bounding box
        GlobalBoundingBox(&scale_tmp, shift_x);

        //std::cout<<"scale_tmp: "<<scale_tmp<<"\n";
        //std::cout<<"shift_x: "<<shift_x[0]<<", "<<shift_x[1]<<", "<<shift_x[2]<<"\n";

        { // scale_x, pt_tree_depth, leaf_size
            ASSERT(scale_tmp!=0, "invalid scale");
            
            Real_t domain_length = 1.0/scale_tmp + 4*r_near;
            //COUT("domain_length_tmp: "<<domain_length);
            Real_t leaf_size = r_near;
            scale_x = 1.0/leaf_size;
            while(domain_length*scale_x>1.0 && tree_depth<MAX_DEPTH-1){
                scale_x *= 0.5;
                tree_depth++;
            }

            leaf_size_ = leaf_size;
        }

        for(size_t j=0;j<DIM;j++){ // Update shift_x
            shift_x[j]=((shift_x[j]/scale_tmp)+2*r_near)*scale_x;
        }

        bbox[0]=shift_x[0];
        bbox[1]=shift_x[1];
        bbox[2]=shift_x[2];
        bbox[3]=scale_x;
    }else{
        bbox[0] = 0;
        bbox[1] = 0;
        bbox[2] = 0;
        bbox[3] = 1.0/box_size_;

        // determine the tree depth
        Real_t leaf_size = box_size_;
        // r_near/2 < leaf_size <= r_near
        while(leaf_size>r_near && tree_depth<MAX_DEPTH-1){
            leaf_size *= 0.5;
            tree_depth++;
        }
        leaf_size_ = leaf_size;
        //COUT("leaf_size: "<<leaf_size);
    }
    pvfmm::Profile::Toc();
}

template<typename Real_t>
void BoundingBoxGrid<Real_t>::ConstructLocalTree(TREEGRID &BB_let)
{
    pvfmm::Profile::Tic("LocalTree",&comm,true);
    PVFMMVec_t           & box_min =BB_let.box_min;
    PVFMMVec_t           & box_max =BB_let.box_max;
    pvfmm::Vector<size_t>& pt_id   =BB_let.pt_id   ;
    pvfmm::Vector<size_t>& box_id  =BB_let.box_id   ;

    pvfmm::Vector<pvfmm::MortonId>& let_mid   =BB_let.mid;
    pvfmm::Vector<size_t>&          let_pt_cnt=BB_let.pt_cnt;
    pvfmm::Vector<size_t>&          let_pt_dsp=BB_let.pt_dsp;
    pvfmm::Vector<pvfmm::MortonId>& let_mins  =BB_let.mins;

    { // build scatter-indices (pt_id) and tree (let_mid, let_pt_cnt, let_pt_dsp)
        pvfmm::Vector<pvfmm::MortonId> pt_mid(N_pts_);
        { // build pt_mid
            Real_t scale_x, shift_x[DIM];
            { // set scale_x, shift_x
                shift_x[0]=bbox_[0];
                shift_x[1]=bbox_[1];
                shift_x[2]=bbox_[2];
                scale_x=bbox_[3];
            }

            #pragma omp parallel for
            for(size_t tid=0;tid<omp_p_;tid++){
                Real_t c[DIM];
                size_t a=((tid+0)*N_pts_)/omp_p_;
                size_t b=((tid+1)*N_pts_)/omp_p_;
                for(size_t i=a;i<b;i++){
                    for(size_t k=0;k<DIM;k++){
                        c[k]=BB_pts_[i*DIM+k]*scale_x+shift_x[k];
                        while(c[k]< 0.0) c[k]+=1.0;
                        while(c[k]>=1.0) c[k]-=1.0;
                    }
                    pt_mid[i]=pvfmm::MortonId(c,tree_depth_);
                }
            }
        }

        pt_id   .ReInit(N_pts_);
        pvfmm::Profile::Tic("SortPoints",&comm,true);
        pvfmm::par::SortScatterIndex(pt_mid, pt_id, comm);
        pvfmm::Profile::Toc();
        pvfmm::Profile::Tic("ScatterPtMid",&comm,true);
        pvfmm::par::ScatterForward  (pt_mid, pt_id, comm);
        pvfmm::Profile::Toc();
        { // build let_mins
            pvfmm::Profile::Tic("LetMins",&comm,true);
            let_mins.ReInit(np_);
            MPI_Allgather(&  pt_mid[0], 1, pvfmm::par::Mpi_datatype<pvfmm::MortonId>::value(),
                          &let_mins[0], 1, pvfmm::par::Mpi_datatype<pvfmm::MortonId>::value(), comm);
            if(rank_) assert(let_mins[rank_]!=let_mins[rank_-1]);
            let_mins[0]=pvfmm::MortonId(0,0,0,tree_depth_);
            pvfmm::Profile::Toc();
        }
        { // Exchange shared octant with neighbour
            pvfmm::Profile::Tic("NeighbourComm",&comm,true);
            int send_size=0;
            int recv_size=0;
            if(rank_<np_-1){ // send_size
                send_size=pt_mid.Dim()-(std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), let_mins[rank_+1])-&pt_mid[0]);
            }
            { // recv_size
                MPI_Status status;
                MPI_Sendrecv(&send_size,1,MPI_INT,(rank_<np_-1?rank_+1:0),0,
                             &recv_size,1,MPI_INT,(rank_>0?rank_-1:np_-1),0,comm,&status);
            }

            { // Set new pt_id
                pvfmm::Vector<size_t> pt_id_new(pt_id.Dim()+recv_size-send_size);
                memcpy(&pt_id_new[0]+recv_size, &pt_id[0], (pt_id.Dim()-send_size)*sizeof(size_t));

                MPI_Status status;
                MPI_Sendrecv(&pt_id[0]+pt_id.Dim()-send_size,send_size,pvfmm::par::Mpi_datatype<size_t>::value(),(rank_<np_-1?rank_+1:0),0,
                             &pt_id_new[0]                  ,recv_size,pvfmm::par::Mpi_datatype<size_t>::value(),(rank_>0?rank_-1:np_-1),0,comm,&status);
                pt_id.Swap(pt_id_new);
            }
            { // Set new pt_mid
                pvfmm::Vector<pvfmm::MortonId> pt_mid_new(pt_mid.Dim()+recv_size-send_size);
                memcpy(&pt_mid_new[0]+recv_size, &pt_mid[0], (pt_mid.Dim()-send_size)*sizeof(pvfmm::MortonId));
                for(size_t i=0;i<recv_size;i++) pt_mid_new[i]=let_mins[rank_];
                pt_mid.Swap(pt_mid_new);
            }
            pvfmm::Profile::Toc();
        }
        { // Sort points by pt_id in each octant
            #pragma omp parallel for
            for(size_t tid=0;tid<omp_p_;tid++){
                size_t a=(pt_mid.Dim()*(tid+0))/omp_p_;
                size_t b=(pt_mid.Dim()*(tid+1))/omp_p_;
                if(a>0           ) a=std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), pt_mid[a])-&pt_mid[0];
                if(b<pt_mid.Dim()) b=std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), pt_mid[b])-&pt_mid[0];
                for(size_t i=a;i<b;){
                    size_t j=i; while(j<b && pt_mid[j]==pt_mid[i]) j++;
                    std::sort(&pt_id[0]+i, &pt_id[0]+j);
                    i=j;
                }
            }
        }

        { // set let_mid
            std::vector<std::vector<pvfmm::MortonId> > mid_omp(omp_p_);
            #pragma omp parallel for
            for(size_t tid=0;tid<omp_p_;tid++){
                std::vector<pvfmm::MortonId>& mid=mid_omp[tid];
                size_t a=(pt_mid.Dim()*(tid+0))/omp_p_;
                size_t b=(pt_mid.Dim()*(tid+1))/omp_p_;
                if(a>0           ) a=std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), pt_mid[a])-&pt_mid[0];
                if(b<pt_mid.Dim()) b=std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), pt_mid[b])-&pt_mid[0];
                if(a<b) mid.push_back(pt_mid[a]);
                for(size_t i=a;i<b;i++){
                    if(mid.back()!=pt_mid[i]) mid.push_back(pt_mid[i]);
                }
            }
            { // Resize let_mid
                size_t size=0;
                for(size_t tid=0;tid<omp_p_;tid++){
                    size+=mid_omp[tid].size();
                }
                let_mid.ReInit(size);
            }
            #pragma omp parallel for
            for(size_t tid=0;tid<omp_p_;tid++){ // Set let_mid
                size_t offset=0;
                for(size_t i=0;i<tid;i++){
                    offset+=mid_omp[i].size();
                }

                std::vector<pvfmm::MortonId>& mid=mid_omp[tid];
                for(size_t i=0;i<mid.size();i++){
                    let_mid[i+offset]=mid[i];
                }
            }
        }
        { // set let_pt_dsp
            let_pt_dsp.ReInit(let_mid.Dim());
            #pragma omp parallel for
            for(size_t i=0;i<let_mid.Dim();i++){
                let_pt_dsp[i]=std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), let_mid[i])-&pt_mid[0];
            }
        }
        { // set let_pt_cnt
            let_pt_cnt.ReInit(let_mid.Dim());
            #pragma omp parallel for
            for(size_t i=1;i<let_mid.Dim();i++){
                let_pt_cnt[i-1]=let_pt_dsp[i]-let_pt_dsp[i-1];
            }
            if(let_mid.Dim()) let_pt_cnt[let_mid.Dim()-1]=pt_mid.Dim()-let_pt_dsp[let_mid.Dim()-1];
        }
    }

    { // scatter pt_coord, box_id, box_min, box_max
        pvfmm::Profile::Tic("ScatterRemain",&comm,true);
        box_min=BB_pts_min_;
        pvfmm::par::ScatterForward(box_min, pt_id, comm);
        box_max=BB_pts_max_;
        pvfmm::par::ScatterForward(box_max, pt_id, comm);
        box_id = BB_id_;
        pvfmm::par::ScatterForward(box_id, pt_id, comm);
        pvfmm::Profile::Toc();
    }
    pvfmm::Profile::Toc();
}

template<typename Real_t>
void BoundingBoxGrid<Real_t>::FindNearPair(std::vector< std::pair<size_t, size_t> > &BBIPairs)
{
    pvfmm::Profile::Tic("NearPair",&comm,true);
    PVFMMVec_t &box_min = BB_let_.box_min;
    PVFMMVec_t &box_max = BB_let_.box_max;
    pvfmm::Vector<size_t> &pt_id = BB_let_.pt_id;
    pvfmm::Vector<size_t> &box_id = BB_let_.box_id;
    pvfmm::Vector<pvfmm::MortonId> &tree_mid = BB_let_.mid;
    pvfmm::Vector<size_t> &tree_pt_cnt = BB_let_.pt_cnt;
    pvfmm::Vector<size_t> &tree_pt_dsp = BB_let_.pt_dsp;

    std::vector<std::vector<std::pair<size_t, size_t> > > near_pair_omp(omp_p_); // (box_id, box_id)
    #pragma omp parallel num_threads(omp_p_)
    {
        size_t tid=omp_get_thread_num();
        std::vector<std::pair<size_t, size_t> >& near_pair=near_pair_omp[tid];

        size_t FLOP=0;
        size_t a=((tid+0)*tree_mid.Dim())/omp_p_;
        size_t b=((tid+1)*tree_mid.Dim())/omp_p_;
        for(size_t i=a;i<b;i++){
            size_t tcnt=tree_pt_cnt[i];
            size_t tdsp=tree_pt_dsp[i];
            PVFMMVec_t tmin, tmax;
            pvfmm::Vector<size_t> tbox_id;
            if(tcnt){ // Set t_coord
                tmin.ReInit(tcnt*DIM,  &box_min[tdsp*DIM],false);
                tmax.ReInit(tcnt*DIM,  &box_max[tdsp*DIM],false);
                tbox_id.ReInit(tcnt,         &box_id[tdsp]           ,false);
            }

            { // Find near pairs
                for(size_t t1=0;t1<tcnt;t1++){
                    for(size_t t2=0;t2<tcnt;t2++){
                        if( (tbox_id[t1]!=tbox_id[t2]) &&
                            CheckBBCollision(&tmin[t1*DIM], &tmax[t1*DIM], &tmin[t2*DIM], &tmax[t2*DIM])
                          )
                        {
                            std::pair<size_t, size_t> new_pair;
                            new_pair.first =  tbox_id[t1];
                            new_pair.second = tbox_id[t2];
                            near_pair.push_back(new_pair);
                        }
                    }

                    std::sort(near_pair.begin(),near_pair.end());
                    near_pair.erase(std::unique(near_pair.begin(), near_pair.end()), near_pair.end());
                
                    // TODO: add FLOP
                }
                FLOP+=tcnt*tcnt*10;
            }
        }
        pvfmm::Profile::Add_FLOP(FLOP);
    }

    size_t near_size=0;
    pvfmm::Vector<size_t> near_cnt(omp_p_);
    pvfmm::Vector<size_t> near_dsp(omp_p_); near_dsp[0]=0;
    for(size_t tid=0;tid<omp_p_;tid++){
        if(tid)
            near_dsp[tid]=near_pair_omp[tid-1].size()+near_dsp[tid-1];
        near_cnt[tid]=near_pair_omp[tid  ].size();
        near_size   +=near_pair_omp[tid  ].size();
    }

    pvfmm::Vector<size_t> near_trg_box_id;
    pvfmm::Vector<size_t> near_src_box_id;
    near_trg_box_id.ReInit(near_size);
    near_src_box_id.ReInit(near_size);

    #pragma omp parallel for
    for(size_t tid=0;tid<omp_p_;tid++){
        size_t dsp=near_dsp[tid];
        size_t cnt=near_cnt[tid];
        std::vector<std::pair<size_t, size_t> >& near_pair=near_pair_omp[tid];
        for(size_t i=0;i<cnt;i++){
            near_trg_box_id[(dsp+i)]=near_pair[i].first;
            near_src_box_id[(dsp+i)]=near_pair[i].second;
        }
    }

    // scatter to box_ids
    size_t box_id_offset;
    {
        long long disp = 0;
        long long size = N_bbox_;
        MPI_Scan(&size, &disp, 1, MPI_LONG_LONG, MPI_SUM, comm);
        box_id_offset = disp - size;
    }
    pvfmm::Vector<size_t> scatter_idx;
    pvfmm::par::SortScatterIndex(near_src_box_id, scatter_idx, comm, &box_id_offset);
    pvfmm::par::ScatterForward(near_src_box_id, scatter_idx, comm);
    pvfmm::par::ScatterForward(near_trg_box_id, scatter_idx, comm);

    // construct unique pairs
    near_size = scatter_idx.Dim();
    BBIPairs.resize(near_size);
    #pragma omp parallel for
    for(size_t i=0; i<near_size;i++){
        BBIPairs[i].first = near_src_box_id[i];
        BBIPairs[i].second = near_trg_box_id[i];
    }
    std::sort(BBIPairs.begin(), BBIPairs.end());
    BBIPairs.erase(std::unique(BBIPairs.begin(), BBIPairs.end()), BBIPairs.end());

    pvfmm::Profile::Toc();
}

// TODO: assumption! targets are all in the global box defined by bbox_
template<typename Real_t>
void BoundingBoxGrid<Real_t>::FindTargetNearPair(pvfmm::Vector<size_t> &near_box_id, pvfmm::Vector<size_t> &near_trg_id, 
        PVFMMVec_t &near_trg_coord, pvfmm::Vector<size_t> &near_trg_scatter)
//void BoundingBoxGrid<Real_t>::FindTargetNearPair(std::vector< std::pair<size_t, size_t> > &BBIPairs)
{
    pvfmm::Vector<size_t> pt_id;
    PVFMMVec_t pt_coord = T_;
    size_t trg_id_offset;
    // get trg_id_offset
    {
        long long disp = 0;
        long long size = pt_coord.Dim()/DIM;
        MPI_Scan(&size, &disp, 1, MPI_LONG_LONG, MPI_SUM, comm);
        trg_id_offset = disp-size;
    }

    pvfmm::Vector<pvfmm::MortonId> tree_mid;
    pvfmm::Vector<size_t> tree_pt_cnt;
    pvfmm::Vector<size_t> tree_pt_dsp;
    // construct target point tree
    {
        pvfmm::Profile::Tic("BoxTrgTree", &comm, true);
        pvfmm::Vector<pvfmm::MortonId> pt_mid(pt_coord.Dim()/DIM);

        // set pt_mid
        {
            Real_t scale_x, shift_x[DIM];
            { // set scale_x, shift_x
                shift_x[0]=bbox_[0];
                shift_x[1]=bbox_[1];
                shift_x[2]=bbox_[2];
                scale_x=bbox_[3];
            }
            assert(BB_let_.mid.Dim());
            size_t tree_depth=BB_let_.mid[0].GetDepth();
            #pragma omp parallel for
            for(size_t tid=0;tid<omp_p_;tid++){
                Real_t c[DIM];
                size_t a=((tid+0)*pt_coord.Dim()/DIM)/omp_p_;
                size_t b=((tid+1)*pt_coord.Dim()/DIM)/omp_p_;
                for(size_t i=a;i<b;i++){
                    for(size_t k=0;k<DIM;k++){
                        c[k]=pt_coord[i*DIM+k]*scale_x+shift_x[k];
                        while(c[k]< 0.0) c[k]+=1.0;
                        while(c[k]>=1.0) c[k]-=1.0;
                    }
                    pt_mid[i]=pvfmm::MortonId(c,tree_depth);
                }
            }
        }

        // sort pt data (pt_mid, pt_coord, pt_id);
        {
            pvfmm::par::SortScatterIndex(pt_mid, pt_id, comm, &BB_let_.mins[rank_]);
            pvfmm::par::ScatterForward(pt_mid, pt_id, comm);
            pvfmm::par::ScatterForward(pt_coord, pt_id, comm);
        }

        // build tree_mid, tree_pt_cnt and tree_pt_dsp
        {
            std::vector<std::vector<pvfmm::MortonId> > mid_(omp_p_);
            #pragma omp parallel for
            for(size_t tid=0;tid<omp_p_;tid++){
                std::vector<pvfmm::MortonId>& mid=mid_[tid];
                size_t a=(pt_mid.Dim()*(tid+0))/omp_p_;
                size_t b=(pt_mid.Dim()*(tid+1))/omp_p_;
                if(a>0           ) a=std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), pt_mid[a])-&pt_mid[0];
                if(b<pt_mid.Dim()) b=std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), pt_mid[b])-&pt_mid[0];
                if(a<b) mid.push_back(pt_mid[a]);
                for(size_t i=a;i<b;i++){
                    if(mid.back()!=pt_mid[i]) mid.push_back(pt_mid[i]);
                }
            }
            { // Resize tree_mid
                size_t size=0;
                for(size_t tid=0;tid<omp_p_;tid++){
                    size+=mid_[tid].size();
                }
                tree_mid.ReInit(size);
            }
            #pragma omp parallel for
            for(size_t tid=0;tid<omp_p_;tid++){ // Set tree_mid
                size_t offset=0;
                for(size_t i=0;i<tid;i++){
                    offset+=mid_[i].size();
                }

                std::vector<pvfmm::MortonId>& mid=mid_[tid];
                for(size_t i=0;i<mid.size();i++){
                    tree_mid[i+offset]=mid[i];
                }
            }

            // Resize tree_pt_cnt, tree_pt_dsp
            tree_pt_cnt.ReInit(tree_mid.Dim());
            tree_pt_dsp.ReInit(tree_mid.Dim());

            #pragma omp parallel for
            for(size_t i=0;i<tree_mid.Dim();i++){
                tree_pt_dsp[i]=std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), tree_mid[i])-&pt_mid[0];
            }

            #pragma omp parallel for
            for(size_t i=1;i<tree_mid.Dim();i++){
                tree_pt_cnt[i-1]=tree_pt_dsp[i]-tree_pt_dsp[i-1];
            }
            if(tree_mid.Dim())
                tree_pt_cnt[tree_mid.Dim()-1]=pt_mid.Dim()-tree_pt_dsp[tree_mid.Dim()-1];
        }

        pvfmm::Profile::Toc();
    }

    // find near box_id - loc_trg_id pairs
    {
        pvfmm::Profile::Tic("PTPair", &comm, true);
        std::vector< std::vector<std::pair<size_t, size_t> > > near_pair_(omp_p_);
        #pragma omp parallel num_threads(omp_p_)
        {
            size_t tid = omp_get_thread_num();
            std::vector< std::pair<size_t, size_t> >& near_pair = near_pair_[tid];

            size_t FLOP = 0;
            size_t a = ((tid+0)*tree_mid.Dim())/omp_p_;
            size_t b = ((tid+1)*tree_mid.Dim())/omp_p_;
            for(size_t i=a; i<b; i++)
            {
                PVFMMVec_t s_box_min;
                PVFMMVec_t s_box_max;
                pvfmm::Vector<size_t> s_box_id;
                pvfmm::Vector<pvfmm::MortonId>& s_mids = BB_let_.mid;
                pvfmm::MortonId t_mid = tree_mid[i];
                size_t tcnt = tree_pt_cnt[i];
                size_t tdsp = tree_pt_dsp[i];
                PVFMMVec_t tcoord;
                if(tcnt)
                {
                    tcoord.ReInit(tcnt*DIM, &pt_coord[tdsp*DIM], false);
                }
                int k = std::lower_bound(&s_mids[0], &s_mids[0]+s_mids.Dim(), t_mid) - &s_mids[0];
                if(k<s_mids.Dim() && s_mids[k]==t_mid)
                {
                    s_box_min.ReInit(BB_let_.pt_cnt[k]*DIM, &BB_let_.box_min[0]+BB_let_.pt_dsp[k]*DIM, false);
                    s_box_max.ReInit(BB_let_.pt_cnt[k]*DIM, &BB_let_.box_max[0]+BB_let_.pt_dsp[k]*DIM, false);
                    s_box_id.ReInit( BB_let_.pt_cnt[k]          , &BB_let_.box_id[0] +BB_let_.pt_dsp[k]          , false);
                }
                
                for(size_t t=0; t<tcnt; t++)
                {
                    for(size_t s=0; s<s_box_id.Dim(); s++)
                    {
                        Real_t t_x = tcoord[t*DIM + 0];
                        Real_t t_y = tcoord[t*DIM + 1];
                        Real_t t_z = tcoord[t*DIM + 2];
                        Real_t min_x = s_box_min[s*DIM+0];
                        Real_t min_y = s_box_min[s*DIM+1];
                        Real_t min_z = s_box_min[s*DIM+2];
                        Real_t max_x = s_box_max[s*DIM+0];
                        Real_t max_y = s_box_max[s*DIM+1];
                        Real_t max_z = s_box_max[s*DIM+2];
                        if( (t_x>=min_x && t_x<=max_x) && (t_y>=min_y && t_y<=max_y) && (t_z>=min_z && t_z<=max_z) )
                        {
                            std::pair<size_t, size_t> new_pair;
                            new_pair.first = s_box_id[s];
                            new_pair.second = tdsp+t;
                            near_pair.push_back(new_pair);
                        }
                    }
                }
            }

            std::sort(near_pair.begin(),near_pair.end());
            near_pair.erase(std::unique(near_pair.begin(), near_pair.end()), near_pair.end());
        }
        
        size_t near_size = 0;
        pvfmm::Vector<size_t> near_cnt(omp_p_);
        pvfmm::Vector<size_t> near_dsp(omp_p_); near_dsp[0]=0;
        for(size_t tid=0; tid<omp_p_; tid++){
            if(tid)
                near_dsp[tid] = near_pair_[tid-1].size()+near_dsp[tid-1];
            near_cnt[tid] = near_pair_[tid].size();
            near_size    += near_pair_[tid].size();
        }

        //pvfmm::Vector<size_t> near_box_id;
        //pvfmm::Vector<size_t> near_trg_id;
        //PVFMMVec_t near_trg_coord;
        near_box_id.ReInit(near_size);
        near_trg_id.ReInit(near_size);
        near_trg_coord.ReInit(near_size*DIM);

        #pragma omp parallel for
        for(size_t tid=0; tid<omp_p_; tid++)
        {
            size_t dsp = near_dsp[tid];
            size_t cnt = near_cnt[tid];
            std::vector< std::pair<size_t, size_t> >& near_pair = near_pair_[tid];
            for(size_t i=0; i<cnt; i++)
            {
                size_t loc_trg_id = near_pair[i].second;
                near_box_id[(dsp+i)] = near_pair[i].first;
                near_trg_id[(dsp+i)] = pt_id[loc_trg_id];
                near_trg_coord[(dsp+i)*DIM+0] = pt_coord[loc_trg_id*DIM+0];
                near_trg_coord[(dsp+i)*DIM+1] = pt_coord[loc_trg_id*DIM+1];
                near_trg_coord[(dsp+i)*DIM+2] = pt_coord[loc_trg_id*DIM+2];
            }
        }

        //scatter target
        size_t box_id_offset;
        {
            long long disp = 0;
            long long size = N_bbox_;
            MPI_Scan(&size, &disp, 1, MPI_LONG_LONG, MPI_SUM, comm);
            box_id_offset = disp - size;
        }
        
        pvfmm::Vector<size_t> scatter_idx;
        pvfmm::par::SortScatterIndex(near_box_id, scatter_idx, comm, &box_id_offset);
        pvfmm::par::ScatterForward(near_box_id, scatter_idx, comm);
        pvfmm::par::ScatterForward(near_trg_id, scatter_idx, comm);
        pvfmm::par::ScatterForward(near_trg_coord, scatter_idx, comm);

        // build trg scatter to scatter back trg to it's resident processor
        //pvfmm::Vector<size_t> near_trg_scatter;
        {
            pvfmm::par::SortScatterIndex(near_trg_id, near_trg_scatter, comm, &trg_id_offset);
            //pvfmm::par::ScatterForward(near_trg_id, near_trg_scatter, comm);
        }

        /*
        near_size = scatter_idx.Dim();
        BBIPairs.resize(near_size);
        #pragma omp parallel for
        for(size_t i=0; i<near_size; i++)
        {
            BBIPairs[i].first = near_box_id[i];
            BBIPairs[i].second = near_trg_id[i];
        }

        // could use __gnu_parallel::sort and __gnu_parallel::unique_copy
        std::sort(BBIPairs.begin(), BBIPairs.end());
        BBIPairs.erase(std::unique(BBIPairs.begin(), BBIPairs.end()), BBIPairs.end());
        */

        pvfmm::Profile::Toc();
    }

}

template<typename Real_t>
void BoundingBoxGrid<Real_t>::GlobalBoundingBox(Real_t *scale_xr, Real_t *shift_xr)
{
    Real_t& scale_x = *scale_xr;
    Real_t* shift_x = shift_xr;

    double loc_min_x[DIM];
    double loc_max_x[DIM];
    
    assert(N_bbox_>0);
    size_t n_src = N_bbox_;

    for(size_t k=0;k<DIM;k++){
        loc_min_x[k] = BB_min_[k];
        loc_max_x[k] = BB_max_[k];
    }
    
    for(size_t i=0;i<n_src;i++){
        const Real_t* x_min_=&BB_min_[i*DIM];
        const Real_t* x_max_=&BB_max_[i*DIM];
        for(size_t k=0;k<DIM;k++){
            if(loc_min_x[k]>x_min_[0]) loc_min_x[k]=x_min_[0];
            if(loc_max_x[k]<x_max_[0]) loc_max_x[k]=x_max_[0];
            ++x_min_; ++x_max_;
        }
    }

    double min_x[DIM];
    double max_x[DIM];
    MPI_Allreduce(loc_min_x, min_x, DIM, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(loc_max_x, max_x, DIM, MPI_DOUBLE, MPI_MAX, comm);

    static Real_t eps=machine_eps<Real_t>()*64; // Points should be well within the box.
    scale_x=1.0/(max_x[0]-min_x[0]+2*eps);
    for(size_t k=0;k<DIM;k++){
        scale_x=std::min(scale_x,(Real_t)(1.0/(max_x[k]-min_x[k]+2*eps)));
    }
    if(scale_x*0.0!=0.0) scale_x=1.0; // fix for scal_x=inf
    for(size_t k=0;k<DIM;k++){
        shift_x[k]=-min_x[k]*scale_x+eps;
    }
}
    
template<typename Real_t>
void BoundingBoxGrid<Real_t>::GenerateBBPoints()
{
    pvfmm::Profile::Tic("GenBBPoints", &comm, true);
    // number of points in each dimension
    std::vector<size_t> bbox_nxyz(DIM*N_bbox_, 0);
    // number of points per box
    std::vector<size_t> bbox_n(N_bbox_, 0);
    // points displacement
    std::vector<size_t> bbox_ndsp(N_bbox_, 0);
    size_t n_sum = 0;

    // set number of pts in each dimension
    // set number of points per box
    #pragma omp parallel for
    for(size_t tid=0; tid<omp_p_; tid++){
        size_t a = ((tid+0)*N_bbox_)/omp_p_;
        size_t b = ((tid+1)*N_bbox_)/omp_p_;

        for(size_t i=a; i<b; i++){
            const Real_t* min_i = &BB_min_[i*DIM];
            const Real_t* max_i = &BB_max_[i*DIM];
            size_t* bbox_nxyzi     = &bbox_nxyz[i*DIM];

            int nx = std::ceil((max_i[0] - min_i[0])/leaf_size_) + 1;
            int ny = std::ceil((max_i[1] - min_i[1])/leaf_size_) + 1;
            int nz = std::ceil((max_i[2] - min_i[2])/leaf_size_) + 1;
            
            ASSERT(nx > 1, "invalid nx");ASSERT(ny > 1, "invalid ny");ASSERT(nz > 1, "invalid nz");
            
            bbox_nxyzi[0] = nx;
            bbox_nxyzi[1] = ny;
            bbox_nxyzi[2] = nz;
            bbox_n[i] = nx*ny*nz;
        }
    }
        
    // set point displacement for each bounding box
    bbox_ndsp[0]=0; pvfmm::omp_par::scan(&bbox_n[0], &bbox_ndsp[0], bbox_n.size());
    // total number of points
    n_sum = pvfmm::omp_par::reduce(&bbox_n[0], bbox_n.size());
    N_pts_ = n_sum;

    // init data for generating points
    BB_pts_.ReInit(DIM*n_sum);
    BB_pts_min_.ReInit(DIM*n_sum);
    BB_pts_max_.ReInit(DIM*n_sum);
    BB_id_.ReInit(n_sum);

    // set global bounding box id offset
    size_t box_id_offset;
    {
        long long disp = 0;
        long long size = N_bbox_;
        MPI_Scan(&size, &disp, 1, MPI_LONG_LONG, MPI_SUM, comm);
        box_id_offset = disp - size;
    }

    // generating points
    #pragma omp parallel for
    for(size_t tid=0; tid<omp_p_; tid++){
        size_t a = ((tid+0)*N_bbox_)/omp_p_;
        size_t b = ((tid+1)*N_bbox_)/omp_p_;

        for(size_t i=a; i<b; i++){
            const Real_t* min_i = &BB_min_[i*DIM];
            const Real_t* max_i = &BB_max_[i*DIM];

            size_t box_id = box_id_offset + i;
            size_t nx = bbox_nxyz[i*DIM + 0];
            size_t ny = bbox_nxyz[i*DIM + 1];
            size_t nz = bbox_nxyz[i*DIM + 2];

            size_t ni = 0;
            size_t disp = bbox_ndsp[i];
            for(size_t nxi=0; nxi<nx; nxi++)
            for(size_t nyi=0; nyi<ny; nyi++)
            for(size_t nzi=0; nzi<nz; nzi++)
            {
                Real_t* BB_pts = &BB_pts_[DIM*(disp+ni)];
                BB_pts[0] = nxi*((max_i[0]-min_i[0])/(nx-1)) + min_i[0];
                BB_pts[1] = nyi*((max_i[1]-min_i[1])/(ny-1)) + min_i[1];
                BB_pts[2] = nzi*((max_i[2]-min_i[2])/(nz-1)) + min_i[2];

                Real_t* BB_pts_min = &BB_pts_min_[DIM*(disp+ni)];
                BB_pts_min[0] = min_i[0]; BB_pts_min[1] = min_i[1]; BB_pts_min[2] = min_i[2];
                Real_t* BB_pts_max = &BB_pts_max_[DIM*(disp+ni)];
                BB_pts_max[0] = max_i[0]; BB_pts_max[1] = max_i[1]; BB_pts_max[2] = max_i[2];
                
                BB_id_[disp+ni] = box_id;
                ni++;
            }
        }
    }
    
    pvfmm::Profile::Toc();
}

template<typename Real_t>
bool BoundingBoxGrid<Real_t>::CheckBBCollision(const Real_t *minA, const Real_t *maxA, const Real_t *minB, const Real_t *maxB)
{
    if(box_size_<=0)
    {
        return ( (minA[0]<=maxB[0] && maxA[0]>= minB[0]) &&
                 (minA[1]<=maxB[1] && maxA[1]>= minB[1]) &&
                 (minA[2]<=maxB[2] && maxA[2]>= minB[2])
               );
    }
    else
    {
        bool flag[DIM];
        for(size_t i=0; i<DIM; i++)
        {
            flag[i] = false;
            if( ((maxA[i]-minA[i]) + (maxB[i]-minB[i])) >= box_size_ )
            {
                flag[i] = true;
                continue;
            }

            double disp = box_size_ * round( ((maxA[i]+minA[i]) - (maxB[i]+minB[i]))/(2*box_size_) );
                
            ASSERT(fabs((maxA[i]+minA[i])/2 - disp - (maxB[i]+minB[i])/2)<=box_size_*0.5,"wrong checkBBcollision");
            flag[i] = ( (minA[i]-disp)<=maxB[i] && (maxA[i]-disp)>= minB[i] );
        }
        return (flag[0] && flag[1] && flag[2]);
    }
}

template class BoundingBoxGrid<double>;

END_EBI_NAMESPACE
