#include "pvfmm_interface.h"
#include "pvfmm.hpp"
#include "pvfmm_stks_kernel.tpp"
#include "pvfmm_laplace_kernel.cpp"
#include "pvfmm_mod_helmholtz_kernel.cpp"
#include "pvfmm_navier_kernel.cpp"
//#include "elastic_kernels.h"

#include <cassert>
#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <vector>
#include <time.h>
#define COUT(str) (std::cout<<str<<std::endl)
#define CERR(str,action) (std::cerr<<"[ERROR]["<< __FUNCTION__ <<"] "<<str<<std::endl, action)
#define ASSERT(expr, msg) ( (expr) ? assert(true) : CERR(msg,abort()))
#define COUTDEBUG(str) (std::cout<<"[DEBUG] "<<str<<std::endl)

double POISS;
double SHEAR;

const static int DIM=3;

//const pvfmm::Kernel<double> ker_mod_helmholtz_sl;
//const pvfmm::Kernel<double> ker_mod_helmholtz_dl;
//const pvfmm::Kernel<double> ker_stokes_dl;
typedef std::vector<real_t> vec;
typedef pvfmm::FMM_Node<pvfmm::MPI_Node<real_t> > Node_t;

class PVFMMContext{
  public:
  typedef pvfmm::FMM_Pts<Node_t> Mat_t;
  typedef pvfmm::FMM_Tree<Mat_t> Tree_t;

  PVFMMContext() : tree(NULL), mat(NULL) {};

  ~PVFMMContext(){
    COUTDEBUG("destructing context object");
    if (this->tree != NULL){
        delete this->tree;
    }
    if (this->mat != NULL){
        delete this->mat;
    }
    //delete this->kifmm_to_pvfmm_ids;
    //delete this->pvfmm_to_kifmm_ids;
    COUTDEBUG("deleted all pointers");
  };

  int mult_order;
  int max_pts;
  int max_depth;
  int sl;
  int dl;
  int periodic;
  int source_dof;
  int target_dof;
  int dense_eval;
  MPI_Comm comm;

  pvfmm::BoundaryType boundary;
  const pvfmm::Kernel<real_t>* ker;
  std::map<int, pvfmm::MortonId>* kifmm_to_pvfmm_ids;
  std::map<pvfmm::MortonId, int>* pvfmm_to_kifmm_ids;

  Tree_t* tree;
  Mat_t* mat;
};


void make_pvfmm_context(const int &mult_order, const int &max_pts,
			const int &max_depth, const int &sl,  const int &dl,
			const int &periodic, const Kernel &kernel, const int &sdof,
            const int &tdof, const int &dense_eval,  void **context, 
            const KernelOptions* kernel_opts)
{
  COUTDEBUG("Initializing a pvfmm context object"
	    <<" mult_order="<<mult_order
	    <<", max_pts="<<max_pts
	    <<", max_depth="<<max_depth
	    <<", sl="<<sl
	    <<", dl="<<dl
	    <<", periodic="<<periodic
	    <<", sdof="<<sdof
	    <<", tdof="<<tdof
	    <<", dense_eval="<<dense_eval
	    );

  pvfmm::BoundaryType bndry(periodic ? pvfmm::Periodic : pvfmm::FreeSpace);

  // Create new context.
  PVFMMContext *ctx = new PVFMMContext;

  ctx->mult_order = mult_order;
  ctx->max_pts	  = max_pts;
  ctx->max_depth  = max_depth;
  ctx->periodic	  = periodic;
  ctx->comm	      = MPI_COMM_WORLD;
  ctx->source_dof = sdof;
  ctx->target_dof = tdof;
  ctx->sl         = sl;
  ctx->dl         = dl;
  ctx->dense_eval = dense_eval;


  ctx->boundary   = bndry;
  switch( kernel ){
      case KNL_STK_S_U:
          ctx->ker = &pvfmm::StokesKernel<double>::velocity();
          break;
      case KNL_STK_S_P:
          ctx->ker = &pvfmm::StokesKernel<double>::pressure();
          break;
      case KNL_STK_D_U:
          ctx->ker = &ker_stokes_dl;
          break;
      case KNL_STK_D_P:
          //ctx->ker = &pvfmm::StokesKernel<double>::pressure();
          ctx->ker = &ker_stokes_pressure_dl;
          break;
      case KNL_LAP_D_U:
          //ctx->ker = &pvfmm::LaplaceKernel<double>::potential();
          ctx->ker = &ker_laplace_dl;
          break;
      case KNL_LAP_S_U:
          //ctx->ker = &pvfmm::LaplaceKernel<double>::potential();
          ctx->ker = &ker_laplace_sl;
          break;
          break;
      case KNL_MODHEL_D_U:
          assert(kernel_opts != NULL);

          helmholtz_frequency(kernel_opts->helmholtz_frequency);
          ctx->ker = &ker_mod_helmholtz_dl;
          const_cast<std::string&> (ctx->ker->ker_name) = helmholtz_name(false); 
          break;
      case KNL_MODHEL_S_U:
          assert(kernel_opts != NULL);

          helmholtz_frequency(kernel_opts->helmholtz_frequency);
          ctx->ker = &ker_mod_helmholtz_sl;
          const_cast<std::string&> (ctx->ker->ker_name) = helmholtz_name(true); 

          break;
      case KNL_NAV_D_U:
          assert(kernel_opts != NULL);

          navier_mu(kernel_opts->navier_mu);
          navier_nu(kernel_opts->navier_nu);
          ctx->ker = &ker_navier_dl;
          POISS= navier_nu();
          SHEAR= navier_mu();
          //ctx->ker = &ElasticKernel<double>::Disp();
          const_cast<std::string&> (ctx->ker->ker_name) = navier_name(false); 
          //ctx->ker->ker_poten = navier_sl;
          break;
      case KNL_NAV_S_U:
          assert(kernel_opts != NULL);
          navier_mu(kernel_opts->navier_mu);
          navier_nu(kernel_opts->navier_nu);
          POISS= navier_nu();
          SHEAR= navier_mu();
          //ctx->ker = &ElasticKernel<double>::Disp();
          ctx->ker = &ker_navier_dl;
          const_cast<std::string&> (ctx->ker->ker_name)= navier_name(false); 
          break;
      default:
        ASSERT(false, "KernelNotImplementedError");
  }
  //pvfmm::mem::MemoryManager mem_mgr(10000000);
  //ctx->mat        = new PVFMMContext::Mat_t(&mem_mgr);
  if(!ctx->dense_eval){
      ctx->mat        = new PVFMMContext::Mat_t();
      
      ctx->mat->Initialize(ctx->mult_order, ctx->comm, ctx->ker);
  }
  COUTDEBUG("Finished");

  *context = (void*) ctx;
}

void clear_pvfmm_context(void **context){

  if(*context == NULL) return;
  COUTDEBUG("deleting context");
  PVFMMContext *ctx = (PVFMMContext*) *context;
  delete ctx;
  *context = NULL;
} 
/* ##########################################################################*/
void laplace_sl_pvfmm(
    const size_t &nsrc, const real_t *src, const real_t *den,
    const size_t &ntrg, const real_t *trg, real_t *pot,
    const int &rebuild_tree, void **context)
{
  ASSERT(*context != NULL, "was called with null context. "
	 "Generate a context by calling make_pvfmm_context.");

  // context object
  PVFMMContext *ctx = (PVFMMContext*) *context;
  COUTDEBUG("pvfmm_laplace_sl with nsrc="<<nsrc<<", ntrg="<<ntrg
	    <<", rebuild_tree="<<rebuild_tree << ", sl=" << ctx->sl << ", dl=" <<ctx->dl);

  //ASSERT(ctx->sl, "context doesn't support single layer");
  //ASSERT(!ctx->dl, "context is build for combined sl and dl kernel");

  COUTDEBUG("Using context with"
	    <<" mult_order="<<ctx->mult_order
	    <<", max_pts="<<ctx->max_pts
	    <<", max_depth="<<ctx->max_depth
	    <<", sl="<<ctx->sl
	    <<", dl="<<ctx->dl
	    <<", periodic="<<ctx->periodic
	    <<", source_dof="<<ctx->source_dof
	    <<", target_dof="<<ctx->target_dof
	    );
  int source_dof = ctx->source_dof;
  int target_dof = ctx->target_dof;
  // copy data
  COUTDEBUG("DIM*nsrc: " << DIM*nsrc);

  vec src_coord(src, src + DIM*nsrc);
  vec src_value(den, den + source_dof*nsrc);
  vec trg_coord(trg, trg + DIM*ntrg);
  vec potv(target_dof*ntrg);

  if(rebuild_tree || ctx->tree == NULL){

    COUTDEBUG("(Re)building pvfmm tree");
    ASSERT(ctx->mat!=NULL, "badly initialized context");

    // FMM Setup
    delete ctx->tree;

    pvfmm::PtFMM_Tree<double>* tree = pvfmm::PtFMM_CreateTree<double>(src_coord, src_value,
           trg_coord, ctx->comm, ctx->max_pts, ctx->boundary);
    ctx->tree = tree;
    ctx->tree->SetupFMM(ctx->mat);
    // Run FMM
    COUT("Evaluating point FMM");
    pvfmm::PtFMM_Evaluate<double>(ctx->tree,potv,ntrg);
  } else {
    // only Run FMM
    ctx->tree->ClearFMMData();
    COUT("Evaluating point FMM");
    pvfmm::PtFMM_Evaluate<double>(ctx->tree,potv,ntrg);
  }
  /*
  // rebuild tree
  if(rebuild_tree || ctx->tree == NULL){

    COUTDEBUG("(Re)building pvfmm tree");
    ASSERT(ctx->mat!=NULL, "badly initialized context");

    // FMM Setup
    delete ctx->tree;
    ctx->tree = PtFMM_CreateTree(srcv, denv, trgv,
                     ctx->comm, ctx->max_pts, ctx->boundary);
    ctx->tree->SetupFMM(ctx->mat);

    // Run FMM
    COUT("Evaluating point FMM");
    PtFMM_Evaluate(ctx->tree,potv,ntrg);
  } else {
    // only Run FMM
    ctx->tree->ClearFMMData();
    COUT("Evaluating point FMM");
    PtFMM_Evaluate(ctx->tree,potv,ntrg,&denv);
  }
  */

  // copy to pot
  for(size_t i=0;i<target_dof*ntrg;++i){
    pot[i]=potv[i];
    //COUTDEBUG("Potential at " << i << "is : " << pot[i]);
    if (std::isnan(pot[i])){
      CERR("fmm potential is NaN",abort());
    }
  }
}


/* ##########################################################################*/






void laplace_sldl_pvfmm(
		       const size_t &sl_nsrc, const real_t *sl_src, const real_t *sl_den,
		       const size_t &dl_nsrc, const real_t *dl_src, const real_t *dl_den_nor,
		       const size_t &ntrg, const real_t *trg, real_t *pot,
		       const int &rebuild_tree, void **context)
{
  ASSERT(*context != NULL, "was called with null context. "
	 "Generate a context by calling make_pvfmm_context.");

  // context object
  PVFMMContext *ctx = (PVFMMContext*) *context;
  COUTDEBUG("pvfmm_laplace_sldl with sl_nsrc="<<sl_nsrc<<", dl_nsrc="<<dl_nsrc
	    <<", ntrg="<<ntrg<<", rebuild_tree="<<rebuild_tree);


  COUTDEBUG("Using context with"
	    <<" mult_order="<<ctx->mult_order
	    <<", max_pts="<<ctx->max_pts
	    <<", max_depth="<<ctx->max_depth
	    <<", sl="<<ctx->sl
	    <<", dl="<<ctx->dl
	    <<", periodic="<<ctx->periodic
	    );

  int source_dof = ctx->source_dof;
  int target_dof = ctx->target_dof;
  // copy data
  vec sl_srcv(sl_src,sl_src + DIM*sl_nsrc);
  vec sl_denv(sl_den,sl_den + source_dof*sl_nsrc);

  vec dl_srcv(dl_src    ,dl_src    +  DIM*dl_nsrc);
  vec dl_denv(dl_den_nor,dl_den_nor + (DIM + source_dof)*dl_nsrc);

  vec trgv(trg,trg + DIM*ntrg);
  vec potv(target_dof*ntrg);
  vec empty;
    std::cout << "dense_eval: " << ctx->dense_eval << ", " << 
     "rebuild_tree: " << rebuild_tree << ", " << int(ctx->tree == NULL) <<std::endl;
  if(ctx->dense_eval){
    int omp_p=omp_get_max_threads();
    //#pragma omp parallel for
    for(int i=0;i<omp_p;i++){
      size_t a=( i   *ntrg)/omp_p;
      size_t b=((i+1)*ntrg)/omp_p;

      if( ctx->sl ){
          /*ctx->ker->ker_poten(sl_srcv.data(), sl_nsrc,
                  sl_denv.data(), 1,
                  trgv.data(), ntrg,
                  potv.data()+ctx->ker->ker_dim[1],NULL);*/
          ctx->ker->ker_poten(sl_srcv.data(), sl_nsrc,
                  sl_denv.data(), 1,
                  trgv.data() + a*3, b-a,
                  potv.data()+a*ctx->ker->ker_dim[1],NULL);
      }else { 
          /*
          ctx->ker->dbl_layer_poten(dl_srcv.data(), dl_nsrc,
                  dl_denv.data(), 1,
                  trgv.data(), ntrg,
                  potv.data(),NULL);
                      */
          ctx->ker->dbl_layer_poten(dl_srcv.data(), dl_nsrc,
                  dl_denv.data(), 1,
                  trgv.data() + a*3, b-a,
                  potv.data()+a*ctx->ker->ker_dim[1],NULL);
      }
    }/*
      ctx->ker->dbl_layer_poten(dl_srcv.data(), dl_nsrc,
              dl_denv.data(), 1,
              trgv.data(), ntrg,
              potv.data(),NULL);
        */
  } else{
      // rebuild tree
      if(rebuild_tree || ctx->tree == NULL){
          COUTDEBUG("(Re)building sldl pvfmm tree");
          ASSERT(ctx->mat != NULL, "badly initialized context");

          // FMM Setup
          delete ctx->tree;
          clock_t t = clock();
          if( ctx->sl ){

              ctx->tree = pvfmm::PtFMM_CreateTree<double>(sl_srcv, sl_denv, empty, empty,
                      trgv,ctx->comm, ctx->max_pts, ctx->boundary);
          } else { //Double layer
              ctx->tree = pvfmm::PtFMM_CreateTree<double>( empty, empty, dl_srcv, dl_denv,
                      trgv,ctx->comm, ctx->max_pts, ctx->boundary);
          }
          t = clock() - t;
          COUTDEBUG("Tree Setup time:" << ((float) t/ CLOCKS_PER_SEC));
          COUTDEBUG("Number of boxes in PvFMM tree: " << ctx->tree->GetNodeList().size());
          t = clock();
          ctx->tree->SetupFMM(ctx->mat);
          t = clock() - t;
          COUTDEBUG("Matrix Setup time:" << ((float) t/ CLOCKS_PER_SEC));

          // Run FMM
          COUT("Evaluating point FMM");
          t = clock();
          pvfmm::PtFMM_Evaluate<double>(ctx->tree, potv, ntrg);
          t = clock() - t;
          COUTDEBUG("FMM eval time:" << ((float) t/ CLOCKS_PER_SEC));

      } else {
          // only Run FMM
          ctx->tree->ClearFMMData();
          if (ctx->sl){
              pvfmm::PtFMM_Evaluate<double>(ctx->tree, potv, ntrg, &sl_denv, NULL);
          } else { 
              pvfmm::PtFMM_Evaluate<double>(ctx->tree, potv, ntrg, NULL, &dl_denv);
          }
      }
  }

int omp_p=omp_get_max_threads();
#pragma omp parallel for
for(int i=0;i<omp_p;i++){
  size_t a=( i   *ntrg*target_dof)/omp_p;
  size_t b=((i+1)*ntrg*target_dof)/omp_p;
  //int r; MPI_Comm_rank(MPI_COMM_WORLD, &r);
  for(size_t iT=a; iT < b; iT++){
    pot[iT]=potv[iT];
    //COUTDEBUG("Potential at " << iT << " on processor " << r << " is : " << pot[iT]);
    if (std::isnan(pot[iT])){
      CERR("fmm potential is NaN",abort());
    }
    
  }
}
/*
  for(size_t iT=0;iT<potv.size();++iT){
    pot[iT]=potv[iT];
        //COUTDEBUG("Potential at " << iT << "is : " << pot[iT]);
    if (std::isnan(pot[iT])){
      CERR("fmm potential is NaN",abort());
    }
  }*/
}


void pvfmm_direct_summation(const size_t &sl_nsrc, const real_t *sl_src, 
                    const real_t *sl_den, const size_t &dl_nsrc, 
                    const real_t *dl_src, const real_t *dl_den_nor,
                    const size_t &ntrg, const real_t *trg, real_t *pot,
                    const int &rebuild_tree, void **context){
    ASSERT(*context != NULL, "was called with null context. "
            "Generate a context by calling make_pvfmm_context.");

    // context object
    PVFMMContext *ctx = (PVFMMContext*) *context;
    COUTDEBUG("pvfmm_direct_summation with sl_nsrc="<<sl_nsrc<<", dl_nsrc="<<dl_nsrc
            <<", ntrg="<<ntrg<<", rebuild_tree="<<rebuild_tree);


    COUTDEBUG("Using context with"
            <<" mult_order="<<ctx->mult_order
            <<", max_pts="<<ctx->max_pts
            <<", max_depth="<<ctx->max_depth
            <<", sl="<<ctx->sl
            <<", dl="<<ctx->dl
            <<", periodic="<<ctx->periodic
            );

    int source_dof = ctx->source_dof;
    int target_dof = ctx->target_dof;
    // copy data
    vec sl_srcv(sl_src,sl_src + DIM*sl_nsrc);
    vec sl_denv(sl_den,sl_den + source_dof*sl_nsrc);

    vec dl_srcv(dl_src    ,dl_src    +  DIM*dl_nsrc);
    vec dl_denv(dl_den_nor,dl_den_nor + (DIM + source_dof)*dl_nsrc);

    vec trgv(trg,trg + DIM*ntrg);
    vec potv(target_dof*ntrg);
    vec empty;

    assert(ctx->dense_eval);
    int omp_p=omp_get_max_threads();
#pragma omp parallel for
    for(int i=0;i<omp_p;i++){
        size_t a=( i   *ntrg)/omp_p;
        size_t b=((i+1)*ntrg)/omp_p;


        if(ctx->sl && !ctx->dl){
            ctx->ker->ker_poten(sl_srcv.data(), sl_nsrc,
                    sl_denv.data(), 1,
                    trgv.data() + a*3, b-a,
                    potv.data()+a*ctx->ker->ker_dim[1],NULL);
        } else if(!ctx->sl && ctx->dl){
            ctx->ker->dbl_layer_poten(dl_srcv.data(), dl_nsrc,
                    dl_denv.data(), 1,
                    trgv.data() + a*3, b-a,
                    potv.data()+a*ctx->ker->ker_dim[1],NULL);
        }
    }

#pragma omp parallel for
    for(int i=0;i<omp_p;i++){
        size_t a=( i   *ntrg*target_dof)/omp_p;
        size_t b=((i+1)*ntrg*target_dof)/omp_p;
        for(size_t iT=a; iT < b; iT++){
            pot[iT]=potv[iT];
            //COUTDEBUG("Potential at " << iT << "is : " << pot[iT]);
            if (std::isnan(pot[iT])){
                CERR("fmm potential is NaN",abort());
            }

        }
    }
}

void pvfmm_interaction_matrix(const size_t &sl_nsrc, const real_t *sl_src, 
                    const real_t *sl_den, const size_t &dl_nsrc, 
                    const real_t *dl_src, const real_t *dl_den_nor,
                    const size_t &ntrg, const real_t *trg, real_t *interaction_matrix,
                    const int &rebuild_tree, void **context){
    ASSERT(*context != NULL, "was called with null context. "
            "Generate a context by calling make_pvfmm_context.");

    // context object
    PVFMMContext *ctx = (PVFMMContext*) *context;
    /*COUTDEBUG("pvfmm_interaction_matrix with sl_nsrc="<<sl_nsrc<<", dl_nsrc="<<dl_nsrc
            <<", ntrg="<<ntrg<<", rebuild_tree="<<rebuild_tree);


    COUTDEBUG("Using context with"
            <<" mult_order="<<ctx->mult_order
            <<", max_pts="<<ctx->max_pts
            <<", max_depth="<<ctx->max_depth
            <<", sl="<<ctx->sl
            <<", dl="<<ctx->dl
            <<", periodic="<<ctx->periodic
            );*/

    int source_dof = ctx->source_dof;
    int target_dof = ctx->target_dof;
    assert(source_dof == ctx->ker->ker_dim[0]);
    assert(target_dof == ctx->ker->ker_dim[1]);
    // copy data
    vec sl_srcv(sl_src,sl_src + DIM*sl_nsrc);
    //vec sl_denv(sl_den,sl_den + source_dof*sl_nsrc);

    vec dl_srcv(dl_src    ,dl_src    +  DIM*dl_nsrc);
    vec dl_denv(dl_den_nor,dl_den_nor + (DIM + source_dof)*dl_nsrc);

    vec trgv(trg,trg + DIM*ntrg);
    vec interaction_matrixv((target_dof*ntrg) * (dl_nsrc*source_dof), 0.);


    //yanked from Kernel<T>::BuildMatrix in kernel.txx
    //(T* r_src, int src_cnt, T* r_trg, int trg_cnt, T* k_out) const{
    //int dim=3; //Only supporting 3D
    //memset(interaction_matrix, 0, src_cnt*ker_dim[0]*trg_cnt*ker_dim[1]*sizeof(T));
    if(ctx->sl && !ctx->dl){
        for(int i = 0; i < sl_nsrc; i++) {
            for(int j = 0; j < source_dof; j++){
                std::vector<real_t> v_src(source_dof,0.);
                v_src[j]=1.0;
                ctx->ker->ker_poten(&sl_srcv[DIM*i], 1, &v_src[0], source_dof, trgv.data(), ntrg,
                        &interaction_matrix[(i*source_dof+j)*ntrg*target_dof], NULL);
            }
        }
    }  else if (ctx->dl && !ctx->sl){
        for(int i = 0; i < dl_nsrc; i++) {
            for(int j = 0; j < source_dof; j++){
                std::vector<real_t> v_src(source_dof+DIM,0.);
                for (int d = 0; d < DIM; d++) {
                    v_src[d] = dl_denv[(DIM+source_dof)*i+d];
                }
                v_src[DIM + j]=1.0;
                ctx->ker->dbl_layer_poten(&dl_srcv[DIM*i], 1, 
                        &v_src[0], source_dof, 
                        trgv.data(), ntrg,
                        &interaction_matrix[(i*source_dof+j)*ntrg*target_dof], NULL);
            }
        }
    }
}


//---------------------------------------------------------------------------
// Begin PvFMM tree wrapper 
//---------------------------------------------------------------------------

void construct_id_maps(void **context){
    
    PVFMMContext *ctx = (PVFMMContext*) *context;

    ASSERT(ctx != NULL, "was called with a null context. Generate a context \
                        by calling make_pvfmm_context");
    ASSERT(ctx->tree != NULL, "PvFMM is  already uninitialized");
    //ASSERT(ctx->kifmm_to_pvfmm_ids == NULL, "kifmm_to_pvfmm_ids already initialized");
    //ASSERT(ctx->pvfmm_to_kifmm_ids == NULL, "pvfmm_to_kifmm_ids already initialized");

    // Build map from tree level i to PvFMM box ID's that occur on level i.
    pvfmm::PtFMM_Tree<double>* tree = ctx->tree;
    Node_t* current_node = tree->PreorderFirst();
    std::map<int , std::vector<pvfmm::MortonId> > level_to_box_ids;
    
    int depth = 0 ;
    while (current_node != NULL){
        int depth = current_node->Depth();
        pvfmm::MortonId mid = current_node->GetMortonId();
        level_to_box_ids[depth].push_back(mid);
        current_node = tree->PreorderNxt(current_node);
    } 

    // Constuct a single vector of MortonId's in order of decreasing box depth 
    // in the tree
    
    int i = 0; //root level
    std::vector<pvfmm::MortonId> top_down_level_order;
    
    // Walk down the levels of the tree and concatenate the box ids
    while(level_to_box_ids.count(i)){
        std::vector<pvfmm::MortonId> ith_level_boxes = level_to_box_ids[i];
        
        // Concatenate the ith level vector to the output vector
        top_down_level_order.insert(top_down_level_order.end(), 
                                    ith_level_boxes.begin(),
                                    ith_level_boxes.end());
        i++;
    }
    std::map<int, pvfmm::MortonId>* kifmm_to_pvfmm_ids =
                                        new std::map<int, pvfmm::MortonId>();
    std::map<pvfmm::MortonId, int>* pvfmm_to_kifmm_ids = 
                                        new std::map<pvfmm::MortonId, int>();

    // Dump these boxes into the final maps to/from Kifmm from/to Pvfmm
    for(int j = 0; j < top_down_level_order.size(); j++){
        pvfmm::MortonId jth_box_in_level_order = top_down_level_order[j];
        (*kifmm_to_pvfmm_ids)[j] = jth_box_in_level_order;
        (*pvfmm_to_kifmm_ids)[jth_box_in_level_order] = j;
    }

    ctx->kifmm_to_pvfmm_ids = kifmm_to_pvfmm_ids;
    ctx->pvfmm_to_kifmm_ids = pvfmm_to_kifmm_ids;
    COUTDEBUG("KIFMM-to-PvFMM ids");
    for(int j = 0; j < top_down_level_order.size(); j++){
        std::vector<double> c=box_center(j, context);
        COUTDEBUG(j << "-  MortonId: " << (*kifmm_to_pvfmm_ids)[j]
                    << "; center: " << c[0] << ", " << c[1] <<  ", " << c[2]);
    }
    *context = (void*) ctx;
}

// Context managing functions
void initialize_pvfmm_tree(vec sl_srcv, vec sl_denv, vec dl_srcv, vec dl_denv,
                           vec trgv, void **context){
    PVFMMContext *ctx = (PVFMMContext*) *context;
    ASSERT(ctx != NULL, "was called with a null context. Generate a context \
                        by calling make_pvfmm_context");
    ASSERT(ctx->mat != NULL, "Matrices uninitialized");
    ASSERT(ctx->tree == NULL, "Tree is already initialized.");

    vec empty;
    if( ctx->sl ){
        ctx->tree = pvfmm::PtFMM_CreateTree<double>(sl_srcv, sl_denv, empty, empty,
                trgv,ctx->comm, ctx->max_pts, ctx->boundary);
    } else { //Double layer 
        // TODO eliminate conditional; tree structure shouldn't depend 
        // on normals + density
        ctx->tree = pvfmm::PtFMM_CreateTree<double>( empty, empty, dl_srcv, dl_denv,
                trgv,ctx->comm, ctx->max_pts, ctx->boundary);
    }
    ctx->tree->SetupFMM(ctx->mat);

    // Setup of mapping from KIFMM ids to PVFMM Morton Ids
    construct_id_maps((void**) &ctx);
    *context = (void*) ctx;
}

void clear_pvfmm_tree(void **context){
  if(*context == NULL)
      return;
  COUTDEBUG("deleting pvfmm tree");
  PVFMMContext *ctx = (PVFMMContext*) *context;
  if(ctx->tree){
      delete ctx->tree;
  }
  if(ctx->kifmm_to_pvfmm_ids){
      delete ctx->kifmm_to_pvfmm_ids;
  }
  if(ctx->pvfmm_to_kifmm_ids){
      delete ctx->pvfmm_to_kifmm_ids;
  } 
}

// Trnaslating between PvFMM Z-order indexing and KIFMM level order indexing
pvfmm::MortonId to_morton_id(int node_id, void **context){
    PVFMMContext *ctx = (PVFMMContext*) *context;
    std::map<int, pvfmm::MortonId> m = *(ctx->kifmm_to_pvfmm_ids);
    return m[node_id];
}

int to_kifmm_id(pvfmm::MortonId mid, void **context){
    PVFMMContext *ctx = (PVFMMContext*) *context;
    std::map<pvfmm::MortonId, int> m = *(ctx->pvfmm_to_kifmm_ids);
    return m[mid];
}

Node_t* find_node(int node_id, void **context){
    PVFMMContext *ctx = (PVFMMContext*) *context;
    pvfmm::MortonId mid = to_morton_id(node_id, context);
    return ctx->tree->FindNode(mid, false);
}

// Necessary FMM box functions on that need to be exposed
bool is_leaf(int node_id, void **context){
    return find_node(node_id, context)->IsLeaf();
}

std::vector<double> box_center(int node_id, void **context){
    Node_t* node = find_node(node_id, context);
    int depth = node->Depth();
    std::vector<double> center(node->Coord(), node->Coord() + 3);
    for(int i=0; i < 3; i++){
        center[i] += pow(.5, depth+1);
    }
    return center;
}

int depth(int node_id, void **context){
    return find_node(node_id, context)->Depth();
}

std::vector<int> child_indices(int node_id, void **context){
    pvfmm::MortonId mid = to_morton_id(node_id, context);
    std::vector<pvfmm::MortonId> children_as_mids = mid.Children();

    std::vector<int> children_as_box_ids;
    for(int i=0; i < children_as_mids.size(); i++){
        int box_id = to_kifmm_id(children_as_mids[i], context);
        children_as_box_ids.push_back(box_id);
    }

    return children_as_box_ids;
}

std::vector<int> target_indices_in_box(int node_id, void **context){
    Node_t* node = find_node(node_id, context);
    pvfmm::Vector<size_t> v = node->trg_scatter;
    std::vector<int> copied_scatter;
    int index;
    for(int i = 0; i < v.Capacity(); i++){
        index = v[i];
        copied_scatter.push_back(index);
    }
    return copied_scatter;
}

std::vector<int> source_indices_in_box(int node_id, void **context){
    Node_t* node = find_node(node_id, context);
    pvfmm::Vector<size_t> v = node->src_scatter;
    std::vector<int> copied_scatter;
    int index;
    for(int i = 0; i < v.Capacity(); i++){
        index = v[i];
        copied_scatter.push_back(index);
    }
    return copied_scatter;
}

// TODO cache these funtions within construct_id_maps. No reason to do more
//      than one tree traversal since this data is built up once already.
std::map<int, std::vector<int> > build_level_to_box_ids_map(void **context){
  
    ASSERT(*context != NULL, "was called with null context. "
	 "Generate a context by calling make_pvfmm_context.");
    
    PVFMMContext *ctx = (PVFMMContext*) *context;

    pvfmm::PtFMM_Tree<double>* tree = ctx->tree;
    
    std::map<int , std::vector<int> > level_to_box_ids;
    int depth = 0;
    Node_t* current_node = tree->PreorderFirst();

    while (current_node != NULL){
        int depth = current_node->Depth();
        pvfmm::MortonId mid = current_node->GetMortonId();
        level_to_box_ids[depth].push_back(to_kifmm_id(mid, context));
        current_node = tree->PreorderNxt(current_node);
    } 
    
    return level_to_box_ids;
}

std::vector<int> top_down_level_order(void **context){
    std::map<int, std::vector<int> > level_to_box_ids = 
            build_level_to_box_ids_map(context);
    int i = 0;
    std::vector<int> top_down_level_order;
    
    // Walk down the levels of the tree and concatenate the box ids
    while(level_to_box_ids.count(i)){
        std::vector<int> ith_level_boxes = level_to_box_ids[i];
        // Concatenate the ith level vector to the output vector
        top_down_level_order.insert(top_down_level_order.end(), 
                                    ith_level_boxes.begin(),
                                    ith_level_boxes.end());
        i++;
    }
    return top_down_level_order;
}

std::vector<int> bottom_up_level_order(void **context){
    std::map<int, std::vector<int> > level_to_box_ids = 
            build_level_to_box_ids_map(context);
    
    // Determine the height of the tree
    int i = 0;
    while(level_to_box_ids.count(i)){
        i++;
    }

    std::vector<int> bottom_up_level_order, ith_level_boxes;
    
    // Walk down the levels of the tree and concatenate the box ids
    // i := tree_depth
    while(level_to_box_ids.count(i)){
        std::vector<int> ith_level_boxes = level_to_box_ids[i];
        // Concatenate the ith level vector to the output vector
        bottom_up_level_order.insert(bottom_up_level_order.end(), 
                                    ith_level_boxes.begin(),
                                    ith_level_boxes.end());
        i--;
    }
    return bottom_up_level_order;
}
//---------------------------------------------------------------------------
// End PvFMM tree wrapper 
//---------------------------------------------------------------------------
void mpi_init_pvfmm(int &rank, int &size){
  int argc(0);
  char **argv;
  MPI_Init(&argc, &argv);
  MPI_Comm comm=MPI_COMM_WORLD;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  COUTDEBUG("mpi_init (rank/size="<<rank<<"/"<<size<<")");
}

void mpi_finalize_pvfmm(){
  int rank,size;
  MPI_Comm comm=MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  COUTDEBUG("finalizing MPI (rank/size="<<rank<<"/"<<size<<")");
  MPI_Finalize();
}
