module fmm_interface

  use, intrinsic :: iso_c_binding
  implicit none

  !!----------------------------------------------------
  !!-- C function interface ----------------------------
  !!----------------------------------------------------
  interface c_interface

     !!-- initializing mpi -----------------------------
     subroutine mpi_init_pvfmm(rank, size) bind(C)
       use, intrinsic :: iso_c_binding
       integer (C_INT), intent(out) :: rank, size
     end subroutine mpi_init_pvfmm

     !!-- finalizing mpi -------------------------------
     subroutine mpi_finalize_pvfmm() bind(C)
     end subroutine mpi_finalize_pvfmm

     !!-- making context -------------------------------
     subroutine make_pvfmm_context(mult_order, max_pts, max_depth, sl, dl, periodic, context) bind(C)
       ! suggested values are: mult_order=10, max_pts=400, max_depth=20
       use, intrinsic :: iso_c_binding

       integer (C_INT), intent(in) :: mult_order, max_pts, max_depth
       logical (C_INT), intent(in) :: sl, dl, periodic
       type (C_PTR), intent(inout) :: context
     end subroutine make_pvfmm_context

     !!-- clearing context -----------------------------
     subroutine clear_pvfmm_context(context) bind(C)
       use, intrinsic :: iso_c_binding

       type (C_PTR), intent(inout) :: context
     end subroutine clear_pvfmm_context

     !!-- evaluating fmm  ------------------------------
     subroutine stokes_sl_pvfmm(nsrc, src, den, ntrg, trg, pot, rebuild_tree, context) bind(C)
       use, intrinsic :: iso_c_binding

       real (C_DOUBLE), intent(in) :: src(*), den(*), trg(*)
       integer (C_SIZE_T), intent(in) :: nsrc, ntrg
       logical (C_INT), intent(in) :: rebuild_tree
       real (C_DOUBLE), intent(out) :: pot(*)
       type (C_PTR), intent(inout) :: context
     end subroutine stokes_sl_pvfmm

     subroutine stokes_sldl_pvfmm(sl_nsrc, sl_src, sl_den, dl_nsrc, dl_src, dl_den_nor, ntrg, trg, pot, rebuild_tree, context) bind(C)
       use, intrinsic :: iso_c_binding

       real (C_DOUBLE), intent(in) :: sl_src(*), sl_den(*), dl_src(*), dl_den_nor(*), trg(*)
       integer (C_SIZE_T), intent(in) :: sl_nsrc, dl_nsrc, ntrg
       logical (C_INT), intent(in) :: rebuild_tree
       real (C_DOUBLE), intent(out) :: pot(*)
       type (C_PTR), intent(inout) :: context
     end subroutine stokes_sldl_pvfmm

  end interface c_interface

contains
  !!----------------------------------------------------
  !!-- Wrapper functions -------------------------------
  !!----------------------------------------------------
  subroutine mpi_init(rank, size)
    integer (C_INT), intent(out) :: rank, size

    call mpi_init_pvfmm(rank, size)
  end subroutine mpi_init

  subroutine mpi_finalize()
    call mpi_finalize_pvfmm()
  end subroutine mpi_finalize

  subroutine make_fmm_context(mult_order, max_pts, max_depth, sl, dl, periodic, context)
    integer (C_INT), intent(in) :: mult_order, max_pts, max_depth
    logical (C_INT), intent(in) :: sl, dl, periodic
    type (C_PTR), intent(inout) :: context

    call make_pvfmm_context(mult_order, max_pts, max_depth, sl, dl, periodic, context)
  end subroutine make_fmm_context

  subroutine clear_fmm_context(context)
    type (C_PTR), intent(inout) :: context

    call clear_pvfmm_context(context)
  end subroutine clear_fmm_context

  subroutine stokes_sl_fmm(nsrc, src, den, ntrg, trg, pot, rebuild_tree, context)
    real (C_DOUBLE), intent(in) :: src(*), den(*), trg(*)
    integer (C_SIZE_T), intent(in) :: nsrc, ntrg
    logical (C_INT), intent(in) :: rebuild_tree
    real (C_DOUBLE), intent(out) :: pot(*)
    type (C_PTR), intent(inout) :: context

    call stokes_sl_pvfmm(nsrc, src, den, ntrg, trg, pot, rebuild_tree, context)
  end subroutine stokes_sl_fmm

  subroutine stokes_sldl_fmm(sl_nsrc, sl_src, sl_den, dl_nsrc, dl_src, dl_den_nor, ntrg, trg, pot, rebuild_tree, context)
    real (C_DOUBLE), intent(in) :: sl_src(*), sl_den(*), dl_src(*), dl_den_nor(*), trg(*)
    integer (C_SIZE_T), intent(in) :: sl_nsrc, dl_nsrc, ntrg
    logical (C_INT), intent(in) :: rebuild_tree
    real (C_DOUBLE), intent(out) :: pot(*)
    type (C_PTR), intent(inout) :: context
  
    call stokes_sldl_pvfmm(sl_nsrc, sl_src, sl_den, dl_nsrc, dl_src, dl_den_nor, ntrg, trg, pot, rebuild_tree, context)
end subroutine stokes_sldl_fmm

end module fmm_interface
