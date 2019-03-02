# plot error, QBX solver convergence
plot_qbx_solver_error_convergence() {
    python utils/plot_results.py \
        --prefix "qbx_gmres_ref_lvl_" --name "QBX" \
        --semilogx --legend \
        --xvar="face-map max patch size" \
        --yvar "max absolute error"  \
        --label "error" \
        --output-file "error.png"  \
        --x-label "\$L\$" \
        --y-label "\$\log_{10}\| u - \\tilde{u} \|_\infty\$" \
        --outdir="output/plots-dev" \
        --indir="$1"
}
# plot error, QBX vs singular quadrature
plot_qbx_vs_singular_quad_error_comparision() {
python utils/plot_results.py \
    --prefix "qbx__ref_lvl_" --name "QBX" \
    --prefix "singular_eval__ref_lvl_" --name "Singular Quad" \
    --semilogx --legend \
    --xvar="coarse spacing" \
    --yvar "max absolute error"  \
    --label "error" \
    --xvar="face-map max patch size" \
    --yvar "max absolute error"  \
    --label "error" \
    --output-file "error.png"  \
    --x-label "\$h\$" \
    --y-label "\$\log_{10}\| u - \\tilde{u} \|_\infty\$" \
    --outdir="output/plots-dev" \
    --indir=$1
    }
plot_qbx_vs_singular_quad_timing_comparision() {
python utils/plot_results.py \
    --prefix "qbx__ref_lvl_" --name "QBX" \
    --prefix "singular_eval__ref_lvl_" --name "Sing. Quad" \
    --semilogx --legend \
    --yvar "max absolute error"  \
    --xvar "total matvec time" \
    --label "time" \
    --yvar "max absolute error"  \
    --xvar "total fft upsample time" \
    --label "fft upsample time" \
    --yvar "max absolute error"  \
    --xvar "total polar quad time" \
    --label "polar quad time" \
    --yvar "max absolute error"  \
    --xvar "total inaccurate subtract time" \
    --label "subtract time" \
    --yvar "max absolute error"  \
    --xvar "total density interp time" \
    --label "density interp time" \
    --yvar "max absolute error"  \
    --xvar "total qbx time" \
    --label "extrap time" \
    --yvar "max absolute error"  \
    --xvar "total fmm time" \
    --label "FMM time" \
    --output-file "timing.png"  \
    --x-label "\$h\$" \
    --y-label "\$\log_{10}\| u - \\tilde{u} \|_\infty\$" \
    --outdir="output/plots-dev" \
    --indir=$1
    }


# plot error, QBX greens identity
plot_qbx_greens_identity_error_convergence(){
python utils/plot_results.py \
    --prefix "qbx_greens_id_ref_lvl_" --name "QBX" \
    --semilogx --legend \
    --xvar="face-map max patch size" \
    --yvar "max absolute error"  \
    --label "error" \
    --output-file "error.png"  \
    --x-label "\$L\$" \
    --y-label "\$\log_{10}\| u - \\tilde{u} \|_\infty\$" \
    --outdir="output/plots-dev" \
    --indir=$1
    #--x-label "\$\max_P \sqrt(|P|)\$" \
    #--y-label "\$\log_{10}\| u - \\tilde{u} \|_\infty\$" \
}
# plot timing, QBX greens identity
plot_qbx_greens_identity_timing_convergence(){
python utils/plot_results.py \
    --prefix "qbx_greens_id_ref_lvl_" --name "QBX" \
    --semilogx --legend \
    --yvar "max absolute error"  \
    --xvar "total matvec time" \
    --label "total time" \
    --yvar "max absolute error"  \
    --xvar "total fmm time" \
    --label "FMM time" \
    --yvar "max absolute error"  \
    --xvar "total qbx time" \
    --label "extrapolation time" \
    --yvar "max absolute error"  \
    --xvar "total density interp time" \
    --label "density interp. time" \
    --output-file "timing.png"  \
    --x-label "OMP Wall time" \
    --y-label "\$\log_{10}\| u - \\tilde{u} \|_\infty\$"  \
    --outdir="output/plots-dev"\
    --indir=$1
}

if [ ] ; then
plot_qbx_solver_error_convergence "data/results3/laplace_qbx_eval_conv_cube"
plot_qbx_solver_error_convergence "data/results3/laplace_qbx_eval_conv_pipe"
plot_qbx_solver_error_convergence "data/results3/laplace_qbx_eval_conv_ttorus"
plot_qbx_solver_error_convergence "data/results3/stokes_qbx_eval_conv_cube"
plot_qbx_solver_error_convergence "data/results3/stokes_qbx_eval_conv_pipe"
plot_qbx_solver_error_convergence "data/results3/stokes_qbx_eval_conv_ttorus"
plot_qbx_solver_error_convergence "data/results3/navier_qbx_eval_conv_cube"
plot_qbx_solver_error_convergence "data/results3/navier_qbx_eval_conv_pipe"
plot_qbx_solver_error_convergence "data/results3/navier_qbx_eval_conv_ttorus"
plot_qbx_vs_singular_quad_error_comparision "data/results3/sph_harm_qbx_vs_singular_quad"
plot_qbx_vs_singular_quad_error_comparision "data/results3/laplace_qbx_vs_singular_quad"
plot_qbx_vs_singular_quad_error_comparision "data/results3/navier_qbx_vs_singular_quad"
plot_qbx_vs_singular_quad_error_comparision "data/results3/face_map_qbx_vs_singular_quad"
    echo "hi"

plot_qbx_vs_singular_quad_timing_comparision "data/results3/sph_harm_qbx_vs_singular_quad"
plot_qbx_vs_singular_quad_timing_comparision "data/results3/laplace_qbx_vs_singular_quad"
plot_qbx_vs_singular_quad_timing_comparision "data/results3/navier_qbx_vs_singular_quad"
plot_qbx_vs_singular_quad_timing_comparision "data/results3/face_map_qbx_vs_singular_quad"
fi
plot_qbx_greens_identity_error_convergence "data/results3/laplace_greens_identity_cube"
plot_qbx_greens_identity_error_convergence "data/results3/laplace_greens_identity_newtorus"
plot_qbx_greens_identity_error_convergence "data/results3/navier_greens_identity_cube"
plot_qbx_greens_identity_error_convergence "data/results3/navier_greens_identity_newtorus"
plot_qbx_greens_identity_error_convergence "data/results3/stokes_greens_identity_cube"
plot_qbx_greens_identity_error_convergence "data/results3/stokes_greens_identity_newtorus"
plot_qbx_greens_identity_timing_convergence "data/results3/laplace_greens_identity_cube"
plot_qbx_greens_identity_timing_convergence "data/results3/laplace_greens_identity_newtorus"
plot_qbx_greens_identity_timing_convergence "data/results3/navier_greens_identity_cube"
plot_qbx_greens_identity_timing_convergence "data/results3/navier_greens_identity_newtorus"
plot_qbx_greens_identity_timing_convergence "data/results3/stokes_greens_identity_cube"
plot_qbx_greens_identity_timing_convergence "data/results3/stokes_greens_identity_newtorus"
# plot timing, QBX vs singular quadrature
if [ ] ; then
python utils/plot_results.py \
    --prefix "qbx__ref_lvl_" --name "QBX" \
    --prefix "singular_eval__ref_lvl_" --name "Sing. Quad" \
    --semilogx --legend \
    --yvar "max absolute error"  \
    --xvar "total matvec time" \
    --label "time" \
    --yvar "max absolute error"  \
    --xvar "total fft upsample time" \
    --label "fft upsample time" \
    --yvar "max absolute error"  \
    --xvar "total polar quad time" \
    --label "polar quad time" \
    --yvar "max absolute error"  \
    --xvar "total inaccurate subtract time" \
    --label "subtract time" \
    --yvar "max absolute error"  \
    --xvar "total density interp time" \
    --label "density interp time" \
    --yvar "max absolute error"  \
    --xvar "total qbx time" \
    --label "extrap time" \
    --yvar "max absolute error"  \
    --xvar "total fmm time" \
    --label "FMM time" \
    --output-file "timing.png"  \
    --x-label "\$h\$" \
    --y-label "\$\log_{10}\| u - \\tilde{u} \|_\infty\$" \
    --outdir="output/plots-dev"
    --indir="data/results3"


# plot error, QBX vs singular quadrature
python utils/plot_results.py \
    --prefix "qbx__ref_lvl_" --name "QBX" \
    --prefix "singular_eval__ref_lvl_" --name "Singular Quad" \
    --semilogx --legend \
    --xvar="coarse spacing" \r
    --yvar "max absolute error"  \
    --label "error" \
    --xvar="face-map max patch size" \
    --yvar "max absolute error"  \
    --label "error" \
    --output-file "error.png"  \
    --x-label "\$h\$" \
    --y-label "\$\log_{10}\| u - \\tilde{u} \|_\infty\$" \
    --outdir="output/plots-dev"
python utils/plot_results.py \
    --prefix "qbx_adaptive_gmres_ref_lvl_" --name "Uniform QBX" \
    --prefix "qbx_gmres_ref_lvl_" --name "Adaptive QBX" \
    --semilogx --legend \
    --yvar "max absolute error"  \
    --xvar "total matvec time" \
    --label "time" \
    --yvar "max absolute error"  \
    --xvar "total density interp time" \
    --label "density interp time" \
    --yvar "max absolute error"  \
    --xvar "total qbx time" \
    --label "extrap time" \
    --yvar "max absolute error"  \
    --xvar "total fmm time" \
    --label "FMM time" \
    --output-file "timing.png"  \
    --x-label "OMP walltime" \
    --y-label "\$\log_{10}\| u - \\tilde{u} \|_\infty\$" \
    --outdir="output/plots-dev"
python utils/plot_results.py \
    --prefix "qbx_adaptive_gmres_ref_lvl_" --name "Adaptive QBX" \
    --prefix "qbx_gmres_ref_lvl_" --name "Uniform QBX" \
    --semilogx --legend \
    --yvar "max absolute error"  \
    --xvar "total coarse samples" \
    --label "dofs" \
    --output-file "dofs.png"  \
    --x-label "# degrees of freedom" \
    --y-label "\$\log_{10}\| u - \\tilde{u} \|_\infty\$" \
    --outdir="output/plots-dev"
# plot error, QBX greens identity
python utils/plot_results.py \
    --prefix "qbx_greens_id_ref_lvl_" --name "QBX" \
    --semilogx --legend \
    --xvar="face-map max patch size" \
    --yvar "max absolute error"  \
    --label "error" \
    --output-file "error.png"  \
    --x-label "\$L\$" \
    --y-label "\$\log_{10}\| u - \\tilde{u} \|_\infty\$" \
    --outdir="output/plots-dev"
    #--x-label "\$\max_P \sqrt(|P|)\$" \
    #--y-label "\$\log_{10}\| u - \\tilde{u} \|_\infty\$" \

# plot timing, QBX greens identity
python utils/plot_results.py \
    --prefix "qbx_greens_id_ref_lvl_" --name "QBX" \
    --semilogx --legend \
    --yvar "max absolute error"  \
    --xvar "total matvec time" \
    --label "total time" \
    --yvar "max absolute error"  \
    --xvar "total fmm time" \
    --label "FMM time" \
    --yvar "max absolute error"  \
    --xvar "total qbx time" \
    --label "extrapolation time" \
    --yvar "max absolute error"  \
    --xvar "total density interp time" \
    --label "density interp. time" \
    --output-file "timing.png"  \
    --x-label "OMP Wall time" \
    --y-label "\$\log_{10}\| u - \\tilde{u} \|_\infty\$"  \
    --outdir="output/plots-dev"
    echo "hi"
fi
