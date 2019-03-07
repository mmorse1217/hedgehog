# Run greens identity test
base_vars=(--optfile=opt/morse_cases.opt 
            --host=prince  
            -Q slurm 
            --nodes 1 
            --symlink-qbx 
            --exec-name build/test_bis 
            --init-basedir ${SCRATCH}/mobo-temp
            --exclusive)
            #--no-stage)
submit_job (){
    python utils/stage_job.py "${base_vars[@]}"  --walltime $3  --memory $4 --cpu $5 --threads $5 --job-name $1  --vars $2 
}


# pushing for 1e-8 error; only reason we need so much memory
submit_job "eye_candy_cube" "[eye-candy][results][cube]" 12h 128 28
submit_job "eye_candy_propeller" "[eye-candy][results][propeller]" 12h 128 28
submit_job "stokes_qbx_eval_conv_pipe" "[results][qbkix-solver-conv][stokes][pipe]" 24h 250 28
submit_job "eye_candy_octopus" "[eye-candy][results][octopus]" 20h 250 28
if [ ]; then


submit_job "laplace_qbx_eval_conv_cube" "[results][qbkix-solver-conv][laplace][cube]" 12h 72 20
submit_job "laplace_greens_identity_cube"  "[results][gid][laplace][cube]" 12h 84 20

submit_job "navier_greens_identity_cube"  "[results][gid][navier][cube]" 12h 100 20 
submit_job "stokes_greens_identity_cube"  "[results][gid][stokes][cube]" 12h 100 20

submit_job "laplace_greens_identity_newtorus2" "[results][gid][laplace][newtorus]" 6h 120 20 
submit_job "navier_greens_identity_newtorus2"  "[results][gid][navier][newtorus]" 6h 130 20 
submit_job "stokes_greens_identity_newtorus2"  "[results][gid][stokes][newtorus]" 6h 130 20

submit_job "adaptive_test_propeller" "[adaptive][propeller][results]" 12h 100 20

submit_job "adaptive_test_propeller2" "[adaptive][propeller][results]" 12h 100 20
submit_job "adaptive_test_comb" "[adaptive][comb][results]" 12h 100 20
submit_job "adaptive_test_pinch" "[adaptive][nearly_touching][results]" 12h 100 20
submit_job "adaptive_test_octopus" "[adaptive][octopus][results]" 12h 120 20
submit_job "adaptive_cube" "[adaptive][cube][results]" 6h 120 20
submit_job "stokes_qbx_eval_conv_ttorus" "[results][qbkix-solver-conv][stokes][ttorus]" 24h 84 20 

submit_job "face_map_qbx_vs_singular_quad_navier" "[results][qbx-vs-singular-quad][mixed-bdry-navier]" 10h 80 20
submit_job "face_map_qbx_vs_singular_quad" "[results][qbx-vs-singular-quad][mixed-bdry]" 10h 80 20
submit_job "navier_qbx_eval_conv_cube" "[results][qbkix-solver-conv][navier][cube]" 12h 84 20 
submit_job "stokes_qbx_eval_conv_cube" "[results][qbkix-solver-conv][stokes][cube]" 12h 84 20 

submit_job "laplace_qbx_eval_conv_pipe" "[results][qbkix-solver-conv][laplace][pipe]" 12h 72 20
submit_job "navier_qbx_eval_conv_pipe" "[results][qbkix-solver-conv][navier][pipe]" 12h 84 20 

submit_job "laplace_qbx_eval_conv_ttorus" "[results][qbkix-solver-conv][laplace][ttorus]" 12h 72 20 
submit_job "navier_qbx_eval_conv_ttorus" "[results][qbkix-solver-conv][navier][ttorus]" 24h 84 20 

# greens identity convergence rates
# solver convergence with QBX evaluation 

submit_job "uniform_vs_adaptive_half_donut" "[adaptive-vs-uniform][half_donut][results]" 36h 100 20

submit_job "stokes_qbx_eval_conv_hourglass" "[results][qbkix-solver-conv][stokes][hourglass]" 12h 72 20

submit_job "laplace_qbx_vs_singular_quad_solve" "[results][qbx-vs-singular-quad][laplace-solve]" 6h 72 20
submit_job "navier_qbx_vs_singular_quad_const_den" "[results][qbx-vs-singular-quad][navier-const-den]" 6h 72 20 







# qbx vs. singular quadrature comparison
submit_job "laplace_qbx_vs_singular_quad" "[results][qbx-vs-singular-quad][laplace]" 6h 72 20
submit_job "navier_qbx_vs_singular_quad" "[results][qbx-vs-singular-quad][navier]" 6h 72 20 
submit_job "sph_harm_qbx_vs_singular_quad" "[results][qbx-vs-singular-quad][density]" 4h 72 20 










# qbx face-map uniform vs adaptive comparison

submit_job "laplace_high_curvature_squished_cube" "[qbkix-curve][results][laplace]" 4h 72

fi
