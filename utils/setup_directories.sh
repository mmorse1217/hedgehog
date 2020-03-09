declare -a test_cases=("test_solver_convergence_qbkix" 
                "test_greens_identity"
                "test_qbkix_eval_convergence"
                "test_qbkix_vs_singular_quad"
                "test_qbkix_vs_singular_quad/face_map_vs_blended"
                "test_adaptive/eye_candy"
                "test_adaptive"
                )
declare -a domains=("solve/cube" "const_density/cube" "cube" "sphere" "pipe" "ttorus2" "new_ppp" "newtorus" 
                "half_donut" "larger_vessel_section2" "nearly_touching" "comb2" "closed_octopus")
                
declare -a kernels=("laplace" "stokes" "navier")

mkdir -p data/
for test_case in "${test_cases[@]}"
do
    for domain in "${domains[@]}"
    do
        for kernel in "${kernels[@]}"
        do
            mkdir -p "output/${test_case}/${domain}/${kernel}"
            # or do whatever with individual element of the array
        done
    done
done
mkdir -p "output/plots"
mkdir -p "output/quadrature_data"
