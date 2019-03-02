#script to submit a test case in build/test_bis to run on prince.
base_vars=(--optfile=opt/morse_cases.opt 
            --host=prince  
            -Q slurm 
            --cpu 1
            --nodes 1 
            --symlink-qbx 
            --exec-name build/test_bis 
            --init-basedir ${SCRATCH}/mobo-temp 
            --job-name $1  
             --vars "${@:2}")
python utils/stage_job.py "${base_vars[@]}"  --walltime 6h  --threads 16 --memory 64  --no-stage 
