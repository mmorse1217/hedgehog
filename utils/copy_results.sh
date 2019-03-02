rsync -rav -e ssh --include '*/' --include='*.pydict' --exclude='*' \
    mjm1030@prince0.hpc.nyu.edu:/scratch/mjm1030/mobo-temp/ data/results14
