#!/bin/bash

# code to execute during CI build
mkdir -p /hedgehog/build/ 
cd /hedgehog/build/
cmake ..  
make -j${nprocs}
cd ../
bin/test_patchwork [critical]
