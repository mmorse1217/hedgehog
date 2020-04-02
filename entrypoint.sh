#!/bin/bash

# code to execute during CI build
mkdir -p build/ 
cd build/
cmake ..  
make 
cd ../
build/tests/test_hedgehog [critical]

