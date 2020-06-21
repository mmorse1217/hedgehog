#!/bin/bash

# code to execute during CI build
mkdir -p build/ 
cd build/
NANOSPLINE_DIR=/libs/nanospline cmake ..  
make 
cd ../
build/tests/test_hedgehog [critical]

