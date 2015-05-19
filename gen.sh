#!/bin/bash

python3 ./generate_twoel.py -l 2 \
                            -g build/generator/generator \
                            -d eri/FO

touch eri/CMakeLists.txt
