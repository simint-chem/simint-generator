#!/bin/bash

python3 ./generate_twoel.py -l 2 \
                            -b FO \
                            -p FO \
                            -g build/generator/generator \
                            -d eri/FO

python3 ./generate_twoel.py -l 2 \
                            -b split \
                            -p split \
                            -g build/generator/generator \
                            -d eri/split

touch eri/CMakeLists.txt
