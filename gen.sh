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

python3 ./generate_twoel.py -l 2 \
                            -b vref \
                            -p vref \
                            -g build/generator/generator \
                            -d eri/vref

touch eri/CMakeLists.txt
