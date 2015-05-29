#!/bin/bash

python3 ../generate_twoel.py -l 2 \
                             -b FO \
                             -p FO \
                             -g generator/generator \
                             -d ../eri/FO

python3 ../generate_twoel.py -l 2 \
                             -f \
                             -b FO \
                             -p FO_flat \
                             -g generator/generator \
                             -d ../eri/FO_flat

#python3 ../generate_twoel.py -l 2 \
#                             -b split \
#                             -p split \
#                             -g generator/generator \
#                             -d ../eri/split

python3 ../generate_twoel.py -l 3 \
                             -b vref \
                             -p vref \
                             -g generator/generator \
                             -d ../eri/vref

python3 ../generate_twoel.py -l 3 \
                             -f \
                             -b vref \
                             -p vref_flat \
                             -g generator/generator \
                             -d ../eri/vref_flat

touch eri/CMakeLists.txt
