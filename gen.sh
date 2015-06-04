#!/bin/bash

VE=1000
HE=1000

python3 ../generate_twoel.py -l 2 -ve ${VE} -he ${HE} \
                             -b FO \
                             -p FO \
                             -g generator/eri_generator \
                             -d ../eri/FO

python3 ../generate_twoel.py -l 2 -ve ${VE} -he ${HE} \
                             -f \
                             -b FO \
                             -p FO_flat \
                             -g generator/eri_generator \
                             -d ../eri/FO_flat

#python3 ../generate_twoel.py -l 2 -ve ${VE} -he ${HE} \
#                             -b split \
#                             -p split \
#                             -g generator/generator -he ${HE} \
#                             -d ../eri/split

python3 ../generate_twoel.py -l 3 -ve ${VE} -he ${HE} \
                             -b vref \
                             -p vref \
                             -g generator/eri_generator \
                             -d ../eri/vref

python3 ../generate_twoel.py -l 3 -ve ${VE} -he ${HE} \
                             -f \
                             -b vref \
                             -p vref_flat \
                             -g generator/eri_generator \
                             -d ../eri/vref_flat

touch eri/CMakeLists.txt
