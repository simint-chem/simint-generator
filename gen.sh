#!/bin/bash


MAXHRR=3
MAXVRR=12
VE=1000
HE=1000
P=1000

mkdir -p ../eri/FO.gen
python3 ../generate_twoel.py -l 2 -ve ${VE} -he ${HE} -P ${P} \
                             -b FO \
                             -p FO \
                             -g generator/eri_generator \
                             -d ../eri/FO.gen


mkdir -p ../eri/FO_flat.gen
python3 ../generate_twoel.py -l 2 -ve ${VE} -he ${HE} -P ${P} \
                             -f \
                             -b FO \
                             -p FO_flat \
                             -g generator/eri_generator \
                             -d ../eri/FO_flat.gen


#mkdir -p ../eri/split.gen
#python3 ../generate_twoel.py -l 2 -ve ${VE} -he ${HE} \
#                             -b split \
#                             -p split \
#                             -g generator/generator -he ${HE} \
#                             -d ../eri/split.gen


mkdir -p ../eri/vref.gen
python3 ../generate_twoel.py -l 3 -ve ${VE} -he ${HE} -P ${P} \
                             -b vref \
                             -p vref \
                             -g generator/eri_generator \
                             -d ../eri/vref.gen


mkdir -p ../eri/vref_flat.gen
python3 ../generate_twoel.py -l 3 -ve ${VE} -he ${HE} -P ${P} \
                             -f \
                             -b vref \
                             -p vref_flat \
                             -g generator/eri_generator \
                             -d ../eri/vref_flat.gen


#mkdir -p ../eri/vrr.gen ../eri/hrr.gen
#generator/vrr_generator -o ../eri/vrr.gen -L ${MAXVRR}
#generator/hrr_generator -o ../eri/hrr.gen -L ${MAXHRR}


touch eri/CMakeLists.txt
