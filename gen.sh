#!/bin/bash


MAXHRR=2
MAXVRR=12
VE=1000
HE=1000
P=1000
STACK=2048

CPUFILE=$1

if [ $# -ne 1 ]
then
  echo "Give CPUInfo file as argument"
  exit 1
fi

if [ ! -f ${CPUFILE} ]
then
  echo "Cannot find CPUFile: ${CPUFILE}"
  exit 1
fi

mkdir -p ../eri/gen
rm -Rf ../eri/gen/*

# SIMD, with intrinsics
python3 ../generate_twoel.py -l 2 -ve ${VE} -he ${HE} -P ${P} \
                             -c ${CPUFILE} \
                             -i \
                             -s ${STACK} \
                             -b FO \
                             -g generator \


# Reference boys function
#mkdir -p ../eri/vref.gen
#python3 ../generate_twoel.py -l 2 -ve ${VE} -he ${HE} -P ${P} \
#                             -s ${STACK} \
#                             -S \
#                             -b vref \
#                             -g generator \


#mkdir -p ../eri/FO_flat.gen
#python3 ../generate_twoel.py -l 2 -ve ${VE} -he ${HE} -P ${P} \
#                             -s ${STACK} \
#                             -f \
#                             -b FO \
#                             -g generator \


#mkdir -p ../eri/split.gen
#python3 ../generate_twoel.py -l 2 -ve ${VE} -he ${HE} \
#                             -b split \
#                             -g generator -he ${HE} \




#mkdir -p ../eri/vref_flat.gen
#python3 ../generate_twoel.py -l 3 -ve ${VE} -he ${HE} -P ${P} \
#                             -s ${STACK} \
#                             -f \
#                             -b vref \
#                             -g generator \
#                             -d ../eri/vref_flat.gen


generator/hrr_generator -o ../eri/gen -c ${CPUFILE} -L ${MAXHRR}

touch ../eri/CMakeLists.txt
