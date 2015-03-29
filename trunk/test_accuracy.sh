#!/bin/bash

#declare -a tharr=("1" "2" "3" "4")
declare -a    arrn=("2" "3" "5" "7")
declare -a arrpmin=("2" "2" "2" "2")
declare -a arrpmax=("11" "7" "4" "3")
declare -a arrtype=("DOUBLE" "FLOAT" "COMPLEX_DOUBLE" "COMPLEX_FLOAT")
declare -a arrplace=("IN_PLACE" "OUT_OF_PLACE")


cd ./build/test
rm -f gfftacc.log

for t in "${arrtype[@]}"
do
echo ">>>>> Value type: ${t}"
for p in "${arrplace[@]}"
do
echo ">>>>> Transform mode: ${p}"
for i in `seq 0 3`
do
    cmake -H../../test -B. -DPNUM:STRING="${arrn[$i]}" -DPMIN:STRING="${arrpmin[$i]}" -DPMAX:STRING="${arrpmax[$i]}" -DPLACE:STRING="${p}" -DTYPE:STRING="${t}" -DFULLOUTPUT:STRING="0" ../../test
    time make
    ./gfft_accuracy >> gfftacc.log
    echo "${arrn[$i]}^[${arrpmin[$i]}-${arrpmax[$i]}]"
    #read -p "Press [Enter] key to continue..."
done
done
done
