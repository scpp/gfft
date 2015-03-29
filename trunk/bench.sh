#!/bin/bash

declare -a tharr=("1" "2" "3" "4")
declare -a parr=("1" "2" "3" "4" "5" "6")


echo "============ Powers of 2 ================"

cd ./build

#for n in "${parr[@]}"
for n in `seq 8 10`
do
    cmake -H.. -B. -DNUM:STRING="(1<<${n})" .; time make
    src/gfft
    read -p "Press [Enter] key to continue..."
done
