#!/bin/bash

N=10

filename=ts_runner

if [[ "$(uname)" == "Windows_NT" ]];
    then
	    EXE=$filename.exe
    else
	    EXE=./$filename
fi

CC=g++
CFLAGS="-Wshadow -Wno-unused-result -Wno-sign-compare -Wno-char-subscripts -std=c++20 -Wall -Wextra -g -DPIZZA"

echo $CC
echo $CFLAGS
echo $EXE

opt_values=("23168" "27953" "3025" "4801" "6785" "28087" "9497" "6568" "38773" "19147" "201052" "222892" "355666" "268150" "201418" "294666" "165556" "153568" "55143" "39258")

for i in $(seq 5 $N); do
   echo "running instance $i"
   echo "/home/lucas/TCC/data/lck_instances/" "inst$i.txt" LKM ${opt_values[$i - 1]}
   
   $CC $CFLAGS -o $EXE $filename.cpp && $EXE "/home/lucas/TCC/data/lck_instances/" "inst$i.txt" LKM ${opt_values[$i - 1]}
   
   echo "instance $i finished"
done
