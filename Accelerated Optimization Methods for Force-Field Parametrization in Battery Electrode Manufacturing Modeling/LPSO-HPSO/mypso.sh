#!/bin/sh

echo "Cleaning"
rm a_analyser.txt
./clean.sh 
sleep 1

echo "Compiling"
./compile.sh
sleep 1

echo "Launching PSO"
./lpso &> a_analyser.txt
