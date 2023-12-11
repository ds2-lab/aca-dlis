#!/bin/bash

#APEX evaluation
#./run.sh longitudes-200M.bin.data binary double 100000000 200000000 search 1 apex 1 0 > apex-longitudes-search.data
#./run.sh longitudes-200M.bin.data binary double 2000000 200000000 insert 1 apex 1 1 > apex-longitudes-insert.data
./run.sh ./datasets/longitudes-200M.bin.data binary double 1000000 200000000 mixed 0 apex 1 0.1 0.5 2
