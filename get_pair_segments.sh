#!/bin/bash

# Usage: ./get_pair_segments.sh {patient1} {patient2} {segments file}

p1=$1
p2=$2
file=$3

head -1 $file > ${p1}_${p2}.segments
grep "$p1" $file | grep "$p2" >> ${p1}_${p2}.segments