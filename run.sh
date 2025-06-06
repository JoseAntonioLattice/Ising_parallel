#/bin/bash
a=$(($1*$2))

{ echo input/parameters.nml; echo $1 $2; } | cafrun -n $a build/ising_parallel   
