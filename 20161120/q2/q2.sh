~/.openmpi/bin/mpic++ main.cpp -o q2
~/.openmpi/bin/mpirun -np $1 q2 < inp2.txt