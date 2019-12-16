~/.openmpi/bin/mpic++ main.cpp -o q1
~/.openmpi/bin/mpirun -np $1 q1 < inp1.txt