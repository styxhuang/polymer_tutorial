gmx_mpi grompp -f npt.mdp -c eql.gro -p init.top -o npt-1 
gmx_mpi mdrun -deffnm npt-1 -v -rdd 0.5
gmx_mpi grompp -f npt-x.mdp -c npt-1.gro -p init.top -o y-x
gmx_mpi mdrun -deffnm y-x -v -rdd 0.5

# KEEP WRITING ANY COMMANDS YOU WANT BASH EXECUTE
