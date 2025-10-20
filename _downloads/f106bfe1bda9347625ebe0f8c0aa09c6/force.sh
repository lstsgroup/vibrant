#!/bin/bash

natom=125

for sign in 0.001 -0.001; do                  # displacement
  for i in $(seq 0 $(($natom - 1))); do       # atoms
    for j in $(seq 0 2); do                   # dimensions
      cd FeTMP.i_atom_${i}.i_coord_${j}.displ_${sign}
      force_x=$(grep -P -A "$natom" "#" force.xyz | tail -n "$natom" | awk '{print $4, $5, $6}')
      printf "%-30s\n" "$force_x" >> ../FeTMP-force.data3
      cd ..
    done
  done
done
