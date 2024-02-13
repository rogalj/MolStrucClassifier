# Input and output symmetry type using LAMMPS

<p align="justify"> 
  
  ## Input
  
  In this directory, we have a new atom style called 'full2'. 'full2' is very similar to 'full' atom style, which stores molecular parameters(bonds, angle, dihedrals, and impropers) and charge. In addition to those parameters, 'full2' stores symmetry function ID. If you want to use this atom_style, you need to have the following data file in the specified order. 

```bash
  Atom_ID  Molecule_ID  Atom_type  Charge  Sym-type  pos image_data(optional)
```

## Output

Dump_custom file can output the dump file with symmetry function ID. In order to output, you need the following line in the run file. traj.dump will be the name of the output file. 

```bash
   dump              d1 all custom 1000 traj.dump id mol element sym_type x y z 
```


</p>

