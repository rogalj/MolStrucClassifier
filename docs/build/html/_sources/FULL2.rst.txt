.. _FULL2:

LAMMPS atom style (USER-ATOM-VEC-FULL2)
=========================================


USER-ATOM-VEC-FULL2
^^^^^^^^^^^^^^^^^^^


Here we introduce a new atom style in LAMMPS called 'full2'. 'full2' is similar to the 'full' atom style, which stores molecular parameters (bonds, angles, dihedrals, and impropers) and charge. Additionally, 'full2' stores symmetry function IDs. If you wish to use this atom_style, you need to have the following data file in the specified order:

.. code-block:: console

     Atom_ID  Molecule_ID  Atom_type  Charge  Sym-type  pos image_data(optional)

Output
""""""

The 'dump_custom' command can output the dump file with symmetry function IDs. To enable this output, include the following line in the run file, where 'traj.dump' will be the name of the output file:

.. code-block:: console

   dump              d1 all custom 1000 traj.dump id mol element sym_type x y z 


