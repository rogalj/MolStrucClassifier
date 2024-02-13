# Calculate symmetry function

<p align="justify"> 
  
  ## Input
  
  In this directory, we have a compute command that can calculate the symmetry function and two more compute commands which is used for symmetry function calculation. 
  In addition to run file, you need INPUT.dat. Check the sample code carefully and run the command. 
  This command currently takes two point-vector respresntations.

```bash
    compute  ID  groupID  chunk/atom molecule
    compute  ID  groupID  vector/chunk cUchunk atomtype1 atomtype2
    compute  ID  groupID  G2G3gen/vecchunk  center_atom_type  ID_of_vector/chunk1  ID_of_vector/chunk2  cutoff_distance total_symmetry_function_number
```
This is the sample line 
```bash
    compute  cUchunk  all chunk/atom molecule
    compute  cCOpos  all  vector/chunk cUchunk 10 9
    compute  cNNpos  all  vector/chunk cUchunk 4 6
    compute  cG2G3  all  G2G3gen/vecchunk  9  cCa2Ca1pos  cCNpos  7.0 24   
```

## Output

After you run the code, you will get g2g3func-0.dat as output. This file contains symmetry function results for each molecule. 



</p>

