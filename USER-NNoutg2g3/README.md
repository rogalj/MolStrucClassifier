# Calculate Neural Network using symmetry function

<p align="justify"> 
  
  ## Input
  
  In this directory, we have a compute command that can calculate the Neural Network using symmetry function. 
  In addition to run file, you need INPUT.dat. Check the sample code carefully and run the command. 
  This command currently takes two point-vector respresntations and works only for urea. 
  For the general use, check <a class="reference external" href="https://github.com/rogalj/MolStrucClassifier/tree/main/USER-NNout">USER-NNout directory</a>. 
  Unlike USER-NNout directory, this compute can check the misclassified molecule and modify the classification result as undefined and save the output in different file. 

```bash
    compute  ID  groupID  nnoutg2g3gen/atom  center_atom_type cutoff_distance total_symmetry_function_number  Prob_of_correct_local_environment_needed
```
This is the sample line 
```bash
    compute CNN  all nnoutg2g3gen/atom 3 24 0.4
```

## Output

After you run the code, you will get g2g3output-0.dat, g2g3good-0.dat, g2g3error-0.dat as output. 
These files contains symmetry function results for each molecule. output contains all information, good contains correctly classified information, error contains misclassified information.



</p>

