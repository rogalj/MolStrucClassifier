# Example files

<p align="justify"> 
  
Here we have urea as sample data. 
In the data file you need to have the following format. 
```
atom style: full2
atom ID, molID, atom_type, charge, sym_type, x, y, z
```
This is the sample line for Atoms information. 
```
1 1 5 -0.326615 1 2.785 7.001 1.945
```


When you make the sym_type, you need to have a unique number for the one you use for the symmetry function calculation. For others, you can use the same number. In this example, I used 0 for atoms that I did not use for the calculation.
