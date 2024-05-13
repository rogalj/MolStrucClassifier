.. _Moltemplate:

How to prepare files using moltemplate
=======================================


How to Use moltemplate
----------------------

Below are the basic steps for working with Moltemplate. I'll explain one of the examples below:

If you haven't worked with Moltemplate recently, get system.lt and urea.lt to the test directory.

Run the following command for lammps run:

.. code-block:: console

   moltemplate.sh system.lt

You can get rid of unnecessary lines in the file by typing the following command:

.. code-block:: console

   cleanup_moltemplate.sh

You will get these files after running the previous commands:



.. table::
   :align: left

   +----------------+----------------------+
   | system.in      |lammps input file     |
   +----------------+----------------------+
   | system.data    |lammps data file      |
   +----------------+----------------------+
   | system.in.init |force field functions |
   +----------------+----------------------+

If you finish the urea sample, you can continue with the actual data.

Step 1
^^^^^^

Get the file from ovito with mass, charges, and bond information (in the format of a .data file). If you have other information, you can delete it now. Then type:

.. code-block:: console

   ltemplyfy.sh ovitofilename > outputfilename     (ltemplyfy.sh urea.data > urea.lt)

For example:

.. code-block:: console

   ltemplyfy.sh urea.data > urea.lt



Step 2
^^^^^^

Check the urea.lt file or the output file you made. You may have to modify several things. If you don't see a part that I will mention, skip that part.

1. Add import "gaff.lt" before write_once("In Init") {.

2. Add Urea inherits GAFF { on the next line after the import "gaff.lt".

3. Add another closed curly bracket } at the end of the file.

4. If you have a unique name for the atom type that is not used in GAFF, make sure you change the name in the Data Atoms section. For example, if you have n1 and n2, you can keep the name in Data Masses, but you have to have the right name such as n in Data atom section. Check useful commands shown below.

5. Remove bond type (such as @bond:type1).

6. Change "Data Bonds" to "Data Bond List".

7. Check the atom information. If it says type##, you have to change the name of the type to the corresponding one. If you change this, you also have to change @atom: to the new name. Also you need to change atom in Data Bond List.

8. In Data Atoms section, if you don't have $mol section or the value for this is the same for the whole system, you have to check the input data. Most likely you need to go to Ovito and assign molecule data. Basically, you need to do clustering and go to compute property > molecule identifier > cluster. More details are in How to use Ovito.

9. If you see angle data or some others after Bond information, remove that information now.

10. module load anaconda3/2020.07

11. Change system.lt's first line to whatever the input file, change the cell box size.

12. Run:

.. code-block::

   moltemplate.sh system.lt

13. Check system.data to make sure it ran.

14. Clean up the code by typing:

.. code-block::

   cleanup_moltemplate.sh


15. If you want to change the atom type, do it now. Make sure you change it to the right name.

16. Change the n2 in the Atom List to the corresponding name.

17. Use this file as a standard.

18. When you use it in LAMMPS, you need to have system.data, system.in.init, and system.in.settings in the running directory.

19. You need to modify urea.in.settings if you change the atom information. For instance, if you added n2, you should have a line for that too.

20. Also, change the atom type from 4 to 5 for urea in system.data.


Useful Commands When Using Moltemplate in Vim:


.. code-block::

   :/@bond:type3           <-- highlight the words type3.

.. code-block::

   :%s///g                 <-- delete the highlighted words.

.. code-block::

   %s/@atom:n1 /@atom:n /g	<-- replace all @atom:n1 in the file to @atom:n



