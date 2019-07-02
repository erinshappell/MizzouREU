--- 
# Micturition BMTK README 
---

This version of the Nair Lab BMTK model for the Rat LUT was developed by Erin Shappell (NSF NeuroREU Summer 2019). 

## Primary Contributions
* All (known) sources for values and connections have been added via inline comments.
* Hoc templates have been added (but are not complete due to a lack of data on some values) for the following neurons:
  * Hypogastric
  * IMG
  * IND
  * MPG
  * PGN
* A shell script (run.sh) has been added to make running the simulation much easier. Simply type the following line to run all Python scripts needed for the simulation:

```bash

$ ./run.sh

```

* In generate_input.py, new code has been added (but not completed) for creating the input spike train for the bladder afferent. This code gives a start to a more realistic rate increase/decrease due to bladder filling/voiding.
