~~in implicit assembly.cpp remove the empty boundary tag check, try to move the check to main~~

~~ in Grid.cpp add a throw error for unsupported .msh file formats~~

~~Write a unique check of .json file for every application~~

* Right now, boundary condition could only support constant boundary case, for variable boundary,
  we have to create an array of boundary values. Take this up later on.

* Right now, Cells and Faces are stored in OOD way, See if it can be refactored to DOD way

