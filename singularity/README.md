# Creating a Singularity image with TMplot

This is not quite ready yet, but to build the Singularity image run

     singularity build --sandbox tmpsandbox Singularity
     singularity build TMplot.img tmpsandbox

The sandbox dir _tmpsandbox_ is just used to circumvent the problem "No space left on device".
It can be removed afterwards.





