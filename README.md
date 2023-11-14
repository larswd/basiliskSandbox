# Lars WD Sandbox
I am a PhD student at the University of Oslo and am currently working on circulation models in fjord environments. Scripts and functions will be added here as the project progress.

Github mirror at:
[https://github.com/larswd/basiliskSandbox/](https://github.com/larswd/basiliskSandbox/)

Contact me at: larswd@math.uio.no

Current files in sandbox: 

## Save multilayer fields to paraview file.

The function in [output_vts.h](output_vts.h) is a function which stores a set of variables of interest to a vts file which can then be viewed using paraview (or other software which can open paraview files). Largely based on a similar function by Ã˜ystein Lande which plots a single layer. Sample function call in a basilisk event:

```C
event output_field (t <= tmax; t += dt)
{
    fprintf(stdout, "field vts output at step: %d, time: %.2f \n", i, t);
    static int j = 0;
    char name[100];
    sprintf(name, "fields/field_%.6i.vts", j++);
    fprintf(stdout, "written to: %s\n", name);
    FILE* fp = fopen(name, "w");
    output_vts_ascii_all_layers(fp, {eta,h,u}, N);
    fclose(fp);
}

```
**Download**:
```bash
wget http://basilisk.fr/sandbox/larswd/output_vts.h?raw -O output_vts.h
```

## Multilayer particle tracker. 
Based on the ```particle.h``` extension of Antoon van Hooft, which can be found [here](http://basilisk.fr/sandbox/Antoonvh/particle.h). The code and algorithms remain largely the same, but they are modified to ensure the particles are placed and tracked correctly in the vertical. This extension is only activated if the multilayer solver is imported, and if this is not the case should be identical to Antoon's ```particle.h```.

Code: 
[particle.h](particle.h).

This library depends on [particle_multilayer.h](particle_multilayer.h) and [particle_classic.h](particle_classic.h) for tracking particles in either a multilayer or non-multilayer setting respectively. 

**Current issues**
- Strange issue when including, but not using, multilayer library. The LAYERS variable is still defined, meaning any uncommented statement importing the multilayer solver, including 

```c
#if 0
  #include "layered/hydro.h"
#endif
```
means the multilayer particle tracker is used instead of the classic particle tracker. Should not cause much trouble except when making comparative test cases. 

- Locate layer algorithm in particle pdf is computationally expensive, and results in a significantly slower execution compared to the non-multilayer case. 

**Test cases and examples**

- [brownian.c](brownian.c) illustrating gaussian drift with multilayer. Based on the test case of Antoon with the same name which can be found [here](http://basilisk.fr/sandbox/Antoonvh/brownian.c).
- [brownian_classic.c](brownian_classic.c) The gaussian drift example of Antoon without multilayer. Code identical [to the original which can be found here.](http://basilisk.fr/sandbox/Antoonvh/brownian.c).

- [tlengt.c](tlengt.c) a test case showing that the placement of particles is implemented correctly when using multilayer.


- [tlengt_classic.c](tlengt_classic.c) a test case showing that the placement of particles is implemented correctly when using multilayer.

- [multilayer_2D_diffusion.c](multilayer_2D_diffusion.c) is an example showing diffusion of particles in both vertical and horizontal direction when using mutilayer.

**Download**
```bash
wget http://basilisk.fr/sandbox/larswd/particle.h?raw -O particle.h
wget http://basilisk.fr/sandbox/larswd/particle_classic.h?raw -O particle_classic.h
wget http://basilisk.fr/sandbox/larswd/particle_multilayer.h?raw -O particle_multilayer.h
```

## Pid.h
I have developed the header file [pid.h](pid.h) to allow for an easy implementation of a PID controller to enforce damping of either the layer thickness $h$ or velocity $u$ at the boundary in the multilayer solver. This is done to enable a sort of "do nothing"-boundary condition in the multilayer solver. 

**Test cases and examples**

- [pid_stokes.c](pid_stokes.c) PID damping of the src example (with some slight modifications) of a [2D stokes wave](http://basilisk.fr/src/test/stokes.c). Notice how the waves are continously dampened until the ocean is nearly at rest. 