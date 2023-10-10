# Lars WD Sandbox
I am a PhD student at the University of Oslo and am currently working on circulation models in fjord environments. Scripts and functions will be added here as the project progress.

Github mirror at:
[https://github.com/larswd/basiliskSandbox/](https://github.com/larswd/basiliskSandbox/)

Current files in sandbox: 

## Save multilayer fields to paraview file.

The function in [output_vts.h](output_vts.h) is a function which stores a set of variables of interest to a vts file which can then be viewed using paraview (or other software which can open paraview files). Sample function call in a basilisk event:

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
### Download:
```bash
wget http://basilisk.fr/sandbox/larswd/output_vts.h?raw -O output_vts.h
```

## Multilayer particle tracker. 
Based on the ```particle.h``` extension of Antoon van Hooft, which can be found [here](http://basilisk.fr/sandbox/Antoonvh/particle.h). The code and algorithms remain largely the same, but they are modified to ensure the particles are placed and tracked correctly in the vertical. This extension should work identically to that of Antoon. 
Code: [particle_multilayer.h](particle_multilayer.h).

### Current issues
Fix vertical position
Check the pdf-function

### Download
```bash
wget http://basilisk.fr/sandbox/larswd/particle_multilayer.h?raw -O particle_multilayer.h
```
