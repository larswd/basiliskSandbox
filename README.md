# Lars WD Sandbox
I am a PhD student at the University of Oslo and am currently working on circulation models in fjord environments. Scripts and functions will be added here as the project progress.

Link to sandbox:
[basilisk.fr/sandbox/larswd/README](basilisk.fr/sandbox/larswd/README) 
## Save multilayer fields to paraview file.

The function in [output_vts.h](basilisk.fr/sandbox/larswd/output_vts.h) is a function which stores a set of variables of interest to a vts file which can then be viewed using paraview (or other software which can open paraview files). Sample function call in a basilisk event:

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
