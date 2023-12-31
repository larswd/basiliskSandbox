
/**
# Brownian motion

A hundred million particles take a random walk. We compare their
behaviour as a group against the diffusion of a scalar field.

~~~gnuplot The statistics of a random walk appear as diffusion
set xlabel 'x'
set ylabel 'pdf'
set grid 
plot 'log' w l lw 2 t 'Initial Particle PDF' ,\
     'log' u 1:3 t 'Scalar (t = 0)' ,\
     'out' w l lw 2 t 'Final particle PDF' ,\
     'out' u 1:3 t 'Scalar (t = t_{end})'
~~~
*/

#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "particle.h"
#include "diffusion.h"

scalar analytical[];
Particles Brownian;

int np = 1e8;
double depth = 0.1;
double Diff = 0.25, tend = 5;

int main() {
  periodic (left);
  L0 = 12.0;
  X0 = -L0/3.;
  N  = 128;
  nl = 1;
  run();
}

event init (t = 0) {
  foreach(){
    zb[] = -depth;
    foreach_layer(){
      h[] = -zb[]/nl;
    }
  }
  DT = 0.05;
  Brownian = new_particles (np);
  analytical = new scalar[nl];
  /** 
Initially, the particles are placed `normally', following the
[Box--Muller
transform](https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform).
  */
  foreach_particle () {
    double a = rand()/(double)(RAND_MAX);
    double b = rand()/(double)(RAND_MAX);
    p().x = sqrt(-2*log(a))*cos(2*pi*b);
  }
  scalar s = new scalar[nl];
  particle_pdf (Brownian, s);
  foreach() {
    analytical[] = exp(-sq(x)/2.)/sqrt(2*pi);
    fprintf (stderr, "%g %g %g\n", x, s[], analytical[]);
  }
}

event integrate (i++) {
  dt = dtnext (DT);
  
  const face vector mu[] = {Diff, Diff, Diff};
  diffusion (analytical, dt, mu);
    /**
       A random walk is performed with Einstein's step size.
    */
  random_step (Brownian, sqrt(2.*Diff*dt));
}

event stop (t = tend) {
  particle_boundary (Brownian); //BCs are applied *after* the simulation
  scalar s = new scalar[nl];
  particle_pdf (Brownian, s);
  foreach()   
    printf ("%g %g %g\n", x, s[], analytical[]);
  return 1;
}
