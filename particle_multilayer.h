typedef struct {
  double x;
  double y;
  double z;
#ifdef ADD_PART_MEM
  ADD_PART_MEM
#endif
} particle;

typedef particle * particles; //omit double pointer types and mistakes
typedef int Particles;        //An index reference
particles * pl = NULL;        
long unsigned int * pn = NULL, * pna =  NULL; // Particles and space
long unsigned int terminate_int = ULONG_MAX;

/* DRAMATIS PERSONAE
  particle - struct containing position of one particle
  particles - array of several particle structs

  pl  - Pointer to each list of particles
  pn  - array of each particle list length
  pna - Not sure

*/


//Creates n new particles, iterates array of particles until end of preexisting array
Particles new_particles (long unsigned int n) {
  int j = 1;
  if (pn != NULL) 
    while (pn[j++] != terminate_int);
  //Reallocating memory to the particle arrays to fit new particles.
  pn = realloc (pn, (j + 1)*sizeof(long unsigned int));
  pna = realloc (pna, (j + 1)*sizeof(long unsigned int));
  pl = realloc (pl, j*sizeof(particles));
  pn[j] = terminate_int;
  pn[--j] = n;
  pna[j] = n;
  particle * pt = calloc (n, sizeof(particle));
  pl[j] = pt;
  return j;
}



Particles * add_Particles_to_list (Particles p, Particles * plist) {
  int l = 0, t = 0;
  if (plist != NULL) {
    while (pn[l] != terminate_int) {
      if (l == plist[t]) 
	t++;
      l++;
    }
  }
  plist = realloc (plist, (t + 1)*sizeof(Particles));
  plist[t] = p;
  printf ("%d\n", t);
  return plist;
}

// Utility function increasing number of slots for particles in list
void change_plist_size (Particles Plist, int n) {
  long unsigned int n_new = n + pn[Plist];
  if (n_new > pna[Plist] || 2.5*(n_new + 1) < pna[Plist]) {
    pna[Plist] = 2*n_new;
    pl[Plist] = realloc (pl[Plist] , pna[Plist]*sizeof(particle));
  }
  pn[Plist] += n;
}


// Adds particle p to particle list with index Plist. 
int add_particle (particle p, Particles Plist) {
#if _MPI
  if (!(locate (p.x, p.y, p.z).level > 0)) //Only place particle where it can be found
    return 0;
#endif
  
  change_plist_size (Plist, 1); //Increases number of particles in Plist
  pl[Plist][pn[Plist] - 1] = p; //append data
  return 1;
}

// Removes n particles from Plist
void remove_particles_index (Particles Plist, int n, long unsigned int ind[n]) {
  int m = 0, j = 0;
  //While j is less than the new number of particles in Plist
  while (j < pn[Plist] - n) {
    //Increment m until m + j == ind[m]. Allows to select certain elements for deletion
    while (m < n ? j + m == ind[m] : 0){
      m++;
  }
    //If m is less than n, then
      // iterate until j >=  the number of elements to be removed from Plist or  j + m is not equal to ind[m]
    //Else
      // Iterate until j < the number of elements to be removed
    while (m < n ? j < pn[Plist] - n && j + m != ind[m] : j < pn[Plist] - n) {
      pl[Plist][j] = pl[Plist][j + m];
      j++;
    }
  }
  // Reduce size of Plist
  change_plist_size(Plist, -n);
}

// Cleanup all particles
void free_p (void) {
  int j = 0;
  while (pn[j++] != terminate_int) {
    free (pl[j - 1]);
    pl[j -1] = NULL;
  }
  free (pl);
  pl = NULL;
  free (pn);
  pn = NULL;
  free (pna);
  pna = NULL;
}

event free_particles (t = end)
  free_p();



// Creating iterators using qcc magic

//Defining p as a substitute for pl[l][j], where l is the list id and j the particle index.
#define p() pl[l][j] 

// Defining the particle variable substitute
@def PARTICLE_VARIABLES
  double x = pl[l][j].x; NOT_UNUSED(x);
  double y = pl[l][j].y; NOT_UNUSED(y);
  double z = pl[l][j].z; NOT_UNUSED(z);
@

// Defining the foreach_particle iterator, which iterates over each particle and defines their x,y,z coordinates
@def foreach_particle() {
  int l = 0;
  while (pn[l] != terminate_int) {
    for (int j = 0; j < pn[l]; j++) {
      PARTICLE_VARIABLES
@
@def end_foreach_particle()
    }
    l++;
  }
}
@

// Iterator looping over each particle in given particle list
@def foreach_particle_in(PARTICLES) {
  int l = PARTICLES;
  for (int j = 0; j < pn[l]; j++) {
    PARTICLE_VARIABLES
@
@def end_foreach_particle_in()
  }
}
@

@def foreach_P_in_list(PARTICLES_LIST) {
  int lll = 0, ll = 0;
  while (pn[lll] != terminate_int) {
    if (lll == PARTICLES_LIST[ll]) {
	Particles P = lll; NOT_UNUSED(P);
@
@def end_foreach_P_in_list()
        ll++;
    }
    lll++;
  }
}
@

// Remove single particle from Plist
int remove_particle (particle p, Particles Plist) {
  long unsigned int * ind = malloc (sizeof(long int));
  int i = 0, a = 1; 
  foreach_particle_in(Plist) {
    foreach_dimension()
      if (p().x != p.x)
	continue;
    if (i + 1 >= a) {
      ind = realloc (ind, a*2*sizeof(long int));
      a *= 2;
    }
    ind[i++] = j;
  }
  remove_particles_index (Plist, i, ind);
  free (ind); ind = NULL;
  return i;
}


#define remove_particles(p,CONDITION) do {		\
  int nr = 0;						\
  foreach_particle_in ((p)) {				\
    if ((CONDITION))					\
      nr++;						\
  }							\
  long unsigned int ind[nr];				\
  int i = 0;						\
  foreach_particle_in ((p)) {				\
    if ((CONDITION))					\
      ind[i++] = j;					\
  }							\
  remove_particles_index ((p), nr, ind);		\
  } while (0)

// Z0 removed, as particle is assumed to not be able to leave ocean
void particle_boundary (Particles P) {
  coord mind = {X0, Y0}; 
  foreach_particle_in(P) { 
    foreach_dimension() {
      if (p().x < mind.x) 
	p().x += L0;
      else if (p().x > (mind.x + L0))
	p().x -= L0;
    }
  }
}
void place_in_cells (Particles P) {
  particles pp = pl[P];                 
  long unsigned int np = 0;
  foreach(serial, reduction(+:np)) {
    double z = zb[];
    foreach_layer(){
      z = z + h[]/2;
      coord loc = {x, y,z};
      foreach_dimension()
        pp[np].x = loc.x;
    
      pp[np].z = loc.z;
      np++;
      z = z + h[]/2;
    }
  }
}

// Structure initializing particles within geometric constraints origin in (xm, ym). length scale l. 
struct Init_P {
  int n; //Num particles
  double xm; 
  double ym; 
  double l;  
  double zp; // Deployment depth
};

void place_in_square (Particles p, struct Init_P inp) {
  //Set number of particles to place to ten unless otherwise specified
  if (!inp.n)
    inp.n = 10;
  // Set width of square to L0 unless specified otherwise
  if (!inp.l)
    inp.l = L0;
  
  //Iterator
  long unsigned int j = 0;
  // Fetching list of particles
  particles pp = pl[p];

  //Distance between particles
  double dx = inp.l/(double)inp.n;

  // Location of first particle
  double x0 = inp.xm - inp.l/2. + dx/2;
  double y0 = inp.ym - inp.l/2. + dx/2;
  for (double xp = x0; xp < inp.xm + inp.l/2.; xp += dx) {
    for (double yp = y0; yp < inp.ym + inp.l/2.; yp += dx) {
      pp[j].x = xp;
      pp[j].y = yp;
      pp[j++].z = inp.zp;
    }
  }
}

void place_in_circle (Particles p, struct Init_P inp) {
  particles pp = pl[p];
  for (int j = 0; j < inp.n; j++) {
    double angle = noise()*pi;
    double R = sqrt(fabs(noise()))*inp.l;
    pp[j].x = R*sin(angle) + inp.xm;
    pp[j].y = R*cos(angle) + inp.ym;
    pp[j].z = inp.zp;
  }
}


// particle statistics
typedef struct {
  coord avg;
  coord min;
  coord max;
  coord stddev;
} pstats;

// Function computing the statistics of the deployed Particles P
pstats statsp (Particles P) {
  // Initializing coordinates
  coord avg = {0 ,0, 0}, min = {0,0,0},
    max = {0, 0, 0}, stddev = {0, 0, 0};
  
  // 
  foreach_dimension(2) {
    avg.x = stddev.x = 0;
    min.x = HUGE;
    max.x = -HUGE;
  }
  avg.z = stddev.z = 0;
  min.z = HUGE;
  max.z = -HUGE;
  long unsigned int np = pn[P];
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, &np, 1, MPI_UNSIGNED_LONG,
		 MPI_SUM, MPI_COMM_WORLD);
#endif
  if (np) {
    foreach_dimension() { //reduction of coord members does not work
      double avgx = 0;
      double minx = HUGE;
      double maxx = -HUGE;
      foreach_particle_in(P, reduction(+:avgx) reduction(max:maxx)
			  reduction (min:minx)) {
	      avgx += p().x;
	      if (p().x < minx)
	        minx = p().x;
	      if (p().x > maxx)
	        maxx = p().x;
      }
      avg.x = avgx/(double)np;
      min.x = minx;
      max.x = maxx;
    }
    double avgz = 0;
    double minz = HUGE; 
    double maxz = -HUGE;
    foreach_particle_in(P, reduction(+:avgz) reduction(max:maxz) reduction (min:minz)) {
        avgz += p().z;
        if (p().z < minz){
          minz = p().z;
        }
        if (p().z > maxz){
          maxz = p().z;
        }
    }
    avg.z = avgz/(double)np;
    min.z = minz;
    max.z = maxz;

    foreach_dimension() {
      double stddevx = 0;
      foreach_particle_in(P, reduction(+:stddevx)) {
	      stddevx += sq(avg.x - p().x);
      }
      stddev.x = sqrt(stddevx/(double)np);
    }
    double stddevz = 0;
    foreach_particle_in(P, reduction(+:stddevz)) {
	      stddevz += sq(avg.z - p().z);
    }
    stddev.z = sqrt(stddevz/(double)np);
  }
  pstats s;
  s.max = max, s.min = min, s.avg = avg, s.stddev = stddev;
  return s;
}


// Does this work without changes, investigate?
void particle_pdf (Particles P, scalar s) {
  foreach()
    s[] = 0;
  particle_boundary (P);
  foreach_particle_in(P) {
    Point point = locate (x, y, z);
    if (point.level > 0)
      s[]++;
  }
  foreach()
    s[] /= (pn[P]*dv());
  boundary ({s});
}


void random_step (Particles P, double step) {
  foreach_particle_in(P) {
    double theta = noise()*pi;
#if (dimension == 1)
    coord d = {sign(theta)};
#elif (dimension < 3)
    coord d = {sin(theta), cos(theta)};
    #if (nl > 1)
      double phi = acos(noise());
      coord d = {sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi)};  
    #endif
#else
    double phi = acos(noise());
    coord d = {sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi)};
#endif
    foreach_dimension()
      p().x += d.x*step;
  }
}


double plength (Particles P) {
  double lr = 0;
  int np = pn[P];
  particles pli = pl[P];
  foreach_particle_in(P) {
    int il = j - 1 < 0 ? np - 1: j - 1;
    int ir = j + 1 >= np ? 0 : j + 1;
    coord pl = (coord){pli[il].x, pli[il].y, pli[il].z};
    coord pm = (coord){pli[j].x, pli[j].y, pli[j].z};
    coord pr = (coord){pli[ir].x, pli[ir].y, pli[ir].z};
    coord bv = {0, 0, 0}, a1 = {0,0,0};
    coord a2 = {0,0,0};
    double bl = 0, ka = 0;
    foreach_dimension() {
      a1.x = pm.x - pl.x;
      bv.x = pr.x - pl.x;
      bl += sq(bv.x);
    }
    a1.z = pm.z - pl.z;
    bv.z = pr.z - pl.z;
    bl += sq(bv.z);
    
    bl = sqrt(bl);
    foreach_dimension(){
      ka += a1.x*bv.x/bl;
    }
    ka += a1.z*bv.z/bl;
    foreach_dimension()
      a2.x = a1.x - bv.x*ka/bl;
    a2.z = a1.z - bv.z*ka/bl;
    normalize (&a2);

    double al = 0, am = 0, ar = 0;
    foreach_dimension() {
      al -= a2.x*pl.x;
      am -= a2.x*pm.x;
      ar -= a2.x*pr.x;
    }
    al -= a2.z*pl.z; am -= a2.z*pm.z; ar -= a2.z*pr.z;

    double C = fabs(am - al);
    double xt = 0;
    double b = 0;
    
    foreach_dimension() {
      xt += sq(pl.x - pm.x);
      b += sq(pl.x - pr.x);
    }
    xt += sq(pl.z - pm.z); b += sq(pl.z - pr.z);
    xt = sqrt (xt - sq(C));
    b = sqrt(b);
    double A =  C/((xt*(b - xt)));
    double xp1 = 0.5*b*(1 - 1/sqrt(3));
    double xp2 = 0.5*b*(1 + 1/sqrt(3));
    double l1 = b*sqrt(1. + sq(-2*A*xp1 + A*b));
    double l2 = b*sqrt(1. + sq(-2*A*xp2 + A*b));
    lr += (l1 + l2)/4;
  }
  return lr;
}

struct DumpP {
  char * fname;   //File name       Default: "pdump"
  Particles * td; //Particle lists  Default: all
  FILE * fp;      //File pointer    Default: Not used
  bool Dmem;      //Member data     Default: false, only x,y,z
  bool print;     //Show dump data  Default: false
};

  
bool pdump (struct DumpP p) {
  // The file:
  FILE * fp = p.fp;
  char def[] = "pdump", * file = p.fname ? p.fname : p.fp ? NULL : def;
  if (file) {
    if ((fp = fopen (file, "wb")) == NULL) {
      perror (file);
      exit (1);
    }
  }
  assert (fp);
  
  // Get particle lists data
  Particles * td = p.td;
  int j = 0, n = 0;        
  if (td) {
    foreach_P_in_list (td) 
      j++;
  } else  {// all
    while (pn[j++] != terminate_int);
    td = malloc (sizeof(int)*j);
    j = 0;
    while (pn[j] != terminate_int) { 
      td[j] = j;
      j++;
    }
  }

  // Nr of particles
  long unsigned int np[j];
  foreach_P_in_list (td) {
    np[n] = pn[P];
    n++;
  }
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, np, j, MPI_UNSIGNED_LONG, MPI_SUM,
		 MPI_COMM_WORLD);
#endif

  // Dump members?
  bool Dmem = true;
  if (!p.Dmem)
    Dmem = false ;

  // Print header
  if (pid() == 0) {
    fwrite (&j, sizeof(int), 1, fp);
    fwrite (np, sizeof(long unsigned int), j, fp);
    fwrite (&Dmem, sizeof(bool), 1, fp);
  }
  int Headersize = sizeof(int) + j*sizeof (long unsigned int) + sizeof(bool); 
  fseek (fp, Headersize, SEEK_SET); //offset

  if (p.print && pid() == 0) 
    fputs ("# P\tamount\n", stderr);

  // Print particle data
  for (int m = 0; m < j; m++) { //Nesting within foreach_P... did not work
    Particles P = td[m];
#if _MPI
    int size = Dmem ? sizeof(particle) : 3*sizeof(double);
    long unsigned int outa[npe()], outat[npe()  + 1];
    outat[0] = 0;
    MPI_Allgather (&pn[td[m]], 1, MPI_UNSIGNED_LONG,
		   &outa[0], 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    //Compute and set displacements
    for (int j = 1; j <= npe(); j++) 
      outat[j] = outa[j - 1] + outat[j - 1];
    fseek (fp, outat[pid()]*size, SEEK_CUR);
#endif
    foreach_particle_in (P) {
      if (Dmem) {
	fwrite (&p(), sizeof(particle), 1, fp);
      } else {
	double c[3] = {x, y, z};
	fwrite (c, sizeof(double), 3, fp);
      }
    }
    if (p.print && pid() == 0) 
      fprintf (stderr, "# %d\t%ld\n", m, np[m]);
#if _MPI
      fseek (fp, (np[m] - outat[pid() + 1])*size, SEEK_CUR);
#endif
  }
  fclose (fp);
  if (!p.td)
    free (td);
  return true; 
}

int prestore (struct DumpP p) {
  // Open file
  FILE * fp = p.fp;
  char * file = p.fname;
  if (file && (fp = fopen (file, "r")) == NULL)
    return 0;
  assert (fp);

  //read nr. of lists and nr. of particles
  int j;
  bool Dmem;
  fread (&j, sizeof(int), 1, fp);
  long unsigned int np[j];
  fread (np, sizeof(long unsigned int), j, fp);
  fread (&Dmem, sizeof(bool), 1, fp);
  
  // Print some reference data
  if (p.print) {
    if (pid() == 0) {
      if (Dmem)
	fputs ("# Restoring members...\n", stderr);
      fputs ("# P\tamount\n", stderr);
      for (int m = 0; m < j; m++)
	fprintf (stderr, "# %d\t%ld\n", m, np[m]);
    }
  }
  
  // Remove existing particles
  if (pl != NULL)
    free_p();

  // Allocate and read
  Particles P;
  for (int m = 0; m < j; m++) {
    if (pid() == 0) { // This could be more balanced
      P = new_particles (np[m]);
      if (Dmem) 
	fread (pl[m], sizeof(particle), np[m], fp);
      else {
	foreach_particle_in (P) {
	  double c[3];
	  fread (&c, sizeof(double), 3, fp);
	  p().x = c[0]; p().y = c[1]; p().z = c[2];
	}
      }
    } else // slaves do not read data
      P = new_particles (0);
     // Apply boundary conditions.
    particle_boundary (P);
  }
  if (file)
    fclose (fp);
  return j;
}
