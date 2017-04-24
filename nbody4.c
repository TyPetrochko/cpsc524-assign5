#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define TILESIZE 16

int numThreads = -1;

struct BodySet {
  float *x, *y, *z;
  float *vx, *vy, *vz;
};

struct centerOfMass {
  float x, y, z;
};

struct BodySet createBodySet(int nBodies){
  struct BodySet toReturn;

  toReturn.x = (float *)calloc(sizeof(float), nBodies);
  toReturn.y = (float *)calloc(sizeof(float), nBodies);
  toReturn.z = (float *)calloc(sizeof(float), nBodies);
  toReturn.vx = (float *)calloc(sizeof(float), nBodies);
  toReturn.vy = (float *)calloc(sizeof(float), nBodies);
  toReturn.vz = (float *)calloc(sizeof(float), nBodies);

  return toReturn;
}

void destroyBodySet(struct BodySet bodySet){
  free(bodySet.x);
  free(bodySet.y);
  free(bodySet.z);
  free(bodySet.vx);
  free(bodySet.vy);
  free(bodySet.vz);
}

struct centerOfMass computeCom(const int nBodies, const struct BodySet bodySet){
  float comx = 0.0f, comy=0.0f, comz=0.0f;
  int i;
  #pragma omp parallel for private(i) reduction(+:comx, comy, comz)
  for (i=0; i<nBodies; i++) {
    comx += bodySet.x[i];
    comy += bodySet.y[i];
    comz += bodySet.z[i];
  }
  comx = comx / nBodies;
  comy = comy / nBodies;
  comz = comz / nBodies;
  
  struct centerOfMass toReturn;
  toReturn.x = comx;
  toReturn.y = comy;
  toReturn.z = comz;
  return toReturn;
}

void MoveBodies(const int nBodies, const struct BodySet bodySet, const float dt) {
  int i, j, x;

  // Avoid singularity and interaction with self
  const float softening = 1e-20;

  // Loop over bodies that experience force
  #pragma omp parallel private (i, j, x)
  {
    #pragma omp for schedule(static) // static fastest!
    for (i = 0; i < nBodies; i += TILESIZE) { 

      // Components of the gravity force on body i
      float Fx[TILESIZE] = {0}, Fy[TILESIZE] = {0}, Fz[TILESIZE] = {0}; 
      // Loop over bodies that exert force: vectorization expected here
      #pragma unroll_and_jam(TILESIZE)
      for (j = 0; j < nBodies; j++) { 
        for(x = i; x < i + TILESIZE; x++){
          // Newton's law of universal gravity
          const float dx = bodySet.x[j] - bodySet.x[x];
          const float dy = bodySet.y[j] - bodySet.y[x];
          const float dz = bodySet.z[j] - bodySet.z[x];
          const float drSquared  = dx*dx + dy*dy + dz*dz + softening;
          const float drPower32  = drSquared * sqrtf(drSquared);
          const float invDrPower32 = 1 / drPower32;
          
          // Calculate the net force
          Fx[x - i] += dx * invDrPower32;
          Fy[x - i] += dy * invDrPower32;
          Fz[x - i] += dz * invDrPower32;
        }
      }

      // Accelerate bodies in response to the gravitational force
      for(x = i; x < i + TILESIZE; x++){
        bodySet.vx[x] += dt*Fx[x - i]; 
        bodySet.vy[x] += dt*Fy[x - i]; 
        bodySet.vz[x] += dt*Fz[x - i];
      }
    }
  }

  // Move bodies according to their velocities
  for (i = 0 ; i < nBodies; i++) { 
    bodySet.x[i] += bodySet.vx[i]*dt;
    bodySet.y[i] += bodySet.vy[i]*dt;
    bodySet.z[i] += bodySet.vz[i]*dt;
  }
}

int main(const int argc, const char** argv) {
  const char *threadsVar = getenv("OMP_NUM_THREADS");
  numThreads = atoi(threadsVar);

  struct centerOfMass COM;

  // Problem size and other parameters
  const int nBodies = (argc > 1 ? atoi(argv[1]) : 16384);
  const int nSteps = 10;  // Duration of test
  const float dt = 0.01f; // Body propagation time step

  // Body data stored as an Array of Structures (AoS)
  struct BodySet bodySet = createBodySet(nBodies);

  // Initialize omp threads
  omp_set_num_threads(numThreads);

  // Initialize random number generator and bodies
  srand(0);
  float randmax;
  randmax = (float) RAND_MAX;
  for(int i = 0; i < nBodies; i++) {
    bodySet.x[i] = ((float) rand())/randmax; 
    bodySet.y[i] = ((float) rand())/randmax; 
    bodySet.z[i] = ((float) rand())/randmax; 
    bodySet.vx[i] = ((float) rand())/randmax; 
    bodySet.vy[i] = ((float) rand())/randmax; 
    bodySet.vz[i] = ((float) rand())/randmax; 
  }

  // Compute initial center of mass  
  float comx = 0.0f, comy=0.0f, comz=0.0f;
  COM = computeCom(nBodies, bodySet);
  comx = COM.x;
  comy = COM.y;
  comz = COM.z;

  printf("Initial center of mass: (%g, %g, %g)\n", comx, comy, comz);

  // Perform benchmark
  printf("\n\033[1mNBODY Version 04\033[0m\n");
  printf("\nPropagating %d bodies using %d thread on %s...\n\n", nBodies, numThreads, "CPU");

  double rate = 0, dRate = 0; // Benchmarking data
  const int skipSteps = 3; // Set this to a positive int to skip warm-up steps

  printf("\033[1m%5s %10s %10s %8s\033[0m\n", "Step", "Time, s", "Interact/s", "GFLOP/s"); fflush(stdout);

  for (int step = 1; step <= nSteps; step++) {

    const double tStart = omp_get_wtime(); // Start timing
    MoveBodies(nBodies, bodySet, dt);
    const double tEnd = omp_get_wtime(); // End timing

    // These are for calculating flop rate. It ignores symmetry and 
    // estimates 20 flops per body-body interaction in MoveBodies
    const float HztoInts   = (float)nBodies * (float)(nBodies-1) ;
    const float HztoGFLOPs = 20.0*1e-9*(float)nBodies*(float)(nBodies-1);

    if (step > skipSteps) { 
      // Collect statistics 
      rate  += HztoGFLOPs / (tEnd - tStart); 
      dRate += HztoGFLOPs * HztoGFLOPs / ((tEnd-tStart)*(tEnd-tStart)); 
    }

    printf("%5d %10.3e %10.3e %8.1f %s\n", 
	   step, (tEnd-tStart), HztoInts/(tEnd-tStart), HztoGFLOPs/(tEnd-tStart), (step<=skipSteps?"*":""));
    fflush(stdout);
  }

  rate/=(double)(nSteps-skipSteps); 
  dRate=sqrt(fabs(dRate/(double)(nSteps-skipSteps)-rate*rate));

  printf("-----------------------------------------------------\n");
  printf("\033[1m%s %4s \033[42m%10.1f +- %.1f GFLOP/s\033[0m\n",
	 "Average performance:", "", rate, dRate);
  printf("-----------------------------------------------------\n");
  printf("* - warm-up, not included in average\n\n");

  // Compute final center of mass
  COM = computeCom(nBodies, bodySet);
  comx = COM.x;
  comy = COM.y;
  comz = COM.z;
  
  printf("Final center of mass: (%g, %g, %g)\n", comx, comy, comz);
  destroyBodySet(bodySet);
}
