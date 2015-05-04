// Program: xdatgen
// Sung Sakong, PhD
// sung.sakong _at_ uni-due.de
// ver 0.5  4. June 2007
// XDATCAR generator from finite difference calculation
// ###.XDATCAR movie for each vibrational mode
// CONTCAR muliplier is implemented
// The characteristic length q_i0 is included
// -> maximum stretching represents one quanta for each 
//    normal mode
// input: OUTCAR output: ###.XDATCAR
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define LINESIZE 1024
#define MAXATOMTYPES 100

int main(void){

  FILE *outcar, *out, *pfile, *pos;
  int eachatom[MAXATOMTYPES];   // number of atoms per atom type
  int cumatom[MAXATOMTYPES];    // number of cumulated atoms
  float atommass[MAXATOMTYPES];
  int i, j, k, numatoms, atomcount, typecount, freq, l;
  char lineptr[LINESIZE];
  char filename[LINESIZE];
  char postitle[LINESIZE];
  int multix, multiy;

  float det;
  float cell[3][3];             // lattice vectors of the unit cell
  float inverse[3][3];
  float crdx, crdy, crdz, dx, dy, dz;
  float dcrdx, dcrdy, dcrdz, ddx, ddy, ddz;
  float omega, ddqi, qi0;

  float hbar, cl, am, sqmass;

  //physical constants in eV and Ang
  hbar=6.58211915*pow(10,-16);
  cl=2.99792458*pow(10,18);
  am=931.494043*pow(10,6);

  outcar=fopen("OUTCAR", "r");

  if (!outcar) {

    fclose(outcar);
    printf("ERROR :: Cannot find OUTCAR!\n");
    return NULL;
  }

  while(fgets(lineptr, LINESIZE, outcar)) {
   if (strstr(lineptr, "POSCAR:") != NULL) {
     sscanf(lineptr, "%*9c  %80c", &postitle);
     break;
    }
  }


  numatoms=0;
  while(fgets(lineptr, LINESIZE, outcar) && numatoms == 0) {
   if (strstr(lineptr, "NIONS =") != NULL) {
      sscanf(lineptr, " %*[ a-zA-Z] = %*d %*[ a-zA-Z] = %d", &numatoms);
      break;
    }
  }


  while (fgets(lineptr, LINESIZE, outcar)) {
     if (strstr(lineptr, "direct lattice vectors") != NULL) {
       int i;
       for (i = 0; i < 3; ++i) {
	 fgets(lineptr, LINESIZE, outcar);
	 if (3 != sscanf(lineptr, "%f %f %f", &cell[i][0], &cell[i][1], &cell[i][2])) {
	   fclose(outcar);
           fprintf(stderr, "ERROR: file 'OUTCAR' does not contain lattice vectors.\n");
           return NULL;
	 }
       }
       break;
     }
  }

  // define the inverse matrix of the lattice vectors
  // to convert the catesian coordinates into direct coordinates
  det = cell[0][0]*cell[1][1]*cell[2][2] 
    +   cell[1][0]*cell[2][1]*cell[0][2]
    +   cell[2][0]*cell[0][1]*cell[1][2] 
    -   cell[0][0]*cell[2][1]*cell[1][2]
    -   cell[1][0]*cell[0][1]*cell[2][2]
    -   cell[2][0]*cell[1][1]*cell[0][2];
  
  inverse[0][0] = ( cell[1][1]*cell[2][2]-cell[1][2]*cell[2][1] )/det;
  inverse[1][0] = ( cell[0][2]*cell[2][1]-cell[0][1]*cell[2][2] )/det;
  inverse[2][0] = ( cell[0][1]*cell[1][2]-cell[0][2]*cell[1][1] )/det;
  inverse[0][1] = ( cell[1][2]*cell[2][0]-cell[1][0]*cell[2][2] )/det;
  inverse[1][1] = ( cell[0][0]*cell[2][2]-cell[0][2]*cell[2][0] )/det;
  inverse[2][1] = ( cell[0][2]*cell[1][0]-cell[0][0]*cell[1][2] )/det;
  inverse[0][2] = ( cell[1][0]*cell[2][1]-cell[1][1]*cell[2][0] )/det;
  inverse[1][2] = ( cell[0][1]*cell[2][0]-cell[0][0]*cell[2][1] )/det;
  inverse[2][2] = ( cell[0][0]*cell[1][1]-cell[0][1]*cell[1][0] )/det;

  rewind(outcar);

  typecount = 0;
  atomcount = 0;
  while (fgets(lineptr, LINESIZE, outcar) && atomcount < numatoms) {
    if (strstr(lineptr, "POMASS") != NULL) sscanf(lineptr, " POMASS = %f;", &atommass[typecount++]);

    if (strstr(lineptr, "ions per type =") != NULL) {
      char const *token = strtok(lineptr, "=");
      for (i = 0; i < typecount; ++i) {
        token = strtok(NULL, " ");
        atomcount +=  eachatom[i] = atoi(token);
	cumatom[i] = atomcount;
      }
    }
  }

  if (atomcount != numatoms) {
    fprintf(stderr, "ERROR: file 'OUTCAR' does not have number of each atom.\n");
    return NULL;
  }

  while (fgets(lineptr, LINESIZE, outcar)) {
    if (strstr(lineptr, "dynamical") != NULL) {
      for(i=0; i<4; i++) {
	fgets(lineptr, LINESIZE, outcar);
      }

      for(freq=1; freq<=3*numatoms; freq++){

	l=sscanf(lineptr, "%*62c %f", &omega);
	if(l<0) break;
	//	printf("%d\n", l);
	qi0=sqrt(pow(hbar*cl,2)/(am*omega/1000.0));

	sprintf(filename, "%03d.XDATCAR", freq);
	out=fopen(filename, "w");
	sprintf(filename, "%03d.POSCAR", freq);
	pos=fopen(filename, "w");
	pfile=tmpfile();

	fprintf(out, "vibgen\n");
	for(j=0; j<5; j++){
	  fprintf(out, "\n");
	}

	fprintf(pos, "%s 1.000000\n", postitle);

	for(j=0; j<3; j++){
	  fprintf(pos, "%12.8f  %12.8f  %12.8f\n", cell[j][0], cell[j][1], cell[j][2]);
	}

	for (j = 0; j < typecount; ++j) {
	  fprintf(pos, "  %d  ", 4*eachatom[j]);
	}
	fprintf(pos, "\n");

	fgets(lineptr, LINESIZE, outcar);

	ddqi=0;
	l=0;
	for(j=0; j<numatoms; j++){
	  if(j+1 > cumatom[l]) l++;
	  sqmass = sqrt(atommass[l]);
	  //	  printf("%f\n", sqmass);

	  fgets(lineptr, LINESIZE, outcar);
	  sscanf(lineptr, "%f %f %f %f %f %f", &crdx, &crdy, &crdz, &dx, &dy, &dz);
	  ddqi+=dx*dx+dy*dy+dz*dz;

	  dcrdx = inverse[0][0]*crdx+inverse[0][1]*crdy+inverse[0][2]*crdz;
	  dcrdy = inverse[1][0]*crdx+inverse[1][1]*crdy+inverse[1][2]*crdz;
	  dcrdz = inverse[2][0]*crdx+inverse[2][1]*crdy+inverse[2][2]*crdz;

	  ddx = (inverse[0][0]*dx+inverse[0][1]*dy+inverse[0][2]*dz)/sqmass;
	  ddy = (inverse[1][0]*dx+inverse[1][1]*dy+inverse[1][2]*dz)/sqmass;
	  ddz = (inverse[2][0]*dx+inverse[2][1]*dy+inverse[2][2]*dz)/sqmass;

	  fprintf(pfile, "%12.8f  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f\n", 
		  dcrdx, dcrdy, dcrdz, ddx, ddy, ddz);
	}

	for(k=-20;k<21;k++){
	  rewind(pfile);
	  
	  for(j=0; j<numatoms; j++){
	    fscanf(pfile, "%f %f %f %f %f %f", &dcrdx, &dcrdy, &dcrdz, &ddx, &ddy, &ddz);
	    for(multix=0; multix<2; multix++){
	      for(multiy=0; multiy<2; multiy++){
		//		fprintf(out, "%12.8f  %12.8f  %12.8f\n", dcrdx+k/20.0*ddx*qi0+multix, 
		//			dcrdy+k/20.0*ddy*qi0+multiy, dcrdz+k/20.0*ddz*qi0);
		fprintf(out, "%12.8f  %12.8f  %12.8f\n", dcrdx+k/20.0*ddx+multix, 
			dcrdy+k/20.0*ddy+multiy, dcrdz+k/20.0*ddz);
	      }
	    }
	  }
	  fprintf(out, "\n");
	}

	fclose(pfile);

	for(j=0; j<2; j++){
	  fgets(lineptr, LINESIZE, outcar);
	}
      }

    }
  }

  fclose(outcar);

  return NULL;
}
