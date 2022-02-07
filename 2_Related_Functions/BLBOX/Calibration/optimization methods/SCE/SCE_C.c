/* C script for Shuffled complex evolution by Duan et al. */

/* Include headers */
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/types.h>
#include <unistd.h>


int main(const char options){ 
    FILE *fp;
    /* Declare variables */
    
    /* Read in options from specified file */
      fp = fopen(options, "r");
      if (fp == NULL) {
         printf("I couldn't open %c for reading.\n",options);
         exit(0);
      }
      
    
    
    
    /* End of main function, return 0 */
    return 0;    
}