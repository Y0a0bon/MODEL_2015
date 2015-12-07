/***************************************
 *              main.c                 *
 *         Prgramme principal          *
 **************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "../include/utils.h"

#define MAX_SIZE 1024


int main(int argc, char **argv){
  printf("Launching...\n");

  /* TEST SOLVE */
  char *path=argv[1];
  printf("path='%s'\n", path);

  mpz_t mod1, mod2;

  mpz_init_set_str(mod1, argv[2], 10);
  mpz_init_set_str(mod2, argv[3], 10);

  solve(path, mod1, mod2);
  printf("\nEnd\n");
  
  return EXIT_SUCCESS;
}

