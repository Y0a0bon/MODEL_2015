#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>
#include <stdarg.h>
#include "../include/utils.h"


/**
 * Print polynomial P
 * @param P Polynomial to print
 * @param degree Polynomial's degree
 */
void print_P(mpz_t *polynome, int degre){
  int i;
  for (i=0; i<degre; i++)
    printf("%ld*X^%d + ", mpz_get_si(polynome[i]), degre-i);
  printf("%ld\n", mpz_get_si(polynome[i]));
}


/**
 * Init a mpz_t vector to 0
 */
void init_mpzs(mpz_t *nb, int inf, int sup){
  int i;
  for(i=inf;i<sup;i++)
    mpz_init(nb[i]);
}


/**
 * Return max of nb_param integers
 * @param resultat Contains the result
 */
void max(int *resultat, ...){
  int param = 1;
  *resultat = 1;
  va_list ap;

  va_start(ap, resultat);

  do{
    if(*resultat < param)
      *resultat=param;
    param=va_arg(ap, int);
  }while(param!=NULL);

  va_end(ap);
}


/**
 * Print a mpz_t number
 * @param mpz Mpz number to print
 */
void print_mpz(mpz_t mpz){
  printf("%ld", mpz_get_si(mpz));
}



/**
 * Print a square matrix
 * @param M Square matrix to print
 * @param size Square matrix M's length
 */
void print_M(mpz_t *M, int size){
  int i;
  char buf[256];
  for(i=0; i<size*size; i++){
    if(!(i%size) && i!=0)
      printf("\n");
    mpz_get_str(buf, 10, M[i]);
    printf("%s\t", buf);
  }
  printf("\n\n");
}

/**
 * Print a matrix
 * @param M Matrix to print
 * @param length Matrix M's length
 * @param width Matrix M's width
 */
void print_M_2(mpz_t *M, int columns, int lines){
  int i, j;
  char buf[256];
  for(i=0; i<lines; i++){
		printf("\n");
		for(j=0;j<columns; j++){
			mpz_get_str(buf, 10, M[i*columns+j]);
			printf("%s\t", buf);
		}
	}
  printf("\n\n");
}

