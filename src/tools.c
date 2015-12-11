#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>
#include <stdarg.h>
#include "../include/utils.h"

#define MAX_SIZE 1024

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
 * Return min of nb_param integers
 * @param resultat Contains the result
 */
void min(int *resultat, ...){
  int param = -1;
  *resultat = -1;
  va_list ap;

  va_start(ap, resultat);

  do{
    if(*resultat > param || *resultat < 0)
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


/**
 *Affiche le polynome bivarie PY
 */
void print_biv(mpz_t **PY, int *degres_PY, int deg_P, mpz_t mod){
  int i, j;
  
  for(j=0; j<deg_P; j++){
    printf("( ");
    for (i=0; i<degres_PY[j]; i++){
      printf("%ld*X^%d + ", mpz_get_si(PY[j][i]), degres_PY[j]-i);
    }
    printf("%ld )", mpz_get_si(PY[j][i]));
    printf(" * Y^%d + ", deg_P-j);
  }

  /*coeff en Y⁰ */
  for (i=0; i<degres_PY[j]; i++){
      printf("%ld*X^%d + ", mpz_get_si(PY[j][i]), degres_PY[j]-i);
    }
    printf("%ld \n", mpz_get_si(PY[j][i]));
}


/**
 * evalue le polynome bivarie PY en (value_X, value_Y)
 */
void eval_bivXY(mpz_t res, mpz_t value_X, mpz_t value_Y, mpz_t **PY, int *degres_PY, int deg_P, mpz_t mod){

  int i;
  mpz_t P[deg_P+1];
  for(i=0; i<deg_P+1; i++)
    mpz_init(P[i]);
  
  for (i=0; i<deg_P+1; i++){
    horner(P[i], PY[i], degres_PY[i], value_X, mod);
  }
  horner(res, P, deg_P, value_Y, mod);
}


/**
 *  Recupere les 2 polynomes specifies dans un fichier
 */
int getData(char *path, char *out[2]){

  FILE *fp;
  out[0]=malloc(MAX_SIZE*sizeof(char));
  out[1]=malloc(MAX_SIZE*sizeof(char));
  
  fp=fopen(path, "r");
  if(fp == NULL) {
      perror("Error opening file");
      return -1;
   }

  if( fgets (out[0], MAX_SIZE, fp)==NULL )
      printf("Error reading P\n");

  if( fgets (out[1], MAX_SIZE, fp)==NULL )
      printf("Error reading Q\n");

  fclose(fp);

  return EXIT_SUCCESS;
}


/**
* Extrait des chaines les donnees
*/
void split(char *str[2], char **outP, char **outQ){

  int i=0;
  char *p=NULL, *p2=NULL;
  /* Polynome P */
  p = strtok(str[0], " ");
  outP[0]=malloc((1+strlen(p))*sizeof(char));
  outP[i++]=p;
  while((p=strtok(NULL, " "))!=NULL){
    outP[i]=malloc((1+strlen(p))*sizeof(char));
    outP[i++]=p;
  }
  /* Polynome Q */
  i=0;
  p2 = strtok(str[1], " ");
  outQ[0]=malloc((1+strlen(p2))*sizeof(char));
  outQ[i++]=p2;
  while((p2=strtok(NULL, " "))!=NULL){
    outQ[i]=malloc((1+strlen(p2))*sizeof(char));
    outQ[i++]=p2;
  }
}


/**
 *  Recupere le degre du polynome en Y
 */
void get_Y(char **str, int *deg_Y){
  *deg_Y=atoi(str[0]);
}


/**
 *  Recupere les degres de chaque polynome en X
 *  formant les coeff du polynome en Y
 */
void get_X(char **str, int deg_Y, int *deg_X){
  int i;
  for(i=1;i<=deg_Y+1;i++)
    deg_X[i-1]=atoi(str[i]);
}


/**
 *  Remplit le polynome de ses coefficients
 */
void get_coeff(char **data, mpz_t **pol, int deg_Y, int *deg_X){

  int i, j;
  int curseur=0;
  for(i=0;i<deg_Y+1;i++){
    for(j=0;j<deg_X[i]+1;j++){
      mpz_init_set_str(pol[i][j], data[2+deg_Y+curseur], 10);
      curseur++;
    }
  }
}
