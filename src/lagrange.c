/*****************************************
 *            lagrange.c                 *
 * Appliquer l interpolation de lagrange *
 ****************************************/


/* http://lumimath.univ-mrs.fr/~jlm/travaux/livretab/node17.html */

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
 * Print polynomial P
 * @param P Polynomial to print
 * @param degree Polynomial's degree
 */
void print_poly(mpz_t *P, int degree){
  int i;
  printf("Polynome = ");
  for (i=degree;i>0; i--)
    printf("%ldX^%d + ", mpz_get_si(P[i]), i);
  printf("%ld\n", mpz_get_si(P[i]));
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
 * Retourne le nombre de coefficient a 0 au debut de P
 */
int nb_zeros(mpz_t *P, int deg_P ){
  int i=0;
  
  while(!mpz_cmp_si(P[i], 0)){
    i++;
  }
  return i;
}

/**
 * Remplit new_P avec les coefficients de P sans coeffs initiaux nuls
 */
void del_zeros(mpz_t *new_P, mpz_t *P, int deg_P, int nb_zeros){
  int i;
  
  for(i=nb_zeros; i<deg_P+1; i++){
    mpz_set(new_P[i-nb_zeros], P[i]); 
  }
  
}
							     


/* P est de plus haut degre que Q */
/* res_int reste intermediaire : polynome de degre < deg_P : on prendra deg_P-1 puis on reduira si necessaire en sortie */

void div_aux(mpz_t *res_int, mpz_t *P, mpz_t *Q, int deg_P, int deg_Q, mpz_t mod){
  printf("entree dans div_aux\n");
  int i;
  
  mpz_t inverse, S[deg_P+1];
  mpz_init(inverse);
  
  for(i=0; i<deg_P+1; i++){
    mpz_init(S[i]);
  }
  
  printf("fin init S\n");
 
  inverse_modulaire(inverse, Q[0], mod);

  /* de deg_Q a deg_P, les coeffs sont a 0 */
  for(i=0; i<=deg_Q; i++){
    mpz_set(S[i], inverse);
    mpz_mul(S[i], Q[i], S[i]);
    mpz_mod(S[i], S[i], mod);
    mpz_mul(S[i], S[i], P[0]);
    mpz_mod(S[i], S[i], mod);
  }
  print_poly(S, deg_P);

  printf("fin remplissage S\n");
  
  /* On remplit res_int */
  for(i=1; i<deg_P+1; i++){
    mpz_sub(res_int[i], P[i], S[i]);
    mpz_mod(res_int[i], res_int[i], mod);
  }

  printf("fin remplissage res_int\n");
}

/**
 * Division de 2 polynomes
 */
void div_poly(mpz_t *res, mpz_t *P, mpz_t *Q, int deg_P, int deg_Q, mpz_t mod){

  int i, nb_zero;
  int deg_Pi=deg_P;
  int deg_Qi=deg_Q;
  mpz_t res_int[deg_P];

  for(i=0; i < deg_P+1; i++){
    mpz_init(res_int[i]);
  }
  
  /* mettre P de plus haut degre */
  div_aux(res_int, P, Q, deg_P, deg_Q, mod);
  nb_zero=nb_zeros(res_int, deg_P);

  /** Ligne declarant new_poly renvoie l'erreur:
   * error: invalid operands to binary - (have ‘int’ and ‘int (*)(struct __mpz_struct (*)[1], int)’) (avec une fleche poitant sur le signe "-"

   mpz_t new_poly[deg_P-nb_zeros +1];
  */
  
  /* del_zeros(new_poly, res_int, deg_P, nb_zero); */

  /* on boucle tant que le reste est de degre superieur au diviseur */
  while(deg_Pi-nb_zero >= deg_Q){

    /* Qi devient le reste precedent */
    mpz_t Qi[deg_Pi-nb_zero+1];
    for(i=0; i<deg_Pi-nb_zero+1; i++){
      mpz_init(Qi[i]);
      /*mpz_set(Qi[i], new_poly[i]);*/
    }

    /* Pi devient le diviseur precedent */
    mpz_t Pi[deg_Qi+1];
    for(i=0; i<deg_Q+1; i++){
      mpz_init(Pi[i]);
      mpz_set(Pi[i], Qi[i]);
    }
    /**
     * Declaration de tmp renvoie l erreur:
     * src/lagrange.c:151:5: error: ‘tmp’ undeclared (first use in this function)
     * tmp=deg_Pi;
     */
    deg_Pi=deg_Qi;
    deg_Qi=deg_Pi-nb_zero;
  }
}


/**
 * Multiply 2 polynomials of different sizes
 * @param res Matrix to complete, initialized at 0
 * @param P First polynomial
 * @param Q Second polynomial
 * @param deg_P P's degree
 * @param deg_Q Q's degree
 */
void mpz_mul_poly(mpz_t *res, mpz_t *P, mpz_t *Q, int deg_P, int deg_Q){
  int i, j;
  int deg_Res = deg_P+deg_Q;
  mpz_t tmp;
  mpz_init(tmp);
  mpz_t copy_tmp[N];

  printf("\n*************Multiplication*************\n");
	/* Init copy_tmp */
  for(i=0;i<=deg_Res;i++)
    mpz_init(copy_tmp[i]);
  /* Calcul */
  for(i=0;i<=deg_P;i++){
    for(j=0;j<=deg_Q;j++){
      mpz_mul(tmp, P[i], Q[j]);
      mpz_add(copy_tmp[i+j], copy_tmp[i+j], tmp);
    }
  }
	/* Affichage */
  print_P(copy_tmp, deg_Res);
  /* Copy de copy_tmp dans res*/
	for(i=0;i<=deg_Res;i++)
		mpz_set(res[i], copy_tmp[i]);
  printf("*************FIN MUL*************\n");
  /* Clear */
  mpz_clear(tmp);
  for(i=0;i<=deg_max;i++){
    mpz_clear(tmp_P[i]);
    mpz_clear(tmp_Q[i]);
  }
  for(i=0;i<=deg_Res;i++)
    mpz_clear(copy_tmp[i]);
}


/**
 * Apply Lagrange algorithm
 * @param M Matrix to complete
 * @param points Point value
 * @param images Result value
 * @param degre Degree
 */
void lagrange(mpz_t *res, mpz_t *points, mpz_t *images, int degre, mpz_t mod){
	
  int i, j;
  mpz_t bin, bin_sub;
  mpz_t tmp_poly[2], tmp_Li[degre+1];
  mpz_init_set_si(bin, 1);
  mpz_init(bin_sub);
  for(i=0;i<degre;i++)
    mpz_init(tmp_Li[i]);
  mpz_init_set_si(tmp_Li[degre], 1);
	
  /* Calcul des Li */
  for(i=0; i<degre+1; i++){
    printf("\n****Calcul de L_%d****\n", i);
    printf("--- L_%d initial: ", i);
    print_P(tmp_Li, degre);
    /* Quotient des produits (aj-ai) */
    for(j=0; j<degre+1; j++){
      if(j!=i){
	mpz_sub(bin_sub, points[j], points[i]);
	mpz_mul(bin, bin, bin_sub);
      }
    }
    
    /* Eleve a la puissance -1 */
    inverse_modulaire(bin, bin, mod);

    /* Multiplication par vi */
    mpz_mul(bin, bin, images[i]);
		
    /** 
     * Multiplication par les polynomes X-aj 
     * et remise a 1 de bin
     */
    for(j=0; j<degre+1; j++){
      
      /* Polynomes (X-aj) */
      if(j!=i){
	/* aj <- (-aj) */
	mpz_neg(bin_sub, points[j]);
	
	/* Construction du polynome (X-aj) */
	mpz_init_set_si(tmp_poly[0], 1);
	mpz_init_set(tmp_poly[1], bin_sub);
	printf("--( X - a_%d)\n", j);
	print_P(tmp_poly, 1);
	
	/* Multiplication de tmp_Li par le polynome */
	mpz_mul_poly(tmp_Li, tmp_Li, tmp_poly, degre, 1);
	printf("--- L_%d intermediaire: ", i);
	print_P(tmp_Li, degre);
      }
    }

    /* Multiplication de bin par tmp_Li */
    for(j=0;j<degre+1;j++){
      mpz_mul(tmp_Li[j], tmp_Li[j], bin);
    }

    printf("L_%d : ", i);
    print_poly(tmp_Li, degre);
    /* Ajout de Li a res */
    for(j=0; j<degre+1; j++)
      mpz_add(res[j], res[j], tmp_Li[i]);

    /* Remise de bin et des tmp_Li à 1 */
    for(j=0;j<degre+1;j++)
      mpz_set_si(tmp_Li[i], 1);
    mpz_set_si(bin, 1);
    
  }

  /* Clear */
  mpz_clear(bin);
  mpz_clear(bin_sub);
  mpz_clear(tmp_poly[0]);
  mpz_clear(tmp_poly[1]);
  for(j=0;j<degre+1; j++)
    mpz_clear(tmp_Li[i]);
}


/**
 * Multiply 2 polynomials of different sizes OLD
 * @param res Matrix to complete, initialized at 0
 */
void mpz_mul_poly_old(mpz_t *res, mpz_t *P, mpz_t *Q, int deg_P, int deg_Q){
  int i, j;
  int deg_max;
  max(&deg_max, deg_P, deg_Q, NULL);
  int N = deg_max+1;
  int R = 2*N;
  mpz_t tmp;
  mpz_init(tmp);
  mpz_t tmp_P[deg_max+1], tmp_Q[deg_max+1];
  mpz_t copy_tmp[R];

  printf("\n*************Multiplication*************\n");

  for(i=0;i<R;i++)
    mpz_init(copy_tmp[i]);

  /* Remplissage des polynomes temporaires */
  for(i=0; i<=deg_max; i++){
    if(i>deg_P){
      mpz_init(tmp_P[deg_max-i]);
      mpz_init_set(tmp_Q[deg_max-i], Q[deg_max-i]);
    }
    if(i>deg_Q){
      mpz_init(tmp_Q[deg_max-i]);
      mpz_init_set(tmp_P[deg_max-i], P[deg_max-i]);
    }
    else{
      mpz_init_set(tmp_P[deg_max-i], P[deg_P-i]);
      mpz_init_set(tmp_Q[deg_max-i], Q[deg_Q-i]);
    }
  }

  printf("Polynome intermediaires prets\n");
  print_P(tmp_P, deg_max);
  print_P(tmp_Q, deg_max);

  /* Calcul */
  for(i=0;i<=deg_max;i++){
    for(j=0;j<=i;j++){
      mpz_mul(tmp, tmp_P[i], tmp_Q[i-j]);
      mpz_add(copy_tmp[i], copy_tmp[i], tmp);
    }
  }
  
  /*printf("Calcul fait\n");*/
  print_P(copy_tmp, R-1);

  /* DEL ZEROS*/
  int nb_z=nb_zeros(copy_tmp, R-1);
  mpz_t tmp_t[R-1-nb_z];
  for(i=0;i<R-1-nb_z;i++)
    mpz_init(tmp_t[i]);
  del_zeros(tmp_t, copy_tmp, R-1, nb_z);
  print_P(tmp_t, R-1-nb_z);
  
  /* Copy de copy_t dans res*/

  printf("*************FIN MUL*************\n");
  
  
  /* Clear */
  mpz_clear(tmp);
  for(i=0;i<=deg_max;i++){
    mpz_clear(tmp_P[i]);
    mpz_clear(tmp_Q[i]);
  }
  for(i=0;i<R;i++)
    mpz_clear(copy_tmp[i]);
}