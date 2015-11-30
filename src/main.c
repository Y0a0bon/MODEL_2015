/***************************************
 *              main.c                 *
 *         Prgramme principal          *
 **************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "../include/utils.h"


int main(int argc, char **argv){
  printf("Launching...\n");
  int i, deg_P, deg_Q, matrix_length;
  mpz_t P[2], Q[3];
  /* mpz_t mod; */

  deg_P=1;
  deg_Q=2;
  matrix_length = deg_P+deg_Q;
  /* Nombre total d elements dans la matrice */
  matrix_length *= matrix_length;
  
  mpz_t M[matrix_length];


  /* TEST SYLVESTER */
  /* Init and fill P */
  mpz_init_set_str(P[0], "1", 10);
  mpz_init_set_str(P[1], "2", 10);
  /* mpz_init_set_str(P[2], "3", 10); */
  /* Init and fill Q */
  mpz_init_set_str(Q[0], "3", 10);
  mpz_init_set_str(Q[1], "4", 10);
  mpz_init_set_str(Q[2], "5", 10);
  /* mpz_init_set_str(Q[3], "7", 10); */
  /* Init M */
  for (i=0; i < matrix_length; i++){
    mpz_init(M[i]);
  }
  /* Fill M */
  printf("Sylvester matrix:\n");
  sylvester(P, Q, deg_P, deg_Q, M);
  /* Print M */
  print_M(M, deg_P+deg_Q);
	
  /**
   * -------------------
   */

  /* TEST EUCLIDE */
  /*printf("test_euclide\n");
  mpz_t resultat[3];

  mpz_inits(resultat[0], resultat[1], resultat[2], (void *)NULL);

  mpz_t nb1, nb2;
  mpz_init_set_si(nb1, 21);
  mpz_init_set_si(nb2, 17);

  printf("resultat[1]=%ld\n",mpz_get_si(resultat[1]));
  
  euclide_etendu(nb1, nb2, resultat);
  printf("apres euclide\n");
  printf("euclide sur %ld et %ld : u=%ld, v=%ld, pgcd=%ld\n",mpz_get_si(nb1),mpz_get_si(nb2), mpz_get_si(resultat[0]),mpz_get_si(resultat[1]),  mpz_get_si(resultat[2]));

  mpz_t inverse;
  mpz_init(inverse);
  inverse_modulaire(inverse, nb2, nb1);
  printf("inverse de %ld modulo %ld : %ld", mpz_get_si(nb2), mpz_get_si(nb1), mpz_get_si(inverse)); 


  del_coeffs(1, deg_P+deg_Q, M);
  print_M(M, deg_P+deg_Q);

  del_coeffs(2, deg_P+deg_Q, M);
  print_M(M, deg_P+deg_Q); */

  /**
   * -------------------
   */

  /* TEST GAUSS + DETERMINANT */
  /*
  mpz_t nb_premier, determinant;
  mpz_init_set_si(nb_premier, 127);
  mpz_init(determinant);
  
  gauss(&determinant, M, deg_P+deg_Q, nb_premier);
  
  print_M(M, deg_P+deg_Q);

  printf("determinant : %ld\n", mpz_get_si(determinant));
  */
  
  /**
   * -------------------
   */
  
  /*
  print_P(P, deg_P);
  print_P(Q, deg_Q);

  printf("test_lagrange\n");

  mpz_t res_int[deg_P+1];
  
  for(i=0; i < deg_P+1; i++){
    mpz_init(res_int[i]);
  }

  printf("appel a div_aux\n");
  div_aux(res_int, P, Q, deg_P, deg_Q, mod);

  print_P(res_int, deg_P);

  printf("appel nb_zeros\n");
  int nb=nb_zeros(res_int, deg_P);

  printf("appel del_zeros\n");
  mpz_t new_P[deg_P-nb+1];
  for(i=0; i<deg_P-nb+1; i++){
    mpz_init(new_P[i]);
  }
  del_zeros(new_P, res_int, deg_P, nb);
  print_P(new_P, deg_P-nb);
  */
  
  /**
   * -------------------
   */

  /* TEST MULTIPLICATION DE POLYNOMES */
  
  mpz_t Pol[3], Pol2[2], Q_res[4];

  for(i=0;i<4;i++)
    mpz_init(Q_res[i]);
  
  mpz_init_set_si(Pol[0], 7);
  mpz_init_set_si(Pol[1], 1);
  mpz_init_set_si(Pol[2], 2);

  mpz_init_set_si(Pol2[0], 4);
  mpz_init_set_si(Pol2[1], 3);
  
  print_poly(Pol, 2);
  print_poly(Pol2, 1);
  
  mpz_mul_poly(Q_res, Pol, Pol2, 2, 1);
  printf("Resultat du produit:\n");
  print_poly(Q_res, 3);
  

  /**
   * -------------------
   */


  
  printf("\nEnd.\n");

  return EXIT_SUCCESS;
}

void trash(){

  /* TEST DE LAGRANGE*/
  mpz_t points[3], images[3], res[3];
  mpz_t modulo;
  int deg_p=2;
  int i;
  /* Init modulo */
  mpz_init_set_si(modulo, 10);
  /* Init points */
  mpz_init_set_si(points[0],0);
  mpz_init_set_si(points[1],-1);
  mpz_init_set_si(points[2],1);
  /* Init images */
  mpz_init_set_si(images[0],-1);
  mpz_init_set_si(images[1],1);
  mpz_init_set_si(images[2],-1);
  /* Init res */
  for(i=0;i<3;i++)
    mpz_init(res[i]);
  /* Apply Lagrange */
  /*lagrange(res, points, images, deg_p, modulo);*/
  /* Print res */
  printf("Resultat = ");
  print_poly(res, 2);
  
}