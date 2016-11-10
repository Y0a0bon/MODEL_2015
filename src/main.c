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
  char *path=argv[1];
  printf("path='%s'\n", path);
  int i, j, deg_P, deg_Q;
  mpz_t P[3], Q[4];
  mpz_t mod;
  mpz_init_set_str(mod, "127", 10);
  
  deg_P=2;
  deg_Q=3;
	int size=deg_P+deg_Q;
  int m_lines = size;
	int m_col = m_lines;
  /* Nombre total d elements dans la matrice */
  int matrix_length = m_lines*m_lines;
  
  mpz_t M[matrix_length];

  /* TEST SYLVESTER */
  /* Init and fill P */
  mpz_init_set_str(P[0], "1", 10);
  mpz_init_set_str(P[1], "2", 10);
  mpz_init_set_str(P[2], "3", 10);
  /* Init and fill Q */
  mpz_init_set_str(Q[0], "4", 10);
  mpz_init_set_str(Q[1], "5", 10);
  mpz_init_set_str(Q[2], "6", 10);
  mpz_init_set_str(Q[3], "7", 10);
  /* Init M */
  for (i=0; i < matrix_length; i++){
    mpz_init(M[i]);
  }
  /* Fill M */
  printf("Sylvester matrix:\n");
  sylvester(P, Q, deg_P, deg_Q, M);
  /* Print M */
  print_M(M, deg_P+deg_Q);

  
  /* TEST SOUS-MATRICE */
  int indice=1;
  m_col=m_col-2*indice;
  mpz_t sub_M[m_col*m_lines];
  /* Init sub_M */
  for(i=0; i < m_lines; i++){
    for(j=0; j< m_col; j++)
      mpz_init(sub_M[i*m_col+j]);
  }
  sub_matrix(sub_M, M, deg_P, deg_Q, indice);
  print_M_2(sub_M, m_col, m_lines);
  /* Matrice Ti */
  mpz_t Ti[m_col*(m_lines-indice)];
  /* Init Ti */
  for(i=0; i<m_col*(m_lines-indice); i++)
    mpz_init(Ti[i]);
  /* Fill Ti */
  del_lines(Ti, sub_M, m_col, m_lines, indice);
  printf("Matrice finale Ti:\n");
  int m_lin_Ti = m_lines-indice;
  print_M_2(Ti, m_col, m_lin_Ti);

  /* Construire Vi */
  mpz_t *Vi[(deg_P+deg_Q-2*indice)*(deg_P+deg_Q-indice)];
  creat_Vi(Vi, deg_P, deg_Q, indice);

  /* TEST LECTURE FICHIER */
  char *data_raw[2], *data_P[MAX_SIZE/2], *data_Q[MAX_SIZE/2];
  /* Recupere les donnes et split selon le separateur */
  getData(path, data_raw);
  printf("data : %s\n%s\n", data_raw[0], data_raw[1]);
  split(data_raw, data_P, data_Q);

  /* Remplit les degres du polynome en Y et ses coefficients en X */
  int deg_P_Y, deg_Q_Y;
  int *deg_PX, *deg_QX;
  
  get_Y(data_P, &deg_P_Y);
  get_Y(data_Q, &deg_Q_Y);
  printf("degPY = %d & degQY = %d\n", deg_P_Y, deg_Q_Y);
  deg_PX=malloc((deg_P_Y+1)*sizeof(mpz_t));
  deg_QX=malloc((deg_Q_Y+1)*sizeof(mpz_t));
  get_X(data_P, deg_P_Y, deg_PX);
  get_X(data_Q, deg_Q_Y, deg_QX);

  /* Remplit P et Q */
  mpz_t *polP[deg_P_Y+1], *polQ[deg_Q_Y+1];
  /* Allocation */
  for (i=0; i<deg_P_Y+1; i++)
    polP[i]=malloc((deg_PX[i]+1)*sizeof(mpz_t));
  for (i=0; i<deg_Q_Y+1; i++)
    polQ[i]=malloc((deg_QX[i]+1)*sizeof(mpz_t));
  
  get_coeff(data_P, polP, deg_P_Y, deg_PX);
  get_coeff(data_Q, polQ, deg_Q_Y, deg_QX);
  
  for(i=0;i<deg_P_Y+1; i++){
    print_P_red(polP[i], deg_PX[i]);
  }
  
  print_biv(polP, deg_PX, deg_P_Y, mod);

  return EXIT_SUCCESS;
}


/**
 *  Fonction trash
 */
void trash(){

  printf("Launching...\n");
  int i, deg_P, deg_Q, matrix_length;
  mpz_t P[2], Q[3];
  mpz_t mod;
  mpz_init_set_str(mod, "127", 10);
  
  deg_P=2;
  deg_Q=2;
  matrix_length = deg_P+deg_Q;
  /* Nombre total d elements dans la matrice */
  matrix_length *= matrix_length;
  
  mpz_t M[matrix_length];

  /* TEST SYLVESTER */
  /* Init and fill P */
  mpz_init_set_str(P[0], "1", 10);
  mpz_init_set_str(P[1], "-2", 10);
  /* init_set_str(P[2], "3", 10); */
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
  printf("test_euclide\n");
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


  /**
   * -------------------
   */

  /* TEST GAUSS + DETERMINANT */
  mpz_t nb_premier, determinant;
  mpz_init_set_si(nb_premier, 127);
  mpz_init(determinant);
  
  gauss(&determinant, M, deg_P+deg_Q, nb_premier);
  
  print_M(M, deg_P+deg_Q);

  printf("determinant : %ld\n", mpz_get_si(determinant));

  
  /**
   * -------------------
   */
  
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


  /* TEST MULTIPLICATION DE POLYNOMES */
  mpz_t Pol[3], Pol2[2], Q_res[4];

  for(i=0;i<4;i++)
    mpz_init(Q_res[i]);
  
  mpz_init_set_si(Pol[0], 7);
  mpz_init_set_si(Pol[1], 1);
  mpz_init_set_si(Pol[2], 2);

  mpz_init_set_si(Pol2[0], 4);
  mpz_init_set_si(Pol2[1], 3);
  
  print_P(Pol, 2);
  print_P(Pol2, 1);
  
  mpz_mul_poly(Q_res, Pol, Pol2, 2, 1);
  printf("Resultat du produit:\n");
  print_P(Q_res, 3);


  /* TEST DE LAGRANGE */
  mpz_t points[3], images[3], res[3], res_mod[4];
  mpz_t modulo;
  int deg_Res = 2;

  /* Init modulo */
  mpz_init_set_si(modulo, 11);
  /* Init points */
  mpz_init_set_si(points[0],0);
  mpz_init_set_si(points[1],-1);
  mpz_init_set_si(points[2],1);
  /* Init images */
  mpz_init_set_si(images[0],-1);
  mpz_init_set_si(images[1],1);
  mpz_init_set_si(images[2],-1);
  /* Init res */
  for(i=0;i<3;i++){
    mpz_init(res[i]);
    mpz_init(res_mod[i]);
  }
  mpz_init(res_mod[3]);
  /* Apply Lagrange */
  lagrange(res, points, images, deg_Res, modulo, res_mod);
  /* Print res */
  printf("\nResultat :\n");
  print_P(res, 2);
  printf("mod ");
  print_P(res_mod, 3);


   /* TEST HORNER */
  printf("\ntest Horner\n");
  mpz_t resultatH, valeur;
  mpz_init(resultatH);
  mpz_init_set_str(valeur, "-2", 10);

  print_P(P, deg_P);
  
  horner(resultatH, P, deg_P, valeur, mod);
  printf("P(%ld)=%ld\n", mpz_get_si(valeur), mpz_get_si(resultatH));

  printf("\ntest racines\n");
  /*
  mpz_t rac[deg_P+1];
  for(i=0; i<deg_P+1; i++){
    mpz_init_set_str(rac[i], "-1", 10);
  }
  */
  
  print_racines(P, deg_P, mod);
  /* print_P(rac, deg_P); */


  
  /* TEST PRINT/EVAL */
  /*
  mpz_t **PY;
  int degres_PY, degres_QY, deg_PY, degQY;
  print_biv(PY, degres_PY, deg_PY, mod);
  mpz_t res_eval, value_X, value_Y;
  mpz_init(res_eval);
  mpz_init_set_str(value_X, "0", 10);
  mpz_init_set_str(value_Y, "2", 10); 
  eval_bivXY(res_eval, value_X, value_Y, PY, degres_PY, deg_PY, mod);
  printf("EVALUATION en (X, Y) = (%ld, %ld) = %ld\n", mpz_get_si(value_X), mpz_get_si(value_Y), mpz_get_si(res_eval));
  
  */
  /* TEST QUESTION 1 */
  /*
  mpz_t value;
  mpz_init(value);
  mpz_set_si(value, 4);

  mpz_t Pi[deg_PY], Qi[deg_QY];
  init_mpzs(Pi, 0, deg_PY+1);
  init_mpzs(Qi, 0, deg_QY+1);
  */
  /* taille : degre max du resultant +1 */
  /*
  int borne=2*deg_PY*deg_QY+1;
  mpz_t res[borne];
  for(i=0; i<borne; i++)
    mpz_init(res[i]);
  resultant(res, PY, QY, deg_PY, deg_QY, degres_PY, degres_QY, mod);
  

  int nombre=nb_zeros(res, borne-1);
  mpz_t new_res[borne-nombre];
  for(i=0; i<borne-nombre; i++)
    mpz_init(new_res[i]);
 
  del_zeros(new_res, res, borne-1, nombre);
  printf("resultant : \n");
  print_P(new_res, borne-1-nombre);
  */
  /*eval_biv(value, PY, QY, degres_PY, degres_QY, Pi, Qi, deg_PY, deg_QY, mod );
    print_P(Pi, deg_PY);
    print_P(Qi, deg_QY);*/

  /* QUESTION 2 */
  /* appeler racines */
  /*
  mpz_t rac[deg_P+1];
  for(i=0; i<deg_P+1; i++)
    mpz_init_set_str(rac[i], "-1", 10);

  mpz_t nb_racines;
  mpz_init(nb_racines);
  racines(rac, new_res, borne-1-nombre, nb_racines, mod);
  */
  /* TEST CRT */
  mpz_t res_crt, a, b, mod_a, mod_b;
  mpz_init(res_crt);
  mpz_init_set_str(a, "8", 10);
  mpz_init_set_str(b, "127", 10);
  mpz_init_set_str(mod_a, "31", 10);
  mpz_init_set_str(mod_b, "13", 10);
  crt(res_crt, a, b, mod_a, mod_b);
  printf("resultat crt : %ld\n", mpz_get_si(res_crt));
  
  
}
