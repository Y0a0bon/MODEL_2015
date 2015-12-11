/*****************************************************************
 *              solution.c                                       *
 *  Calcul des solution du systeme associe au fichier en entree  *
 *****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>
#include "../include/utils.h"


#define MAX_SIZE 1024

void solve(char *path, mpz_t mod1, mpz_t mod2){

  int i, j;
  /* LECTURE du fichier */

   char *data_raw[2], *data_P[MAX_SIZE/2], *data_Q[MAX_SIZE/2];
  /* Recupere les donnes et split selon le separateur */
  getData(path, data_raw);
  split(data_raw, data_P, data_Q);

  /* Remplit les degres du polynome en Y et ses coefficients en X */

  /* degres en Y des 2 polynomes bivaries */
  int deg_P_Y, deg_Q_Y;

  /* tableaux des degres des differents coefficients en Y */
  int *deg_PX, *deg_QX;

  get_Y(data_P, &deg_P_Y);
  get_Y(data_Q, &deg_Q_Y);
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


  printf("polynome P modulo %ld :\n", mpz_get_si(mod1));
  print_biv(polP, deg_PX, deg_P_Y, mod1);

  printf("polynome Q modulo %ld :\n", mpz_get_si(mod1));
  print_biv(polQ, deg_QX, deg_Q_Y, mod1);


  /* CALCULS MODULO MOD1 */
  printf("*****************Calculs modulo %ld*****************\n", mpz_get_si(mod1));

  /* deg_result est la borne sur le degre max du resultant */
  int deg_result=2*deg_P_Y*deg_Q_Y;
  mpz_t result[deg_result+1];
  init_mpzs(result, 0, deg_result+1 );

  /* Calcul du resultant */
  resultant(result, polP, polQ, deg_P_Y, deg_Q_Y, deg_PX, deg_QX, mod1);

  printf("resultant : \n");
  print_P(result, deg_result);

  /* on supprime les coefficients a 0 du resultant */
  int nb=nb_zeros(result, deg_result);
  mpz_t new_result[deg_result-nb+1];
  for(i=0; i<deg_result-nb+1; i++){
    mpz_init(new_result[i]);
  }
  del_zeros(new_result, result, deg_result, nb);
  printf("resultant sans 0 initiaux :\n");
  print_P(new_result, deg_result-nb);

  /* Calcul des racines du resultant */
  printf("Calcul des racines du resultant\n");
  mpz_t rac[deg_P_Y+1];
  for(i=0; i<deg_P_Y+1; i++)
    mpz_init_set_str(rac[i], "-1", 10);

  /* nb_racines_X contient le nb de racines du resultant */
  int nb_racines_X=0, nb_racines_Y=0;
  racines(rac, new_result, deg_result-nb, &nb_racines_X, mod1);

  i=0;
  printf("Racines du resultant:\n");
  while(mpz_cmp_si(rac[i], -1))
    printf("%ld\t", mpz_get_si(rac[i++]));
  printf("\n");


  mpz_t P[deg_P_Y+1], Q[deg_Q_Y+1];
  init_mpzs(P, 0, deg_P_Y+1);
  init_mpzs(Q, 0, deg_Q_Y+1);

  /* on a au max deg_P_Y racines au systeme de 2 equations a 2 inconnues */
  mpz_t rac_Y[deg_P_Y];
  init_mpzs(rac_Y, 0, deg_P_Y);

  /* */
  mpz_t *couples_tmp[nb_racines_X];
  int size_couples=0;
  /* boucle sur les racines du resultant */
  for(i=0; i<nb_racines_X; i++){

    /* evaluation du polynome bivarie en x=result[i] */
    eval_biv(rac[i], polP, polQ, deg_PX, deg_QX, P, Q, deg_P_Y, deg_Q_Y, mod1 );

    printf("polynome P evalue en %ld :\n", mpz_get_si(rac[i]));
    print_P(P, deg_P_Y);
    printf("polynome Q evalue en %ld :\n", mpz_get_si(rac[i]));
    print_P(Q, deg_Q_Y);
    printf("\n");

    /* calcul des racines en Y correspondant a x=rac[i] */

    for(j=0; j<deg_P_Y; j++)
      mpz_set_si(rac_Y[j], -1);

    printf("racines en Y pour X=%ld :\n", mpz_get_si(rac[i]));

    /* affichage des racines dans search_roots */
    search_roots(&nb_racines_Y, rac_Y, P, Q, deg_P_Y, deg_Q_Y, mod1);

    couples_tmp[i]=malloc((nb_racines_Y+2)*sizeof(mpz_t));
    size_couples+=nb_racines_Y;
    printf("nb_racines = %d\n", nb_racines_Y);
    if(nb_racines_Y!=0){
      printf("x[%d] = %ld\t", i, mpz_get_si(rac[i]));
      mpz_init_set(couples_tmp[i][0], rac[i]);
      mpz_init_set_si(couples_tmp[i][1], nb_racines_Y);

      /* Remplit avec les racines de rac_Y */
      for(j=0;j<nb_racines_Y;j++){
	mpz_set(couples_tmp[i][j+2], rac_Y[j]);
	printf("%ld  ", mpz_get_si(rac_Y[j]));
      }
      printf("\n");
    }
    nb_racines_Y=0;
  }
    mpz_t couples[2][size_couples];
    printf("Nb de couples = %d\n", size_couples);
    for(i=0;i<size_couples;i++){
      for(j=0;j<mpz_get_si(couples_tmp[i][1]);j++){
	mpz_init_set(couples[0][i], couples_tmp[i][0]);
	mpz_init_set(couples[1][i], couples_tmp[i][j+2]);
      }
    }
    if(size_couples){
    print_P(couples[0], size_couples-1);
    print_P(couples[1], size_couples-1);
    }
    else{
      printf("Pas de racines en Y !\nFin.\n");
      exit(0);
    }

    /* CALCULS MODULO MOD2 */
    printf("*****************Calculs modulo %ld*****************\n", mpz_get_si(mod2));

    /* deg_result est toujours la borne sur le degre max du resultant */

  mpz_t result2[deg_result+1];
  init_mpzs(result2, 0, deg_result+1 );

  /* Calcul du resultant */
  resultant(result2, polP, polQ, deg_P_Y, deg_Q_Y, deg_PX, deg_QX, mod2);

  printf("resultant2 : \n");
  print_P(result2, deg_result);

  /* on supprime les coefficients a 0 du resultant */
  int nb2=nb_zeros(result2, deg_result);
  mpz_t new_result2[deg_result-nb2+1];
  for(i=0; i<deg_result-nb2+1; i++){
    mpz_init(new_result2[i]);
  }
  del_zeros(new_result2, result2, deg_result, nb2);
  printf("resultant2 sans 0 initiaux :\n");
  print_P(new_result2, deg_result-nb2);

  /* Calcul des racines du resultant */
  printf("Calcul des racines du resultant2\n");
  mpz_t rac2[deg_P_Y+1];
  for(i=0; i<deg_P_Y+1; i++)
    mpz_init_set_str(rac2[i], "-1", 10);

  /* nb_racines_X contient le nb de racines du resultant */
  int nb_racines_X2=0, nb_racines_Y2=0;
  racines(rac2, new_result2, deg_result-nb2, &nb_racines_X2, mod2);

  i=0;
  printf("Racines du resultant2:\n");
  while(mpz_cmp_si(rac2[i], -1))
    printf("%ld\t", mpz_get_si(rac2[i++]));
  printf("\n");


  mpz_t P2[deg_P_Y+1], Q2[deg_Q_Y+1];
  init_mpzs(P2, 0, deg_P_Y+1);
  init_mpzs(Q2, 0, deg_Q_Y+1);

  /* on a au max deg_P_Y racines au systeme de 2 equations a 2 inconnues */
  mpz_t rac_Y2[deg_P_Y];
  init_mpzs(rac_Y2, 0, deg_P_Y);

  /* */
  mpz_t *couples_tmp2[nb_racines_X2];
  int size_couples2=0;
  /* boucle sur les racines du resultant */
  for(i=0; i<nb_racines_X2; i++){

    /* evaluation du polynome bivarie en x=result[i] */
    eval_biv(rac2[i], polP, polQ, deg_PX, deg_QX, P2, Q2, deg_P_Y, deg_Q_Y, mod2 );

    printf("polynome P2 evalue en %ld :\n", mpz_get_si(rac2[i]));
    print_P(P2, deg_P_Y);
    printf("polynome Q2 evalue en %ld :\n", mpz_get_si(rac2[i]));
    print_P(Q2, deg_Q_Y);
    printf("\n");

    /* calcul des racines en Y correspondant a x=rac[i] */

    for(j=0; j<deg_P_Y; j++)
      mpz_set_si(rac_Y2[j], -1);

    printf("racines en Y pour X=%ld :\n", mpz_get_si(rac2[i]));

    /* affichage des racines dans search_roots */
    search_roots(&nb_racines_Y2, rac_Y2, P2, Q2, deg_P_Y, deg_Q_Y, mod2);

    couples_tmp2[i]=malloc((nb_racines_Y2+2)*sizeof(mpz_t));
    size_couples2+=nb_racines_Y2;
    printf("nb_racines2 = %d\n", nb_racines_Y2);
    if(nb_racines_Y2!=0){
      printf("x[%d] = %ld\t", i, mpz_get_si(rac2[i]));
      mpz_init_set(couples_tmp2[i][0], rac2[i]);
      mpz_init_set_si(couples_tmp2[i][1], nb_racines_Y2);

      /* Remplit avec les racines de rac_Y */
      for(j=0;j<nb_racines_Y2;j++){
	mpz_set(couples_tmp2[i][j+2], rac_Y2[j]);
	printf("%ld  ", mpz_get_si(rac_Y2[j]));
      }
      printf("\n");
    }
    nb_racines_Y2=0;
  }
    mpz_t couples2[2][size_couples2];
    printf("Nb de couples 2= %d\n", size_couples2);
    for(i=0;i<size_couples2;i++){
      for(j=0;j<mpz_get_si(couples_tmp2[i][1]);j++){
	mpz_init_set(couples2[0][i], couples_tmp2[i][0]);
	mpz_init_set(couples2[1][i], couples_tmp2[i][j+2]);
      }
    }
    if(size_couples2){
    print_P(couples2[0], size_couples2-1);
    print_P(couples2[1], size_couples2-1);
    }
    else{
      printf("Pas de racines en Y !\nFin.\n");
      exit(0);
    }

    mpz_t mod;
    mpz_init(mod);
    mpz_mul(mod, mod1, mod2);
    /* Recherche des racines */
    find_roots(couples, couples2, size_couples, size_couples2,  mod1, mod2, polP, polQ, deg_PX, deg_QX, deg_P_Y, deg_Q_Y, mod);
		
		printf("\n********************\n");
		int mini;
		min(&mini, deg_P_Y, deg_Q_Y, NULL);
		
		for(i=0; i<mini; i++){
			mpz_t sub_M[(m+n)*(n-i)];
			for(j=0; j<(m+n)*(n-i); j++)
				mpz_init(sub_M[j]);
			sub_matrix(sub_M, );
		}
		
}
