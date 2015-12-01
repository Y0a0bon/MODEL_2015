/*****************************************
 *            horner.c                   *
 * Trouver les racines d un polynome     *
 ****************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>  /* pour le sleep */
#include <gmp.h>
#include "../include/utils.h"

/*
 *Evaluation du polynome P en valeur
 *@param res resultat de l evaluation
 */
void horner(mpz_t res, mpz_t *P, int deg_P, mpz_t valeur, mpz_t mod){
  int i;
  
  mpz_mul(res, P[0], valeur);
  for (i=1; i<deg_P; i++){
    mpz_add(res, res, P[i]);
    mpz_mul(res, res, valeur);
  }
  mpz_add(res, res, P[deg_P]);
  mpz_mod(res, res, mod);
}

/**
 *Affiche les racines de P
 */
void print_racines(mpz_t *P, int deg_P, mpz_t mod){
  int i, nb_racines=0;
  mpz_t rac[deg_P+1];
  for(i=0; i<deg_P+1; i++){
    mpz_init_set_str(rac[i], "-1", 10);
  }

  printf("\n\nEtude du Polynome :  ");
  print_P(P, deg_P);
  
  racines(rac, P, deg_P,&nb_racines, mod);

  for(i=0; i<deg_P+1; i++){
    if (!mpz_cmp_si(rac[i],-1)){
      printf("Polynome :  ");
      print_P(P, deg_P);
      printf("%d racine(s)\n", nb_racines);
      
      return;
    }
    printf("%ld est racine\n", mpz_get_si(rac[i]));
  }
  printf("Polynome :  ");
  print_P(P, deg_P);
  printf("%d racines\n", nb_racines);
  
}

/* rac de taille deg_P+1 rempli de -1*/
/* nb_racines sera rempli du nombre de racines de P */
void racines(mpz_t *rac, mpz_t *P, int deg_P, int *nb_racines, mpz_t mod){
  /* k indice de remplissage de rac */
  long int i, k=0;
  mpz_t res_horner, m, m_plus;
  /* m_plus= m+1 */
  mpz_inits(res_horner, m, NULL);
  mpz_init_set_str(m_plus, "1", 10);
  /* on reduit P modulo mod */
  for (i=0; i<deg_P+1; i++){
    mpz_mod(P[i], P[i], mod);
  }
  
  while( mpz_cmp(mod, m_plus) ){
    /* printf("mod=%ld, m=%ld\n", mpz_get_si(mod),mpz_get_si(m)); */
    /* on teste la valeur m */
    horner(res_horner, P, deg_P, m, mod);

    /* si on a une racine, on l ajoute a rac */
    if (!mpz_cmp_si(res_horner, 0)){
      mpz_set(rac[k], m);
      /* printf("-----------K=%ld---------------------\n", k); */
      k++;
    }
    /* printf("apres if : k=%ld\n", k); */
    /* on ne peut pas avoir plus de racines que le degre du polynome +1 */
    if (k >=deg_P+1){
      printf("deg_P+1 atteint : k=%ld, deg_P+1=%d\n", k, deg_P+1);
      *nb_racines=k;
      mpz_clear(res_horner);
      mpz_clear(m);
      mpz_clear(m_plus);
      return;
    }
    mpz_add_ui(m, m, 1);
    mpz_add_ui(m_plus, m_plus, 1);
    /* printf("fin while :mod=%ld, m=%ld\n", mpz_get_si(mod),mpz_get_si(m)); */
  }
  /* printf("sortie:mod=%ld, m=%ld\n", mpz_get_si(mod),mpz_get_si(m)); */
  *nb_racines=k;
  mpz_clear(res_horner);
  mpz_clear(m);
  mpz_clear(m_plus);
  return;
}
