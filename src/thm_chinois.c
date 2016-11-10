/*****************************************************************
 *            thm_chinois.c                                      *
 * Application du theoreme chinois pour reconstruire des racines *
 *****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "../include/utils.h"

/**
 * Applique le crt et remplit res solution de X = a mod mod_a, X = b mod mod_b
 */
void crt(mpz_t res, mpz_t a, mpz_t b, mpz_t mod_a, mpz_t mod_b){
  
  mpz_t res_euclide[3], tmp1, tmp2, mod_fin;
  mpz_inits(res_euclide[0], res_euclide[1], res_euclide[2], tmp1, tmp2, mod_fin, NULL);
  mpz_mul(mod_fin, mod_a, mod_b);
  mpz_mod(a, a, mod_a);
  mpz_mod(b, b, mod_b);
  euclide_etendu(mod_a, mod_b, res_euclide);

  if (!mpz_cmp_si(res_euclide[2], 0)){
    printf("Choisir des moduli premiers entre eux\n");
    return;
  }

  
  mpz_mul(tmp1, mod_a, res_euclide[0]);
  mpz_mul(tmp1, tmp1, b);

  mpz_mul(tmp2, mod_b, res_euclide[1]);
  mpz_mul(tmp2, tmp2, a);

  mpz_add(res, tmp1, tmp2);
  mpz_mod(res, res, mod_fin);
}

/**
 * Recombine les racines modulo mod1 et mod2 pour trouver les racines finales.
 * rac_modi[0]contient les valeurs en x des racines et rac_modi[1] celles en y. 
 * @param rac_mod1 tableau des couples (x, y) de racines modulo mod1
 * @param rac_mod2 tableau des couples (x, y) de racines modulo mod2
 * @param nb1 nombre de racines (x, y) dans rac_mod1
 * @param nb2 nombre de racines (x, y) dans rac_mod2 
 */
void find_roots(mpz_t *rac_mod1[2], mpz_t *rac_mod2[2], int nb1, int nb2,  mpz_t mod1, mpz_t mod2, mpz_t **PY, mpz_t **QY, int *degres_PY, int *degresQY, int deg_P, int deg_Q, mpz_t mod){
  
  int i, j, nb_racines=0;
    mpz_t rx, ry, isroot;
  mpz_inits(rx, ry, isroot, NULL);

  for(i=0; i<nb1; i++){
    for(j=0; j<nb2; j++){

      /* crt sur la i eme racine mod1 et la j eme mod2 */ 
      crt(rx, rac_mod1[0][i], rac_mod2[0][j], mod1, mod2);
      crt(ry, rac_mod1[1][i], rac_mod2[1][j], mod1, mod2);

      /* on evalue en le resultat (rx, ry) trouve */
      eval_bivXY(isroot, rx, ry, PY, degres_PY, deg_P, mod);

      /* on teste si le resultat est bien racine */
      if (!mpz_cmp_si(isroot, 0)){
	printf("(%ld, %ld) est racine\n", mpz_get_si(rx), mpz_get_si(ry));
	nb_racines++;
      }
    }
  }
  /*
  printf("%d racines au total\n", k);
   */
  
  
}
