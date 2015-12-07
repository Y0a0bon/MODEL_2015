/***************************************
 *            sylvester.c              *
 * Calcule du resultant de 2 polynomes *
 **************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "../include/utils.h"

#define MAX_SIZE 256


/**
 * Renvoie l'indice de la ligne du maximum de la colonne column
 * a partir de son element ligne r
 */
int bound_Max(mpz_t *M, int column, int matrix_length, int r){
  int i, indice=r;
  mpz_t max;
  mpz_init(max);
  mpz_set(max, M[r*matrix_length+column]);
  for(i=r+1; i<matrix_length; i++){
    if(mpz_get_si(M[i*matrix_length+column])>mpz_get_si(max)){
      mpz_set(max, M[i*matrix_length+column]);
      indice=i;
    }
  }
  mpz_clear(max);
  return indice;
}


/**
 * Calcule la transposee d'une matrice
 * @param M Matrice a transposer
 */
void transpose(mpz_t *M, int matrix_length){
  int i, j;
  mpz_t tmp;
  mpz_init(tmp);
  for(i=1; i<matrix_length; i++){
    for(j=0; j<i; j++){
      mpz_set(tmp, M[j*matrix_length+i]);
      mpz_set(M[j*matrix_length+i], M[i*matrix_length+j]);
      mpz_set(M[i*matrix_length+j], tmp);
    }
  }
  mpz_clear(tmp);
}

/**
 * Retourne l inverse modulaire de nb modulo mod
 */
void inverse_modulaire(mpz_t inverse, mpz_t nb, mpz_t mod){

  mpz_t resultat[3];
  mpz_inits(resultat[0], resultat[1], resultat[2], (void *)NULL);

  euclide_etendu(nb, mod,resultat);
  
  /* pas d inverse modulaire */
  if (mpz_get_si(resultat[2]) != 1){
    printf("Pas d'inverse modulo %ld, choisir %ld premier\n", mpz_get_si(mod), mpz_get_si(mod));
    exit(0);
  }
  mpz_set(inverse, resultat[0]);
	
  /* Clear */
  mpz_clear(resultat[0]);
  mpz_clear(resultat[1]);
  mpz_clear(resultat[2]);
}




/**
 * Applique euclide etendu a nb1 et nb2
 * @param resultat tableau de 3 entiers : u, v et r tels que r=pgcd(nb1, nb2) et a*u+b*v=r
 */
void euclide_etendu(mpz_t nb1, mpz_t nb2, mpz_t *resultat){

  mpz_t r, r1, u, u1, v, v1;
  mpz_inits(r, r1, NULL);
  mpz_set(r, nb1);
  mpz_set(r1, nb2);
  mpz_init_set_si(u, 1);
  mpz_init_set_si(u1, 0);
  mpz_init_set_si(v, 0);
  mpz_init_set_si(v1, 1);
  
  /* printf("u=%ld, v=%ld, pgcd=%ld\n",mpz_get_si(u),mpz_get_si(v), mpz_get_si(r)); */

  mpz_t q, r_temp, u_temp, v_temp, tmp;
  mpz_inits(q, r_temp, u_temp, v_temp, tmp, NULL);
  while (mpz_get_si(r1) != 0){
    mpz_cdiv_q(q,r,r1); /* Division entiere */

    mpz_set(r_temp,r);
    mpz_set(u_temp,u);
    mpz_set(v_temp,v);

    mpz_set(r,r1);
    mpz_set(u,u1);
    mpz_set(v,v1);

    mpz_mul(tmp, q, r1);
    mpz_sub(r1, r_temp, tmp);
    /* r1=r_temp-q*r1; */
    mpz_mul(tmp, q, u1);
    mpz_sub(u1, u_temp, tmp);
    mpz_mul(tmp, q, v1);
    mpz_sub(v1, v_temp, tmp);
  }
  if(mpz_get_si(r)<0){
    mpz_neg(u,u);
    mpz_neg(v,v);
    mpz_neg(r,r);
  }
  
  mpz_set(resultat[0],u);
  mpz_set(resultat[1],v);
  mpz_set(resultat[2],r);
	
  /* Creer une fonction */
	
  /* Clear */
  mpz_clear(r);
  mpz_clear(r1);
  mpz_clear(u);
  mpz_clear(u1);
  mpz_clear(v);
  mpz_clear(v1);
}


/**
 * Echange les lignes d indices k et r de la matrice M
 */
void switch_lines(int k, int r, int matrix_length, mpz_t *M){
  int i;
  mpz_t tmp;
  mpz_init(tmp);
  for (i=0; i<matrix_length; i++){
    mpz_set(tmp, M[k*matrix_length+i]);
    mpz_set(M[k*matrix_length+i], M[r*matrix_length+i]);
    mpz_set(M[r*matrix_length+i], tmp);
  }
  mpz_clear(tmp);
}


/**
 * Set pivot column to 0
 * @param r indice de la premiere ligne a annuler
 * on recalcule les coefficients modulo mod a chaque fois
 */
void del_coeffs(int r, int matrix_length, mpz_t *M, mpz_t mod){
  int i, j;
  mpz_t M_tmp[matrix_length*matrix_length];
  for (i=0; i < matrix_length*matrix_length; i++){
    mpz_init(M_tmp[i]);
    mpz_set(M_tmp[i], M[i]);
  }
    
  for (i=r; i<matrix_length; i++){

    for (j=r-1; j<matrix_length; j++){
      /* a modifier si on enleve la normalisation */
      mpz_submul(M_tmp[i*matrix_length+j], M[i*matrix_length+r-1], M[(r-1)*matrix_length+j]);
    }
  }
  for (i=0; i < matrix_length*matrix_length; i++){
    mpz_mod(M[i], M_tmp[i], mod);
  }
	/* Clear */
	for(i=0;i<matrix_length*matrix_length;i++)
		mpz_clear(M_tmp[i]);
}


/**
 * Apply Gauss algorithm
 * @param determinant le determinant de la matrice M
 */
void gauss(mpz_t *determinant, mpz_t *matrice, int matrix_length, mpz_t mod){
    
  int i, j, k=-1;
  int r=0;
  mpz_t inverse;
  mpz_init(inverse);

  /* changement de signe du determinant si signe == 1 */
  long int signe=1;
  /* coefficient par lequel multiplier le determinant */
  mpz_t coeff_det;
  mpz_init_set_si(coeff_det,1);


  /* on copie matrice dans M */
  mpz_t M[matrix_length*matrix_length];
  for(i=0; i<matrix_length*matrix_length; i++){
    mpz_init_set(M[i], matrice[i]);
  }

  /* Reduit la matrice M modulo mod */
  for (i=0; i < matrix_length*matrix_length; i++){
    mpz_mod(M[i], M[i], mod);
  }

  /* Algorithme de Gauss */
  for(j=0; j<matrix_length; j++){
    k=bound_Max(M, j, matrix_length, r);
    if (k != r)
      signe =-signe;

    if(mpz_get_si(M[k*matrix_length+j])!=0){
      /* r indice de la premiere ligne que l on va annuler */
      r++;
      /* On normalise la ligne du pivot (supprimer pour determinant ?) */
      inverse_modulaire(inverse, M[k*matrix_length+j], mod);
      /* printf("inverse : %ld\n", mpz_get_si(inverse)); */

      for (i=j; i<matrix_length; i++){
	mpz_mul(M[k*matrix_length+i], M[k*matrix_length+i],inverse);
      }
      /* Actualisation du coefficient pour le determinant */
      /* printf("coeff_det : %ld\n", mpz_get_si(coeff_det)); */
      mpz_mul(coeff_det, coeff_det, inverse );
    
      /* On place le pivot en haut des lignes a traiter */
     
      switch_lines(k, r-1,matrix_length, M);
    
      del_coeffs(r, matrix_length, M, mod);
     
    }
  }
  /* printf("coeff_det : %ld, signe : %ld\n", mpz_get_si(coeff_det), signe); */
  mpz_mul_si(*determinant, coeff_det, signe);
  mpz_mod(*determinant, *determinant, mod);
  inverse_modulaire(*determinant, *determinant, mod);
	
  mpz_clear(inverse);
  mpz_clear(coeff_det);
}


/**
 * calcule le determinant d une matrice a coefficients entiers donnee en entree
 *
 */
void determinant (mpz_t det, mpz_t *M, int matrix_length, mpz_t mod){

}


/**
 * Remplit la matrice de Sylvester M en parametre associee au systeme de 2 equations S donne en entree
 * @param P Equation du polynome P
 * @param Q Equation du polynome Q
 * @param deg_P Degré du polynome P
 * @param deg_Q Degré du polynome Q
 * @param M: Matrice de Sylvester transposee, carrée de taille degre_P+degre_Q, implementee en ligne, a remplir
 */
void sylvester(mpz_t *P, mpz_t *Q, int deg_P, int deg_Q, mpz_t *M){

  int i, k;
  int size=deg_P+deg_Q;

  for(i=0; i<size; i++){
    /* Si on remplit "la partie P" */
    if(i<deg_Q){
      for(k=i; k<i+deg_P+1; k++)
	mpz_set(M[i*size+k], P[k-i]);
    }
    /* Si on remplit "la partie Q" */
    else{
      for(k=i; k<i+deg_Q+1; k++)
	mpz_set(M[(i)*size+k-deg_Q], Q[k-i]);
    }
  }
  /* Transpose la matrice obtenue */
  transpose(M, size);
}


/**
 * Calcul la sous-matrice Ti de la matrice de Sylvester,
 * compose des colonnes de 1 a n-indice et de n+1 a m+n+indice
 */
void sub_matrix(mpz_t *res, mpz_t *input, int m, int n, int indice){
  int col, lin;
	int size=m+n;
	/* On ne copie que les colonnes voulues dans la nouvelle matrice */
	for(lin=0; lin<size; lin++){
		/* Colonnes 1 a n-indice */
		for(col=0; col<n-indice; col++)
			mpz_set(res[lin*(size-2*indice)+col], input[lin*size+col]);
		/* Colonnes n a m+n-indice */
		for(col=n; col<m+n-indice; col++)
			mpz_set(res[lin*(size-2*indice)+(col-indice)], input[lin*size+col]);
	}
}

/*
* Delete i_th lines of zeros
*/
void del_lines(mpz_t *res, mpz_t *input, int indice, int m_col, int m_lines){
	int col, lin;
	/* Init & fill */
	for(lin=0; lin<m_lines-indice; lin++){
		for(col=0; col<m_col; col++)
			mpz_init_set(res[lin*(m_lines-indice)+col], input[lin*m_lines+col]);
	}
}


/**
 *  Evalue PY et QY en value et remplit les polynomes en Y P et Q
 */
void eval_biv(mpz_t value, mpz_t **PY, mpz_t **QY, int *degres_PY, int *degres_QY, mpz_t *P, mpz_t *Q, int deg_P, int deg_Q, mpz_t mod ){
  int i;

  for (i=0; i<deg_P+1; i++){
    horner(P[i], PY[i], degres_PY[i], value, mod);
 
  }

  for(i=0; i<deg_Q+1; i++ ){
   
    horner(Q[i], QY[i], degres_QY[i], value, mod);
 
  }

}


/**
 * Calcule le resultant des polynomes bivaries PY et QY
 * Appeler nb_zeros et del_zeros en sortie pour avoir le resultant sans les premiers coeffs nuls
 * @param resultant le resultant de PY et QY de taille 2*deg_P*deg_Q+1
 * @param deg_P le degre du polynome PY en Y (nb de colonnes)
 * @param deg_Q degre de QY en Y
 * @param degres_PY liste des degres en X des coefficients de PY
 * @param degres_QY liste des degres en X des corfficients de QY
 */
void resultant(mpz_t *resultant, mpz_t **PY, mpz_t **QY, int deg_P, int deg_Q, int *degres_PY, int *degres_QY, mpz_t mod){
  printf("resultant\n");
  int i,j, borne;
  int matrix_length=deg_P+deg_Q;
  int matrix_size= matrix_length*matrix_length;
  
  /* allocation de la matrice de sylvester */
  mpz_t M[matrix_size];
  for(i=0; i<matrix_size; i++)
    mpz_init(M[i]);

  /* Calcul de la borne superieure sur le degre du resultant */
  borne=2*deg_P*deg_Q+1;

  /* si le modulo est trop petit, on ne peut pas choisir assez de valeurs pour l interpolation */
  if(mpz_cmp_si(mod, borne)<=0){
    printf("Modulo trop petit !\n");
    exit(0);
  }
  printf("borne = %d\n", borne);

  /* Choix des valeurs d interpolation*/
  mpz_t values[borne];
  for(i=0; i<borne; i++){
    mpz_init_set_si(values[i], i);
    /*printf("values[%d] = %ld", i, mpz_get_si(values[i]));*/
  }
  /* initialisation de P, Q et du tableau des determinant */
  /*printf("\nfin init values\n");*/
  mpz_t P[deg_P+1], Q[deg_Q+1];
  for(i=0; i<deg_P+1; i++){
    mpz_init(P[i]);
  }
  for(i=0; i<deg_Q+1; i++){
    mpz_init(Q[i]);
  }
  /* determinant contiendra les images des points de values pour l interpolation */
  mpz_t determinant[borne];
  for(i=0;i<borne;i++){
    mpz_init(determinant[i]);
  }
  /* A chaque iteration, on traite une valeur de values et on calcule le determinant associe */
  for(i=0; i<borne; i++){
    /* evaluation du polynome en values[i] */
    eval_biv(values[i], PY, QY, degres_PY, degres_QY, P, Q, deg_P, deg_Q, mod );
    /*print_P(P, deg_P);
    print_P(Q, deg_Q);
    
    printf("appel sylvester\n");
      */
    /* Calcul de la matrice de Sylvester dans M*/
    sylvester(P, Q, deg_P, deg_Q, M);
    /*printf("matrice de sylvester:\n");
    print_M(M, deg_P+deg_Q);

    printf("appel Gauss\n");
      */
    /* Appel de Gauss */
    /* remplit le determinant associe a values[i] */
    gauss(&determinant[i], M, matrix_length, mod);

    /* on remet M a 0 pour la prochaine iteration */
    for(j=0; j<matrix_size; j++)
      mpz_set_si(M[j], 0);
  }
  /* affichage des determinants (avec la fonction d affichage de ploynomes) */
  /*printf("determinants : ");
  print_P(determinant, borne-1);*/

  /* appel a Lagrange */
  
  mpz_t res_mod[borne+1];
  for(j=0; j<borne+1; j++){
    mpz_init(res_mod[j]);
  }
  lagrange(resultant, values, determinant, borne-1, mod, res_mod);
  /*printf("fin lagrange\n");*/

  for(j=0; j<borne; j++){
    mpz_mod(resultant[j],resultant[j], mod);
  }
  
}
