#ifndef _UTILS_H_
#define _UTILS_H_

/* sylvester.c */

int bound_Max(mpz_t *M, int column, int matrix_length, int r);

void transpose(mpz_t *M, int matrix_length);

void inverse_modulaire(mpz_t inverse, mpz_t nb, mpz_t mod);

void sub_matrix(mpz_t *res, mpz_t input, int deg_m, int deg_n, int indice);

void euclide_etendu(mpz_t nb1, mpz_t nb2, mpz_t *resultat);

void switch_lines(int k, int r, int matrix_length, mpz_t *M);

void del_coeffs(int r, int matrix_length, mpz_t *M, mpz_t mod);
  
void gauss(mpz_t *determinant, mpz_t *M, int matrix_length, mpz_t mod);

void sylvester(mpz_t *P, mpz_t *Q, int degre_P, int degre_Q, mpz_t *M);

void sub_matrix(mpz_t *res, mpz_t input, int m, int n, int indice);

void eval_biv(mpz_t value, mpz_t **PY, mpz_t **QY, int *degres_PY, int *degres_QY, mpz_t *P, mpz_t *Q, int deg_P, int deg_Q, mpz_t mod );

void resultant(mpz_t *resultant, mpz_t **PY, mpz_t **QY, int deg_P, int deg_Q, int *degres_PY, int *degres_QY, mpz_t mod);



/* lagrange.c */

void div_aux(mpz_t *res_int, mpz_t *P, mpz_t *Q, int deg_P, int deg_Q, mpz_t mod);

void div_poly(mpz_t *res, mpz_t *P, mpz_t *Q, int deg_P, int deg_Q, mpz_t mod);

int nb_zeros(mpz_t *P, int deg_P );

void del_zeros(mpz_t *new_P, mpz_t *P, int deg_P, int nb_zeros);

void mpz_mul_poly(mpz_t *res, mpz_t *P, mpz_t *Q, int deg_P, int deg_Q);

void lagrange(mpz_t *res, mpz_t *points, mpz_t *images, int degre, mpz_t mod, mpz_t *res_mod);


/* horner.c */

void horner(mpz_t res, mpz_t *P, int deg_P, mpz_t valeur, mpz_t mod);

void print_racines(mpz_t *P, int deg_P, mpz_t mod);

void racines(mpz_t *rac, mpz_t *P, int deg_P, int *nb_racines, mpz_t mod);

/* tools.c */

void print_P(mpz_t *polynome, int degre);

void init_mpzs(mpz_t *nb, int inf, int sup);

void max(int *resultat, ...);

void print_mpz(mpz_t mpz);

void print_M(mpz_t *M, int size);

#endif

