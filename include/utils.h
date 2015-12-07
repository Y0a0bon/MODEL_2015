#ifndef _UTILS_H_
#define _UTILS_H_

/* sylvester.c */

int bound_Max(mpz_t *M, int column, int matrix_length, int r);

void transpose(mpz_t *M, int matrix_length);

void inverse_modulaire(mpz_t inverse, mpz_t nb, mpz_t mod);

void euclide_etendu(mpz_t nb1, mpz_t nb2, mpz_t *resultat);

void switch_lines(int k, int r, int matrix_length, mpz_t *M);

void del_coeffs(int r, int matrix_length, mpz_t *M, mpz_t mod);
  
void gauss(mpz_t *determinant, mpz_t *M, int matrix_length, mpz_t mod);

void sylvester(mpz_t *P, mpz_t *Q, int degre_P, int degre_Q, mpz_t *M);

void sub_matrix(mpz_t *res, mpz_t *input, int m, int n, int indice);

void del_lines(mpz_t *res, mpz_t *input, int indice, int m_col, int m_lines);

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

void search_roots(int *nb_racines, mpz_t *racines, mpz_t *P, mpz_t *Q, int deg_P, int deg_Q, mpz_t mod);

/* tools.c */

void print_P(mpz_t *polynome, int degre);

void init_mpzs(mpz_t *nb, int inf, int sup);

void max(int *resultat, ...);

void print_mpz(mpz_t mpz);

void print_M(mpz_t *M, int size);

void print_M_2(mpz_t *M, int length, int width);

int getData(char *path, char *out[2]);

void split(char *str[2], char **outP, char **outQ);

void get_Y(char **str, int *deg_Y);

void get_X(char **str, int deg_Y, int *deg_X);

void get_coeff(char **data, mpz_t **pol, int deg_Y, int *deg_X);

void print_biv(mpz_t **PY, int *degres_PY, int deg_P, mpz_t mod);

void eval_bivXY(mpz_t res, mpz_t value_X, mpz_t value_Y, mpz_t **PY, int *degres_PY, int deg_P, mpz_t mod);


/* thm_chinois.c */

void crt(mpz_t res, mpz_t a, mpz_t b, mpz_t mod_a, mpz_t mod_b);

void find_roots(mpz_t *rac_mod1[2], mpz_t *rac_mod2[2], int nb1, int nb2,  mpz_t mod1, mpz_t mod2, mpz_t **PY, mpz_t **QY, int *degres_PY, int *degresQY, int deg_P, int deg_Q, mpz_t mod);

/* solution.c */

void solve(char *path, mpz_t mod1, mpz_t mod2);

#endif

