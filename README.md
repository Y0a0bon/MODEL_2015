#######################
#    MODEL project    #
#######################

************************
Pour instancier 2 polynomes P et Q appartenant a |K[X, Y], creer le fichier comme suit:
* deg_PY deg_PX[deg_PY] coeffs
* deg_QY deg_QX[deg_QY] coeffs

Tous les separateurs sont des espaces.

* deg_PY: degre du polynome P en Y
* deg_PX[deg_PY]: degre de chaque polynome en X formant les coefficients du polynome en Y
* coeffs Coefficients des polynoomes en Y, selon X, du plus grand au plus petit

************************
Exemple:

* P=(3*X + 1)*Y^3 + (X^2)*Y^2 + (2*X)*Y + (X^3)

s'ecrit:

* 3 1 2 1 3 3 1 1 0 0 2 0 1 0 0 0
* 3: P de degre 2 en Y;
* 1 2 1 3: polynomes selon X de degre respectivement 1, 2, 1 et 3;
* 1 1: coeff du premier polynome en X, coeff de Y^3;
* 1 0 0: coeff du 2e polynome;
* 2 0: coeff du 3e;
* 1 0 0 0: coeff du 4e

************************
COMPILATION:

  $ make


************************
EXECUTION

  $ bin/main votrefichier.txt modulo1 modulo2
(modulo1 et modulo2 premiers)

************************
NETTOYAGE

  $ make clean
