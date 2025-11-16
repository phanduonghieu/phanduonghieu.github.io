#include "gmp.h"
#include"ecm.h"

#ifndef PARFAIT_H
#define PARFAIT_H


typedef struct premier_exposant premier_exposant;

struct premier_exposant
{
    mpz_t z_premier;
	int cmpt_expo_list;
    unsigned long int exposant;
};
extern int test_sur_composante(mpz_t z_p,int special);
extern void mpz_sigma(mpz_t z_somme, mpz_t z_premier_N, int exposant);
extern void mpz_recherche_maximum(mpz_t z_max, int n, mpz_t z_tab_sigma []);
extern int factoriser(mpz_t z_n,mpz_t tableau_10000_premier[10000],premier_exposant tableau_diviseur[45][45],mpz_t z_max,int j) ;
extern int contrainte_b2(mpz_t z_new,unsigned long int expo_new);
extern void initialisation_matrice(premier_exposant tableau_diviseur[45][45]);
extern void initialisation_tableau_fixe(premier_exposant tableau_fixe[45]);
extern void initialisation_tableau_ajuste(premier_exposant tableau_ajuste[2025]);
extern void tableau_exposant(unsigned long int exposant[40]);
extern void remplissage_tableau_10000_premier(mpz_t tableau_10000_premier[10000]);
extern void tri_tableau(premier_exposant tableau_trie[2025],premier_exposant tableau_ajuste[2025]);
extern int contrainte_b1(premier_exposant tableau_trie[2025],premier_exposant tableau_fixe[8],int j,mpz_t z_b1, int special);
extern int contrainte_b3(mpz_t z_b3);
extern int contrainte_s(premier_exposant tableau_trie[2025],premier_exposant tableau_fixe[8],int j,mpz_t z_s,mpz_t z_b1, int special);
#endif
