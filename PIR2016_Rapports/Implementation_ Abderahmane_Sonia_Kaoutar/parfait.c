#include <math.h>
#include "parfait.h"


//fonction de factorisation test sur 10 000 nombre premiers
int factoriser(mpz_t z_n,mpz_t tableau_10000_premier[10000],premier_exposant tableau_diviseur[45][45], mpz_t z_max,int j) 
{
    int i,k,comparaison_expo,comparaison,nb_facteur,cmp_max ;
    mpz_t z_r;
    mpz_init(z_r);
    k = 0;
   
	
for (i =0 ; i < 10000; i++)
{ 
	mpz_mod (z_r,z_n, tableau_10000_premier[i]);
	comparaison = mpz_cmp_ui (z_r, 0);
	if ( comparaison == 0 )
        {
		if (k<1)	
		 {
		  	mpz_set(tableau_diviseur[k][j].z_premier,tableau_10000_premier[i]);             
		  	gmp_printf("----------->[ %Zu ]\n",tableau_10000_premier[i]);
			cmp_max=mpz_cmp(z_max,tableau_10000_premier[i]);
				if(cmp_max == -1)
					{
						mpz_set(z_max,tableau_10000_premier[i]);
						//gmp_printf("\n Le max ici est : %Zu",z_max);					
					}	
		mpz_divexact (z_n,z_n,tableau_10000_premier[i]);
                k++;
		}
		else
		{
		comparaison_expo=mpz_cmp(tableau_diviseur[k-1][j].z_premier,tableau_10000_premier[i]);
			if (comparaison_expo==0)
			{
				 tableau_diviseur[k-1][j].exposant++;
				mpz_divexact (z_n,z_n,tableau_10000_premier[i]);
				
			}			
			else
			{	
			
			mpz_set(tableau_diviseur[k][j].z_premier,tableau_10000_premier[i]);
			mpz_divexact (z_n,z_n,tableau_10000_premier[i]);k++;
			gmp_printf("----------->[ %Zu ]\n",tableau_10000_premier[i]);
			cmp_max=mpz_cmp(z_max,tableau_10000_premier[i]);
				if(cmp_max == -1)
					{mpz_set(z_max,tableau_10000_premier[i]);
						//Zgmp_printf("\n Le max ici est : %Zu\n",z_max);
					}	
			}
		
		}
	i=i-1;
	 
        }
      
      
}
//printf("\n");
nb_facteur = k;   
return nb_facteur;
}


// fonctions sur la composante premiere modulo 4
 int test_sur_composante(mpz_t z_premier, int special)
{
	mpz_t z_r,z_d,z_c;
	mpz_init (z_r);
	int resultat;
	mpz_init_set_ui	(z_d,4);
	mpz_init_set_ui	(z_c,1);
	mpz_mod (z_r, z_premier, z_d);
	resultat = mpz_cmp (z_r, z_c);
	if ((resultat == 0)&&(special==0))
	{ 
		printf("\n------------------------------------\n > Il faut prendre l'exposant '1', ce nombre peut être le premier spécial\n");
		return 1;
	}
	else
	{
		printf("\n------------------------------------\n > Il faut prendre un exposant pair\n");		
		return 0;
	}
	
}

// fonction de recherche du max dans un tableau
void mpz_recherche_maximum (mpz_t z_max, int n, mpz_t z_tab_sigma []) 
{
	int i,resultat;
	unsigned long int m;
	m = 0;


 for (i = 0; i < n; i++)
     {
	resultat = mpz_cmp(z_tab_sigma[i] , z_tab_sigma[m]);
        if(resultat > 0)  
	{ m = i ; }

     } 
mpz_set(z_max,z_tab_sigma[m]);
return ;
}

// fonction sigma
void mpz_sigma(mpz_t z_somme,mpz_t z_premier_N, int exposant)
{
    int i;
	i = 0;
	mpz_t resultat_expo;
	mpz_init_set_ui(resultat_expo,0);
	mpz_init_set_ui(z_somme,0);
	
    for (i = 0 ; i <= exposant ; i++)
    {
     mpz_pow_ui(resultat_expo, z_premier_N, i);
     mpz_add(z_somme,resultat_expo,z_somme);
   }

return;
}


// fonction contraintes 
//fonction B2

int contrainte_b2(mpz_t z_new,unsigned long int expo_new)
{	
	//int taille;
	mpz_t z_b2;
	mpz_init_set(z_b2,z_new);
	mpz_t z_racine_borne;
	mpz_init_set_ui(z_racine_borne,10);
	int cmp_b2;
		mpz_pow_ui(z_racine_borne, z_racine_borne, 80);
		mpz_pow_ui(z_b2, z_new, expo_new);
		cmp_b2=mpz_cmp(z_b2,z_racine_borne);//mpz_sizeinbase (z_b2, 10);//printf("\n taille de b2 = %i \n",taille);
if ((cmp_b2>=0)) { return 0; }
else { return 1;  }
}	

void initialisation_matrice(premier_exposant tableau_diviseur[45][45])
{
	int j,i;
	for (j = 0; j < 45 ; j++)
	{
		for (i = 0; i < 45 ; i++)
		{
		mpz_init_set_ui(tableau_diviseur[i][j].z_premier,1);
		tableau_diviseur[i][j].exposant = 1;
		}

	}
}

void initialisation_tableau_fixe(premier_exposant tableau_fixe[])
{
	int i;
	for (i = 0; i < 45 ; i++)
	{
	mpz_init_set_ui(tableau_fixe[i].z_premier,1);
	tableau_fixe[i].exposant = 1;
	tableau_fixe[i].cmpt_expo_list=0;
	}
}

void initialisation_tableau_ajuste(premier_exposant tableau_ajuste[2025])
{
	int i;
	for (i = 0; i < 2025 ; i++)
	{
	mpz_init_set_ui(tableau_ajuste[i].z_premier,1);
	tableau_ajuste[i].exposant = 1;
	}
}

void tableau_exposant (unsigned long int exposant[])
{
	mpz_t z_e,z_p;
	mpz_init_set_ui(z_e,2);
	mpz_init(z_p);
	int entier_p,cmpt_expo_list;
	for (cmpt_expo_list = 0;cmpt_expo_list < 40; cmpt_expo_list++)
	{
 	 mpz_nextprime (z_p,z_e);
 	 entier_p = mpz_get_ui(z_p);
 	 exposant[cmpt_expo_list] = entier_p-1;
 	 mpz_set(z_e,z_p);
	 	
	}
}

void remplissage_tableau_10000_premier(mpz_t tableau_10000_premier[])
{
	int i;
	unsigned long int variable_temporaire;
	
	FILE* fichier = NULL;
	fichier = fopen("premier.txt", "r");
	if (fichier != NULL)
	{
		for(i = 0; i < 10000 ;i++)
		{
		mpz_init(tableau_10000_premier[i]);
		fscanf(fichier, "%ld", &variable_temporaire); 	
		mpz_set_ui(tableau_10000_premier[i],variable_temporaire); 
		
        	}
	}	
	fclose(fichier);
}

void tri_tableau(premier_exposant tableau_trie[2025],premier_exposant tableau_ajuste[2025])
{

int i,j,variable_comparaison;
unsigned long int variable_temporaire_exposant;
mpz_t variable_temporaire;
mpz_init_set_ui(variable_temporaire,1);

for (i = 0; i < 2025; i++)
     {
      for (j = i+1;j < 2025; j++)
	{
	variable_comparaison = mpz_cmp(tableau_ajuste[j].z_premier,tableau_ajuste[i].z_premier);
          
		if (variable_comparaison > 0) 
		{		
		
      		mpz_set(variable_temporaire,tableau_ajuste[i].z_premier);
      		mpz_set(tableau_ajuste[i].z_premier,tableau_ajuste[j].z_premier);
      		mpz_set(tableau_ajuste[j].z_premier,variable_temporaire);
		variable_temporaire_exposant = tableau_ajuste[i].exposant;
		tableau_ajuste[i].exposant = tableau_ajuste[j].exposant ;
		tableau_ajuste[j].exposant = variable_temporaire_exposant;
      		}
      
        }
     }
for(i = 0;i<2025;i++)
{
	mpz_set(tableau_trie[i].z_premier,tableau_ajuste[i].z_premier);
	tableau_trie[i].exposant = tableau_ajuste[i].exposant;
}
}


int contrainte_b1(premier_exposant tableau_trie[2025],premier_exposant tableau_fixe[8],int j,mpz_t z_b1, int special)
{

	int i,k,comp,comp2;
	mpz_t tableau_temporaire[j];
	mpz_t z_m;
	
 	mpz_init_set_ui(z_m,1);
	for(i=0;i<=j;i++)
	{ mpz_init_set_ui(tableau_temporaire[i],1);
	}
	
	for(i=0;i<=j;i++)
	{
	mpz_pow_ui(z_m,tableau_fixe[i].z_premier,tableau_fixe[i].exposant);
	//gmp_printf("\n*****Le calcul de l'exponentiation est : %Zu ^% ld \n",tableau_fixe[i].z_premier,tableau_fixe[i].exposant);
	mpz_mul(z_b1,z_b1,z_m);
	//gmp_printf("z_b1 à l'étape %i = %Zu \n",i,z_b1);
	mpz_set(tableau_temporaire[i],tableau_fixe[i].z_premier);
	mpz_init_set_ui(z_m,1);	
	}
	k=0;
	while(k<2025)
	{
	 for(i=0;i<=j;i++)
		{
			comp=mpz_cmp(tableau_trie[k].z_premier,tableau_temporaire[i]);
			if(comp==0)
			{
			break;}
		}
			
	comp2=mpz_cmp_ui(tableau_trie[k].z_premier,1);	
	if((comp==0)||(comp2==0))
	{ k++;}
	else
	{ mpz_pow_ui(z_m,tableau_trie[k].z_premier,tableau_trie[k].exposant);
	//gmp_printf("\n*****Le calcul de l'exponentiation est : %Zu ^% ld \n",tableau_trie[k].z_premier,tableau_trie[k].exposant);
	mpz_mul(z_b1,z_b1,z_m);k++;//gmp_printf("z_b1 à l'étape %i = %Zu \n",k,z_b1);
	}
	
	}
mpz_set_ui(z_m,10);
mpz_pow_ui(z_m,z_m,160);
comp=mpz_cmp(z_m,z_b1);
return comp;
} 

// test sur b3
int contrainte_b3(mpz_t z_b3)
{
	int comp;
	mpz_t z_m;

  mpz_init(z_m);
  mpz_set_ui(z_m,10);

  mpz_pow_ui(z_m,z_m,160);

comp = mpz_cmp(z_m,z_b3);
return comp;
}

// containte s
int contrainte_s(premier_exposant tableau_trie[2025],premier_exposant tableau_fixe[8],int j,mpz_t z_s,mpz_t z_b1, int special)
{

	int i,k,comp,comp2;
	mpz_t tableau_temporaire[j];
	mpz_t z_m;
	mpz_t div_s;
	mpz_init_set_ui(div_s,1);
 	mpz_init_set_ui(z_m,1);
	for(i = 0;i <= j;i++)
	{ 
	mpz_init_set_ui(tableau_temporaire[i],1);
	}
	
	for(i = 0;i <= j;i++)
	{
	//mpz_pow_ui(z_m,tableau_fixe[i].z_premier,tableau_fixe[i].exposant);
	//gmp_printf("\n ON CALCULE LA FONCTION SIGMA POUR S DE : %Zu à exposant %ld \n",tableau_fixe[i].z_premier,tableau_fixe[i].exposant);
	mpz_sigma(z_m,tableau_fixe[i].z_premier,tableau_fixe[i].exposant);
	mpz_mul(z_s,z_s,z_m);
	//gmp_printf("z_s à l'étape %i = %Zu \n",i,z_s);
	mpz_set(tableau_temporaire[i],tableau_fixe[i].z_premier);
	}
	k = 0;
	while(k < 2025)
	{
	 for(i=0;i<=j;i++)
		{
			comp=mpz_cmp(tableau_trie[k].z_premier,tableau_temporaire[i]);
			if(comp==0)
			{
			break;}
		}
			
	comp2=mpz_cmp_ui(tableau_trie[k].z_premier,1);	
	if((comp==0)||(comp2==0))
	{ k++;}
	else
	{ //mpz_pow_ui(z_m,tableau_trie[k].z_premier,tableau_trie[k].exposant);
	//gmp_printf("\n *ON CALCULE LA FONCTION SIGMA POUR S DE : %Zu à exposant %ld \n",tableau_trie[k].z_premier,tableau_trie[k].exposant);
	mpz_sigma(z_m,tableau_trie[k].z_premier,tableau_trie[k].exposant);
	mpz_mul(z_s,z_s,z_m);k++;
	//gmp_printf("z_s à l'étape %i = %Zu \n",k,z_s);
	}
	
	}


//gmp_printf("sigma rend %Zu et b1 rend %Zu\n",z_s,z_b1);
mpz_mul_ui(z_b1,z_b1,2);
comp = mpz_cmp(z_s,z_b1);

return comp;
} 



