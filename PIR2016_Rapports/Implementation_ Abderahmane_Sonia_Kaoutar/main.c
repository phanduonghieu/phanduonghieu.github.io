#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

#include"gmp.h"
#include"ecm.h"
#include"parfait.h"


int main()
{

// variable du tableau des 8 premiers
 mpz_t tab_premiers[8];

 // variable du tableau des 8 premiers de la contrainte D
mpz_t tab_premiers_D[8];

// initialisation et declaration des 8 nombres premiers
 mpz_init_set_ui(tab_premiers[0],127);
 mpz_init_set_ui(tab_premiers[1],19);
 mpz_init_set_ui(tab_premiers[2],7);
 mpz_init_set_ui(tab_premiers[3],11); 
 mpz_init_set_ui(tab_premiers[4],31);
 mpz_init_set_ui(tab_premiers[5],13);
 mpz_init_set_ui(tab_premiers[6],3);
 mpz_init_set_ui(tab_premiers[7],5);

// initialisation et declaration des 8 nombres premiers de la contrainte D
 mpz_init_set_ui(tab_premiers_D[0],0);
 mpz_init_set_ui(tab_premiers_D[1],0);
 mpz_init_set_ui(tab_premiers_D[2],0);
 mpz_init_set_ui(tab_premiers_D[3],0); 
 mpz_init_set_ui(tab_premiers_D[4],0);
 mpz_init_set_ui(tab_premiers_D[5],0);
 mpz_init_set_ui(tab_premiers_D[6],0);
 mpz_init_set_ui(tab_premiers_D[7],0);

//Ouverture du fichier contenant les 10000 premiers et les mettre dans un tableau de taille [10000]
int nb_facteur,test_composante,test_2;
mpz_t tableau_10000_premier[10000];
remplissage_tableau_10000_premier(tableau_10000_premier);

//Initialisation des elements de la structure à 1 et declaration de la matrice contenant les diviseurs ainsi que du vecteur "fixe"
premier_exposant tableau_diviseur[45][45];
premier_exposant tableau_fixe[45];
premier_exposant tableau_ajuste[2025]; premier_exposant tableau_ajuste_pour_b3[2025];
initialisation_matrice(tableau_diviseur);
initialisation_tableau_fixe(tableau_fixe);

//Creation du tableau contenant les exposants pairs jusqu a 178
unsigned long int exposant[40];
tableau_exposant (exposant);
int j,i,resultat_b1,contrainte_d,id;
j=0;

// Debut du programme
mpz_t z_res_mod,z_quatre;
mpz_init(z_res_mod);
mpz_init_set_ui(z_quatre,4);
double B1;
unsigned long int expo_max;
int resultat_b2;
mpz_t z_max; 

//test sur les composant pour l'eriture d'une nombre parfait impair
mpz_t z_prime;
mpz_init(z_prime);
mpz_set(z_prime,tab_premiers[j]);
int special;
special=0;
test_composante=test_sur_composante(z_prime,special);

// Passage au tableau fixe
if(test_composante==0)
{
mpz_set(tableau_fixe[j].z_premier,tab_premiers[j]);
tableau_fixe[j].exposant=exposant[tableau_fixe[j].cmpt_expo_list];
}
else
{
mpz_set(tableau_fixe[j].z_premier,tab_premiers[j]);
tableau_fixe[j].exposant=1;
if(special ==0){special=1;}
}
// variables
mpz_t z_b3; mpz_t z_x,z_y;
mpz_init(z_x); mpz_init_set_ui(z_y,1);
mpz_t z_b3_initiale;
mpz_t z_b1;
mpz_t z_somme;
mpz_t z_diviseur;
int resultat_b3;
mpz_t z_s;
int resultat_s;
int variable_temporaire;
int resultat_facto_ecm;
int courbe_elliptique;
int cmp_max;
int compteur_ajuste,ia,ja;
int test_expo,test_expo2;
int test_0_fixe,reponse_premier;
resultat_s=0;

mpz_t z_special;
mpz_init(z_special); int test_special;

int k;
//la boucle
do
{nb_facteur=0;
//initialisation du vecteur de la matrice
for(k=j;k<45;k++)
{
	for(i=0;i<45;i++)
	{ mpz_set_ui(tableau_diviseur[i][k].z_premier,1);
 	  tableau_diviseur[i][k].exposant=1;
	}
}
premier_exposant tableau_trie[2025];
premier_exposant tableau_trie_pour_b3[2025];
	for (ia = 0; ia < 2025; ia++)
	{
	mpz_init_set_ui(tableau_trie[ia].z_premier,1);
	tableau_trie[ia].exposant = 0;
	}
	


	for (ia = 0; ia < 2025; ia++)
	{
	mpz_init_set_ui(tableau_ajuste[ia].z_premier,1);
	tableau_ajuste[ia].exposant = 0;
	}

	for (ia = 0; ia < 2025; ia++)
	{
	mpz_init_set_ui(tableau_ajuste_pour_b3[ia].z_premier,1);
	tableau_ajuste_pour_b3[ia].exposant = 0;
	}

printf("Nous sommes à l'étape J=%i \n",j);
mpz_init_set_ui(z_max,1);

// application de la contrainte D
for (id=0;id<8;id++)
 {	contrainte_d=mpz_cmp(tab_premiers_D[id],tableau_fixe[j].z_premier);
 		if (contrainte_d == 0)
		 {
			printf("contrainte D détectée\n");
		        break;
		 }
 } 

if (contrainte_d != 0)
{
	{
	 printf("contrainte D passée\n");
	}
//application de la contrainte B2

gmp_printf("----->[%Zu ^ %li ] \n",tableau_fixe[j].z_premier,tableau_fixe[j].exposant);
resultat_b2= contrainte_b2(tableau_fixe[j].z_premier,tableau_fixe[j].exposant);

if (resultat_b2!=0)
{	
	printf("contrainte B2 passée\n");

// test sigma
mpz_init(z_somme);
mpz_sigma(z_somme,tableau_fixe[j].z_premier,tableau_fixe[j].exposant);
gmp_printf("Sigma de %Zu^%ld est %Zu\n",tableau_fixe[j].z_premier,tableau_fixe[j].exposant, z_somme);

//test primalité
reponse_premier = mpz_probab_prime_p (z_somme, 5);
 	//printf("\nLA REPONSE DE mpz_probab prime est : %d \n",reponse_premier);
	if (reponse_premier != 0)
	{
		// le nombre est premier, on l'ecrit directement dans la matrice
	test_2=mpz_cmp_ui(z_somme,2);
	if(test_2!=0)
	{		
	mpz_set(tableau_diviseur[nb_facteur][j].z_premier,z_somme);
	cmp_max=mpz_cmp(z_max,tableau_diviseur[nb_facteur][j].z_premier);
	gmp_printf("---------->[ %Zu ]\n",z_somme);
				if(cmp_max == -1)
					{mpz_set(z_max,tableau_diviseur[nb_facteur][j].z_premier);
					//gmp_printf("\n Le max ici est: %Zu",z_max);
					}	
	mpz_set(z_max,tableau_diviseur[nb_facteur][j].z_premier);nb_facteur++;
	}
	}
else
{
//test facto
nb_facteur = factoriser(z_somme,tableau_10000_premier,tableau_diviseur,z_max,j);
//printf("\nLe nombre de facteur est :%d \n",nb_facteur);
//gmp_printf("\n Z_somme est égal à %Zu \n",z_somme);
variable_temporaire = mpz_cmp_ui(z_somme,3);
//gmp_printf("\nZ_somme = %Zu\n",z_somme);
//printf("variable temporaire = %i \n",variable_temporaire);

//second test de primalité
reponse_premier = mpz_probab_prime_p (z_somme, 5);
 	//printf("\nLA REPONSE DE mpz_probab prime est : %d \n",reponse_premier);
	if (reponse_premier != 0)
	{test_2=mpz_cmp_ui(z_somme,2);
	if(test_2!=0)
	{
		// le nombre est premier, on l'ecrit directement dans la matrice
		mpz_set(tableau_diviseur[nb_facteur][j].z_premier,z_somme);
		cmp_max=mpz_cmp(z_max,tableau_diviseur[nb_facteur][j].z_premier);
		gmp_printf("----------->[ %Zu ]\n",z_somme);
				if(cmp_max == -1)
					{mpz_set(z_max,tableau_diviseur[nb_facteur][j].z_premier);
					//gmp_printf("\n Le max ici est: %Zu",z_max);
					}
	mpz_set(z_max,tableau_diviseur[nb_facteur][j].z_premier);nb_facteur++;
	}}
else
{
while (variable_temporaire > 0)//Grand while !
{
	//test ecm
	
	mpz_init(z_diviseur);
	
	courbe_elliptique = 0;
	
	B1 =75000;
 	do
  	{
    		courbe_elliptique++;
    		printf(".\n");
		resultat_facto_ecm = ecm_factor(z_diviseur,z_somme,B1,NULL);	
		
		//printf("nb_elliptique = %u : ecm_factor a retournée %d\n", courbe_elliptique, resultat_facto_ecm);  
  		if ((courbe_elliptique %20)==0)
		{B1=B1*2;printf("+");
		}
		}
  	while (resultat_facto_ecm == 0);
  	if (resultat_facto_ecm < 0)
   		printf("\n Une erreur s'est poduite lors de l'éxecution\n");
 	else
   	{ 	reponse_premier = mpz_probab_prime_p (z_diviseur, 2);
if (reponse_premier!=0)		
		{gmp_printf("---------->ecm[ %Zu ]\n",z_diviseur);
			
		mpz_set(tableau_diviseur[nb_facteur][j].z_premier,z_diviseur);
		
		cmp_max=mpz_cmp(z_max,z_diviseur);
gmp_printf("----------->[ %Zu ]\n",z_somme);
				if(cmp_max == -1)
					{mpz_set(z_max,z_diviseur);
					//gmp_printf("\n Le max ici est: %Zu",z_max);
					}	
		mpz_divexact(z_somme,z_somme,z_diviseur);
		nb_facteur++;
		variable_temporaire = mpz_cmp_ui(z_somme,1);
		// dernier test de primalité
		reponse_premier = mpz_probab_prime_p (z_somme, 5);
 	//printf("\nLA REPONSE DE mpz_probab prime est : %d \n",reponse_premier);
	if (reponse_premier != 0)
	{test_2=mpz_cmp_ui(z_somme,2);
	if(test_2!=0)
	{
		// le nombre est premier, on l'ecrit directement dans la matrice
		mpz_set(tableau_diviseur[nb_facteur][j].z_premier,z_somme);
		cmp_max=mpz_cmp(z_max,tableau_diviseur[nb_facteur][j].z_premier);

				if(cmp_max == -1)
					{mpz_set(z_max,tableau_diviseur[nb_facteur][j].z_premier);
					//gmp_printf("\n Le max ici est: %Zu",z_max);
					}
	;nb_facteur++;variable_temporaire=0;
	}
	else
	{variable_temporaire=0;}
	}
		} 

}}

}
}
//printf("\n");
// affichage de la colone des diviseurs trouvés
expo_max=0;

for(i = 0; i < 7; i++)
{
	//gmp_printf("\nle vecteur à la case %d est: %Zu ",i,tableau_diviseur[i][j].z_premier);
	//printf("\t\tl'exposant à la case %d est:%li ",i,tableau_diviseur[i][j].exposant);
	cmp_max=mpz_cmp(z_max,tableau_diviseur[i][j].z_premier);
	if (cmp_max == 0)
	{//printf("\n expo_max du tableau = %li va être attribué à expo_max = %li \n",tableau_diviseur[i][j].exposant,expo_max);
	expo_max= expo_max + tableau_diviseur[i][j].exposant;	
	}
}	
//gmp_printf("\nLe max est: %Zu \b son exposant est %li\n",z_max,expo_max);

//Remplissage du vecteur qui servira pour le calcul des exposants ajustés et les contraintes
 compteur_ajuste=0,ia=0;


	for(ia=0;ia<=j;ia++)
	{
  		for(ja=0;ja<45;ja++)
  		{
  		mpz_set(tableau_ajuste[compteur_ajuste].z_premier,tableau_diviseur[ja][ia].z_premier);
		(tableau_ajuste[compteur_ajuste].exposant)=(tableau_diviseur[ja][ia].exposant);
   		//gmp_printf("\n la case %i eme du tableau à ajuster est: %Zu \n",compteur_ajuste,tableau_ajuste[compteur_ajuste].z_premier);
		compteur_ajuste++;
   		}
	}

// remplissage du vecteur qui servira à calculer b3
compteur_ajuste=0;
for(ia=0;ia<j;ia++)
	{
  		for(ja=0;ja<j;ja++)
  		{
  		mpz_set(tableau_ajuste_pour_b3[compteur_ajuste].z_premier,tableau_diviseur[ja][ia].z_premier);
		(tableau_ajuste_pour_b3[compteur_ajuste].exposant)=(tableau_diviseur[ja][ia].exposant);
   		//gmp_printf("\n la case %i eme du tableau à ajuster est: %Zu \n",compteur_ajuste,tableau_ajuste[compteur_ajuste].z_premier);
		compteur_ajuste++;
   		}
	}


//printf("\n coucou\n");
tri_tableau(tableau_trie,tableau_ajuste);
	/*for(ia = 0;ia < 25; ia++)
  		{
   		gmp_printf("\n la case %i eme du tableau trié est: %Zu \n",ia,tableau_trie[ia].z_premier);
		printf("\n la case %i eme du tableau trié a pour exposant: %ld \n",ia,tableau_trie[ia].exposant);
		}*/
//printf("bonjour\n");	


tri_tableau(tableau_trie_pour_b3,tableau_ajuste_pour_b3);
// tasser les éléments du tableau trié
compteur_ajuste = 0;
for(ia = 1;ia < 2025;ia++)
{	while(tableau_trie[compteur_ajuste].exposant==0)
	{compteur_ajuste++;}
  test_expo=mpz_cmp(tableau_trie[compteur_ajuste].z_premier,tableau_trie[ia].z_premier);
if(test_expo == 0)
{tableau_trie[compteur_ajuste].exposant=(tableau_trie[compteur_ajuste].exposant)+(tableau_trie[ia].exposant);
 tableau_trie[ia].exposant=0;//à changer peut être?
mpz_set_ui(tableau_trie[ia].z_premier,1);}
else
{compteur_ajuste++;}
}
tri_tableau(tableau_trie,tableau_trie);


/*for(ia = 0;ia < 25; ia++)
  		{
   		gmp_printf("|%Zu ^ %ld |",tableau_trie[ia].z_premier,tableau_trie[ia].exposant);
		//printf("\n la case %i eme du tableau trié a pour exposant: %ld \n",ia,tableau_trie[ia].exposant);
		} */

//tasser les éléments du tableau trié de b3
compteur_ajuste = 0;
for(ia = 1;ia < 2025;ia++)
{	
  test_expo=mpz_cmp(tableau_trie_pour_b3[compteur_ajuste].z_premier,tableau_trie_pour_b3[ia].z_premier);
if(test_expo==0)
{tableau_trie_pour_b3[compteur_ajuste].exposant=(tableau_trie_pour_b3[compteur_ajuste].exposant)+(tableau_trie_pour_b3[ia].exposant);
 tableau_trie_pour_b3[ia].exposant=0;//à changer peut être?
mpz_set_ui(tableau_trie_pour_b3[ia].z_premier,1);ia++;}
else
{compteur_ajuste++;}
}
tri_tableau(tableau_trie_pour_b3,tableau_trie_pour_b3);

for(ia = 0;ia < 25; ia++)
  		{
   		//gmp_printf("|%Zu ^ %ld |",tableau_trie_pour_b3[ia].z_premier,tableau_trie_pour_b3[ia].exposant);
		//printf("\n la case %i eme du tableau trié a pour exposant: %ld \n",ia,tableau_trie_pour_b3[ia].exposant);
		} 

//ajuster l'exposant
ia=0;ja=0; 
if (special ==0)
{while(ia<45)
{
 mpz_mod(z_res_mod,tableau_trie[ia].z_premier,z_quatre);
 test_expo2=mpz_cmp_ui(z_res_mod,1); 
//gmp_printf("\n le test sur %Zu = %i \n ",tableau_trie[ia].z_premier,test_expo2);
 if( (test_expo2!=0))
   {//printf("\n ici");
	if(tableau_trie[ia].exposant <= exposant[ja])
	  {tableau_trie[ia].exposant=exposant[ja];
		//gmp_printf("%i l'exposant de %Zu est %li \n",ia,tableau_trie[ia].z_premier,tableau_trie[ia].exposant);
	 ia++;ja=0;}
	else
        {ja++;}
   }
   else
   {ia++;}
}
}
else
{
while(ia<2025)
{
test_expo=mpz_cmp(tableau_trie[ia].z_premier,z_special); 
//printf("\n %i là",ia);
//gmp_printf("\n le test entre %Zu et %Zu (fixe) = %i \n",tableau_trie[ia].z_premier,tableau_fixe[j].z_premier, test_expo);
   if (test_expo!=0)
   {
	if(tableau_trie[ia].exposant <= exposant[ja])
	  {tableau_trie[ia].exposant=exposant[ja];
	//gmp_printf("l'exposant de %Zu est %li \n",tableau_trie[ia].z_premier,tableau_trie[ia].exposant);
	ia++;ja=0;}
	else
        {ja++;}
   }
   else
   {ia++;}
}
}



//ajuster l'exposant pour b3
ia=0;ja=0; 
if (special ==0)
{while(ia<45)
{
 mpz_mod(z_res_mod,tableau_trie_pour_b3[ia].z_premier,z_quatre);
 test_expo2=mpz_cmp_ui(z_res_mod,1); 
//gmp_printf("\n le test sur %Zu = %i \n ",tableau_trie_pour_b3[ia].z_premier,test_expo2);
 if( (test_expo2!=0))
   {//printf("\n ici");
	if(tableau_trie_pour_b3[ia].exposant <= exposant[ja])
	  {tableau_trie_pour_b3[ia].exposant=exposant[ja];
		//gmp_printf("%i l'exposant de %Zu est %li \n",ia,tableau_trie_pour_b3[ia].z_premier,tableau_trie_pour_b3[ia].exposant);
	 ia++;ja=0;}
	else
        {ja++;}
   }
   else
   {ia++;}
}
}
else
{
while(ia<2025)
{
test_expo=mpz_cmp(tableau_trie_pour_b3[ia].z_premier,z_special); 
//printf("\n %i là",ia);
//gmp_printf("\n le test entre %Zu et %Zu (fixe) = %i \n",tableau_trie_pour_b3[ia].z_premier,tableau_fixe[j].z_premier, test_expo);
   if (test_expo!=0)
   {
	if(tableau_trie_pour_b3[ia].exposant <= exposant[ja])
	  {tableau_trie_pour_b3[ia].exposant=exposant[ja];
	//gmp_printf("l'exposant de %Zu est %li \n",tableau_trie_pour_b3[ia].z_premier,tableau_trie_pour_b3[ia].exposant);
	ia++;ja=0;}
	else
        {ja++;}
   }
   else
   {ia++;}
}
}


//printf("\n");
 printf("tableau avec exposants ajustés\n");
for(ia = 0;ia < 20; ia++)
  		{
   		gmp_printf(" [%Zu ^ %ld] \t",tableau_trie_pour_b3[ia].z_premier,tableau_trie_pour_b3[ia].exposant);
		if(tableau_trie_pour_b3[ia+1].exposant>70)
		{break;}		
		}
printf("\n");
//Contrainte B3

	mpz_init_set_ui(z_b1,1);

	mpz_init_set_ui(z_b3_initiale,1);
if (j>0)
{
contrainte_b1(tableau_trie_pour_b3,tableau_fixe,j,z_b3_initiale,special);
mpz_set(z_b3,z_b3_initiale);
resultat_b3=contrainte_b3(z_b3);
//gmp_printf("b3 = %Zu\n",z_b3);
}
else
{resultat_b3=1;}
if(resultat_b3>0)
{	//printf("contrainte B3 passée\n");    
		// Contrainte B1
	resultat_b1=contrainte_b1(tableau_trie,tableau_fixe,j,z_b1,special);
	//printf("résultat b1 = %i \n",resultat_b1);
	//gmp_printf("b1 = %Zu\n",z_b1);
int taille;	
	taille=mpz_sizeinbase (z_b1, 10);
printf("taille de b1 = %i \n",taille);
	if(resultat_b1>0)
	{//printf("contrainte B1 passée\n");
		//contrainte S
		mpz_init_set_ui(z_s,1);
		resultat_s=contrainte_s(tableau_trie,tableau_fixe,j,z_s,z_b1,special);
		//printf("\n le résultat S est :%i \n ",resultat_s);
		if(resultat_s<0)
		{//printf("contrainte s passée\n");
		mpz_set(z_prime,z_max);
			j++;
		tableau_fixe[j].cmpt_expo_list=0;
		printf("j = %i \n",j);
		test_composante=test_sur_composante(z_prime, special);
	if (test_composante==1)
	{mpz_set(z_special,z_prime);}
				
		// Passage au tableau fixe
		if((test_composante==0) || (special==1))
		{
		mpz_set(tableau_fixe[j].z_premier,z_prime);
		tableau_fixe[j].exposant=exposant[tableau_fixe[j].cmpt_expo_list];
		}
		else
		{
		mpz_set(tableau_fixe[j].z_premier,z_prime);
		tableau_fixe[j].exposant=1;tableau_fixe[j].cmpt_expo_list--;
		if(special ==0){special=1;printf("compteur special = %i \n",special);}
		}
		//gmp_printf("le nouvel élément du tableau fixe est %Zu\n",tableau_fixe[j].z_premier);
		
		}
		else 
		{
		if(resultat_s ==0)
		{//printf("le nombre est parfait\n");
		}
		else
		{printf("contrainte S détectée !!!\n");
		printf("______________________________\n");
		tableau_fixe[j].cmpt_expo_list++;
		tableau_fixe[j].exposant=exposant[tableau_fixe[j].cmpt_expo_list];
		// special
		test_special=mpz_cmp(tableau_fixe[j].z_premier,z_special);
		//printf("test special = %i \n",test_special);
		if((test_special==0)&&(tableau_fixe[j].exposant>=2))
		{special=0; //printf("special = %i",special);
		}	
		}
	}
	}
	else
	{printf("contrainte B1 détectée !!!\n");printf("______________________________\n");
	tableau_fixe[j].cmpt_expo_list++;
		tableau_fixe[j].exposant=exposant[tableau_fixe[j].cmpt_expo_list];
	//special
	test_special=mpz_cmp(tableau_fixe[j].z_premier,z_special);
	//printf("test special = %i \n",test_special);
	if((test_special==0)&&(tableau_fixe[j].exposant>=2))
	{//special=0; printf("special = %i",special);
	}		
	}               
}
else
{
printf("contrainte B3 détectée !!!\n");printf("______________________________\n");

	mpz_set_ui(tableau_fixe[j].z_premier,1);
	j--; 	tableau_fixe[j].cmpt_expo_list++;
	if(j>=0)
	{tableau_fixe[j].exposant=exposant[tableau_fixe[j].cmpt_expo_list];}
	test_special=mpz_cmp(tableau_fixe[j].z_premier,z_special);
printf("j = %i \n",j);
	//gmp_printf("fixe[%Zu] - special[%Zu]",tableau_fixe[j].z_premier,z_special);	
	//printf("test special = %i \n",test_special);
	if((test_special==0)&&(tableau_fixe[j].exposant>=2))
	{//special=0; printf("special = %i",special);
	}	
}//ici ici ici


}// accolade B2
else
{	
	printf("contrainte B2 detectée !!!\n");printf("______________________________\n");
	mpz_set_ui(tableau_fixe[j].z_premier,1);
	j--; 	tableau_fixe[j].cmpt_expo_list++; 
	if(j!=0)
	{tableau_fixe[j].exposant=exposant[tableau_fixe[j].cmpt_expo_list];}
	test_special=mpz_cmp(tableau_fixe[j].z_premier,z_special);
	//printf("test special = %i \n",test_special);
	if((test_special==0)&&(tableau_fixe[j].exposant>=2))
	{special=0; printf("special = %i",special);
	}	
}

}//accolade D
else
{
	printf("contrainte D non passée \n");
	mpz_set_ui(tableau_fixe[j].z_premier,1);
	j--; tableau_fixe[j].cmpt_expo_list++;
	if(j!=0)
	{ tableau_fixe[j].exposant=exposant[tableau_fixe[j].cmpt_expo_list]; }	
}
test_0_fixe=mpz_cmp_ui(tableau_fixe[j].z_premier,1);
if(test_0_fixe==0)
{
j=-1;
}
}while(j>=0);//end while
//fin de l'algorithme
}
