from math import sqrt
import matplotlib.pyplot as plt
plt.ion()
#Le but est de calculer les periodes des Un+3= aUn+2 + bUn+1 + cUn modulo un grand nombre de nombres premiers, l'idée est de creer ensuite un histogramme pour mieux visualier les donnée
#Il faut faire attention aux valeurs de a,b,c en effet, si par exemple c=0 modulo p alors la suite avec comme etat initial u0=1 u1=0 u2=0 aura pour valeurs 1 0 0 0 0.... 
#"inserer une étude sur l'etat initial et les valeurs de a,b,c "interdites"-> Pour que la suite ne soit pas constance égale a 0 acr( a partir d'un certaing rang ) il faut que c != o modulo p
def Periode (a,b,c,p): #On calcule la periode modulo p de la suite Un+3= aUn+2 + bUn+1 + cUn
    EtatInitial=[0,0,1]; #U0 U1 U2
    periode=1;
    Etatn=[EtatInitial[1],EtatInitial[2],(c*EtatInitial[0]+b*EtatInitial[1]+a*EtatInitial[2])%p]; #U1 U2 U3
    if c%p == 0:
        return(0);
    while Etatn != EtatInitial:
        Etatn=[Etatn[1],Etatn[2],(c*Etatn[0]+b*Etatn[1]+a*Etatn[2])%p]; #ici on passe de l'etat n, a l'etat n+1 CAD la matrice Etatn etait egale a [Un,Un+1,Un+2] elle vaut, apres cette etape [Un+1,Un+2,Un+3]
        periode=periode+1;
    return(periode);

def Premier(n): #ici je vais regarder si n est premier, ensuite je vais faire une fonction qui liste tout les nombres premiers jusqu'a un certain nombre
    if n==1:
        return(False);
    else :
        for i in range(2,int(sqrt(n)+1)) :
            if ( (n%i)== 0):
                return(False);
    return(True);
    
    
def ListePremier(n): #Ici je vais faire la liste des nombres premiers plus petits que n
    a=[];
    for i in range(1,n+1) :
        if Premier(i):
            a.append(i);
    return(a);
    
def ListePeriodes(a,b,c,n): #Ici on calcule les periodes de la suite Un+3= aUn+2 + bUn+1 + cUn pour tout les premiers plus petits que n
    ListeP=ListePremier(n);
    periodes={};
    for i in ListeP:
        periodes[i]=Periode(a,b,c,i);
    return(periodes);
#Ici on renvoi un dico, c'est a dire une liste de couple nombre premier - periode modulo ce nombre, Il nous suffit d'appeller la methode .get(p) et cela nous renvoi la periode de la suite modulo p


def main(): #interface avec l'utilisateur
    
    print("On va calculer la periode modulo p de la suite Un+3= aUn+2 + bUn+1 + cUn pour tout les p plus petit que n, si la suite est constante, par convention la periode sera 0");
    print("Rentrez a");
    coefx2=int(input());
    print("Rentrez b");
    coefx1=int(input());
    print("Rentrez c");
    coefunit=int(input());
    print("Rentrez n");
    maximum=int(input());

    print(ListePeriodes(coefx2,coefx1,coefunit,maximum));

    LP=ListePeriodes(coefx2,coefx1,coefunit,maximum);
    Liste=ListePremier(maximum);

    Image=[];
    for i in Liste:
        Image.append(LP[i]);
    
    X=[];
    for i in range(1,maximum):
        X.append(i); 
        
    Y0=[];
    Y1=[];
    Y2=[];
    Y3=[];
    #ici je crée y=x^3 -1 je crée aussi y=x^1 -1 et y=x(x-1)
    for i in X:
        Y0.append(i-1);
        Y1.append(i*(i-1));
        Y2.append(i*i -1 );
        Y3.append(i*i*i -1 );
    
    #Ici on choisit celui qu'on veut afficher !
    
    plt.plot(X, Y0); # P-1
    plt.plot(X, Y1); # P(P-1)
    plt.plot(X, Y2); # P²-1
    plt.plot(X, Y3); # P^3 -1
    plt.plot(Liste, Image);#Ici on trace le graphe des periodes de la suite modulo P Nous allons aussi tracer les courbes des majorants naturels p²-1 p^3 -1 et p(p-1), pour voir quand ils sont atteints

