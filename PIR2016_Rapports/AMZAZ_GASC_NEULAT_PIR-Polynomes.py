def ordre(n,p): #calcul de l'ordre mult de n modulo p
    N=n;
    o=1;
    while N!=1:
        N=(N*n)%p;
        o+=1;
    return(o);
    

#ici on va ecrire tout les poly de deg d modulo p, pour faire toutes les SRL possibles
def Poly(d,p):
    A=[];
    if d==0:
        for i in range(0,p):
            A.append([i]);
    else:
        for i in range(0,p):
            for j in Poly(d-1,p):
                A.append([i]+j);
                
    return(A);

#ainsi on peut faire tout les poly de degres 3 et regarder la SRL qui leur est associée


def Periode (a,b,c,p):
    EtatInitial=[0,0,1]; #U0 U1 U2
    periode=1;
    Etatn=[EtatInitial[1],EtatInitial[2],(c*EtatInitial[0]+b*EtatInitial[1]+a*EtatInitial[2])%p]; #U1 U2 U3
    if c%p == 0:
        return(0);
    while Etatn != EtatInitial:
        Etatn=[Etatn[1],Etatn[2],(c*Etatn[0]+b*Etatn[1]+a*Etatn[2])%p]; #ici on passe de l'etat n, a l'etat n+1 CAD la matrice Etatn etait egale a [Un,Un+1,Un+2] elle vaut, apres cette etape [Un+1,Un+2,Un+3]
        periode=periode+1;
    return(periode);

def PeriodeSRLPoly(p):
    A=Poly(2,p); #On genere tout les poly de deg 3 modulo p
    B=[];
    for i in A:
        B.append("X^3 = "+str(i[0])+"X² + "+str(i[1])+"X + "+str(i[2])+" : "+str(Periode(i[0],i[1],i[2],p))); # Ca affiche X^3=aX²+bX+C : Tp , ou Tp est la periode de la suite avec ce poly comme poly carracteristique
    return(B);
    
    
