import math
import random
from matplotlib.pyplot import plot, pause, clf, scatter, show, axis

#On modélise l'espace par un carré de taille L dans lequel les électrons évoluent, partant du bord supérieur.
#v0y est négatif

def ListeElectrons(L, n, v0x, v0y): #L est la longueur du segment, n le nombre d'électrons envoyés. Correspond au bord supérieur. i*L/n est l'abscisse, L l'ordonnée, v0x est la vitesse en x et v0y la vitesse en y
    T = []
    for i in range (n): #On veut des électrons aux deux extrémités
        T.append([i*L/(n-1), L, v0x, v0y]) #où v0x et v0y sont à définir.
    return T

def Obstacle(L, n): #n est le nombre d'obstacles souhaités. On tire au hasard les deux coordonnées de chaque extrémité du segment.
    P = []
    for i in range (n):
        a, c = L*random.random(), L*random.random()
        b = a + L*random.random()/5
        d = c + L*random.random()/5
        if a > c:
            (a,b,c,d) = (c,d,a,b)
        P.append([a, b, c, d])
        plot([a, b], [c, d], 'r', linewidth = 2) #trace le segment obstacle
    return P

def ObstacleGlaceClassique(L, n): #modélise des échantillons figés dans de la glace classique
    P = []
    for i in range (n):
        a, c = L*random.random(), L*random.random()
        b = a + L*random.random()/5
        P.append([a, b, c, c])
        plot([a, b], [c, c], 'r', linewidth = 2) #trace le segment obstacle
    return P

def ListversTab(L, n): #convertit la liste en tableau
    P = Obstacle(L, n)
    tab1 = []
    tab2 = []
    for i in range(n):
        a = P[i][0]
        b = P[i][1]
        c = P[i][2]
        d = P[i][3]
        tab1.append(a)
        tab1.append(b)
        tab2.append(c)
        tab2.append(d)
    tab = [tab1, tab2]
    return tab
            
def plusBas(tab, n): #renvoie l’ordonnée du point le plus bas. En cas d’égalité, renvoie la plus petite des abscisses des deux points d'égale ordonnée
    j = 0
    for i in range(1, 2*n):
        if (tab[1][i] < tab[1][j]) or ((tab[1][i] == tab[1][j]) and (tab[0][i] < tab[0][j])):
            j = i
    return j
    
def orient(tab, i, j, k): #renvoie l'orientation du triangle formé par les trois points
    det = (tab[0][j] - tab[0][i]) * (tab[1][k] - tab[1][i]) - (tab[1][j] - tab[1][i]) * (tab[0][k] - tab[0][i])
    if det > 0:
        return 1
    elif det < 0:
        return -1
    else:
        return 0
    
def prochainPoint(tab, n, i): #renvoie l'indice du prochain point à insérer dans l'enveloppe convexe
    j = 0 if i else 1
    for k in range(2*n):
        if k != i and k != j and orient(tab, i, j, k) < 0:
            j = k
    return j
    
def enveloppe(tab, n): #renvoie la liste des sommets de l'enveloppe
    i = plusBas(tab, n)
    j = prochainPoint(tab, n, i)
    env = [i]
    while j != i:
        env.append(j)
        j = prochainPoint(tab, n, j)
    return env
        
def TabversList(tab, n): #convertit tab en liste lisible par Trajectoire
    P = enveloppe(tab, n)
    Q = []
    for i in (P) :
        j = prochainPoint(tab, n, i)
        Q.append([tab[0][i], tab[0][j], tab[1][i], tab[1][j]])
        plot([tab[0][i], tab[0][j]], [tab[1][i], tab[1][j]], 'b', linewidth = 2) #trace le segment obstacle
    return Q
    
def ObstacleNegativeStaining(L, n): #modélise des échantillons figés par negative staining par une enveloppe convexe    
    tab = ListversTab(L, n)
    P = TabversList(tab, n)
    return P

def ObstacleSort(O, P, nObstacle, i):  # range les obstacles de la liste O en fonction du plus proche au plus lointain de l'électron i
    H = []
    nO = len(O)
    for j in range(nO):
        M = (O[j][2] - O[j][3])/(O[j][0] - O[j][1]) #coefficient directeur de l'obstacle
        V = O[j][2] - M*O[j][0] #ordonnée à l'origine de l'obstacle
        D = math.fabs(P[i][1] - M*P[i][0] - V)/(math.sqrt(1 + M**2)) #distance de l'électron à l'obstacle
        H.append([D, O[j][0], O[j][1], O[j][2], O[j][3]])
    H.sort() #tri de H selon D (c'est un tri lexicographique : or on veut justement trier les obstacles en fonction de leurs distances à l'électron)
    for i in range(nO):
        O[i] = [H[i][1], H[i][2], H[i][3], H[i][4]] #compilation de la nouvelle liste obstacle
    return O
    
#e est l'intervalle de temps entre chaque actualisation de ListeElectron. Trajectoire renvoie la liste actualisée des emplacements et vitesse des électrons.
#Les trois fonctions suivantes calculent et affichent la trajectoire des électrons pour les différents types d'obstacles générés

def Afficher_Obstacle(O):
    for i in range(len(O)):
        a,b,c,d = (O[i][0],O[i][1],O[i][2],O[i][3])
        plot([a, b], [c, d], 'r', linewidth = 2) #trace le segment obstacle

def Intersection(I,J,O,K):  #renvoie un booléen qui indique si [I,J] intercepte [O,K]
    #l'équation de (IJ) est appelée y1 = C x + B
    #l'équation de (OK) est appelée y2 = M x + V
    def y1():
            C = (I[1]-J[1])/(I[0]-J[0])#coefficient directeur de la trajectoire
            B = I[1] - I[0]*C #ordonnée à l'origine de la trajectoire
            return (C,B)
    def y2():  
            M = (O[1]-K[1])/(O[0]-K[0])#coefficient directeur de la trajectoire
            V = O[1] - O[0]*M #ordonnée à l'origine de la trajectoire
            return (M,V)
    (C,B) = y1()
    (M,V) = y2()
    xeq = (V-B)/(C-M) #potentiel point d'intersection
    return xeq <= max(J[0],I[0]) and xeq <= max(K[0],O[0]) and min(K[0],O[0]) < xeq and min(J[0],I[0]) < xeq


Couleurs = ['b','g','r','c','m','y','k','w']

"""def Ranger_croiss(O,j): #fait en sorte que le jème obstacle de O soit représenté par un couple de points d'abscisses croissantes
    if O[j][0] < O[j][1]:
        tmpx = O[j][1]
        tmpy = O[j][3]
        O[j][1] = O[j][0]
        O[j][0] = tmpx
        O[j][3] = O[j][2]
        O[j][2] = tmpy"""

def TrajectoireGlaceClassique(L, e, n, nObstacle, v0x, v0y):
    
    pas_moyens = L/(v0x**2 + v0y**2)    #donne un nombre moyen de points utilisés pour tracer la trajectoire de chaque électron : utile pour l'affichage
    indic_affichage = pas_moyens//50  #on veut à peu près 50 points par électron 
    
    
    P = ListeElectrons(L, n, v0x, v0y)
    O = ObstacleGlaceClassique(L, nObstacle)
    
    
    for i in range (n):
        def Equation(P,i):  #renvoie le coefficient directeur et l'ordonnée à l'origine de la droite qui prolonge le segment trajectoire du ième électron
            C = P[i][3]/P[i][2]#coefficient directeur de la trajectoire
            B = P[i][1] - P[i][0]*C #ordonnée à l'origine de la trajectoire
            return (C,B)
        def EquationObstacle(O,j):  #renvoie le coefficient directeur et l'ordonnée à l'origine de la droite qui prolonge le j-ème obstacle
            M = (O[j][2]-O[j][3])/(O[j][0]-O[j][1])#coefficient directeur de l'obstacle
            V = O[j][2]-M*O[j][0]#ordonnée à l'origine de l'obstacle
            return (M,V)
        
        t_affichage = -1
        while (P[i][0] >= 0) and (P[i][0] <= L) and (P[i][1] >= 0) and (P[i][1] <= L): #tant qu'on n'a pas rencontré de paroi
            t_affichage += 1   
            if t_affichage % indic_affichage == 0:  
                 X = [P[i][0]]
                 Y = [P[i][1]]
                 scatter(X, Y, c = Couleurs[i % 8])    #chaque électron a une couleur
                 axis([0,5/4*L,0,5/4*L])
                 show()
                 pause(0.1)
                 if t_affichage % 5 == 0:   #on efface la mémoire de la trajectoire tous les 5 points
                     clf()
                     Afficher_Obstacle(O)
            O = ObstacleSort(O, P, nObstacle, i)
            j = 0
            while j != nObstacle + 1 :  # j = nObstacle va produire l'arrêt de la boucle par la première condition if
                if j == nObstacle:  # si aucun obstacle n'est rencontré par le ième électron
                    P[i][0] += P[i][2]*e
                    P[i][1] += P[i][3]*e
                    j = nObstacle + 1   # on sort de la boucle
                else:
                    (C,B) = Equation(P,i)
                    (M,V) = EquationObstacle(O,j)
                    if M-C != 0:
                        x = (B-V)/(M-C) #abscisse de l'éventuel point d'intersection de la trajectoire et de l'obstacle
                        I = P[i][0],P[i][1] #extremité début de trajectoire
                        mi = min(P[i][0],P[i][0]+e*P[i][2])
                        ma = max(P[i][0],P[i][0]+e*P[i][2]) # [mi,ma] décrit les valeurs du segment sur l'axe des abscisses
                        if  x > mi and x < ma and x < O[j][1] and x > O[j][0]: #condition de rencontre (intersection des segments). On suppose e suffisamment petit pour que le rebond se fasse après le déplacement et non à la collision.
                            y = C*x + B#ordonnée du point d'intersection
                            a2 = (O[j][0]-x)**2+(O[j][2]-y)**2
                            b2 = (I[0]-x)**2+(I[1]-y)**2
                            c2 = (I[1]-O[j][2])**2+(O[j][0]-I[0])**2
                        
                            
                            alpha = math.acos((a2+b2-c2)/(2*math.sqrt(a2*b2)))  #théorème d'Al-Kashi -> alpha est l'angle entre la trajectoire initiale et l'obstacle
                            d = abs(O[j][3]-y)
                            f = abs(O[j][1]-x)
                            beta = math.atan(d/f)
                            
                            v = e*math.sqrt((P[i][2])**2+(P[i][3])**2)
                            vy = math.sin(alpha+beta)*v
                            vx = math.cos(alpha+beta)*v
                        
                            #choix des signes : rustine au problème géomètrique
                            
                            choix = [[1,1],[1,-1],[-1,1],[-1,-1]]
                            boo = True
                            c = 0
                            while boo and c < 4:
                                vxn = choix[c][0]*vx
                                vyn = choix[c][1]*vy
                                
                                if vyn == -1*P[i][3] and vxn == -1*P[i][2]:   #cela correspondrait à un retour en arrière
                                    c+=1
                                else:
                                    P[i][0] = x + e*vxn #l'électron rebondit
                                    P[i][1] = y + e*vyn #""
                                    boo = Intersection(I,[P[i][0],P[i][1]],[O[j][0],O[j][2]],[O[j][1],O[j][3]])
                                    c+=1

                                    
                            P[i][2] = vxn #et a des vitesses conséquentes au rebond
                            P[i][3] = vyn
                                
                                
                            j = nObstacle + 1 #il y a eu rebond : l'électron s'est déplacé donc on peut sortir de la boucle
                        else :
                            j += 1  
                    else:
                        j += 1
    return P





def TrajectoireCryo(L, e, n, nObstacle, v0x, v0y):
    
    pas_moyens = L/(v0x**2 + v0y**2)    #donne un nombre moyen de points utilisés pour tracer la trajectoire de chaque électron : utile pour l'affichage
    indic_affichage = pas_moyens//50  #on veut à peu près 50 points par électron 
    
    
    P = ListeElectrons(L, n, v0x, v0y)
    O = Obstacle(L, nObstacle)
    
    
    for i in range (n):
        def Equation(P,i):  #renvoie le coefficient directeur et l'ordonnée à l'origine de la droite qui prolonge le segment trajectoire du ième électron
            C = P[i][3]/P[i][2]#coefficient directeur de la trajectoire
            B = P[i][1] - P[i][0]*C #ordonnée à l'origine de la trajectoire
            return (C,B)
        def EquationObstacle(O,j):  #renvoie le coefficient directeur et l'ordonnée à l'origine de la droite qui prolonge le j-ème obstacle
            M = (O[j][2]-O[j][3])/(O[j][0]-O[j][1])#coefficient directeur de l'obstacle
            V = O[j][2]-M*O[j][0]#ordonnée à l'origine de l'obstacle
            return (M,V)
        
        t_affichage = -1
        while (P[i][0] >= 0) and (P[i][0] <= L) and (P[i][1] >= 0) and (P[i][1] <= L): #tant qu'on n'a pas rencontré de paroi
            t_affichage += 1   
            if t_affichage % indic_affichage == 0:  
                 X = [P[i][0]]
                 Y = [P[i][1]]
                 scatter(X, Y, c = Couleurs[i % 8])    #chaque électron a une couleur
                 axis([0,5/4*L,0,5/4*L])
                 show()
                 pause(0.1)
                 if t_affichage % 5 == 0:   #on efface la mémoire de la trajectoire tous les 5 points
                     clf()
                     Afficher_Obstacle(O)
            O = ObstacleSort(O, P, nObstacle, i)
            j = 0
            while j != nObstacle + 1 :  # j = nObstacle va produire l'arrêt de la boucle par la première condition if
                if j == nObstacle:  # si aucun obstacle n'est rencontré par le ième électron
                    P[i][0] += P[i][2]*e
                    P[i][1] += P[i][3]*e
                    j = nObstacle + 1   # on sort de la boucle
                else:
                    (C,B) = Equation(P,i)
                    (M,V) = EquationObstacle(O,j)
                    if M-C != 0:
                        x = (B-V)/(M-C) #abscisse de l'éventuel point d'intersection de la trajectoire et de l'obstacle
                        I = P[i][0],P[i][1] #extremité début de trajectoire
                        mi = min(P[i][0],P[i][0]+e*P[i][2])
                        ma = max(P[i][0],P[i][0]+e*P[i][2]) # [mi,ma] décrit les valeurs du segment sur l'axe des abscisses
                        if  x > mi and x < ma and x < O[j][1] and x > O[j][0]: #condition de rencontre (intersection des segments). On suppose e suffisamment petit pour que le rebond se fasse après le déplacement et non à la collision.
                            y = C*x + B#ordonnée du point d'intersection
                            a2 = (O[j][0]-x)**2+(O[j][2]-y)**2
                            b2 = (I[0]-x)**2+(I[1]-y)**2
                            c2 = (I[1]-O[j][2])**2+(O[j][0]-I[0])**2
                        
                            
                            alpha = math.acos((a2+b2-c2)/(2*math.sqrt(a2*b2)))  #théorème d'Al-Kashi -> alpha est l'angle entre la trajectoire initiale et l'obstacle
                            d = abs(O[j][3]-y)
                            f = abs(O[j][1]-x)
                            beta = math.atan(d/f)
                            
                            v = e*math.sqrt((P[i][2])**2+(P[i][3])**2)
                            vy = math.sin(alpha+beta)*v
                            vx = math.cos(alpha+beta)*v
                        
                            #choix des signes : rustine au problème géomètrique
                            
                            choix = [[1,1],[1,-1],[-1,1],[-1,-1]]
                            boo = True
                            c = 0
                            while boo and c < 4:
                                vxn = choix[c][0]*vx
                                vyn = choix[c][1]*vy
                                
                                if vyn == -1*P[i][3] and vxn == -1*P[i][2]:   #cela correspondrait à un retour en arrière
                                    c+=1
                                else:
                                    P[i][0] = x + e*vxn #l'électron rebondit
                                    P[i][1] = y + e*vyn #""
                                    boo = Intersection(I,[P[i][0],P[i][1]],[O[j][0],O[j][2]],[O[j][1],O[j][3]])
                                    c+=1

                                    
                            P[i][2] = vxn #et a des vitesses conséquentes au rebond
                            P[i][3] = vyn
                                
                                
                            j = nObstacle + 1 #il y a eu rebond : l'électron s'est déplacé donc on peut sortir de la boucle
                        else :
                            j += 1  
                    else:
                        j += 1
    return P


def TrajectoireCryoNS(L, e, n, nObstacle, v0x, v0y):
    
    pas_moyens = L/(v0x**2 + v0y**2)    #donne un nombre moyen de points utilisés pour tracer la trajectoire de chaque électron : utile pour l'affichage
    indic_affichage = pas_moyens//50  #on veut à peu près 50 points par électron 
    
    
    P = ListeElectrons(L, n, v0x, v0y)
    O = ObstacleNegativeStaining(L, nObstacle)
    nObstacle = len(O)  # Attention : il y a potentiellement moins d'obstacles représentés par l'enveloppe des obstacles que d'obstacles eux-mêmes
    
    for i in range (n):
        def Equation(P,i):  #renvoie le coefficient directeur et l'ordonnée à l'origine de la droite qui prolonge le segment trajectoire du ième électron
            C = P[i][3]/P[i][2]#coefficient directeur de la trajectoire
            B = P[i][1] - P[i][0]*C #ordonnée à l'origine de la trajectoire
            return (C,B)
        def EquationObstacle(O,j):  #renvoie le coefficient directeur et l'ordonnée à l'origine de la droite qui prolonge le j-ème obstacle
            M = (O[j][2]-O[j][3])/(O[j][0]-O[j][1])#coefficient directeur de l'obstacle
            V = O[j][2]-M*O[j][0]#ordonnée à l'origine de l'obstacle
            return (M,V)
        
        t_affichage = -1
        while (P[i][0] >= 0) and (P[i][0] <= L) and (P[i][1] >= 0) and (P[i][1] <= L): #tant qu'on n'a pas rencontré de paroi
            t_affichage += 1   
            if t_affichage % indic_affichage == 0:  
                 X = [P[i][0]]
                 Y = [P[i][1]]
                 scatter(X, Y, c = Couleurs[i % 8])    #chaque électron a une couleur
                 axis([0,5/4*L,0,5/4*L])
                 show()
                 pause(0.1)
                 if t_affichage % 5 == 0:   #on efface la mémoire de la trajectoire tous les 5 points
                     clf()
                     Afficher_Obstacle(O)
            O = ObstacleSort(O, P, nObstacle, i)
            j = 0
            while j != nObstacle + 1 :  # j = nObstacle va produire l'arrêt de la boucle par la première condition if
                if j == nObstacle:  # si aucun obstacle n'est rencontré par le ième électron
                    P[i][0] += P[i][2]*e
                    P[i][1] += P[i][3]*e
                    j = nObstacle + 1   # on sort de la boucle
                else:
                    (C,B) = Equation(P,i)
                    (M,V) = EquationObstacle(O,j)
                    if M-C != 0:
                        x = (B-V)/(M-C) #abscisse de l'éventuel point d'intersection de la trajectoire et de l'obstacle
                        I = P[i][0],P[i][1] #extremité début de trajectoire
                        mi = min(P[i][0],P[i][0]+e*P[i][2])
                        ma = max(P[i][0],P[i][0]+e*P[i][2]) # [mi,ma] décrit les valeurs du segment sur l'axe des abscisses
                        if  x > mi and x < ma and x < O[j][1] and x > O[j][0]: #condition de rencontre (intersection des segments). On suppose e suffisamment petit pour que le rebond se fasse après le déplacement et non à la collision.
                            y = C*x + B#ordonnée du point d'intersection
                            a2 = (O[j][0]-x)**2+(O[j][2]-y)**2
                            b2 = (I[0]-x)**2+(I[1]-y)**2
                            c2 = (I[1]-O[j][2])**2+(O[j][0]-I[0])**2
                        
                            
                            alpha = math.acos((a2+b2-c2)/(2*math.sqrt(a2*b2)))  #théorème d'Al-Kashi -> alpha est l'angle entre la trajectoire initiale et l'obstacle
                            d = abs(O[j][3]-y)
                            f = abs(O[j][1]-x)
                            beta = math.atan(d/f)
                            
                            v = e*math.sqrt((P[i][2])**2+(P[i][3])**2)
                            vy = math.sin(alpha+beta)*v
                            vx = math.cos(alpha+beta)*v
                        
                            #choix des signes : rustine au problème géomètrique
                            
                            choix = [[1,1],[1,-1],[-1,1],[-1,-1]]
                            boo = True
                            c = 0
                            while boo and c < 4:
                                vxn = choix[c][0]*vx
                                vyn = choix[c][1]*vy
                                
                                if vyn == -1*P[i][3] and vxn == -1*P[i][2]:   #cela correspondrait à un retour en arrière
                                    c+=1
                                else:
                                    P[i][0] = x + e*vxn #l'électron rebondit
                                    P[i][1] = y + e*vyn #""
                                    boo = Intersection(I,[P[i][0],P[i][1]],[O[j][0],O[j][2]],[O[j][1],O[j][3]])
                                    c+=1

                                    
                            P[i][2] = vxn #et a des vitesses conséquentes au rebond
                            P[i][3] = vyn
                                
                                
                            j = nObstacle + 1 #il y a eu rebond : l'électron s'est déplacé donc on peut sortir de la boucle
                        else :
                            j += 1  
                    else:
                        j += 1
    return P