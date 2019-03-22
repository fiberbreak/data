#/usr/local/bin/python3
# -*- coding: utf-8 -*-

import os.path
import subprocess
import time
import string
import math
import numpy as np
import warnings
import sys


def normalise(v, tolerance=1e-18):
    ''' Normalise un quaternion'''
    mag2 = np.sum(np.square(v))
    if abs(mag2 - 1.0) > tolerance:
        mag = np.sqrt(mag2)
        v/=(mag+tolerance)
    return v

def axeangle_to_q(v, theta):
    ''' Renvoie le quaternion associé à une rotation d'un angle theta et autour de l'axe v '''
    v = normalise(v)
    x, y, z = v
    theta /= 2.
    w = np.cos(theta)
    x = x * np.sin(theta)
    y = y * np.sin(theta)
    z = z * np.sin(theta)
    return w, x, y, z



# fichier de depart
file = open("16-Mar-2017_export_4_classification_moy_exp.txt","r")
newfile = open("16-Mar-2017_export_4_classification_moy_exp_quaternion_corriges_canonique_shortestpath.txt","w")


# format fichier initial
# colonne 1 : id de la courbe
# colonne 2 : timestamps
# colonne 3 : S1 (X)
# colonne 4 : S2 (Y)
# colonne 5 : S3 (Z)


k=0
id_old=float(-2)
id=float(-1)
t=float(0)
s1=float(0)
s2=float(0)
s3=float(0)
ps=float(0.000000000)
xinit=float(0)
yinit=float(0)
zinit=float(1)
debutserie=0
notsp=0;

for line in file:
    k=k+1

    temp=line.rstrip('\n\r')


    if (k==1):
        newfile.write("id" + '\t' + "t" + '\t' +  "a" + '\t' + "X" + '\t' + "Y" + '\t' + "Z" + '\n')

    # id, t, S1, S2, S3
    donnees = temp.rstrip('\n\r').split("\t")
    #print k


    if (k>1):
        id=float(donnees[0]);
        t=float(donnees[1]);
        s1=float(donnees[2]);
        s2=float(donnees[3]);
        s3=float(donnees[4]);

    #si c'est une nouvelle time series pour le premier point on ne fait rien
    #sinon
    if ((id!=id_old) and (k>1)):
        #print id
        sys.stdout.write (str('\n') + str(id) )
        # on memorise les coordonnees initiale du point
        xinit=s1
        yinit=s2
        zinit=s3

        #### on calcul le quaternion permettant de se ramener au pole nord
        x1=float(0)
        y1=float(0)
        z1=float(1)
        x2=s1
        y2=s2
        z2=s3
        # on  calcul le prodcuit scalaire
        a=(x1*x2)+(y1*y2)+(z1*z2)
		# on se fait attention pour le arccos après
        if a>1.0:
			a=1.0
			sys.stdout.write (str(' a. '))
		
			
        # on calcul le produit vectoriel
        b=(y1*z2)-(z1*y2)
        c=(z1*x2)-(x1*z2)
        d=(x1*y2)-(y1*x2)

        # correction SR
        (a,b,c,d) = axeangle_to_q(-np.array((b,c,d)),np.arccos(a))
		
		# indicateur debut serie remis a zero
        debutserie=0
        
        if 0:
            t2 =   a*b
            t3 =   a*c
            t4 =   a*d
            t5 =  -b*b
            t6 =   b*c
            t7 =   b*d
            t8 =  -c*c
            t9 =   c*d
            t10 = -d*d
            v1new = 2*( (t8 + t10)*s1 + (t6 -  t4)*s2 + (t3 + t7)*s3 ) + s1
            v2new = 2*( (t4 +  t6)*s1 + (t5 + t10)*s2 + (t9 - t2)*s3 ) + s2
            v3new = 2*( (t7 -  t3)*s1 + (t2 +  t9)*s2 + (t5 + t8)*s3 ) + s3
            print(v1new,v2new,v3new)

    if ((id==id_old) and (k>1)):

        debutserie=debutserie+1
        # rotation pour avoir la reference au pole
        # voir wikipedia pour cette optimsation du code
        t2 =   a*b
        t3 =   a*c
        t4 =   a*d
        t5 =  -b*b
        t6 =   b*c
        t7 =   b*d
        t8 =  -c*c
        t9 =   c*d
        t10 = -d*d
        v1new = 2*( (t8 + t10)*s1 + (t6 -  t4)*s2 + (t3 + t7)*s3 ) + s1
        v2new = 2*( (t4 +  t6)*s1 + (t5 + t10)*s2 + (t9 - t2)*s3 ) + s2
        v3new = 2*( (t7 -  t3)*s1 + (t2 +  t9)*s2 + (t5 + t8)*s3 ) + s3
        s1=v1new
        s2=v2new
        s3=v3new


        # calcul du quaternion via le point en cours et le precedent
        # on passe en vecteurs unitaires au cas ou
        n=float(math.sqrt((s1*s1)+(s2*s2)+(s3*s3)))
        #print n
        s1=(s1/n)
        s2=(s2/n)
        s3=(s3/n)

        x1=s1
        y1=s2
        z1=s3
        x2=s1_old
        y2=s2_old
        z2=s3_old
        # on  calcul le produit scalaire
        ps=(x1*x2)+(y1*y2)+(z1*z2)
        # on fait attention pour le arccos après
        if ps>1.0:
			sys.stdout.write (str(' ps. '))
			#sys.stdout.write (str(' ps. ') + str(t) + str('t= ') + str('ps=%0.24f'%ps) + str('\n') )
			ps=1.0
			
			
        # on calcul le produit vectoriel
        a1=(y1*z2)-(z1*y2)
        a2=(z1*x2)-(x1*z2)
        a3=(x1*y2)-(y1*x2)

        # correction SR
        (ps,a1,a2,a3) = axeangle_to_q(-np.array((a1,a2,a3)),np.arccos(ps))
		
        # test canonique
		# autre possibilite signe(max(|ps|,|a1|,|a2|,|a3|))>0
		# q = np.asarray((ps,a1,a2,a3))
		# i = np.argmax(np.abs(q))
		# if q|i]<0:
		
        if (ps<0):
                ps=-ps
                a1=-a1
                a2=-a2
                a3=-a3
				#sys.stdout.write (str('notcano ps. ') + str(t) + str('t= ') + str('ps=%0.24f'%ps) + str('\n') )
                #sys.stdout.write (str(' cano=') + str(debutserie) + str(' '))
		
        # test du shortest path
        if debutserie>1:
                q1=(ps,a1,a2,a3)
                q2=quatold
                d=np.dot(q1,q2)  # = cos(angle)
                if (d<0):
					notsp=notsp+1
					sys.stdout.write (str(' sp=') + str(debutserie) + str(' ') + str('notsp') + str(notsp))		
					ps=-ps
					a1=-a1
					a2=-a2
					a3=-a3
			

        # ecriture id
        if (id<10): newfile.write('0000' + str(int(id)) + '\t')
        if ((id<100) and (id>=10)): newfile.write('000' + str(int(id)) + '\t')
        if ((id<1000) and (id>=100)): newfile.write('00' + str(int(id)) + '\t')
        if ((id<10000) and (id>=1000)): newfile.write('0' + str(int(id)) + '\t')
        if ((id<100000) and (id>=10000)): newfile.write(  str(int(id)) + '\t')

        # ecriture t
        newfile.write(str(t) + '\t')
		
        # ecriture quaternion
        newfile.write(str(ps) + '\t')
        newfile.write(str(a1) + '\t')
        newfile.write(str(a2) + '\t')
        newfile.write(str(a3) + '\n')
		
        # on stoque quaternion actuel
        quatold=(ps,a1,a2,a3)

    #fin du if


    #on memorise
    id_old=id
    t_old=t
    s1_old=s1
    s2_old=s2
    s3_old=s3
