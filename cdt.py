from __future__ import division
import numpy as np
nb1 = raw_input("Uniform? y or n ")
if nb1=="y":
    a = True
else:
    a = False

nb = raw_input('The k-mesh:')
N = int(nb)

if a == False:
    P=np.array([[0.5 ,1 ,0],#   40  ! W
    [0.5 ,0.5, 0.5], #  40    ! L
    [0, 0, 0], #    40  ! G
    [0, 1, 0], #    40  ! X
    [0.5, 1, 0], #   40  ! W
    [0.75, 0.75, 0]]) # 1 ! K
    n = 10
    nk = P.shape[0]
    print P.shape[0]
    for i in range(0,nk-1):
        p_temp1 = np.linspace(P[i][0],P[i+1][0],num=n,endpoint=False)
        p_temp2 = np.linspace(P[i][1],P[i+1][1],num=n,endpoint=False)
        p_temp3 = np.linspace(P[i][2],P[i+1][2],num=n,endpoint=False)
        for j in range(0,n):
                print p_temp1[j],p_temp2[j],p_temp3[j],0
    print P[nk-1][0], P[nk-1][1], P[nk-1][2],0
else:
    number = str(N)
    f = open(number+".dat", "wb")
    f.write(str(N**3)+'\n')
    for i in range(N):
	for j in range(N):
	    for k in range(N):
		    f.write( ('{0:12.8f} {1:12.8f} {2:12.8f} {3:12.6e}'.format(i/N, j/N, k/N, 1/(N**3)))+'\n')
#a=np.linspace(0.5,0.5,n)
#b=np.linspace(1,0.5,n)
#c=np.linspace
#print
