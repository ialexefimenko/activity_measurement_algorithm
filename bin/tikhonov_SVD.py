'''
Программа находит приближенное решение СЛАУ
применяя метод Тихонова для общей формулировки 
наименьших квадратов
'''

import numpy
import coefficients 



K=numpy.matrix([[coefficients.k11,coefficients.k12],[coefficients.k21,coefficients.k22],
				[coefficients.k11,coefficients.k12],[coefficients.k21,coefficients.k22]])
K_square=numpy.matrix([[coefficients.k11,coefficients.k12],[coefficients.k21,coefficients.k22]])
R=numpy.matrix([[3.022],[0.430],[2.998],[0.441]])
R_square=numpy.matrix([[2.988],[0.430]])


Q_0=numpy.linalg.solve(K_square,R_square)

KR=numpy.append(K,R[:,0],axis=1)
U,S,V=numpy.linalg.svd(KR)
U_improve=numpy.delete(U,numpy.s_[3:4],axis=1) 
S_improve=numpy.matrix([[S[0], 0,0], [0, S[1],0],[0,0,S[2]]])

# Матрица SVD
KR_expansion=U_improve*S_improve*V

E=numpy.eye(2)
L=numpy.matrix([[1,-1],[0,1]])

sigma_3=lambda_E=(S_improve[2,2]**2)

VT=numpy.transpose(V)

# Стандартная форма
Q_1=-1/VT[2,2]
Q_2=numpy.matrix([[VT[0,2],VT[1,2]]])
Q=numpy.dot(Q_1,Q_2)

LT=numpy.transpose(L)
ll=LT*L

q=numpy.dot(numpy.transpose(K),K)
sigma_E=lambda_E*E
w=q+sigma_E
er=numpy.linalg.inv(w)
f=numpy.dot(numpy.transpose(K),R)
qq=er*f



mu_1=numpy.dot(Q,numpy.transpose(K_square))
r=R_square-K_square*numpy.transpose(Q)
mu_2=numpy.dot(mu_1,r)
mu_3=Q*ll*numpy.transpose(Q)
mu=mu_2/mu_3

lambda_L=mu*(1+numpy.linalg.norm(numpy.transpose(Q)))

# Общая форма с L
ww=q+sigma_E*lambda_L[0,0]*ll
inv_w=numpy.linalg.inv(ww)
qqq=inv_w*f


if __name__=='__main__':
	print(Q_0,'номинальная')
	print(qqq,'общая форма')
	print(Q,'стандартная форма')
	print()