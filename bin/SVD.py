'''
Сингулярное разложение матрицы для решения
переопределенной СЛАУ:
KQ=R, m>n (число уравнений больше числа неизвестных)
'''
import numpy 
import coefficients 
from pprint import pprint 

K=numpy.matrix([[coefficients.k11,coefficients.k12],[coefficients.k21,coefficients.k22],
				[coefficients.k11,coefficients.k12],[coefficients.k21,coefficients.k22]])

R=numpy.matrix([[3.022],[0.430],[2.998],[0.441]])

# Применение модуля numpy для нахлждения сингулярной матрицы
U,S,Vh=numpy.linalg.svd(K)
S_add=numpy.matrix([[S[0], 0], [0, S[1]]])

# Правка матрицы U
U_add=numpy.delete(U,numpy.s_[2:4],axis=1)

# Матрица SVD
K_expansion=U_add*S_add*Vh

# Приведение к псевдообратной матрице
Vh_transpose=numpy.transpose(Vh)	# Транспонирование транспонированной матрицы
S_inverse=numpy.linalg.inv(S_add)
U_transpose=numpy.transpose(U_add)

# Псевдообратная матрица
K_pseudomatrix=Vh_transpose*S_inverse*U_transpose

# Решение уравнения
Q=K_pseudomatrix*R

if __name__=='__main__':
	pprint(Q)