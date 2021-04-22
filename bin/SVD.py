'''
Сингулярное разложение матрицы для решения
переопределенной СЛАУ:
KQ=R, m>n (число уравнений больше числа неизвестных)
'''
import numpy as np
import coefficients as cof

K=np.matrix([[0.278,0.082],[0,0.042],[0.278,0.082],[0,0.042]])
R=np.matrix([[3.022],[0.430],[2.998],[0.441]])

# Применение модуля numpy для нахлждения сингулярной матрицы
U,S,Vh=np.linalg.svd(K)
S_add=np.matrix([[S[0], 0], [0, S[1]]])

# Правка матрицы U
U_add=np.delete(U,np.s_[2:4],axis=1)

# Матрица SVD
K_expansion=U_add*S_add*Vh

# Приведение к псевдообратной матрице
Vh_transpose=np.transpose(Vh)	# Транспонирование транспонированной матрицы
S_inverse=np.linalg.inv(S_add)
U_transpose=np.transpose(U_add)

# Псевдообратная матрица
K_pseudomatrix=Vh_transpose*S_inverse*U_transpose

# Решение уравнения
Q=K_pseudomatrix*R

print(Q)