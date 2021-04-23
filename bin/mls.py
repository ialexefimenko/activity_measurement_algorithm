'''
Метод наименьших квадратов.
Исходные данные идентифицировать как список списков
'''

import functions_LA as function 
import coefficients 
from pprint import pprint 


def least_squares(A,b, tol=2):
	'''
	Функция возвращает переменные x уравнения Ax=b 
	методом наименьших квадратов
	'''
	A_transpose=function.transpose(A)
	A_transpose_A=function.matrix_multiply(A_transpose, A)
	A_transpose_b=function.matrix_multiply(A_transpose, b)
	x=function.solve_equations(A_transpose_A,A_transpose_b,tol=tol)
    
	return x

# Массив значений переменных матрицы коэффициентов преобразования
K=[[coefficients.k11,coefficients.k12],[coefficients.k21,coefficients.k22],
				[coefficients.k11,coefficients.k12],[coefficients.k21,coefficients.k22]]
# Массив значений переменных скорости счета детектора
R=[[3.022],[0.430],[2.998],[0.441]]
# Решение уравнения LS
Q=least_squares(K,R)



if __name__=='__main__':
	pprint(Q)


