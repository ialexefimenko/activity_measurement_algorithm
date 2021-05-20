
import numpy 
import coefficients 
import csv

results=[]
with open("count_compendium.csv") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) 
    for row in reader: 
        results.append(row)
    
K=numpy.matrix([[coefficients.k11,coefficients.k12],
				[coefficients.k21,coefficients.k22],
				[coefficients.k11,coefficients.k12],
				[coefficients.k21,coefficients.k22],
				[coefficients.k11,coefficients.k12],
				[coefficients.k21,coefficients.k22],
				[coefficients.k11,coefficients.k12],
				[coefficients.k21,coefficients.k22],
				[coefficients.k11,coefficients.k12],
				[coefficients.k21,coefficients.k22]])

R=numpy.matrix([results[0],
				results[1],
				results[2],
				results[3],
				results[4],
				results[5],
				results[6],
				results[7],
				results[8],
				results[9]])
	
square_K = numpy.matrix([[coefficients.k11,coefficients.k12],
				[coefficients.k21,coefficients.k22]
				])
R_1=numpy.matrix([results[0],
				results[1]])
R_2=numpy.matrix([results[2],
				results[3]])
R_3=numpy.matrix([results[4],
				results[5]])
R_4=numpy.matrix([results[6],
				results[7]])
R_5=numpy.matrix([results[8],
				results[9]])
		
		
		
def SVD(A,b):
	#Применение модуля numpy для нахлждения сингулярной матрицы
	U,S,Vh=numpy.linalg.svd(A)
	S_add=numpy.matrix([[S[0], 0], [0, S[1]]])
	if len(K)==10:
		U_add=numpy.delete(U,numpy.s_[2:10],axis=1)
	elif len(K)==6:
		U_add=numpy.delete(U,numpy.s_[2:6],axis=1)
	elif len(K)==8:
		U_add=numpy.delete(U,numpy.s_[2:8],axis=1)
	elif len(K)==4:
		U_add=numpy.delete(U,numpy.s_[2:4],axis=1)
	elif len(K)==2:
		U_add=U
	K_expansion=U_add*S_add*Vh
	
	# Приведение к псевдообратной матрице
	Vh_transpose=numpy.transpose(Vh)	# Транспонирование транспонированной матрицы
	S_inverse=numpy.linalg.inv(S_add)
	U_transpose=numpy.transpose(U_add)

	# Псевдообратная матрица
	pseudomatrix_A=Vh_transpose*S_inverse*U_transpose
	
	# Решение уравнения SVD
	x=pseudomatrix_A*b

	return x
	
def least_squares(A,b):
	'''
	Функция возвращает переменные x уравнения Ax=b 
	методом наименьших квадратов
	'''
	A_transpose=numpy.transpose(A)
	A_transpose_A=numpy.dot(A_transpose, A)
	A_transpose_b=numpy.dot(A_transpose, b)
	x=numpy.linalg.solve(A_transpose_A,A_transpose_b)
    
	return x

if __name__=='__main__':
	print('SVD')
	print(SVD(K,R))
	print('Sum:',numpy.sum(SVD(K,R)))
	print('LS')
	print(least_squares(K,R))
	print('Sum:',numpy.sum(least_squares(K,R)))

	
	
