import numpy
import coefficients 
import csv
import SVD_LS
import matplotlib.pyplot as plt

plt.style.use('seaborn-ticks') 

x=[1,2,3,4,5]
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

R_1=numpy.matrix([
				results[0],
				results[1],
				results[8],
				results[9],
				results[8],
				results[9],
				results[4],
				results[5],
				results[4],
				results[5]
				])
R_2=numpy.matrix([
				results[0],
				results[1],
				results[8],
				results[9],
				results[2],
				results[3],
				results[4],
				results[5],
				results[4],
				results[5]
				])
R_3=numpy.matrix([
				results[0],
				results[1],
				results[8],
				results[9],
				results[2],
				results[3],
				results[4],
				results[5],
				results[6],
				results[7]
				])
R_4=numpy.matrix([
				results[6],
				results[7],
				results[8],
				results[9],
				results[2],
				results[3],
				results[4],
				results[5],
				results[6],
				results[7]
				])
R_5=numpy.matrix([
				results[8],
				results[9],
				results[0],
				results[1],
				results[4],
				results[5],
				results[4],
				results[5],
				results[6],
				results[3]
				])		


def tikhonov_regularization(A,b):
	E=numpy.eye(2)
	Ab=numpy.append(A,b[:,0],axis=1)
	U,S,V=numpy.linalg.svd(Ab)
	
	V_transpose=numpy.transpose(V)

	if len(Ab)==10:
		U_improve=numpy.delete(U,numpy.s_[3:10],axis=1) 
	elif len(Ab)==8:
		U_improve=numpy.delete(U,numpy.s_[3:8],axis=1)
	elif len(Ab)==6:
		U_improve=numpy.delete(U,numpy.s_[3:6],axis=1)
	elif len(Ab)==4:
		U_improve=numpy.delete(U,numpy.s_[3:4],axis=1)
	elif len(Ab)==2:
		U_improve=U
	
	if len(S)==3:
		S_improve=numpy.matrix([ [S[0], 0,0], [0, S[1],0],[0,0,S[2]] ])
	elif len(S)==2:
		S_improve=numpy.matrix([ [S[0], 0], [0, S[1]] ])
	
	if len(Ab)==2:
		V_i=numpy.delete(V,numpy.s_[2:3],axis=0) 
		V_improve=numpy.delete(V_i,numpy.s_[2:3],axis=1)
		SVD_Ab=U_improve*S_improve*V_improve # SVD матрица (A,b)=USV^T
		lambda_EL=(S_improve[1,1]**2)# Параметр регуляризации

		left_equation=numpy.transpose(A)*A+lambda_EL*E
		right_equation=numpy.transpose(A)*b
		x=numpy.linalg.solve(left_equation,right_equation)
	
	else:
		SVD_Ab=U_improve*S_improve*V # SVD матрица (A,b)=USV^T
		lambda_EL=(S_improve[2,2]**2) # Параметр регуляризации
		left_equation=numpy.transpose(A)*A+lambda_EL*E
		right_equation=numpy.transpose(A)*b

		x=numpy.linalg.solve(left_equation,right_equation)
	
	equal_first=-1/V_transpose[2,2]
	equal_second=numpy.matrix([[V_transpose[0,2],V_transpose[1,2]]])
	Q=numpy.dot(equal_first,equal_second)

	return x,Q


norm_solve_q=SVD_LS.SVD(K,R)
list_norm_solve_q=[	norm_solve_q[0,0],
					norm_solve_q[0,0],
					norm_solve_q[0,0],
					norm_solve_q[0,0],
					norm_solve_q[0,0]]
list_norm_solve_q2=[norm_solve_q[1,0],
					norm_solve_q[1,0],
					norm_solve_q[1,0],
					norm_solve_q[1,0],
					norm_solve_q[1,0]]

list_solve_activity_1=[	7.82754979,
						7.78213663,
						7.80635069,
						7.83326356,
						7.79648012]	
list_solve_activity_2=[	10.50451338,
						10.31787599,
						10.48437687,
						10.21737795,
						10.41787648]	


					

									
direct_method_Q1=[	7.64132205,
					7.86632811,
					7.65481583,
					7.77588686,
					8.09338928]
direct_method_Q2=[	10.55238219,
					9.66703493,
					10.38488406,
					10.21738593,
					10.6002388]



figure, axis = plt.subplots( )

plt.plot(x,list_norm_solve_q2,'r--',label='Нормальное решение')
plt.plot(x,list_solve_activity_2,color='k',linestyle='',marker='.',label='$\lambda_{EL}=9.092*10^{-9}$')
plt.plot(x,direct_method_Q2,color='k',linestyle='',marker='+',label='Обратная матрица')

plt.grid(which='major',
        color = 'gray',
        linestyle = ':')
plt.axis('equal')

plt.legend()

plt.xlabel('Номер наблюдения')
plt.ylabel('Активность $Q_2$')


if __name__=='__main__':
	print('LS')
	print(SVD_LS.SVD(K,R))
	print('TLS')
	print(tikhonov_regularization(K,R))
# 	print('Sum:',numpy.sum(tikhonov_regularization(K,R)))
	print(norm_solve_q[1,0])
	print	(
			'1.',tikhonov_regularization(K,R_1),
			'\n2.',tikhonov_regularization(K,R_2),
			'\n3.',tikhonov_regularization(K,R_3),
			'\n4.',tikhonov_regularization(K,R_4),
			'\n5.',tikhonov_regularization(K,R_5),
			)
# 	plt.show()
# 	figure.savefig('parametr.png',dpi=1000)

	

