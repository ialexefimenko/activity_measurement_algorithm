from scipy.linalg import toeplitz
import coefficients 
import numpy 

f=numpy.matrix([[coefficients.k11,coefficients.k12],[coefficients.k21,coefficients.k22],
				[coefficients.k11,coefficients.k12],[coefficients.k21,coefficients.k22]])



toeplitz=toeplitz(f)
toeplitz_transpose=numpy.transpose(toeplitz)

tt=numpy.dot(toeplitz_transpose,toeplitz)

tt_inverse=numpy.linalg.inv(tt)




ff=numpy.transpose(f)
j=numpy.dot(ff,f)
j_inverse=numpy.linalg.inv(j)


E=numpy.array([[1,0],[0,1]])
E_T=numpy.transpose(E)
EE=numpy.dot(E,E_T)



print(EE)