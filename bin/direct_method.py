import numpy 
import coefficients 
import csv

results=[]
with open("count_compendium.csv") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) 
    for row in reader: 
        results.append(row)
    
K=numpy.matrix([[coefficients.k11,coefficients.k12],
				[coefficients.k21,coefficients.k22]])

R=numpy.matrix([[3.022],
				[0.430]])
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

vector_Q=numpy.linalg.solve(K,R)
		
if __name__=='__main__':
	print('Direct:')
	print(vector_Q)
	print('Sum:',numpy.sum(vector_Q))
	print(	'1.',numpy.linalg.solve(K,R_1),
			'2.', numpy.linalg.solve(K,R_2),
			'3.', numpy.linalg.solve(K,R_3),
			'4.', numpy.linalg.solve(K,R_4),
			'5.', numpy.linalg.solve(K,R_5)		)
