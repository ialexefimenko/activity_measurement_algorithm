import scipy.linalg
from numpy import arange, dot, concatenate,linalg,random,linspace,array, concatenate,zeros
import matplotlib.pyplot as plt 
import coefficients 
import activity_source 
from scipy.stats import linregress


# Скалярное произведение матриц K (коэффициентов чувствительности 
# по полученным образцовым мерам) и Q (истинная активность (deltaQ=6%))
scalar_AQ_varitable=dot(coefficients.K,activity_source.Q)
# Скалярное произведение матриц K (коэффициентов чувствительности 
# по полученным образцовым мерам) и Q (активность, полученная путем решения СЛАУ)
scalar_AQ_experiment=dot(coefficients.K,coefficients.matrix_Q_11_21)

# Единичная матрица
E=array([[1,0],[0,1]])


# L = scipy.linalg.toeplitz(coefficients.K)

# Tikhonov regularization

L_0=array([[0],[0]])
separator=[0.5,0.4,0.35,0.3]
R_squared=[]

C=concatenate([coefficients.matrix_R_11_21,L_0])

figure, axis = plt.subplots()
plt.style.use('ggplot')
plt.rcParams['figure.figsize']=(7,5)

for parametr_regularization in separator:
	B=concatenate([coefficients.K,parametr_regularization*E])
	R_lstsq=linalg.lstsq(B,C, rcond=None)
	R_1=R_lstsq[0][0]
	R_2=R_lstsq[3][0]
	R=array([[R_1],[R_2]])
	axis.plot(R,label="$\lambda=$" + str(parametr_regularization),linewidth=0.5)
	

plt.title('Selection regularization parameter $\lambda$')
plt.ylabel('value $Q$')
axis.plot(scalar_AQ_varitable,color='black',linewidth=2,label='$Q$')
axis.plot(scalar_AQ_experiment,'b--',label='$\delta{Q}$',linewidth=1)

axis.legend()


if __name__=='__main__':
	
# 	plt.show()
 	print()
# 	figure.savefig('parametr.png',dpi=1000)
