import numpy
import activity_source
import math

list_N_beta_Sr90=numpy.array([3.51,3.47,3.52,3.59,3.67])
list_N_beta_Am241=numpy.array([0.787,0.814,0.782,0.839,0.788])
list_background_beta=numpy.array([0.806,0.807,0.798,0.800,0.796])
list_N_alpha_Am241=numpy.array([0.374,0.440,0.400,0.421,0.403])

counts_sr90=list_N_beta_Sr90-list_background_beta
counts_am2411=list_N_beta_Am241
counts_am2411_alpha=list_N_alpha_Am241

coefficients_11=counts_sr90/activity_source.A_Srt
coefficients_21=numpy.array([0,0,0,0,0])

coefficients_12=counts_am2411/activity_source.A_Amt
coefficients_21=counts_am2411/activity_source.A_Amt
coefficients_22=counts_am2411_alpha/activity_source.A_Amt

norm_Kw=0.294
norm_sigma_k=9.1*10**(-3)						



def cov(x, y):
	with numpy.errstate(divide='ignore'):
		numpy.float64(1.0) / 0.0
		xbar, ybar = x.mean(), y.mean()
		h=numpy.sum((x - xbar)*(y - ybar))/(len(x) - 1)
		return h

# Дисперсионная матрица вектора Kw1
Vk11=numpy.array([[	cov(coefficients_11,coefficients_11),cov(coefficients_11,coefficients_12)],
					[cov(coefficients_12,coefficients_11),cov(coefficients_12,coefficients_12)]])

# Дисперсионная матрица вектора Kw2
Vk22=numpy.array([[	cov(coefficients_21,coefficients_21),cov(coefficients_21,coefficients_22)],
					[cov(coefficients_22,coefficients_21),cov(coefficients_22,coefficients_22)]])

# Кросс-ковариационная матрица векторов Kw1 and Kw2
Vk12=numpy.array([[	cov(coefficients_11,coefficients_21),cov(coefficients_11,coefficients_22)],
					[cov(coefficients_12,coefficients_21),cov(coefficients_12,coefficients_22)]])

# Кросс-ковариационная матрица векторов Kw2 and Kw1
Vk21=numpy.array([[	cov(coefficients_21,coefficients_11),cov(coefficients_21,coefficients_12)],
					[cov(coefficients_22,coefficients_11),cov(coefficients_22,coefficients_12)]])

# Совместная ковариационная матрица
Vk=numpy.array([[Vk11,Vk12],[numpy.transpose(Vk21),Vk22]])

# Норма Совместная ковариационная матрица
norm_Vk=numpy.linalg.norm(Vk)

norm_Vk11=numpy.linalg.norm(Vk11)
norm_Vk22=numpy.linalg.norm(Vk22)

cof_var=numpy.sqrt(norm_Vk)/norm_Kw
cof_var11=numpy.sqrt(norm_Vk11)/norm_Kw
cof_var22=numpy.sqrt(norm_Vk22)/norm_Kw



delta_p=-((numpy.sqrt(norm_Vk)-norm_sigma_k)/norm_Kw)
delta_p11=-((numpy.sqrt(norm_Vk11)-norm_sigma_k)/norm_Kw)
delta_p22=-((numpy.sqrt(norm_Vk22)-norm_sigma_k)/norm_Kw)



if __name__=='__main__':

	print(delta_p,delta_p11,delta_p22)

	



