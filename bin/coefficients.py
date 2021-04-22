# Алгоритм определения активности [Q_1] and [Q_2] путем прямого решения СЛАУ 
# (должно давать состоятельные оценки активности (!))
# для смеси Sr-90+Y-90 + Am-241 по равноточным кожффициентам преобразования

import math
import matplotlib.pyplot as plt
from numpy import mean 
import numpy as np
import activity_source as ac

# Идентификация переменных (скорость счета N) по каналу beta
list_N_beta_Sr90=np.array([3.51,3.47,3.52,3.59,3.67])
list_N_beta_Am241=np.array([0.787,0.814,0.782,0.839,0.788])
list_background_beta=np.array([0.806,0.807,0.798,0.800,0.796])
# N' -- скорость счета средняя
av_background_beta=mean(list_background_beta) # N_b
av_list_N_beta_Sr90=mean(list_N_beta_Sr90) # N' (beta) Sr90
av_list_N_beta_Am241=mean(list_N_beta_Am241) # N' (beta) Am241
# Разница N-N_b
av_list_N_beta_Am241_mback=av_list_N_beta_Am241-av_background_beta
av_list_N_beta_Sr90_mback=av_list_N_beta_Sr90-av_background_beta
# Идентификация переменных (скорость счета N) по каналу alpha
list_N_alpha_Sr90=0 # Скорость счета для Sr-90 при стат.погрешности >100% = 0 (1/(Бк*с))
list_N_alpha_Am241=np.array([0.374,0.440,0.400,0.421,0.403]) #скорость счета для Am-241 
# при стат.погрешности 10%
list_background_alpha=np.array([0])
# N' -- скорость счета средняя
av_list_N_alpha_Am241=mean(list_N_alpha_Am241) # N' (alpha) Am241
# Определение коэффициентов чувствительности (получены с равноточными набдлюдениями)
# Равноточные измерения -- измерения, выполненые одним и тем же прибором, одним и тем 
# же методом, одинаковым числом приемов и в одинаковых условиях
k11=av_list_N_beta_Sr90_mback/ac.A_Srt
k21=0
k12=av_background_beta/ac.A_Amt
k22=av_list_N_alpha_Am241/ac.A_Amt
# Идентификация матрицы коэффициентов преобразования [K_w]
K=np.array([[k11,k12],[k21,k22]])

# Идентификация переменных для смеси (Am-241 + Sr-90)
# Идентификация массива скорости счета от смеси (Am-241 + Sr-90)
list_beta_Sr_Am=np.array([3.8,3.79,3.79,3.81,3.93])
mean_beta_Sr_Am=mean(list_beta_Sr_Am)
list_alpha_Sr_Am=np.array([0.441,0.404,0.434,0.427,0.443])
mean_alpha_Sr_Am=mean(list_alpha_Sr_Am)


#Решение системы уравнений методом Крамера
determinant_K=k11*k22-k21*k12 						#определитель матрица [К_w]

R_1=mean_beta_Sr_Am-av_background_beta
R_2=mean_alpha_Sr_Am

matrix_R_11_21=np.array([[mean_beta_Sr_Am-av_background_beta],[mean_alpha_Sr_Am]])
matrix_Q_11_21=np.linalg.solve(K,matrix_R_11_21)

############################
# Идентификация переменных для нахождения относительной погрешности решения СЛАУ 
delta_k_a11=((list_beta_Sr_Am[0]-mean_beta_Sr_Am)/ac.A_Sr_compendum_t+
(list_beta_Sr_Am[1]-mean_beta_Sr_Am)/ac.A_Sr_compendum_t+(list_beta_Sr_Am[2]-
mean_beta_Sr_Am)/ac.A_Sr_compendum_t+(list_beta_Sr_Am[3]-mean_beta_Sr_Am)/ac.A_Sr_compendum_t+
(list_beta_Sr_Am[4]-mean_beta_Sr_Am)/ac.A_Sr_compendum_t)

delta_k_a21=0

delta_k_a12=((list_beta_Sr_Am[0]-mean_beta_Sr_Am)/ac.A_Am_compendum_t+
(list_beta_Sr_Am[1]-mean_beta_Sr_Am)/ac.A_Am_compendum_t+(list_beta_Sr_Am[2]-
mean_beta_Sr_Am)/ac.A_Am_compendum_t+(list_beta_Sr_Am[3]-mean_beta_Sr_Am)/ac.A_Am_compendum_t+
(list_beta_Sr_Am[4]-mean_beta_Sr_Am)/ac.A_Am_compendum_t)

delta_k_a22=((list_alpha_Sr_Am[0]-mean_alpha_Sr_Am)/ac.A_Am_compendum_t+
(list_alpha_Sr_Am[1]-mean_alpha_Sr_Am)/ac.A_Am_compendum_t+(list_alpha_Sr_Am[2]-
mean_alpha_Sr_Am)/ac.A_Am_compendum_t+(list_alpha_Sr_Am[3]-mean_alpha_Sr_Am)/ac.A_Am_compendum_t+
(list_alpha_Sr_Am[4]-mean_alpha_Sr_Am)/ac.A_Am_compendum_t)

mat_delta_aj_K=np.array([[delta_k_a11,delta_k_a12],[delta_k_a21,delta_k_a22]])

# Применения модуля numpy
inv_K=np.linalg.inv(K)
Q=np.dot(inv_K,matrix_R_11_21)

# Проверка
KQ=np.dot(K,matrix_Q_11_21)
result_res=KQ-matrix_R_11_21

# Невязка 
residual=np.linalg.norm(result_res)

if __name__=='__main__':
	print('Активность с linalg.solve:\n',matrix_Q_11_21)
	print('Равные активности, решая различным способом?',matrix_Q_11_21==Q)
	print('Матрица коэффициентов чувствительности:\n',K*1000)
	print('Определитель матрицы преобразования K:\n',determinant_K)
	






