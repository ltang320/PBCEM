#!/usr/share/anaconda2/bin/python
from __future__ import division
from numpy import matrix
import numpy as np
import math as ma
from itertools import combinations,permutations, product
import copy
import pp
import time

 #A Matrix
A = []
A_tmp = list(product([0, 1,-1], repeat=4))
for i in range(len(A_tmp)):
		tmp_ind = np.where(np.array(A_tmp[i])!=0)[0]
		if len(tmp_ind)>0:
				if A_tmp[i][tmp_ind[0]]>0:
						A.append(A_tmp[i])
## print len(A)

## Data obtained from hematite02_RT_MEM75%
hema = input("Which temperature I want to calculate? hema_ _?")
mem = input("the mem scale is ? ")
E = list(np.loadtxt(('./Data_hema02_09/hema0%d_mem%d_energy.txt')%(hema,mem))) #Peak Position
I = list(np.loadtxt(('./Data_hema02_09/hema0%d_mem%d_intensity.txt')%(hema,mem))) #Intensity
number = input("How many numbers you want to calculate? ")
Energy = E[:number]
Intensity = I[:number]
##Nomalize all the Intensity Line
## 1st Step: Add up all the intensity line
Intensity_total = 0
for n in xrange(len(Intensity)):
	Intensity_total += Intensity[n]
## 2nd Step: Every Intensity Line occupied the percent of total
for x in xrange(len(Intensity)):
	Intensity[x] = Intensity[x]/Intensity_total

Energy_solution = []
Intensity_solution = []

for l in xrange(len(Energy)):
## Fine the Energy position in the possible cursor
		if (1.5E-07<=Energy[l]<3.4E-07):
				Energy_solution.append(Energy[l])
				Intensity_solution.append(Intensity[l])
## Caculate the possible comination of bases
Ener = list(permutations(Energy_solution,4))
Inten = list(permutations(Intensity_solution,4))

Ener_new = []
Inten_new = []
## Reduce the caculation number using some conditions
for i in xrange(len(Ener)):
	if (Ener[i][0]>Ener[i][1] and Ener[i][0]>Ener[i][2] and Ener[i][0]>Ener[i][3]):
		if (Ener[i][1]<Ener[i][2]<Ener[i][3]) or (Ener[i][1]>Ener[i][2]>Ener[i][3]):
			if 1.5<((Ener[i][1]+Ener[i][2]+Ener[i][3])/Ener[i][0])<=1.8:
				if 0.85<(Ener[i][1]+Ener[i][3])/2/Ener[i][2] < 1.15:
					Ener_new.append(Ener[i])
					Inten_new.append(Inten[i])
del Ener
del Inten
## Give file name

f = open('findsolution_hema0%d_mem%d_thres_%d_0128.txt'%(hema,mem,number),'w')
f1 = open('process_hema0%d_mem%d_thres_%d_0128.txt'%(hema,mem,number),'w')

### A function of finding the solution of all those lines who is the real base
def calculation_tang(A,Energy,Intensity,Ener_new,Inten_new,i):
	res = [[0 for j in xrange(len(Energy))] for k in xrange(len(A))]
	res_abs = [[0 for j in xrange(len(Energy))] for k in xrange(len(A))]
	mean_res = []
	res_min = []
	Energy_combi = []
	Inten_combi = []
## Step 1: Give the caculation results using 4 bases combination (minus or add)
	Energy_reference = matrix(A) * matrix(Ener_new[i]).T

## Step 2: Obtain the residue using different 4 bases
	for x in xrange(len(A)):
		for j in xrange(len(Energy)):
			res[x][j] = (abs(Energy_reference[x]) - Energy[j])/Energy[j]
			res_abs[x][j] = abs(res[x][j])
##	Energy_reference = np.array(Energy_reference)
##	print Energy_reference.shape  (40,1)

## Step 3: Find the best 40(the number decided by the equation number) minimum residue and the corresponding Energy and Intensity
	index_A = []
	res_abs = np.array(res_abs)
	while len(Energy_combi)<len(A):
##		res_abs.shape = len(A),len(Energy)
##		print res_abs.shape (40,122,1,1)
		res_min_id = np.where(res_abs == np.min(res_abs))

		for x in xrange(len(res_min_id[0])):
			Energy_combi.append(Energy[res_min_id[1][x]])
			Inten_combi.append(Intensity[res_min_id[1][x]])
			res_min.append(res[res_min_id[0][x]][res_min_id[1][x]])
			index_A.append(res_min_id[0][x])
##			Modified here res_min_id[0][x]
##			shape[1] the length of hang shape[2] the length of lie id[0] stands for x id[1] stands for y
##			print res_abs.shape[0]
##			print res_abs.shape[1]
			
			for m in xrange(res_abs.shape[0]):
				res[m][res_min_id[1][x]] = 10
				res_abs[m][res_min_id[1][x]] = 10
			for m2 in xrange(res_abs.shape[1]):
				res[res_min_id[0][x]][m2] = 10
				res_abs[res_min_id[0][x]][m2] = 10
## Step 4: Caculate the mean and standard diviation of the residue between 4 bases in 40 equations and the actual peak position
	mean_res = np.mean(res_min)
	std_res = np.std(res_min)
##	print len(index_A)
	return mean_res,std_res,Energy_combi,Inten_combi,index_A

start = time.time()
job_server = pp.Server()
mean_res_1 = []
mean_res_abs = []
std_res_1 = []
num_res = []

Inten_nomal = []
Itotal_nomalize = []
Inten_score = []

Solution_id = []
Solution_mean = []
Solution_std = []
Solution_Inten_base = []
Solution_Inten_equ = []
Solution_score = []
Solution_Energy = []
Solution_bi = []

Solution_mean_nomalize = 0
Solution_std_nomalize = 0
Solution_Inten_base_nomalize = 0
Solution_Inten_equ_nomalize = 0

Energy_combi_s = [[0 for j in xrange(len(A))] for k in xrange(len(Ener_new))]
Intensity_combi_single = [[0 for j in xrange(len(A))] for k in xrange(len(Ener_new))]
Solution_A = [[0 for j in xrange(len(A))] for k in xrange(len(Ener_new))]
Inten_combi_nomal = [0 for j in xrange(len(A))]
Inten_combi_s = [0 for k in xrange(len(Ener_new))]

f.write("The Intensity Max group is: \nEa Eb Ec Ed mean abs_mean std I_4peak I_all I_Sum\n")
f1.write("The group is: \nEa_min Eb_min Ec_min Ed_min mean_min I_4peak_nomal I_40_nomal\n")

## The main part of this process to find the 4 bases
for i in xrange(len(Ener_new)):
	jobs = [i,job_server.submit(calculation_tang,(A,Energy,Intensity,Ener_new,Inten_new,i),(),())]
	job = calculation_tang(A,Energy,Intensity,Ener_new,Inten_new,i)
	if job is None: pass
	else:
		Ea,Eb,Ec,Ed = Ener_new[i]
		Ia,Ib,Ic,Id = Inten_new[i]
		mean_res_1.append(job[0])
		mean_res_abs.append(abs(job[0]))
		std_res_1.append(job[1])
		Energy_combi_s[i] = job[2]
		Intensity_combi_single[i] = job[3]
		Solution_A[i] = job[4]

## Nomalize the 40 Intensity Line to know the weight of the 4 base line
		## 1st Step: Find the position of 4 bases' Intensity position in input files
		Ia_id = job[3].index(Inten_new[i][0])
		Ib_id = job[3].index(Inten_new[i][1])
		Ic_id = job[3].index(Inten_new[i][2])
		Id_id = job[3].index(Inten_new[i][3])
		## 2nd Step: Add up 40 Intensity lines
		Inten_combi_s[i] = sum(job[3][:])

		## 3rd Step: Nomalize
		for xx in xrange(len(job[3])):
			Inten_nomal.append(job[3][xx]/Inten_combi_s[i])
		## 4th Step: Get the information of 4 bases Nomalization Result
		Ia_nomalize = job[3][Ia_id]
		Ib_nomalize = job[3][Ib_id]
		Ic_nomalize = job[3][Ic_id]
		Id_nomalize = job[3][Id_id]
		## 5th Step: Know the weight of 4 bases Intensity in all the 40 lines
		Itotal_nomalize.append(Ia_nomalize+Ib_nomalize+Ic_nomalize+Id_nomalize)
		Inten_score.append((Itotal_nomalize[i]+Inten_combi_s[i]))
##		Score_all.append((abs(ma.log(abs(job[0]))+ ma.log(job[1])))+(Itotal_nomalize[i]+Inten_score[i]))
## Find the 4 bases who satisfy the condition
##				if ( abs(job[0]) < 5E-9 and job[1] < 5E-9 ):
##  and Itotal_nomalize[i]>0.1 and Inten_combi_s[i]>0.33
		##if (Ea+Eb+Ec+Ed)/4.79E-8 <17.2:
		f.write("%5.3e %5.3e %5.3e %5.3e %e %e %e %f %f %f \n" % (Ea,Eb,Ec,Ed,job[0],abs(job[0]),job[1],(Itotal_nomalize[i]),(Inten_combi_s[i]),Inten_score[i]))
		Solution_id.append(i)
		Solution_mean.append(abs(job[0]))
		Solution_std.append(abs(job[1]))
		Solution_Inten_base.append(Itotal_nomalize[i])
		Solution_Inten_equ.append(Inten_combi_s[i])
		Solution_bi.append(abs((Eb+Ec+Ed)/Ea-1.716))

## Find all the minimum or maxmium value of all the choice
for mm in xrange(len(Solution_id)):
	Solution_Energy.append(Ener_new[Solution_id[mm]][0]+ Ener_new[Solution_id[mm]][1]+ Ener_new[Solution_id[mm]][2]+ Ener_new[Solution_id[mm]][3])
##	Solution_score.append(1/(Solution_std[mm]+Solution_mean[mm]))
	Solution_score.append((Solution_Inten_equ[mm]*Solution_Inten_base[mm])/(Solution_mean[mm])/(Solution_std[mm])/(Solution_bi[mm]))
	
##	Solution_score.append(abs(ma.log(Solution_Inten_equ[mm])+ma.log(Solution_Inten_base[mm]))-ma.log(Solution_mean[mm])-ma.log(Solution_std[mm]))
##	Solution_score.append((Solution_mean[mm]/Solution_std[mm]))
##print len(Solution_score)
## Give all the group 4 bases that satisfy the mean < 1E-11 and sort by the intensity weight
for ii in xrange(len(Solution_id)):
	Ea_max,Eb_max,Ec_max,Ed_max = Ener_new[Solution_id[Solution_score.index(max(Solution_score))]]
	f1.write("%d %5.3e %5.3e %5.3e %5.3e %e %e %f %fmm/s %f\n" % (ii,Ea_max,Eb_max,Ec_max,Ed_max,mean_res_1[Solution_id[Solution_score.index(max(Solution_score))]],std_res_1[Solution_id[Solution_score.index(max(Solution_score))]],max(Solution_score),((Ea_max+Eb_max+Ec_max+Ed_max)/4.79E-8),((Eb_max+Ec_max+Ed_max)/Ea_max)))
	f1.write("The Energy and Intensity which we can choose %d group: index combination         choose_energy combi_energy intensity\n" % ii)
	for m in xrange(len(A)):
##		print Solution_A[Solution_score.index(max(Solution_score))][m]
##		print Energy_combi_s[Solution_id[Solution_score.index(max(Solution_score))]][m]
		f1.write("%4d %4d %4d %4d %4d %e %e %e %e\n" % (m,A[Solution_A[Solution_score.index(max(Solution_score))][m]][0],A[Solution_A[Solution_score.index(max(Solution_score))][m]][1],A[Solution_A[Solution_score.index(max(Solution_score))][m]][2],A[Solution_A[Solution_score.index(max(Solution_score))][m]][3],Energy_combi_s[Solution_id[Solution_score.index(max(Solution_score))]][m],abs(Ea_max*A[Solution_A[Solution_score.index(max(Solution_score))][m]][0]+Eb_max*A[Solution_A[Solution_score.index(max(Solution_score))][m]][1]+Ec_max*A[Solution_A[Solution_score.index(max(Solution_score))][m]][2]+Ed_max*A[Solution_A[Solution_score.index(max(Solution_score))][m]][3]),(Energy_combi_s[Solution_id[Solution_score.index(max(Solution_score))]][m]-abs(Ea_max*A[Solution_A[Solution_score.index(max(Solution_score))][m]][0]+Eb_max*A[Solution_A[Solution_score.index(max(Solution_score))][m]][1]+Ec_max*A[Solution_A[Solution_score.index(max(Solution_score))][m]][2]+Ed_max*A[Solution_A[Solution_score.index(max(Solution_score))][m]][3])),Intensity_combi_single[Solution_id[Solution_score.index(max(Solution_score))]][m]))
	Solution_score[Solution_score.index(max(Solution_score))] = -10

elapsed = (time.time()-start)
print ("Time used: %5.5f \n" % elapsed)
print('Finished.')
f.close
f1.close
