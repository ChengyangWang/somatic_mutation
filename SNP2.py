### generate mpileup
### read data for each position for treatment and control
### recursively calculation
### weighted fisher.test

import numpy
import subprocess
import time
import math
import os
import random
import scipy
import scipy.stats
from collections import Counter
from optparse import OptionParser

#parser = OptionParser()
#parser.add_option("-o","--output",dest="output",type="str")##output
#parser.add_option("-m","--mixture",dest="mixture",type="float")##mixture probability

#(options,args)=parser.parse_args()

def sp(cmd):
	a=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell='TRUE')
	ac = a.communicate()
	return ac

def generate_mpileup_for_each_chrom_on_specific_region(bam,chrom,start,reading_length): ###
	cmd="samtools mpileup -r %s:%s-%s -s "%c(chrom,str(start),str(start+reading_length))
	return sp(cmd)[0]

def sequencing_depth_filter(sequencing_depth,sequencing_depth_cutoff):
	if sequencing_depth>sequencing_depth_cutoff:
		return True
	else:
		return False

def mapping_quality_filter(map_quality,map_quality_cutoff):
	if map_quality>map_quality_cutoff:
		return True
	else:
		return False
		
def base_quality_filter(base_quality,base_quality_cutoff):
	if base_quality>base_quality_cutoff:
		return True
	else:
		return False

def mutation_filter(base_quality,base_quality_cutoff,map_quality,map_quality_cutoff):
	if base_quality_filter(base_quality,base_quality_cutoff):
		base_p=quality_to_probability(base_quality)
	else:
		return False
	if mapping_quality_filter(map_quality,map_quality_cutoff):
		map_p=quality_to_probability(map_quality)
	else:
		return False
	return [base_p,map_p]
	
def insertion_filter(map_quality,map_quality_cutoff):
	if mapping_quality_filter(map_quality,map_quality_cutoff):
		map_p=quality_to_probability(map_quality)
	else:
		return False
	return map_p
	
def char_to_quality(char):
	return ord(char)-33
	
def quality_to_probability(quality):
	return math.pow(10,-quality/10.0)


def analyze_one_position(line,base_quality_cutoff,map_quality_cutoff,sequencing_depth_cutoff):  ###analyze a list of line in mpileup 
	start=time.time()
	mutation_result=[]
	deletion_result=[]
	insertion_result=[]
	base_list=[]
	
	
	base_char=line[5]
	map_char=line[6]
	bases=line[4]
	
	i=-1
	read_counts=-1
	
	while True:
		i+=1
		
		if i==len(bases):
			break
		
		if bases[i]=="^":
			i+=2
			base=bases[i]
			read_counts+=1
		elif bases[i]=="$":
			i+=1
			base=bases[i]
			read_counts+=1
		elif bases[i]=="*":
			
			read_counts+=1
			base_p_map_p=mutation_filter(char_to_quality(base_char[read_counts]),base_quality_cutoff,char_to_quality(base_char[read_counts]),map_quality_cutoff)
			if not base_p_map_p:
				base_p=base_p_map_p[0]
				map_p=base_p_map_p[1]
			else:
				continue
			
			deletion_result.append([base_p,map_p])### deletion deletion quality and mapping probability 
			continue
		elif bases[i]=="+":
			step=[]
			for j in range(1,4):
				if bases[i+j] in [str(v) for v in range(10)]:
					step.append(int(bases[i+j]))
				else:
					break
			if len(step)==1:  ##default the length is negative than 1000 
				number=step[0]
			elif len(step)==2:
				number=step[0]*10+step[1]
			else:
				number=step[0]*100+step[1]*10+step[2]	
			
			map_p=insertion_filter(char_to_quality(base_char[read_counts]),map_quality_cutoff)
			if not map_p:
				continue
			base_string=bases[i+len(step)+1:i+number+len(step)+1]
			
			insertion_result.append([base_string,map_p])
			
			i+=number+len(step)
			continue
			
		elif bases[i]=="-" :
			step=[]
			for j in range(1,3):
				if bases[i+j] in [str(v) for v in range(10)]:
					step.append(int(bases[i+j]))
				else:
					break
			if len(step)==1:  ##default the length is negative than 1000 
				number=step[0]
			elif len(step)==2:
				number=step[0]*10+step[1]
			else:
				number=step[0]*100+step[1]*10+step[2]	
			
			i+=number+len(step)
			continue
		
		else:
			base=bases[i]
			read_counts+=1
			
		
		base_p_map_p=mutation_filter(char_to_quality(base_char[read_counts]),base_quality_cutoff,char_to_quality(base_char[read_counts]),map_quality_cutoff)
		if base_p_map_p:
			
			base_p=base_p_map_p[0]
			map_p=base_p_map_p[1]
		else:
			continue
		
		if base==base.upper():
			mutation_result.append([base.upper(),base_p,map_p,'+'])### base,base_probability, mapping_probability, strand
		else:
			mutation_result.append([base.upper(),base_p,map_p,'-'])### base,base_probability, mapping_probability, strand
		
		base_list.append(base.upper())
	
	if sequencing_depth_filter(len(mutation_result),sequencing_depth_cutoff):
		pass
	else:
		return False
	print time.time()-start,"dd",len(mutation_result)
	return [mutation_result,deletion_result,insertion_result,base_list] ## extract information on specific region if not information is available, return [[],[],[]]


def probability_calculation_for_mutation_test(mutation_result,base_list,mutation):

	depth=len(mutation_result)-1
	
	try:
		frequency=dict(Counter(base_list).most_common())[mutation]/float(depth)
	except KeyError:
		frequency=1/float(depth)
	a=time.time()
	lower=scipy.stats.binom.ppf(0.05,depth,frequency)
	upper=scipy.stats.binom.ppf(0.95,depth,frequency)-1
	
	lower_quantile_limit=int(depth-lower)
	
	upper_quantile_limit=int(upper)

	
	pre_update=[1.0]
	post_update=[]
	
	pre_update_plus=[1.0]
	post_update_plus=[]
	
	pre_update_minus=[1.0]
	post_update_minus=[]
	
	for i in range(min(lower_quantile_limit,upper_quantile_limit)+1):
		base=mutation_result[i][0]
		base_p=mutation_result[i][1]
		map_p=mutation_result[i][2]
		mutation_p=(1-map_p)*base_p+map_p/4.0
		strand=mutation_result[i][3]
		
		if base==mutation:
			post_update.append(mutation_p*pre_update[0])
			for j in range(1,len(pre_update)):
				post_update.append(pre_update[j]*mutation_p+pre_update[j-1]*(1-mutation_p))
			post_update.append((1-mutation_p)*pre_update[-1])
			pre_update=post_update
			post_update=[]
		
		else:
			post_update.append((1-mutation_p/3)*pre_update[0])
			for j in range(1,len(pre_update)):
				post_update.append(pre_update[j]*(1-mutation_p/3)+pre_update[j-1]*mutation_p/3)
			post_update.append(mutation_p/3*pre_update[-1])
			pre_update=post_update
			post_update=[]
			
	if lower_quantile_limit<=upper_quantile_limit:
		for i in range(lower_quantile_limit+1,upper_quantile_limit+1):
			base=mutation_result[i][0]
			base_p=mutation_result[i][1]
			map_p=mutation_result[i][2]
			mutation_p=(1-map_p)*base_p+map_p/4.0
			strand=mutation_result[i][3]
		
			if base==mutation:
				
				for j in range(1,len(pre_update)):
					post_update.append(pre_update[j]*mutation_p+pre_update[j-1]*(1-mutation_p))
				post_update.append((1-mutation_p)*pre_update[-1])
				pre_update=post_update
				post_update=[]
		
			else:
				
				for j in range(1,len(pre_update)):
					post_update.append(pre_update[j]*(1-mutation_p/3)+pre_update[j-1]*mutation_p/3)
				post_update.append(mutation_p/3*pre_update[-1])
				pre_update=post_update
				post_update=[]
	else:
		for i in range(upper_quantile_limit+1,lower_quantile_limit+1):
			base=mutation_result[i][0]
			base_p=mutation_result[i][1]
			map_p=mutation_result[i][2]
			mutation_p=(1-map_p)*base_p+map_p/4.0
			strand=mutation_result[i][3]
		
			if base==mutation:
				post_update.append(mutation_p*pre_update[0])
				for j in range(1,len(pre_update)):
					post_update.append(pre_update[j]*mutation_p+pre_update[j-1]*(1-mutation_p))
				pre_update=post_update
				post_update=[]
		
			else:
				post_update.append((1-mutation_p/3)*pre_update[0])
				for j in range(1,len(pre_update)):
					post_update.append(pre_update[j]*(1-mutation_p/3)+pre_update[j-1]*mutation_p/3)
				pre_update=post_update
				post_update=[]
		
	for i in range(max(lower_quantile_limit,upper_quantile_limit)+1,depth+1):
		base=mutation_result[i][0]
		base_p=mutation_result[i][1]
		map_p=mutation_result[i][2]
		mutation_p=(1-map_p)*base_p+map_p/4.0
		strand=mutation_result[i][3]
		
		if base==mutation:
			for j in range(1,len(pre_update)):
				post_update.append(pre_update[j]*mutation_p+pre_update[j-1]*(1-mutation_p))
			pre_update=post_update
			post_update=[]
		
		else:
			for j in range(1,len(pre_update)):
				post_update.append(pre_update[j]*(1-mutation_p/3)+pre_update[j-1]*mutation_p/3)
			pre_update=post_update
			post_update=[]
				
	
	print 'f',time.time()-a
	return [pre_update,pre_update_plus,pre_update_minus] ###[[p,p,p,..],[p,p,p],[p,p,p,p,p]] 
	
	


def probability_calculation_for_mutation(mutation_result,mutation):
	a=time.time()
	pre_update=[1.0]
	post_update=[]
	
	pre_update_plus=[1.0]
	post_update_plus=[]
	
	pre_update_minus=[1.0]
	post_update_minus=[]
	

	
	for i in range(len(mutation_result)):
		
		base=mutation_result[i][0]
		base_p=mutation_result[i][1]
		map_p=mutation_result[i][2]
		mutation_p=(1-map_p)*base_p+map_p/4.0
		strand=mutation_result[i][3]
		
		if base==mutation:
			
			if strand=="+":
				post_update_plus.append(mutation_p*pre_update_plus[0])
				for j in range(1,len(pre_update_plus)):
					post_update_plus.append(pre_update_plus[j]*mutation_p+pre_update_plus[j-1]*(1-mutation_p))
				post_update_plus.append((1-mutation_p)*pre_update_plus[-1])
				pre_update_plus=post_update_plus
				post_update_plus=[]
			else:
				
				post_update_minus.append(mutation_p*pre_update_minus[0])
				for j in range(1,len(pre_update_minus)):
					post_update_minus.append(pre_update_minus[j]*mutation_p+pre_update_minus[j-1]*(1-mutation_p))
				post_update_minus.append((1-mutation_p)*pre_update_minus[-1])
				pre_update_minus=post_update_minus
				post_update_minus=[]
				
			post_update.append(mutation_p*pre_update[0])
			for j in range(1,len(pre_update)):
				post_update.append(pre_update[j]*mutation_p+pre_update[j-1]*(1-mutation_p))
			post_update.append((1-mutation_p)*pre_update[-1])
			pre_update=post_update
			post_update=[]
		
		else:
			if strand=="+":
				post_update_plus.append((1-mutation_p/3)*pre_update_plus[0])
				for j in range(1,len(pre_update_plus)):
					post_update_plus.append(pre_update_plus[j]*(1-mutation_p/3)+pre_update_plus[j-1]*mutation_p/3)
				post_update_plus.append(mutation_p/3*pre_update_plus[-1])
				pre_update_plus=post_update_plus
				post_update_plus=[]
				
			else:
				post_update_minus.append((1-mutation_p/3)*pre_update_minus[0])
				for j in range(1,len(pre_update_minus)):
					post_update_minus.append(pre_update_minus[j]*(1-mutation_p/3)+pre_update_minus[j-1]*mutation_p/3)
				post_update_minus.append(mutation_p/3*pre_update_minus[-1])
				pre_update_minus=post_update_minus
				post_update_minus=[]
			
			
			post_update.append((1-mutation_p/3)*pre_update[0])
			for j in range(1,len(pre_update)):
				post_update.append(pre_update[j]*(1-mutation_p/3)+pre_update[j-1]*mutation_p/3)
			post_update.append(mutation_p/3*pre_update[-1])
			pre_update=post_update
			post_update=[]
			
		
	
	print 'f',time.time()-a
	return [pre_update,pre_update_plus,pre_update_minus] ###[[p,p,p,..],[p,p,p],[p,p,p,p,p]] 
	
def probability_calculation_for_deletion(deletion_result):
	pre_update=[1.0]
	post_update=[]
	
	for i in range(len(deletion_result)):
		
		base_p=deletion_result[i][0]
		map_p=deletion_result[i][1]
		
		post_update.append(base_p*pre_update[0])
		for j in range(1,len(pre_update)):
			post_update.append(pre_update[j]*base_p+pre_update[j-1]*(1-base_p))	
		post_update.append((1-base_p)*pre_update[-1])
		
		pre_update=post_update
		
		post_update=[]
			
		
	return pre_update ###[p,p,p,..] 

def probability_calculation_for_insertion(insertion_result):
	if insertion_result==[]:
		return [[[]],[[]]]
	max_insertion_length=max([len(v[0]) for v in insertion_result])
	pre_update_list=[]
	post_update_list=[]
	base_list=[]
	for i in range(max_insertion_length):
		pre_update_list.append([1.0])
		post_update_list.append([])
		base_list.append([])
		
	
	for i in range(len(insertion_result)):
		
		base_string=insertion_result[i][0]
		map_p=insertion_result[i][1]
		
		for v in range(len(base_string)):
			base=base_string[v]
			post_update_list[v].append(map_p*pre_update_list[v][0])
			for j in range(1,len(pre_update_list[v])):
				post_update_list[v].append(pre_update_list[v][j]*map_p+pre_update_list[v][j-1]*(1-map_p))	
			post_update_list[v].append((1-map_p)*pre_update_list[v][-1])
			
			
			base_list[v].append(base.upper())
			pre_update_list[v]=post_update_list[v]
			post_update_list[v]=[]
	
	
	for i in range(len(base_list)):
		base_list[i]=Counter(base_list[i]).most_common(1)
	
	return [pre_update_list,base_list] ###for each position[[p,p,p,p,p],[p,p,p,p],[p]],['A','G'...]


def weighted_avg_and_std(values, weights,depth):
    """
    Returns the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    V1=weights.sum()
    norm_weights=weights/V1
    average=numpy.dot(values, norm_weights)
    p_estimate=average/float(depth)
    
    V2=(norm_weights**2).sum()
    D2=average*(1-p_estimate)
    D22=numpy.dot(weights, (values-average)**2)/(1-V2)
    return (p_estimate,V2*D22/(depth**2))
    
    

def weighted_fisher_test_for_mutation_test(treatment,control):
	a=time.time()
	treatment_read_depth=len(treatment[0])-1
	control_read_depth=len(control[0])-1
	treatment_mean_variance=weighted_avg_and_std(numpy.array(range(treatment_read_depth+1)), numpy.array(treatment[0]),1490)
	control_mean_variance=weighted_avg_and_std(numpy.array(range(control_read_depth+1)), numpy.array(control[0]),717)
	b=(treatment_mean_variance[0]-control_mean_variance[0])/math.sqrt(treatment_mean_variance[1]+control_mean_variance[1])
	print treatment_mean_variance[0]*1490,control_mean_variance[0]*717
	print 1-scipy.stats.norm.cdf(b)
	
	
	print "l",time.time()-a
	
	

def weighted_fisher_test_for_mutation(treatment,control,accuracy):
	a=time.time()
	treatment_read_depth=len(treatment[0])-1
	control_read_depth=len(control[0])-1
	
	approximate_cutoff=math.sqrt(accuracy/(treatment_read_depth*control_read_depth))
	
	treatment_index=[]
	for i in range(len(treatment[0])):
		if treatment[0][i]>=approximate_cutoff:
			treatment_index.append(i)
	
	control_index=[]
	for i in range(len(control[0])):
		if control[0][i]>=approximate_cutoff:
			control_index.append(i)
	
	
	p_value=0.0
	for i in treatment_index:
		for j in control_index:
			treatment_match=i
			control_match=j
			treatment_mismatch=treatment_read_depth-i
			control_mismatch=control_read_depth-j
			fisher_p_value=scipy.stats.fisher_exact([[treatment_mismatch,treatment_match],[control_mismatch,control_match]],alternative="less")[1]
			p_value+=fisher_p_value*treatment[0][i]*control[0][j]
			
	print "l",time.time()-a
	return p_value
	
def weighted_fisher_test_for_strand(treatment,accuracy):
	treatment_plus_read_depth=len(treatment[1])-1
	treatment_minus_read_depth=len(treatment[2])-1
	
	
	approximate_cutoff=math.sqrt(accuracy/(treatment_plus_read_depth*treatment_minus_read_depth))
	
	plus_index=[]
	for i in range(len(treatment[1])):
		if treatment[1][i]>=approximate_cutoff:
			plus_index.append(i)
			
	
	
	minus_index=[]
	for i in range(len(treatment[2])):
		if treatment[2][i]>=approximate_cutoff:
			minus_index.append(i)
			
	
			
		
	p_value_strand=0.0
	for i in plus_index:
		for j in minus_index:
			plus_match=i
			minus_match=j
			plus_mismatch=treatment_plus_read_depth-i
			minus_mismatch=treatment_minus_read_depth-j
			fisher_p_value=scipy.stats.fisher_exact([[plus_mismatch,plus_match],[minus_mismatch,minus_match]])[1]
			p_value_strand+=fisher_p_value*treatment[1][i]*treatment[2][j]
			
	return p_value_strand	
	
def weighted_fisher_test_for_deletion(treatment,control,treatment_mutation_depth,control_mutation_depth,accuracy):
	
	
	treatment_deletion_depth=len(treatment)-1
	control_deletion_depth=len(control)-1
	
	
	
	if control_deletion_depth==0:
	
		p_value=0.0
		approximate_cutoff=math.sqrt(accuracy/treatment_deletion_depth)
		
		treatment_index=[]
		for i in range(len(treatment)):
			if treatment[i]>=approximate_cutoff:
				treatment_index.append(i)
		
		for i in treatment_index:
			
			treatment_deletion=i
			control_deletion=0
			treatment_not_deletion=treatment_mutation_depth
			control_not_deletion=control_mutation_depth
			print treatment_not_deletion,treatment_deletion,control_not_deletion,control_deletion
			fisher_p_value=scipy.stats.fisher_exact([[treatment_not_deletion,treatment_deletion],[control_not_deletion,control_deletion]],alternative="less")[1]
			p_value+=fisher_p_value*treatment[i]
			
		return p_value 
			
	else:
	
		p_value=0.0
		approximate_cutoff=math.sqrt(accuracy/(treatment_deletion_depth*control_deletion_depth))
	
		treatment_index=[]
		for i in range(len(treatment)):
			if treatment[i]>=approximate_cutoff:
				treatment_index.append(i)
	
		control_index=[]
		for i in range(len(control)):
			if control[i]>=approximate_cutoff:
				control_index.append(i)
				
		for i in treatment_index:
			for j in control_index:
				
				treatment_deletion=i
				control_deletion=j
				treatment_not_deletion=treatment_mutation_depth
				control_not_deletion=control_mutation_depth
				fisher_p_value=scipy.stats.fisher_exact([[treatment_not_deletion,treatment_deletion],[control_not_deletion,control_deletion]],alternative="less")[1]
				p_value+=fisher_p_value*treatment[i]*control[j]
				
		return p_value
				
def weighted_fisher_test_for_insertion(treatment,control,treatment_mutation_depth,control_mutation_depth,accuracy):
	max_treatment_insertion_length=len(treatment[0])
	max_control_insertion_length=len(control[0])
	p_value_list=[]
	for position in range(max_treatment_insertion_length):
		if position<=max_control_insertion_length-1:
		
			treatment_insertion_depth=len(treatment[0][position])-1
			control_insertion_depth=len(control[0][position])-1
			p_value=0.0
			approximate_cutoff=math.sqrt(accuracy/(treatment_insertion_depth*control_insertion_depth))
			
			treatment_index=[]
			for i in range(len(treatment[0][position])):
				if treatment[0][position][i]>=approximate_cutoff:
					treatment_index.append(i)
	
			control_index=[]
			for i in range(len(control[0][position])):
				if control[0][position][i]>=approximate_cutoff:
					control_index.append(i)
				
			for i in treatment_index:
				for j in control_index:
					
					treatment_insertion=i
					control_insertion=j
					treatment_not_insertion=treatment_mutation_depth
					control_not_insertion=control_mutation_depth
					fisher_p_value=scipy.stats.fisher_exact([[treatment_not_insertion,treatment_insertion],[control_not_insertion,control_insertion]],alternative="less")[1]
					p_value+=fisher_p_value*treatment[0][position][i]*control[0][position][j]
					
		else:
			
			treatment_insertion_depth=len(treatment[0][position])-1
			
			p_value=0.0
			approximate_cutoff=math.sqrt(accuracy/treatment_insertion_depth)
			
			treatment_index=[]
			for i in range(len(treatment[0][position])):
				if treatment[0][position][i]>=approximate_cutoff:
					treatment_index.append(i)
	
			
			for i in treatment_index:
				
				treatment_insertion=i
				control_insertion=0
				treatment_not_insertion=treatment_mutation_depth
				control_not_insertion=control_mutation_depth
				fisher_p_value=scipy.stats.fisher_exact([[treatment_not_insertion,treatment_insertion],[control_not_insertion,control_insertion]],alternative="less")[1]
				p_value+=fisher_p_value*treatment[0][position][i]
				
		p_value_list.append(p_value)
	
	return p_value_list

def main():
	t="""1	10017	N	1500	CCCCCCcCC+2CTCC*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC+2CTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcCCCCCCCCCCCCCCCCCCCCCCCCCCCCC+2CTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCCCCC+2CTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCccccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcccccccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCccccc^!C^!C^!C^5C^!C^"C^]C^!C^!C^8C^!C^!C^!C^!C^!C^OC^!C^!C^!C^!C^!C^!C^"C^!C^!C^!C^!C^!C^!C^"C^!C^*C^!C^8C^!C^*C^!C^OC^!C^!C^8C^!C^!C^FC^!C^!C^2C^!C^!C^!C^!C^%C^!C^!C^!C^!C^!C^>C^#C^!C^!C^!C^*C^!C^!C^!C^!C^!C^!c^!c^5c^Fc	91#9DE7#5988##99999999999999999999999999999:999999999999999999999999999?9??>?????????DEEEEDEEDEEEEEEEEEEEFEEFFEFE@EEEEE?>??999:9999999999999999999999999999999999999999999999999999999999999???????EDDDD@DDEDEEDEEEEEEEFEEE=E?9999999999999999999999999999999899999:9999999999>;??#;?>?DDDDDDDDEED@DDCDDEEEE>EDDDEEEDDDDDDEEEEEE999999999999:999999999999999999999:9999999999999999999999:9999::9>>?>??>??DDEEEEDDDBE;EEEDDD;EDDDDDDDDDDDDEEEDEEDC?>7::9::99:9:9::9999::9::99:99:99::9::9:::999::99:9::::99:::9::99::::::99:::::9:9:::9:::9>>>>>>?>@>??>9DDDDDDDDDDDDDDDDDDDDEEEEEEEEEDDDDEEEEDFDCDDDDDEDDDFDD>::::::::::::::::9:9:99:9::9::::999:::::::::::::::::::9:::::9::9::>1?>??>>DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD>DDDDCDDDDEDDDDDED@DDDD?::::::::::::::::::::::::::9:::::::>?>?><DCCDDDDDDDDDDDDDDDDDDDDEDDDDBD:::::::::6::::::::::::::::::::':::::::::??9DDDDCCDDDDDDDDD@DDDDDDDDDDDDC:;::::;:;;;;;;;:::;:::;;:;;:;;;;;9;:;::;;::;:;;;;:;:6:;=?:==>CD<CDDDDCD2DCDCCDCDDDA=DDDDDDDD>?455;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;<=<<C>DCACCDCDDDDA;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;<;;;;=6?;<<>DCCCDDBBBB:DDDDDBDBC<><;;<;;<;;;;<;;<;<;;<<<;;;<<<<;<;<<;<<>;<>>;DD>D?ADADD>ADAADDBD;>3A@=<=<<====<=<====<<<<<<=<<<<<=====<<<=<=<=<<====<=<?;<<?;;@D@AA@@D@@@@D@@@@D;>@=============================<=======<==?<?BCCCC=CC?CCABCBCC;;00<<>>>>>>>>>=>=>>>>>>>>>>>>><>>==>>>>><><>>>>>>?>>>>>CBC6C<CCBB?CC5ACCBCC?CC>//50BEAA?BB<BB>BBBBBBBBAB>BB?B@@@BBBBBBBBBB=<==BC@?CBCBD<1/9D7@@:@9>9@@99@999999999@9999@999999999@=9A;A;?;?;C=@C@@@@CCC@AC>AA7C?B#.##	!!+!!!+!!!*!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"!!!!!!!!!!!!!"!!!!!!!!!-$!-*!$!-0!>!!!!!"!>!0!%!!!>!!!!!!!%!!!!!!!!>!%!!**-!!!5!!!!!!!!!"!!!!!!!!!!!!!!!!!!!!$"!""!!!!!!!"!!!!"!!!!!!!!!!!!!!!!>*+!!!%!-!!$>*!>!!!!!!!!*!!!!!!!!!!!!!!!!$*8!!!!!!/!!!!!!!!!!!!$!!!!!!!!!!*!!!!+2!2!!!*!$>!%!!*!!!!$%$$*3!!>$!!,$$!5$!>$!!!$!*$!!!$$!!!!!*!*/.!!!!!'*!*!!$!!**!!!!!!!!/5!!!:**!!!*$!&$!!*!$$!!*/!">!>5!!!!2!!!!*$$!!!*$$$$!$!$!$!!!*$!$3$$0*2!!0!0!0-*!3/!'/!!!!!!!*!*!*!/!$!!$*!$!'!!/!!!!!!!/"$!$!$*!$!!!!!/!!!"!!!&!*!//"!$!!!!!!!!!!>!$/!$!!""!!+*!!&!*3>!!!*$%3!!$!*$*!!!$(!!!!$$!!8!!$!$!!*!*$****7!$>!1!!/!/./*!!!>!!!!>!!!!$!*!*$)%>!!!!/!*!$$!!*!*!*!*>*//*!/!!!*$>!!!**!!>!20!-$$*!$$$3!$!$2>>$1!!$!$!$4$>*!**!!>!2!$>>*>*>!!!$10!*!$$1$5!!"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!$!!!"!!>!!$!!!!!!!!!!!!!!>!!!!!%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%!!!!!%!!!!!!!!!""%!!-!!!!!!1!!!!!!!!!!F!5!!!!!!!!!!!!!!!!!!"!!!!!!!!!!!!!!!!$!!#!!!!%!!!!(!!!!!!!!!!!!!!!!!!%!!!+!!!!!!!!!1!!!!!*.!!!!"*'!1!!!!!!!!!#!!!!!!!!!!%!!!!!!!!!!!!!!!!!!"!!!L!!!!!!!!!!!!!!!!!!!!!"!!!!!!$!!!!!!!!!!!%!2!!!!!F!>!%%F!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!L!!!!--!!!!!!!!!!!!!!%"!!!4!!FF!!!!!!!!!,!!F!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"!!!!!!!2$!!!2"*!!!!!!!!!!!T!!!!!!!!!!!!!!!!!!!0!!!!!!!!!!!!!!!!!!!!!!!!!!!!>!2%!!!!!!!!!!!!!!!!!!8FF!!!!!!8!!!!!!!!!!!3!!!!!!O!!!!!!!!!!!!!!!"!!!!!#>!!!!!!!!!F!#!!!!!!!!!!!!!!88!83!F!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!$!!!!!0!>!$!58FF!!!!5!"]!!8!!!!!O!!!!!!"!!!!!!"!*!8!*!O!!8!!F!!2!!!!%!!!!!>#!!!*!!!!!!!5F"""
	c="""1	10006	N	721	CCCCCCcCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC+1TCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-1NCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC^/C^!C^/C^.C^/C^*C^!C^!C^!C^>C^!C^!C^!C^!C^>C^!C^!C^!C^!C^$C^!C^*C^!C^*C^$C^)C^%C^>C^!C^!C^!C^!C^/C^!C^*C^!C^$C^$C^!C^!C^*C^!C^*C^!C^*C^!C^*C^>C^*C^/C^/C^*C^!C^/C^!C^!C^!C^*C^$C^>C^!C^!C^!C^*C^*C^!C^!C^>C^!C^2C^0C^!C^-C^$C^$C^*C^!C^$C^$C^$C^3C^!C^$C^!C^$C^2C^>C^>C^$C^1C^!C^!C^$C^!C^$C^!C^$C^4C^$C^>C^*C^!C^*C^*C^!C^!C^>C^!C^2C^!C^$C^>C^>C^*C^>C^*C^>C^!C^!C^!C^$C^1C^0C^!C^*C^!C^$C^$C^1C^$C	::#;CD6#8<;2##<<;<<;<<<<<<<<<<<<<<<;<<<;;;;<<<;<<<<<;<;;<<<;;<;<;<<;<;<?>><>>:>;?;;<>ADDDDADAADDDDDDDDDDDABDABDAD:BCDDD>=<?<<<<<<<=<==<====<=<=<====<<<===<=<=<<==<===<<<<===<<<==<=<=<=<==<??;=??;D@@AAA@@D@D@CDCCDADD@DDDAD;=================================================?@=6?===BCBB=BCBCCA+AACAACCD244A?=7CCABB8CCCCCAAA>>>>>>>>>>>>=<>>>>>>>>>>>>>>>>>>>>=>>>>>>>>>>>>>>>=>=>>>>=>>>>==>>B?>??>>;BACCCCCCBBC>CCCBC;@C>C?@CC@CCCCC??CCBDCC?>0BBBBBBBABABBBBBBBBABBABBBBBBBBBBBBBBBBBBBBABBBBBBBBABBBBBBBBBBBBBBABBBBBBB?BBBBBBBBBBB==><====>=>>=CCCBCCCCCBCCCCCCCCBCBCCCCCCCCCCCBBCCCC@DBAAACCCBAAAD9C=9@@@9?@9@99@@9=@9@9@99@9@@9@@?@9999=A;A9A999@AA99AAA?9AA@AA9AA9>AA;A;?A:;CCCCCCCCC@@@@@@CC@CCCCCC@C@C@CCCCACCCCAACCAA7CBBBAC=CCCCA	!!+!!!+!!!*!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"!!!!!!!!!!!!!"!!!!!!!!!-$!-*!$!-0!>!!!!!"!>!0!%!!!>!!!!!!!%!!!!!!!!>!%!!**-!!!5!!!!!!!!!"!!!!!!!!!!!!!!!!!!!!$"!""!!!!!!!"!!!!"!!!!!!!!!!!!!!!!>*+!!!%!-!!$>*!>!!!!!!!!*!!!!!!!!!!!!!!!!$*8!!!!!!/!!!!!!!!!!!!$!!!!!!!!!!*!!!!+2!2!!!*!$>!%!!*!!!!$%$$*3!!>$!!,$$!5$!>$!!!$!*$!!!$$!!!!!*!*/.!!!!!'*!*!!$!!**!!!!!!!!/5!!!:**!!!*$!&$!!*!$$!!*/!">!>5!!!!2!!!!*$$!!!*$$$$!$!$!$!!!*$!$3$$0*2!!0!0!0-*!3/!'/!!!!!!!*!*!*!/!$!!$*!$!'!!/!!!!!!!/"$!$!$*!$!!!!!/!!!"!!!&!*!//"!$!!!!!!!!!!>!$/!$!!""!!+*!!&!*3>!!!*$%3!!$!*$*!!!$(!!!!$$!!8!!$!$!!*!*$****7!$>!1!!/!/./*!!!>!!!!>!!!!$!*!*$)%>!!!!/!*!$$!!*!*!*!*>*//*!/!!!*$>!!!**!!>!20!-$$*!$$$3!$!$2>>$1!!$!$!$4$>*!**!!>!2!$>>*>*>!!!$10!*!$$1$"""
	t=t.strip().split()
	c=c.strip().split()
	treatment=analyze_one_position(t,5,0,30)
	control=analyze_one_position(c,5,0,30)
	
	treatment_mutation=treatment[0]
	control_mutation=control[0]
	treatment_mutation_p=probability_calculation_for_mutation_test(treatment[0],treatment[3],'A')
	print treatment_mutation_p
	control_mutation_p=probability_calculation_for_mutation_test(control[0],control[3],'A')
	print control_mutation_p
	
	weighted_fisher_test_for_mutation_test(treatment_mutation_p,control_mutation_p)
	#print weighted_fisher_test_for_deletion(treatment_mutation_p,control_mutation_p,1500,720,0.01)
	#print weighted_fisher_test_for_mutation(treatment_mutation_p,control_mutation_p,0.5)
	#print scipy.stats.fisher_exact([[1000,10],[1000,1]],alternative="less")[1]
	#print weighted_fisher_test_for_strand(treatment_mutation_p,0.5)
	
	
def initiate(workspace,chrom):
	file=open(workspace+'/'+chrom+'_'+'variants.txt','w')
	file.close()
	return workspace+'/'+chrom+'_'+'variants.txt'
	
def handle_files_and_output(treatment_bam,control_bam,output_path,base_quality_cutoff,map_quality_cutoff,sequencing_depth_cutoff,start,reading_length=100000,accuracy=1e-2):
	
	chrom=1
	treatment_list=generate_mpileup_for_each_chrom_on_specific_region(treatment_bam,chrom,start,reading_length).strip().split("\n")
	control_list=generate_mpileup_for_each_chrom_on_specific_region(control_bam,chrom,start,reading_length).strip().split("\n")
	
	output_file=open(output_path,'w+')  ###
	
	
	treatment_length=len(treatment_list)
	control_length=len(control_list)
	control_pivot=0
	
	for i in range(treatment_length):
		t_line=treatment_list[i].split()
		t_position=int(t_line[1])
		for j in range(control_pivot,control_length):
			c_line=control_list[j].split()
			c_position=int(c_line[1])
			if c_position<t_position:
				pass
			elif c_position==t_position:
				treatment_information=analyze_one_position(t_line,base_quality_cutoff,map_quality_cutoff,sequencing_depth_cutoff)
				if not treatment_information:
					print "treatment on %s:%s is not qualified"%c(chrom,str(t_position))
					break
				
				control_information=analyze_one_position(c_line,base_quality_cutoff,map_quality_cutoff,sequencing_depth_cutoff)
				if not treatment_information:
					print "control on %s:%s is not qualified"%c(chrom,str(t_position))
					break
				
				treatment_mutation=treatment_information[0]
				control_mutation=control_information[0]
				
				#######mutation
				treatment_base_list=treatment_information[3]
				control_base_list=control_information[3]
				treatment_potential_mutations=Counter(treatment_base_list).most_common()[1:]
				
				
				for mutation in treatment_potential_mutations:
					treatment_mutation_p=probability_calculation_for_mutation(treatment_mutation,mutation)
					control_mutation_p=probability_calculation_for_mutation(control_mutation,mutation)
					print weighted_fisher_test_for_mutation(treatment_mutation_p,control_mutation_p,accuracy)
					print weighted_fisher_test_for_strand(treatment_mutation,accuracy)
					
				####deletion
				treatment_deletion=treatment_information[1]
				treatment_deletion_p=probability_calculation_for_deletion(treatment_deletion)
				if treatment_deletion_p==[1.0]:
					print "no deletion on %s:%s"%c(chrom,str(t_position))
				else:
					pass
				
				control_deletion=control_information[1]
				control_deletion_p=probability_calculation_for_deletion(control_deletion)
				print weighted_fisher_test_for_deletion(treatment_deletion_p,control_deletion_p,len(treatment_base_list),len(control_base_list),accuracy)
				
				###insertion
				treatment_insertion=treatment_information[2]
				treatment_insertion_p=probability_calculation_for_insertion(treatment_insertion)
				
				if treatment_insertion_p==[[[]],[[]]]:
					print "no insertion on %s:%s"%c(chrom,str(t_position))
				else:
					pass
					
				control_insertion=control_information[2]
				control_insertion_p=probability_calculation_for_insertion(control_insertion)
				print weighted_fisher_test_for_insertion(treatment_insertion_p,control_insertion_p,len(treatment_base_list),len(control_base_list),accuracy)
				
				control_pivot=j+1
			
			else:
				print "no available control on %s:%s"%c(chrom,str(t_position)) 
				control_pivot=j
				break
	
def read_bed(bed_path):
	file=open(bed_path)
	a=[]
	while True:
		line=file.readline()
		if not line:break
		line=line.strip().split()
		a.append([line[0],line[1],line[2]]) ###chrom,start,end
	return a 

main()