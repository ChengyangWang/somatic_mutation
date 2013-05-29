### generate mpileup
### read data for each position for treatment and control
### recursively calculation
### weighted fisher.test


import sys
import subprocess
import argparse
import time
import math
import os
import random

import numpy
import scipy
import scipy.stats
import scipy.misc
from scipy.special import beta
from scipy.special import betainc
from scipy.misc import comb
from scipy.stats import chi2_contingency
import argparse
import optparse

def prepare_argparser ():
    """Prepare optparser object. New options will be added in this
    function first.
    
    """
    description = "DETECTION ---  Detect mutation in paired sequnencing sample"
    
    argparser = argparse.ArgumentParser(description = description)
      
    DETECTION_VERSION = '1.0.0'
    argparser.add_argument("--version", action="version", version="%(prog)s "+DETECTION_VERSION)
    argparser.add_argument("-t","--treatment",dest="treatment",type=str,required = True,help="Treatment alighment in bam format with header. The bam must be sorted and indexed properly")
    argparser.add_argument("-c","--control",dest="control",type=str,required = True,help="Treatment alighment in bam format with header. The bam must be sorted and indexed properly")
    argparser.add_argument("--t_name",dest = "treatment_name",type = str, help = "Input the name of treatment sample",default = "tumor")
    argparser.add_argument("--c_name",dest = "control_name",type = str, help = "Input the name of control sample",default = "control")
    argparser.add_argument("--species",dest = "species",choices = ("hg", "mm"), required = True,help = "specify species between human and mouse, it can be hg for human and mm for mouse")
    argparser.add_argument("--base_q",dest = "base_quality",type = int, help = "Input a number from 0 to 50 as a threshold of phred-like score to filter out low quality base. DEFAULT = 10, keeping bases > 10",default = 10)
    argparser.add_argument("--map_q",dest = "mapping_quality",type = int, help = "Input a number from 0 to 50 as a threshold of phred-like score to filter out low mapping quality read. DEFAULT = 0, keeping reads > 0",default = 0)
    argparser.add_argument("--seq_depth",dest = "sequencing_depth",type = int, help = "Input a number greater than 10 as a threshold to filter out low coverage sites. Only qualified reads and bases instructed by --base_q and --map_q are used to calculate coverage. \
    DEFAULT = 10, keeping sites with coverage > 10.",default = 10)
    argparser.add_argument("-p","--p_value",dest = "p_value",type = float, help = "Input a number between 0 and 1 as a significance threshold for SNV detection",default = 0.01)
    argparser.add_argument("-k","--kind", dest = "kind",  choices = ("whole","chrom","region","region_append","position","position_append"), required = True, help="specify the range where DETECTION detects SNV. \
    whole for scanning all positions on treatment sample; chrom for scanning certain chromsomes to call SNVs. You would use --chrom to indicate your interested chromsome afterwards.\
    region and region_append for scanning specified regions to call SNVs. You would assign a file in bed format to --region or specify regions in command line using --region_append afterwards. position and position_append for scanning your interested positions. \
    You would assign a bed file to --position or specify positions in command line using --position_append afterwards")
    argparser.add_argument("--chrom",dest="chrom",action='append',help="Once specify -k as chrom, please input your interested chromsome. For example, --chrom chr1 --chrom chr2.")
    argparser.add_argument("--region",dest="region",type=str,help="Once specify -k as region, please indicate interested regions using file in bed format with three columns, in which the format of each line is\
    chromsome\tstart\tend\n. For example, chr1\t100\t300\nchrX\t400\t1000\n. DETECTION will ignore other columns.")
    argparser.add_argument("--region_append",dest="region_append",action='append',help="Once specify -k as region_append, please type chromosome:start-end in command line. For example, \
    --region_append chr1:100-5000 --region_append chrX:1000-5000")
    argparser.add_argument("--position",dest="position",type=str,help="Once specify -k as position, please indicate interested positions using file in bed format as chromosome\tposition\tposition\n. For example, the file is \
    chr1\t400\t400\nchrX\t600\t600\n")
    argparser.add_argument("--position_append",dest="position_append",action='append',help="Once specify -k as position_append, please type chromosome:position in command line. For example, \
    --position_append chr1:100 --position_append chrX:300")
    
    return argparser

def validate_argparser(argparser):
	options = optparser.parse_args()
	return options
    
def sp(cmd):
	a=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell='TRUE')
	ac = a.communicate()
	return ac

def generate_mpileup_for_each_chrom_on_specific_region(bam,chrom,start,reading_length): ###
	cmd="samtools mpileup -r %s:%s-%s -s %s"%(chrom,str(start),str(start+reading_length),bam)
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

def minimum_mutation_counts(mutation_observation_counts,mutation_observation_counts_cutoff):
	if mutation_observation_counts>mutation_observation_counts_cutoff:
		return True
	else:
		return False
		
def char_to_quality(char):
	return ord(char)-33
	
def quality_to_probability(quality):
	return math.pow(10,-quality/10.0)

def weighted_avg(values, weights):
    """
    Returns the weighted average
    """
    V1=sum(weights)
    return sum(map(lambda x,y:x*y,values,weights))/V1
    #norm_weights=weights/V1
    #average=numpy.dot(values, norm_weights)
    #p_estimate=average/float(depth)
    #
    #V2=(norm_weights**2).sum()
    #D2=average*(1-p_estimate)
    #D22=numpy.dot(weights, (values-average)**2)/(1-V2)
    #return (p_estimate,V2*D22/(depth**2))
    
def parser_one_base_dic(base_dic):
	base_list=base_dic.items()
	base_list.sort(key=lambda x:x[1],reverse=True)
	return [(base_list[i][0],base_list[i][1]) for i in range(4) if base_list[i][1]!=0]  ###include the most frenquency values
	
def parser_two_base_dic(first_base_dic,second_base_dic):
	merge_list=[(key,first_base_dic[key]+second_base_dic[key]) for key in ['A','T','C','G']]
	merge_list.sort(key=lambda x:x[1],reverse=True)
	return [(merge_list[i][0],merge_list[i][1]) for i in range(4) if merge_list[i][1]!=0]	###include the most frenquency values

def analyze_one_position(line,base_quality_cutoff,map_quality_cutoff,sequencing_depth_cutoff):  ###analyze a list of line in mpileup
	
	if int(line[3])>sequencing_depth_cutoff:
		pass
	else:
		return False
	
	start=time.time()
	mutation_result=[]
	deletion_result=[]
	insertion_result=[]
	mutation_plus={'A':0,'T':0,'C':0,'G':0}
	mutation_minus={'A':0,'T':0,'C':0,'G':0}
	
	base_char=line[5]
	map_char=line[6]
	bases=line[4]
	
	i=-1
	read_counts=-1
	plus_depth=0
	minus_depth=0
	
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
			if base_p_map_p:
				base_p=base_p_map_p[0]
				map_p=base_p_map_p[1]
			else:
				continue
			
			deletion_result.append([base_p,map_p])### deletion probability and mapping probability 
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
			base_string=[v.upper() for v in bases[i+len(step)+1:i+number+len(step)+1]]
			
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
		
		if base==".":
			plus_depth+=1
			mutation_result.append([line[2],base_p,map_p])### base,base_probability, mapping_probability
			##mutation_plus[line[2]]+=1
		elif base==",":
			minus_depth+=1
			mutation_result.append([line[2],base_p,map_p])### base,base_probability, mapping_probability
			##mutation_minus[line[2]]+=1
		elif base==base.upper():
			plus_depth+=1
			mutation_result.append([base.upper(),base_p,map_p])### base,base_probability, mapping_probability
			mutation_plus[base.upper()]+=1
		else:
			minus_depth+=1
			mutation_result.append([base.upper(),base_p,map_p])### base,base_probability, mapping_probability
			mutation_minus[base.upper()]+=1 
			
		
	
	if sequencing_depth_filter(len(mutation_result),sequencing_depth_cutoff):
		pass
	else:
		return False
	
	plus_m=0
	minus_m=0
	for allele in ["A","T","C","G"]:
		if allele!=line[2]:
			plus_m+=mutation_plus[allele]
			minus_m+=mutation_minus[allele]
	mutation_plus[line[2]]=plus_depth-plus_m
	mutation_minus[line[2]]=minus_depth-minus_m
	
	print "read data for certain site",time.time()-start
	return [mutation_result,deletion_result,insertion_result,mutation_plus,mutation_minus,len(mutation_result),plus_depth,minus_depth,line[2]] 
	### extract information on specific region if neither deletion nor insertion is information available, return [] for certain content
	### mutation_result [[base,base_p,map_p].....]
	### deletion_result [[base_p,map_p].....]
	### insertion_result
	### mutation_plus {}
	### mutation_minus {}
	### valid sequencing depth
	### valid plus sequencing depth
	### valid minus sequencing depth
	### ref

def probability_calculation_for_mutation(analyse_result,mutation):
	
	a=time.time()
	depth=analyse_result[5]
	frequency=(analyse_result[3][mutation]+analyse_result[4][mutation])/float(depth)
	if frequency==0.0:
		frequency=1/float(depth)
	
	
	lower=scipy.stats.binom.ppf(0.01,depth,frequency)
	upper=scipy.stats.binom.ppf(0.995,depth,frequency)-1
	
	lower_quantile_limit=int(depth-lower)
	upper_quantile_limit=int(upper)
	
	pre_update=[1.0]
	post_update=[]
	
	for i in range(min(lower_quantile_limit,upper_quantile_limit)):
		base=analyse_result[0][i][0]
		base_p=analyse_result[0][i][1]
		map_p=analyse_result[0][i][2]
		mutation_p=(1-map_p)*base_p+map_p/4.0
		
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
		for i in range(lower_quantile_limit,upper_quantile_limit):
			base=analyse_result[0][i][0]
			base_p=analyse_result[0][i][1]
			map_p=analyse_result[0][i][2]
			mutation_p=(1-map_p)*base_p+map_p/4.0
		
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
		for i in range(upper_quantile_limit,lower_quantile_limit):
			#print i,analyse_result[0][i],depth,lower_quantile_limit
			base=analyse_result[0][i][0]
			base_p=analyse_result[0][i][1]
			map_p=analyse_result[0][i][2]
			mutation_p=(1-map_p)*base_p+map_p/4.0
		
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
		
	for i in range(max(lower_quantile_limit,upper_quantile_limit),depth):
		base=analyse_result[0][i][0]
		base_p=analyse_result[0][i][1]
		map_p=analyse_result[0][i][2]
		mutation_p=(1-map_p)*base_p+map_p/4.0
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
				
	print 'probability calculation time',time.time()-a
	return [pre_update,int(lower),int(upper),depth] ###[[p,p,p,..],absolute_lower_bound,absolute_upper_bound] 
	
#def probability_calculation_for_mutation_previous(mutation_result,mutation):
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
	return [pre_update,pre_update_plus,pre_update_minus] ###[[p,p,p,..],absolute_lower_bound,absolute_upper_bound] 
	
def probability_calculation_for_deletion(analyse_result):
	pre_update=[1.0]
	post_update=[]
	deletion_result=analyse_result[1]
	depth=analyse_result[5]
	
	for i in range(len(deletion_result)):
		base_p=deletion_result[i][0]
		map_p=deletion_result[i][1]
		
		post_update.append(base_p*pre_update[0])
		for j in range(1,len(pre_update)):
			post_update.append(pre_update[j]*base_p+pre_update[j-1]*(1-base_p))	
		post_update.append((1-base_p)*pre_update[-1])
		
		pre_update=post_update
		
		post_update=[]
			
		
	return [pre_update,depth] ###[[p,p,p,..],sequencing_depth], when deletion is null, return [1.0]

def probability_calculation_for_insertion(analyse_result):
	insertion_result=analyse_result[2]
	depth=analyse_result[5]
	
	if insertion_result==[]:
		return [[],[],depth]
	
	max_insertion_length=max([len(v[0]) for v in insertion_result])
	pre_update_list=[]
	post_update_list=[]
	base_dic_list=[]
	for i in range(max_insertion_length):
		pre_update_list.append([1.0])
		post_update_list.append([])
		base_dic_list.append({'A':0,'T':0,'C':0,'G':0})
		
	
	for i in range(len(insertion_result)):
		
		base_string=insertion_result[i][0]
		map_p=insertion_result[i][1]
		
		for v in range(len(base_string)):
			base=base_string[v]
			post_update_list[v].append(map_p*pre_update_list[v][0])
			for j in range(1,len(pre_update_list[v])):
				post_update_list[v].append(pre_update_list[v][j]*map_p+pre_update_list[v][j-1]*(1-map_p))	
			post_update_list[v].append((1-map_p)*pre_update_list[v][-1])
			
			
			base_dic_list[v][base.upper()]+=1
			pre_update_list[v]=post_update_list[v]
			post_update_list[v]=[]
	
	return [pre_update_list,base_dic_list,depth] ###[[first_position_p_list,second_position_p_list,...],[{},{},{}....],depth]


def SNP_dicision(treatment_mutation_fraction,control_mutation_fraction,SNP_cutoff):  ##return observation by the way
	
	if treatment_mutation_fraction>=SNP_cutoff and control_mutation_fraction>=SNP_cutoff:
		return True
	else:
		return False

def strand_bias_test(analyse_result,mutation):
	treatment_plus_mutation=analyse_result[3][mutation]
	treatment_minus_mutation=analyse_result[4][mutation]
	treatment_plus_not_mutation=analyse_result[6]-treatment_plus_mutation
	treatment_minus_not_mutation=analyse_result[7]-treatment_minus_mutation
	fisher_p_value=scipy.stats.fisher_exact([[treatment_plus_not_mutation,treatment_plus_mutation],[treatment_minus_not_mutation,treatment_minus_mutation]],alternative="two-sided")[1]
	return fisher_p_value

#def somatic_mutation_test(treatment_mutation_probability,control_mutation_probability):
	a=time.time()
	treatment_p=treatment_mutation_probability[0]
	treatment_p_lower=treatment_mutation_probability[1]
	treatment_p_upper=treatment_mutation_probability[2]
	treatment_depth=treatment_mutation_probability[3]
	
	treatment_observation=weighted_avg(range(treatment_p_lower,treatment_p_upper+1), treatment_p)
	
	control_p=control_mutation_probability[0]
	control_p_lower=control_mutation_probability[1]
	control_p_upper=control_mutation_probability[2]
	control_depth=control_mutation_probability[3]
	
	p_value=0.0
	
	B_denominator=map(lambda x,y:betainc(x,y,0.5),range(control_p_lower+1,control_p_upper+2),range(control_depth-control_p_lower+1,control_depth-control_p_upper,-1))
	for extreme in range(int(treatment_observation),-1,-1):
		print extreme
		B_numerator=map(beta,range(extreme+control_p_lower+1,extreme+control_p_upper+2),range(treatment_depth+control_depth-extreme-control_p_lower+1,treatment_depth+control_depth-extreme-control_p_upper,-1))
		med=sum(map(lambda x,y,z:x*y/z,control_p,B_numerator,B_denominator))
		additional_p_value=med*int(comb(treatment_depth,extreme))
		#print additional_p_value
		if additional_p_value<=1e-10:
			p_value+=additional_p_value
			return 1-p_value
		else:
			p_value+=additional_p_value
	return 1-p_value

def somatic_mutation_test(treatment_mutation_probability,control_mutation_probability):
	a=time.time()
	treatment_p=treatment_mutation_probability[0]
	treatment_p_lower=treatment_mutation_probability[1]
	treatment_p_upper=treatment_mutation_probability[2]
	treatment_depth=treatment_mutation_probability[3]
	
	treatment_observation=weighted_avg(range(treatment_p_lower,treatment_p_upper+1), treatment_p)
	
	control_p=control_mutation_probability[0]
	control_p_lower=control_mutation_probability[1]
	control_p_upper=control_mutation_probability[2]
	control_depth=control_mutation_probability[3]
	
	control_observation=weighted_avg(range(control_p_lower,control_p_upper+1), control_p)
	
	p_value=0.0
	
	B_denominator=map(lambda x,y:beta(x,y),range(control_p_lower+1,control_p_upper+2),range(control_depth-control_p_lower+1,control_depth-control_p_upper,-1))
	for extreme in range(int(treatment_observation),-1,-1):
		#print extreme
		B_numerator=map(beta,range(extreme+control_p_lower+1,extreme+control_p_upper+2),range(treatment_depth+control_depth-extreme-control_p_lower+1,treatment_depth+control_depth-extreme-control_p_upper,-1))
		med=sum(map(lambda x,y,z:x*y/z,control_p,B_numerator,B_denominator))
		additional_p_value=med*int(comb(treatment_depth,extreme))
		
		if additional_p_value<=1e-10:
			p_value+=additional_p_value
			return (1-p_value,treatment_observation/float(treatment_depth),control_observation/float(control_depth),treatment_depth,control_depth)
		else:
			p_value+=additional_p_value
	
	return (1-p_value,treatment_observation/float(treatment_depth),control_observation/float(control_depth),treatment_depth,control_depth) ###p_value,treatment_fraction,control_fraction,treatment_depth,control_depth

def deletion_test(treatment_deletion_probability,control_deletion_probability):
	treatment_p=treatment_deletion_probability[0]
	treatment_deletion_counts=len(treatment_p)-1
	treatment_depth=treatment_deletion_probability[1]
	
	treatment_observation=weighted_avg(range(treatment_deletion_counts+1), treatment_p)
	
	control_p=control_deletion_probability[0]
	control_deletion_counts=len(control_p)-1
	control_depth=control_deletion_probability[1]
	
	p_value=0.0
	
	if control_deletion_counts==0:
		for extreme in range(int(treatment_observation),-1,-1):
			additional_p_value=int(comb(treatment_depth,extreme))*beta(extreme+1,treatment_depth+control_depth-extreme+1)/beta(1,control_depth+1)
			if additional_p_value<=1e-10:
				p_value+=additional_p_value
				return (1-p_value,treatment_observation/float(treatment_depth),0.0,treatment_depth,control_depth)
			else:
				p_value+=additional_p_value
		return (1-p_value,treatment_observation/float(treatment_depth),0.0,treatment_depth,control_depth)
	
	control_observation=weighted_avg(range(control_deletion_counts+1), control_p)

	B_denominator=map(lambda x,y:beta(x,y),range(1,control_deletion_counts+2),range(control_depth+1,control_depth-control_deletion_counts,-1))
	for extreme in range(int(treatment_observation),-1,-1):
		
		B_numerator=map(beta,range(extreme+1,extreme+control_deletion_counts+2),range(treatment_depth+control_depth-extreme+1,treatment_depth+control_depth-extreme-control_deletion_counts,-1))
		med=sum(map(lambda x,y,z:x*y/z,control_p,B_numerator,B_denominator))
		additional_p_value=med*int(comb(treatment_depth,extreme))
		#print additional_p_value
		if additional_p_value<=1e-10:
			p_value+=additional_p_value
			return (1-p_value,treatment_observation/float(treatment_depth),control_observation/float(control_depth),treatment_depth,control_depth)
		else:
			p_value+=additional_p_value
	return (1-p_value,treatment_observation/float(treatment_depth),control_observation/float(control_depth),treatment_depth,control_depth)  ##pvalue treatment_fraction control_fraction treatment_depth control_depth

def insertion_test(treatment_insertion_probability,control_insertion_probability):
	###[[first_position_p_list,second_position_p_list,...],[{},{},{}....],depth]
	
	max_treatment_insertion_length=len(treatment_insertion_probability[0])
	max_control_insertion_length=len(control_insertion_probability[0])
	
	treatment_depth=treatment_insertion_probability[2]
	control_depth=control_insertion_probability[2]
	
	p_value_list=[]
	allele_list=[]
	treatment_fraction_list=[]
	control_fraction_list=[]
	
	for position in range(max_treatment_insertion_length):
		
		allele_list.append(parser_one_base_dic(treatment_insertion_probability[1][position])[0][0])
		if position<=max_control_insertion_length-1:
		
			treatment_p=treatment_insertion_probability[0][position]
			treatment_insertion_counts=len(treatment_p)-1
			
			treatment_observation=weighted_avg(range(treatment_insertion_counts+1), treatment_p)
			
			control_p=control_insertion_probability[0][position]
			control_insertion_counts=len(control_p)-1
			
			control_observation=weighted_avg(range(control_insertion_counts+1), control_p)
			
			p_value=0.0
			
			B_denominator=map(lambda x,y:beta(x,y),range(1,control_insertion_counts+2),range(control_depth+1,control_depth-control_insertion_counts,-1))
			for extreme in range(int(treatment_observation),-1,-1):
		
				B_numerator=map(beta,range(extreme+1,extreme+control_insertion_counts+2),range(treatment_depth+control_depth-extreme+1,treatment_depth+control_depth-extreme-control_insertion_counts,-1))
				med=sum(map(lambda x,y,z:x*y/z,control_p,B_numerator,B_denominator))
				additional_p_value=med*int(comb(treatment_depth,extreme))
				#print additional_p_value
				if additional_p_value<=1e-10:
					p_value+=additional_p_value
					p_value_list.append(1-p_value)
				else:
					p_value+=additional_p_value
			p_value_list.append(1-p_value)
			treatment_fraction_list.append(treatment_observation/float(treatment_depth))
			control_fraction_list.append(control_observation/float(control_depth))
			
		else:
			
			treatment_p=treatment_insertion_probability[0][position]
			treatment_insertion_counts=len(treatment_p)-1
			
			treatment_observation=weighted_avg(range(treatment_insertion_counts+1), treatment_p)
			control_observation=0.0
			
			p_value=0.0
			
			for extreme in range(int(treatment_observation),-1,-1):
				additional_p_value=int(comb(treatment_depth,extreme))*beta(extreme+1,treatment_depth+control_depth-extreme+1)/beta(1,control_depth+1)
				if additional_p_value<=1e-10:
					p_value+=additional_p_value
					p_value_list.append(1-p_value)
				else:
					p_value+=additional_p_value
				
			p_value_list.append(p_value)
			treatment_fraction_list.append(treatment_observation/float(treatment_depth))
			control_fraction_list.append(control_observation/float(control_depth))
	
	
	
	return [p_value_list,allele_list,treatment_fraction_list,control_fraction_list]  #####[[pvalue...],[allele...],[treatment_fraction...],[control_fraction...]]

def read_mpileup(file):
	file=open(file)
	a=[]
	while True:
		line=file.readline()
		if not line:break
		line=line.strip()
		a.append(line)
		
	return a 
def test():
	a=time.time()
	t="""1	10017	N	1500	CCCCCCcCC+2CTCC*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC+2CTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcCCCCCCCCCCCCCCCCCCCCCCCCCCCCC+2CTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCCCCC+2CTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCccccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcccccccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCccccc^!C^!C^!C^5C^!C^"C^]C^!C^!C^8C^!C^!C^!C^!C^!C^OC^!C^!C^!C^!C^!C^!C^"C^!C^!C^!C^!C^!C^!C^"C^!C^*C^!C^8C^!C^*C^!C^OC^!C^!C^8C^!C^!C^FC^!C^!C^2C^!C^!C^!C^!C^%C^!C^!C^!C^!C^!C^>C^#C^!C^!C^!C^*C^!C^!C^!C^!C^!C^!c^!c^5c^Fc	91#9DE7#5988##99999999999999999999999999999:999999999999999999999999999?9??>?????????DEEEEDEEDEEEEEEEEEEEFEEFFEFE@EEEEE?>??999:9999999999999999999999999999999999999999999999999999999999999???????EDDDD@DDEDEEDEEEEEEEFEEE=E?9999999999999999999999999999999899999:9999999999>;??#;?>?DDDDDDDDEED@DDCDDEEEE>EDDDEEEDDDDDDEEEEEE999999999999:999999999999999999999:9999999999999999999999:9999::9>>?>??>??DDEEEEDDDBE;EEEDDD;EDDDDDDDDDDDDEEEDEEDC?>7::9::99:9:9::9999::9::99:99:99::9::9:::999::99:9::::99:::9::99::::::99:::::9:9:::9:::9>>>>>>?>@>??>9DDDDDDDDDDDDDDDDDDDDEEEEEEEEEDDDDEEEEDFDCDDDDDEDDDFDD>::::::::::::::::9:9:99:9::9::::999:::::::::::::::::::9:::::9::9::>1?>??>>DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD>DDDDCDDDDEDDDDDED@DDDD?::::::::::::::::::::::::::9:::::::>?>?><DCCDDDDDDDDDDDDDDDDDDDDEDDDDBD:::::::::6::::::::::::::::::::':::::::::??9DDDDCCDDDDDDDDD@DDDDDDDDDDDDC:;::::;:;;;;;;;:::;:::;;:;;:;;;;;9;:;::;;::;:;;;;:;:6:;=?:==>CD<CDDDDCD2DCDCCDCDDDA=DDDDDDDD>?455;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;<=<<C>DCACCDCDDDDA;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;<;;;;=6?;<<>DCCCDDBBBB:DDDDDBDBC<><;;<;;<;;;;<;;<;<;;<<<;;;<<<<;<;<<;<<>;<>>;DD>D?ADADD>ADAADDBD;>3A@=<=<<====<=<====<<<<<<=<<<<<=====<<<=<=<=<<====<=<?;<<?;;@D@AA@@D@@@@D@@@@D;>@=============================<=======<==?<?BCCCC=CC?CCABCBCC;;00<<>>>>>>>>>=>=>>>>>>>>>>>>><>>==>>>>><><>>>>>>?>>>>>CBC6C<CCBB?CC5ACCBCC?CC>//50BEAA?BB<BB>BBBBBBBBAB>BB?B@@@BBBBBBBBBB=<==BC@?CBCBD<1/9D7@@:@9>9@@99@999999999@9999@999999999@=9A;A;?;?;C=@C@@@@CCC@AC>AA7C?B#.##	!!+!!!+!!!*!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"!!!!!!!!!!!!!"!!!!!!!!!-$!-*!$!-0!>!!!!!"!>!0!%!!!>!!!!!!!%!!!!!!!!>!%!!**-!!!5!!!!!!!!!"!!!!!!!!!!!!!!!!!!!!$"!""!!!!!!!"!!!!"!!!!!!!!!!!!!!!!>*+!!!%!-!!$>*!>!!!!!!!!*!!!!!!!!!!!!!!!!$*8!!!!!!/!!!!!!!!!!!!$!!!!!!!!!!*!!!!+2!2!!!*!$>!%!!*!!!!$%$$*3!!>$!!,$$!5$!>$!!!$!*$!!!$$!!!!!*!*/.!!!!!'*!*!!$!!**!!!!!!!!/5!!!:**!!!*$!&$!!*!$$!!*/!">!>5!!!!2!!!!*$$!!!*$$$$!$!$!$!!!*$!$3$$0*2!!0!0!0-*!3/!'/!!!!!!!*!*!*!/!$!!$*!$!'!!/!!!!!!!/"$!$!$*!$!!!!!/!!!"!!!&!*!//"!$!!!!!!!!!!>!$/!$!!""!!+*!!&!*3>!!!*$%3!!$!*$*!!!$(!!!!$$!!8!!$!$!!*!*$****7!$>!1!!/!/./*!!!>!!!!>!!!!$!*!*$)%>!!!!/!*!$$!!*!*!*!*>*//*!/!!!*$>!!!**!!>!20!-$$*!$$$3!$!$2>>$1!!$!$!$4$>*!**!!>!2!$>>*>*>!!!$10!*!$$1$5!!"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!$!!!"!!>!!$!!!!!!!!!!!!!!>!!!!!%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%!!!!!%!!!!!!!!!""%!!-!!!!!!1!!!!!!!!!!F!5!!!!!!!!!!!!!!!!!!"!!!!!!!!!!!!!!!!$!!#!!!!%!!!!(!!!!!!!!!!!!!!!!!!%!!!+!!!!!!!!!1!!!!!*.!!!!"*'!1!!!!!!!!!#!!!!!!!!!!%!!!!!!!!!!!!!!!!!!"!!!L!!!!!!!!!!!!!!!!!!!!!"!!!!!!$!!!!!!!!!!!%!2!!!!!F!>!%%F!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!L!!!!--!!!!!!!!!!!!!!%"!!!4!!FF!!!!!!!!!,!!F!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"!!!!!!!2$!!!2"*!!!!!!!!!!!T!!!!!!!!!!!!!!!!!!!0!!!!!!!!!!!!!!!!!!!!!!!!!!!!>!2%!!!!!!!!!!!!!!!!!!8FF!!!!!!8!!!!!!!!!!!3!!!!!!O!!!!!!!!!!!!!!!"!!!!!#>!!!!!!!!!F!#!!!!!!!!!!!!!!88!83!F!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!$!!!!!0!>!$!58FF!!!!5!"]!!8!!!!!O!!!!!!"!!!!!!"!*!8!*!O!!8!!F!!2!!!!%!!!!!>#!!!*!!!!!!!5F"""
	c="""1	10006	N	721	CCCCCCcCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC+1TCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-1NCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC^/C^!C^/C^.C^/C^*C^!C^!C^!C^>C^!C^!C^!C^!C^>C^!C^!C^!C^!C^$C^!C^*C^!C^*C^$C^)C^%C^>C^!C^!C^!C^!C^/C^!C^*C^!C^$C^$C^!C^!C^*C^!C^*C^!C^*C^!C^*C^>C^*C^/C^/C^*C^!C^/C^!C^!C^!C^*C^$C^>C^!C^!C^!C^*C^*C^!C^!C^>C^!C^2C^0C^!C^-C^$C^$C^*C^!C^$C^$C^$C^3C^!C^$C^!C^$C^2C^>C^>C^$C^1C^!C^!C^$C^!C^$C^!C^$C^4C^$C^>C^*C^!C^*C^*C^!C^!C^>C^!C^2C^!C^$C^>C^>C^*C^>C^*C^>C^!C^!C^!C^$C^1C^0C^!C^*C^!C^$C^$C^1C^$C	::#;CD6#8<;2##<<;<<;<<<<<<<<<<<<<<<;<<<;;;;<<<;<<<<<;<;;<<<;;<;<;<<;<;<?>><>>:>;?;;<>ADDDDADAADDDDDDDDDDDABDABDAD:BCDDD>=<?<<<<<<<=<==<====<=<=<====<<<===<=<=<<==<===<<<<===<<<==<=<=<=<==<??;=??;D@@AAA@@D@D@CDCCDADD@DDDAD;=================================================?@=6?===BCBB=BCBCCA+AACAACCD244A?=7CCABB8CCCCCAAA>>>>>>>>>>>>=<>>>>>>>>>>>>>>>>>>>>=>>>>>>>>>>>>>>>=>=>>>>=>>>>==>>B?>??>>;BACCCCCCBBC>CCCBC;@C>C?@CC@CCCCC??CCBDCC?>0BBBBBBBABABBBBBBBBABBABBBBBBBBBBBBBBBBBBBBABBBBBBBBABBBBBBBBBBBBBBABBBBBBB?BBBBBBBBBBB==><====>=>>=CCCBCCCCCBCCCCCCCCBCBCCCCCCCCCCCBBCCCC@DBAAACCCBAAAD9C=9@@@9?@9@99@@9=@9@9@99@9@@9@@?@9999=A;A9A999@AA99AAA?9AA@AA9AA9>AA;A;?A:;CCCCCCCCC@@@@@@CC@CCCCCC@C@C@CCCCACCCCAACCAA7CBBBAC=CCCCA	!!+!!!+!!!*!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"!!!!!!!!!!!!!"!!!!!!!!!-$!-*!$!-0!>!!!!!"!>!0!%!!!>!!!!!!!%!!!!!!!!>!%!!**-!!!5!!!!!!!!!"!!!!!!!!!!!!!!!!!!!!$"!""!!!!!!!"!!!!"!!!!!!!!!!!!!!!!>*+!!!%!-!!$>*!>!!!!!!!!*!!!!!!!!!!!!!!!!$*8!!!!!!/!!!!!!!!!!!!$!!!!!!!!!!*!!!!+2!2!!!*!$>!%!!*!!!!$%$$*3!!>$!!,$$!5$!>$!!!$!*$!!!$$!!!!!*!*/.!!!!!'*!*!!$!!**!!!!!!!!/5!!!:**!!!*$!&$!!*!$$!!*/!">!>5!!!!2!!!!*$$!!!*$$$$!$!$!$!!!*$!$3$$0*2!!0!0!0-*!3/!'/!!!!!!!*!*!*!/!$!!$*!$!'!!/!!!!!!!/"$!$!$*!$!!!!!/!!!"!!!&!*!//"!$!!!!!!!!!!>!$/!$!!""!!+*!!&!*3>!!!*$%3!!$!*$*!!!$(!!!!$$!!8!!$!$!!*!*$****7!$>!1!!/!/./*!!!>!!!!>!!!!$!*!*$)%>!!!!/!*!$$!!*!*!*!*>*//*!/!!!*$>!!!**!!>!20!-$$*!$$$3!$!$2>>$1!!$!$!$4$>*!**!!>!2!$>>*>*>!!!$10!*!$$1$"""
	t=t.strip().split()
	c=c.strip().split()
	treatment=analyze_one_position(t,5,0,30)
	control=analyze_one_position(c,5,0,30)
	print treatment[1],treatment[2]
	treatment_mutation_p=probability_calculation_for_mutation(treatment,'A')
	control_mutation_p=probability_calculation_for_mutation(control,'A')
	print mutation_test(treatment_mutation_p,control_mutation_p)
	print "strand bias",strand_bias_test(treatment,'A')
	print time.time()-a,"L"
		
def additional_string(chrom,position,ref_allele,tumor_sample,tumor_depth,control_sample,control_depth,SNV_kind,main_allele,mutation_allele,fraction_in_tumor,fraction_in_normal,p_value,strand_p_value,fdr,tag):
	a="%s\t"*15+"%s\n"
	if strand_p_value!="*":
		return a%(chrom,str(position),ref_allele,tumor_sample,str(tumor_depth),control_sample,str(control_depth),SNV_kind,main_allele,mutation_allele,str(round(fraction_in_tumor,3)),str(round(fraction_in_normal,3)),str(round(p_value,3)),str(round(strand_p_value,3)),fdr,tag)
	else:
		return a%(chrom,str(position),ref_allele,tumor_sample,str(tumor_depth),control_sample,str(control_depth),SNV_kind,main_allele,mutation_allele,str(round(fraction_in_tumor,3)),str(round(fraction_in_normal,3)),str(round(p_value,3)),"*",fdr,tag)


def rapid_mutation_test(treatment_information,control_information,minimum_allele_fraction,p_value_cutoff):
	allele_list=[]
	for allele in ["A","T","C","G"]:
		#print treatment_information[3],treatment_information[8]
		if allele!=treatment_information[8] and (treatment_information[3][allele]+treatment_information[4][allele])!=0:
			
			
			treatment_mutation=treatment_information[3][allele]+treatment_information[4][allele]
			treatment_not_mutation=treatment_information[5]-treatment_mutation
			if treatment_mutation/float(treatment_information[5])>minimum_allele_fraction:
				pass
			else:
				continue
			
			p_value=0.0
			control_mutation=control_information[3][allele]+control_information[4][allele]
			control_not_mutation=control_information[5]-control_mutation
			
			denominator=beta(control_mutation+1,control_not_mutation+1)
			for extreme in range(treatment_mutation,-1,-1):
				#print treatment_information[5],extreme
				numerator=beta(extreme+control_mutation+1,treatment_not_mutation+control_not_mutation+1)
				additional_p_value=comb(treatment_information[5],extreme)*numerator/float(denominator)
				if additional_p_value<=1e-10:
					p_value+=additional_p_value
					break
				else:
					p_value+=additional_p_value
			
			if p_value<p_value_cutoff:
				allele_list.append(allele)
			else:
				pass
		
	return allele_list


def handle_certain_segment_on_specific_chrom(chrom,treatment_chrom,control_chrom,treatment_bam,control_bam,start,tumor_name,control_name,base_quality_cutoff,map_quality_cutoff,sequencing_depth_cutoff,SNP_cutoff,significant_cutoff,reading_length):
	
	output_string=""
	#treatment_list=generate_mpileup_for_each_chrom_on_specific_region(treatment_bam,treatment_chrom,start,reading_length).strip().split("\n")
	#control_list=generate_mpileup_for_each_chrom_on_specific_region(control_bam,control_chrom,start,reading_length).strip().split("\n")
	
	treatment_list=read_mpileup("/Users/chenfei/Desktop/test_down_sample1.txt")
	control_list=read_mpileup("/Users/chenfei/Desktop/test_down_sample2.txt")
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
				##sequencing depth cutoff
				treatment_information=analyze_one_position(t_line,base_quality_cutoff,map_quality_cutoff,sequencing_depth_cutoff)
				
				if not treatment_information:
					print "%s sample on %s:%s is not deep enough."%(tumor_name,chrom,str(t_position))
					break
				
				control_information=analyze_one_position(c_line,base_quality_cutoff,map_quality_cutoff,sequencing_depth_cutoff)
				if not control_information:
					print "%s sample on %s:%s is not deep enough."%(control_name,chrom,str(t_position))
					break
				
				
				#######mutation
				allele_list=rapid_mutation_test(treatment_information,control_information,0.002,0.1)
				treatment_depth=treatment_information[5]
				control_depth=control_information[5]
				for allele in allele_list:
					treatment_mutation_probability=probability_calculation_for_mutation(treatment_information,allele[0])
					control_mutation_probability=probability_calculation_for_mutation(control_information,allele[0])
					[p_value,treatment_fraction,control_fraction,treatment_depth,control_depth]=somatic_mutation_test(treatment_mutation_probability,control_mutation_probability)
					if  p_value<significant_cutoff:
						strand_p_value=strand_bias_test(treatment_information,allele[0])  ###roughly estimate strand bias
						output_string+=additional_string(chrom,c_position,"*",tumor_name,treatment_depth,control_name,control_depth,"somatic",allele_ranking[0][0],allele[0],treatment_fraction,control_fraction,p_value,strand_p_value,"*",'1/2')
						
				#check_main_allele=False
				#allele_ranking=parser_two_base_dic(treatment_information[3],treatment_information[4])
				#control_allele_ranking=parser_two_base_dic(control_information[3],control_information[4])
				#print allele_ranking
				#print parser_two_base_dic(control_information[3],control_information[4])
				#for allele in allele_ranking[1:]:
				#	treatment_mutation_probability=probability_calculation_for_mutation(treatment_information,allele[0])
				#	control_mutation_probability=probability_calculation_for_mutation(control_information,allele[0])
				#	[p_value,treatment_fraction,control_fraction,treatment_depth,control_depth]=somatic_mutation_test(treatment_mutation_probability,control_mutation_probability)
				#	if SNP_dicision(treatment_fraction,control_fraction,SNP_cutoff):
				#		output_string+=additional_string(chrom,c_position,"*",tumor_name,treatment_depth,control_name,control_depth,"SNP",allele_ranking[0][0],allele[0],treatment_fraction,control_fraction,p_value,"*","*",'3/4')
				#	elif p_value<significant_cutoff:
				#		strand_p_value=strand_bias_test(treatment_information,allele[0])  ###roughly estimate strand bias
				#		output_string+=additional_string(chrom,c_position,"*",tumor_name,treatment_depth,control_name,control_depth,"somatic",allele_ranking[0][0],allele[0],treatment_fraction,control_fraction,p_value,strand_p_value,"*",'1/2')
				#		
				#	else:
				#		check_main_allele=True
				#if check_main_allele:
				#	treatment_mutation_probability=probability_calculation_for_mutation(treatment_information,allele_ranking[0][0])
				#	control_mutation_probability=probability_calculation_for_mutation(control_information,allele_ranking[0][0])
				#	[p_value,treatment_fraction,control_fraction,treatment_depth,control_depth]=somatic_mutation_test(treatment_mutation_probability,control_mutation_probability)
				#	if p_value<significant_cutoff:
				#		strand_p_value=strand_bias_test(treatment_information,allele_ranking[0][0])
				#		output_string+=additional_string(chrom,c_position,"*",tumor_name,treatment_depth,control_name,control_depth,"somatic",allele_ranking[0][0],allele_ranking[0][0],treatment_fraction,control_fraction,p_value,strand_p_value,"*",'6')
						
					
				####deletion
				
				treatment_deletion_probability=probability_calculation_for_deletion(treatment_information)
				if treatment_deletion_probability[0]==[1.0]:
					print "no deletion on %s:%s from %s sample"%(chrom,str(t_position),tumor_name)
				else:
					control_deletion_probability=probability_calculation_for_deletion(control_information)
					#print treatment_deletion_probability,control_deletion_probability
					(p_value,treatment_fraction,control_fraction,treatment_depth,control_depth)=deletion_test(treatment_deletion_probability,control_deletion_probability)
					if p_value<significant_cutoff:
						output_string+=additional_string(chrom,c_position,"*",tumor_name,treatment_depth,control_name,control_depth,"del","*","-",treatment_fraction,control_fraction,p_value,"*","*",'*')
				
				
				###insertion
				treatment_insertion_probability=probability_calculation_for_insertion(treatment_information)
				
				if treatment_insertion_probability[0]==[]:
					print "no insertion on %s:%s from %s sample"%(chrom,str(t_position),tumor_name)
					control_pivot=j+1
					break
				else:
					control_insertion_probability=probability_calculation_for_insertion(control_information)
					print treatment_insertion_probability,control_insertion_probability
					(p_value_list,allele_list,treatment_fraction_list,control_fraction_list)=insertion_test(treatment_insertion_probability,control_insertion_probability)
					for v in range(len(p_value_list)):
						p_value=p_value_list[v]
						if p_value>=significant_cutoff:
							break
						allele=allele_list[v]
						treatment_fraction=treatment_fraction_list[v]
						control_fraction=control_fraction_list[v]
						output_string+=additional_string(chrom,c_position,"*",tumor_name,treatment_depth,control_name,control_depth,"ins","*",allele,treatment_fraction,control_fraction,p_value,"*","*",'*')
	
					control_pivot=j+1
					break
			
			else:
				

def read_bed(bed_path):
	file=open(bed_path)
	a=[]
	while True:
		line=file.readline()
		if not line:break
		line=line.strip().split()
		a.append([line[0],int(line[1]),int(line[2])]) ###chrom,start,end
	return a
	
def read_postion(position_path):
	file=open(position_path)
	a=[]
	while True:
		line=file.readline()
		if not line:break
		line=line.strip().split()
		a.append([line[0],int(line[1]),int(line[2])]) ###chrom,position,position
	return a

def mm9_parser_bam_header(bam_file):
	chrom_list=[str(i) for i in range(1,20)]+["chr"+str(i) for i in range(1,23)]+["X","chrX","Y","chrY"]
	bam_chrom_dic={} ###{chromosome:length}
	cmd="samtools view -H %s"%bam_file
	header=sp(cmd)[0].strip().split("\n")
	for value in header:
		value=value.split()
		if value[0]=="@SQ":
			if value[1].split(":")[1] in chrom_list:
				bam_chrom_dic[value[1].split(":")[1]]=int(value[2].split(":")[1])			
		else:
			break
	return bam_chrom_dic

def hg19_parser_bam_header(bam_file):
	chrom_list=[str(i) for i in range(1,23)]+["chr"+str(i) for i in range(1,23)]+["X","chrX","Y","chrY"]
	bam_chrom_dic={} ###{chromosome:length}
	cmd="samtools view -H %s"%bam_file
	header=sp(cmd)[0].strip().split("\n")
	for value in header:
		value=value.split()
		if value[0]=="@SQ":
			if value[1].split(":")[1] in chrom_list:
				bam_chrom_dic[value[1].split(":")[1]]=int(value[2].split(":")[1])			
		else:
			break
	return bam_chrom_dic

def chromsome_information(treatment_bam,control_bam,species):
	if options.species=="hg":
		chrom_complete_list=[("chr"+str(i),str(i)) for i in range(1,23)]+[("chrX","X"),("chrY","Y")]
		chrom_list=[]
		print_chrom_list=[]
		
		treatment_chrom_dic=hg19_parser_bam_header(treatment_bam)
		control_chrom_dic=hg19_parser_bam_header(control_bam)
		for value in chrom_complete_list:
			if treatment_chrom_dic.has_key(value[0]):
				if control_chrom_dic.has_key(value[0]):
					chrom_list.append((value[0],value[0],treatment_chrom_dic[value[0]]))
					print_chrom_list.append(value[0])
				elif control_chrom_dic.has_key(value[1]):
					chrom_list.append((value[0],value[1],treatment_chrom_dic[value[0]]))
					print_chrom_list.append(value[0])
				else:
					pass
			elif treatment_chrom_dic.has_key(value[1]):
				if control_chrom_dic.has_key(value[0]):
					chrom_list.append((value[1],value[0],treatment_chrom_dic[value[1]]))
					print_chrom_list.append(value[0])
				elif control_chrom_dic.has_key(value[1]):
					chrom_list.append((value[1],value[1],treatment_chrom_dic[value[1]]))
					print_chrom_list.append(value[0])
				else:
					pass
			else:
				pass
		a="%s "*len(print_chrom_list)
		print "Alignment on "+a%print_chrom_list+"are availble for Both tumor sample and control sample~~"
		return (print_chrom_list,chrom_list) ####([chr1,chr2....],[(1,chr1,length).....])
					
	elif options.species=="mm":
		chrom_complete_list=[("chr"+str(i),str(i)) for i in range(1,20)]+[("chrX","X"),("chrY","Y")]
		chrom_list=[]
		print_chrom_list=[]
		treatment_chrom_dic=mm9_parser_bam_header(treatment_bam)
		control_chrom_dic=mm9_parser_bam_header(control_bam)
		for value in chrom_complete_list:
			if treatment_chrom_dic.has_key(value[0]):
				if control_chrom_dic.has_key(value[0]):
					chrom_list.append((value[0],value[0],treatment_chrom_dic[value[0]]))
					print_chrom_list.append(value[0])
				elif control_chrom_dic.has_key(value[1]):
					chrom_list.append((value[0],value[1],treatment_chrom_dic[value[0]]))
					print_chrom_list.append(value[0])
				else:
					pass
			elif treatment_chrom_dic.has_key(value[1]):
				if control_chrom_dic.has_key(value[0]):
					chrom_list.append((value[1],value[0],treatment_chrom_dic[value[1]]))
					print_chrom_list.append(value[0])
				elif control_chrom_dic.has_key(value[1]):
					chrom_list.append((value[1],value[1],treatment_chrom_dic[value[1]]))
					print_chrom_list.append(value[0])
				else:
					pass
			else:
				pass
		a="%s "*len(print_chrom_list)
		print "Alignment on "+a%print_chrom_list+"are availble for Both tumor sample and control sample~~"
		return (print_chrom_list,chrom_list) ####([chr1,chr2....],[(1,chr1,length).....])
	else:
		print "please refer to assign your interested chromsome."

def handle_large_interval(start,end,output_file,chrom,treatment_chrom,control_chrom,treatment_bam,control_bam,tumor_name,control_name,base_quality_cutoff,map_quality_cutoff,sequencing_depth_cutoff,SNP_cutoff,significant_cutoff,reading_length):
	while True:
		if start+reading_length<end:
			start_time=time.time()
			output_string=handle_certain_segment_on_specific_chrom(chrom,treatment_chrom,control_chrom,treatment_bam,control_bam,start,tumor_name,control_name,base_quality_cutoff,map_quality_cutoff,sequencing_depth_cutoff,SNP_cutoff,significant_cutoff,reading_length)
			output_file.writelines(output_string)
			end_time=time.time()
			start+=reading_length
			print "handling position %s to %s on %s consumes %s minutes."%(str(start),str(start+reading_length),chrom,str((end_time-start_time)/60))
		else:
			start_time=time.time()
			output_string=handle_certain_segment_on_specific_chrom(chrom,treatment_chrom,control_chrom,treatment_bam,control_bam,start,tumor_name,control_name,base_quality_cutoff,map_quality_cutoff,sequencing_depth_cutoff,SNP_cutoff,significant_cutoff,reading_length,end-start)
			output_file.writelines(output_string)
			end_time=time.time()
			print "handling position %s to %s on %s consumes %s minutes."%(str(start),str(end),chrom,str((end_time-start_time)/60))
			break


	
	
def handle_one_position(chrom,treatment_chrom,control_chrom,treatment_bam,control_bam,position,tumor_name="tumor",control_name="control",base_quality_cutoff=10,map_quality_cutoff=0,sequencing_depth_cutoff=20,SNP_cutoff=0.4,significant_cutoff=0.01):
	
	output_string=""
	treatment=generate_mpileup_for_each_chrom_on_specific_region(treatment_bam,treatment_chrom,start,0).strip().split("\n").split()
	control=generate_mpileup_for_each_chrom_on_specific_region(control_bam,control_chrom,start,0).strip().split("\n").split()
	position=int(treatment[1])
	
	treatment_information=analyze_one_position(treatment,base_quality_cutoff,map_quality_cutoff,sequencing_depth_cutoff)
	if treatment_information:
		print "%s sample on %s:%s is not deep enough."%(tumor_name,chrom,str(t_position))
		break
				
	control_information=analyze_one_position(control,base_quality_cutoff,map_quality_cutoff,sequencing_depth_cutoff)
	if control_information:
		print "%s sample on %s:%s is not deep enough."%(control_name,chrom,str(t_position))
		break
	
	#######mutation
	check_main_allele=False
	allele_ranking=parser_two_base_dic(treatment_information[3],treatment_information[4])
	for allele in allele_ranking[1:]:
		treatment_mutation_probability=probability_calculation_for_mutation(treatment_information,allele)
		control_mutation_probability=probability_calculation_for_mutation(control_information,allele)
		(p_value,treatment_fraction,control_fraction,treatment_depth,control_depth)=somatic_mutation_test(treatment_mutation_probability,control_mutation_probability,allele)
		if SNP_dicision(treatment_mutation_fraction,control_mutation_fraction,SNP_cutoff):
			output_string+=additional_string(chrom,position,"*",tumor_name,treatment_depth,control_name,control_depth,"SNP",allele_ranking[0],allele,treatment_fraction,control_fraction,p_value,"*","*",'3/4')
		elif p_value<significant_cutoff:
			strand_p_value=strand_bias_test(treatment_information,allele)
			output_string+=additional_string(chrom,position,"*",tumor_name,treatment_depth,control_name,control_depth,"somatic",allele_ranking[0],allele,treatment_fraction,control_fraction,p_value,strand_p_value,"*",'1/2')
		else:
			check_main_allele=True
		
		if check_main_allele:
			(p_value,treatment_fraction,control_fraction,treatment_depth,control_depth)=somatic_mutation_test(treatment_mutation_probability,control_mutation_probability,allele_ranking[0])
			if p_value<significant_cutoff:
				strand_p_value=strand_bias_test(treatment_information,allele_ranking[0])
				output_string+=additional_string(chrom,position,"*",tumor_name,treatment_depth,control_name,control_depth,"somatic",allele_ranking[0],allele_ranking[0],treatment_fraction,control_fraction,p_value,strand_p_value,"*",'6')
						
			
	####deletion
				
	treatment_deletion_probability=probability_calculation_for_deletion(treatment_information)
	if treatment_deletion_probability[0]==[1.0]:
		print "no deletion on %s:%s from %s sample"%(chrom,str(position),tumor_name)
	else:
		pass
				
	control_deletion_probability=probability_calculation_for_deletion(control_information)
				
	(p_value,treatment_fraction,control_fraction,treatment_depth,control_depth)=deletion_test(treatment_deletion_probability,control_deletion_probability)
	if p_value<significant_cutoff:
		output_string+=additional_string(chrom,position,"*",tumor_name,treatment_depth,control_name,control_depth,"del",allele_ranking[0],"-",treatment_fraction,control_fraction,p_value,"*","*",'*')
				
				
	##insertion
	treatment_insertion_probability=probability_calculation_for_insertion(treatment_information)
				
	if treatment_insertion_p[0]==[]:
		print "no insertion on %s:%s"%(chrom,str(position))
	else:
		pass
					
	control_insertion_probability=probability_calculation_for_insertion(control_information)
	(p_value_list,allele_list,treatment_fraction_list,control_fraction_list]=insertion_test(treatment_insertion_probability,control_insertion_probability)
	for v in range(len(p_value_list):
		p_value=p_value_list[v]
		if p_value>=significant_cutoff:
			break
		allele=allele_list[v]
		treatment_fraction=treatment_fraction_list[v]
		control_fraction=control_fraction_list[v]
		output_string+=additional_string(chrom,position,"*",tumor_name,treatment_depth,control_name,control_depth,"ins",allele_ranking[0],allele,treatment_fraction,control_fraction,p_value,"*","*",'*')
	
				
	return output_string

def main():
	chrom_information=chromsome_information(options.treatment_bam,options.control_bam,species)
	output_string="chrom\tposition\tref_allele\ttumor_sample\ttumor_depth\tcontrol_sample\tcontrol_depth\tSNV_kind\tmain_allele\talter_allele\tfraction_in_tumor\tfraction_in_normal\tp_value\tstrand_p_value\tfdr\ttag\n"
	output_file=open(options.output_file,'w')
	
	if ***:
		for i in range(len(chrom_information[0])):
			chrom=chrom_information[0][i]
			treatment_chrom=chrom_information[1][i][0]
			control_chrom=chrom_information[1][i][1]
			chrom_length=chrom_information[1][i][2]
		
			start=0
			end=chrom_length
			handle_large_interval(start,end,output_file,chrom,treatment_chrom,control_chrom,options.treatment_bam,options.control_bam,options.tumor_name,options.control_name,options.base_quality_cutoff,options.map_quality_cutoff,options.sequencing_depth_cutoff,options.SNP_cutoff,options.significant_cutoff,options.reading_length):
	elif ***:
		for interval in interval_list:
			if interval[0] in chrom_information[0]:
				chrom=interval[0]
				index=chrom_information[0].index(chrom)
				treatment_chrom=chrom_information[1][index][0]
				control_chrom=chrom_information[1][index][1]
				start=interval[1]
				end=interval[2]
				handle_large_interval(start,end,output_file,chrom,treatment_chrom,control_chrom,options.treatment_bam,options.control_bam,options.tumor_name,options.control_name,options.base_quality_cutoff,options.map_quality_cutoff,options.sequencing_depth_cutoff,options.SNP_cutoff,options.significant_cutoff,options.reading_length):
			else:
				print "can not find %s alignment on %s bam file!"%(chrom,tumor_name)
	else:
		for position in position_list:
			if position[0] in chrom_information[0]:
				chrom=position[0]
				index=chrom_information[0].index(chrom)
				treatment_chrom=chrom_information[1][index][0]
				control_chrom=chrom_information[1][index][1]
				string=handle_one_position(chrom,treatment_chrom,control_chrom,treatment_bam,control_bam,position,tumor_name,control_name=,base_quality_cutoff,map_quality_cutoff,sequencing_depth_cutoff,SNP_cutoff,significant_cutoff)
				output_file.writelines(string)
			else:
				print "Can not find %s alignment on %s bam file!"%(chrom,tumor_name)

	output_file.close()
				print "no available control on %s:%s"%(chrom,str(t_position))
				control_pivot=j
				break
			
				
	return output_string
	
string=handle_certain_segment_on_specific_chrom("chr1",'1','1',2,3,10027,"tumor","control",10,0,30,0.4,1,10)
output=open("/Users/chenfei/Desktop/test_result.txt",'w')
output.writelines(string)
output.close()


	
		