#!/usr/bin/env python

__author__ = "Mohammad Hadigol"
__version__ = "1.2"
__date__ = "Date: 01-2018"

import argparse
import time
import shlex
import sys 
import re
import subprocess
import os
import itertools
from operator import itemgetter 
import math
from distutils import spawn
import csv
import numpy as np

frw_charsm = ['.', 'A', 'T', 'C', 'G', 'N']	
rev_charsm = [',', 'a', 't', 'c', 'g', 'n']	

def main():
	"""Main block"""
	print "#################################################"
	print "################# Start MERIT ###################"
	print "#################################################"
	
	### Get arguments 
	print args

	### Make required directories
	res_dir = args.RES_dir + args.Name + "/" 	### main res directory
	if not os.path.exists(res_dir):
		os.makedirs(res_dir)
	
	for i in args.steps:
		### get the class from the classes in the current script
		myclass = globals()['Step' + i]
		
		### instantiate a step obj for the class
		mystep = myclass(args, i)
		
		### execute the run method
		mystep.run()
	
	#### delete temporary intermediate files	
	if not args.noclean: 
	    	if (args.ann): 
		    
	    		numann = len(args.ann.split(",")) 
			intermediate_files = [".var", ".var.qual", ".mp.vcf.gz", ".ann." + str(numann) + ".vcf.gz", ".mp.vcf.gz", ".mpileup"]
			
		elif not "2" in args.steps:
		
			intermediate_files = [".var", ".var.qual", ".mp.vcf.gz", ".mp.vcf.gz", ".mpileup"]
			
		else:
			    
			intermediate_files = [".var", ".var.qual", ".mp.vcf.gz", ".ann.0.vcf.gz", ".mp.vcf.gz", ".mpileup"]
			
		for j in intermediate_files:
			junkfile =  args.RES_dir + "/" + args.Name + "/" + args.Name + j
			if os.path.exists(junkfile):
				cmd_rm = "rm {}".format(junkfile)
				os.system(cmd_rm)
				
###---------------------------------------------------------------------	

def get_arg():
	"""Get Arguments"""
	
	MERIT_description = """This is MERIT - a Mutation Error Rate Identification Toolkit - that includes the following steps:  
	
	 - Step 1: Run samtools mpileup to detect variants and generate Pileup and BCF/VCF
	 - Step 2: Optional annotation via snpEff and SnpSift
	 - Step 3: Convert VCf to human readable file via bcftools query and generate SampleName.var
	 - Step 4: Extract Phred scores and more from Pileup and store in SampleName.var.qual
	 - Step 5: Adding context to variants and generating SampleName.var.ctx
	 - Step 6: Estimate the error rates
	 
	 For a complete description of how to use this software, including dependencies and usage examples, see www.software.khiabanian-lab.org
	"""
	
	### MERIT directory  			
	MERITdir = os.path.dirname(os.path.realpath(__file__))
	
	# Current working directory 
	cwdir = os.getcwd()
	
	parser = argparse.ArgumentParser(description = MERIT_description, prog='MERIT_v1.py', usage='python %(prog)s [OPTIONS] ...')

	parser.add_argument("-B", "--BAM_dir", default = MERITdir + '/BAM/' , help="Directory of input BAM file")
	parser.add_argument("-R", "--RES_dir", default = MERITdir + '/RES/' , help="Directory of results")
	parser.add_argument("-I", "--positions", help="BED file containing a list of regions or sites where pileup or BCF should be generated")
	parser.add_argument("-H", "--REF", default = MERITdir + '/bin/ref_genome/hg19.fa', help="Reference genome")
	parser.add_argument("-N", "--Name", help="Sample name")
	parser.add_argument("--max_hp_no", type=int, default = 4, help="Maximum number of homopolymer repeats to be considered in indel error rate analysis")
	parser.add_argument("--cut_off_freq", default = 10.0, help='Cut off frequency (percentage) to identify errors')
	parser.add_argument("--min_depth", default = 100, help="Minimum total depth of variant to be considered in error rate anlysis")
	parser.add_argument("--qual-offset", default = 33, help="Quality offset") 
	parser.add_argument("--ann", help="Comma-separated list of annotating VCFs with which to provide additional annotation in SnpSift. Absolute path should be provided.")
	parser.add_argument("--steps", default = '123456', help="steps to run")
	parser.add_argument("--memory", default= "4", help="the memory for the (SnpEff) Java virtual machine in gigabytes")
	parser.add_argument("-d", "--max_depth", default = 100000000, help="Max depth option -d for samtools mpileup")
	parser.add_argument("-q", "--min_MQ", default = 0, help="Minimum mapping quality for an alignment to be used")
	parser.add_argument("-Q", "--min_BQ", default = 0, help="Minimum base quality for a base to be used")
	parser.add_argument("-e", "--ext_prob", default = 20, help="Phred-scaled gap extension sequencing error probability")
	parser.add_argument("-F", "--gap_frac", default = 0.000001, help="Minimum fraction of gapped reads")
	parser.add_argument("--tandem_qual", default = 100, help="Coefficient for modeling homopolymer errors")
	parser.add_argument("-L", "--max_idepth", type=int, default = 100000000, help="Skip INDEL calling if the average per-input-file depth is above INT")
	parser.add_argument("-o", "--open_prob", type=int, default = 40, help="Phred-scaled gap open sequencing error probability")
	parser.add_argument("--verbose","-v", action="store_true", help="Print commands used in each step (default: off)")
	parser.add_argument("--two_nt_indel", action="store_true", help="Consider two nt indels in error rate analysis (default: off)")
	parser.add_argument("--noclean", action="store_true", help="do not delete temporary intermediate files (default: off)")

	args = parser.parse_args()
	
	return args

###---------------------------------------------------------------------		
def terminal_size():
    ### Source: https://stackoverflow.com/questions/566746/how-to-get-linux-console-window-width-in-python
    import fcntl, termios, struct
    hh, ww, hhp, wwp = struct.unpack('HHHH',
        fcntl.ioctl(0, termios.TIOCGWINSZ,
        struct.pack('HHHH', 0, 0, 0, 0)))
    return ww, hh


###---------------------------------------------------------------------	

def check_file_exists_and_nonzero(myfile):
	"""Check for the existence and nonzero-ness of a file"""

	# loop through comma-delimited list of files
	for i in myfile.split(","):
		if (os.path.isfile(os.path.expanduser(i))):
			if (os.path.getsize(os.path.expanduser(i)) == 0):
				print(i + " is empty. Exiting")
				sys.exit(1)
		else:
			print("Can't find " + i + ". Exiting.")
			sys.exit(1)
	
###---------------------------------------------------------------------	

def mytimer(myfunc):
	"""Decorator for timing a Step's run function"""
	# http://stackoverflow.com/questions/5478351/python-time-measure-function

	def mynewfunc(*args, **kwargs):
		startTime = time.time()
		myfunc(*args, **kwargs)
		print('Step execution time: {} sec'.format(int(time.time() - startTime)))
		print "#################################################"

	return mynewfunc
	
###---------------------------------------------------------------------	

def getFrwRevTotalDepth(items):
	'''get total depth on forward and reverse strands'''
	frw_chars = ['.', 'A', 'T', 'C', 'G', 'N', 'F']	
	rev_chars = [',', 'a', 't', 'c', 'g', 'n', 'R']	
	
	t_d_frw = 0
	t_d_rev = 0
	
	for i in range(len(frw_chars)):
		t_d_frw += items.count(frw_chars[i])
		t_d_rev += items.count(rev_chars[i])
		
	return t_d_frw, t_d_rev
		
###---------------------------------------------------------------------	

def getPrefix( s1, s2 ):
    '''get common prefix of strings s1 and s2.'''
    n = min( len(s1), len(s2))
    for x in range( n ):
	if s1[x] != s2[x]: return s1[:x]
    return s1[:n]
    
###---------------------------------------------------------------------	

def getSuffix( s1, s2 ):
    '''get common sufix of strings s1 and s2.'''
    n = min( len( s1), len( s2 ) )
    if s1[-1] != s2[-1]: return ""
    for x in range( -2, -n - 1, -1 ):
	if s1[x] != s2[x]: return s1[x+1:]
    return s1[-n:]
    
###---------------------------------------------------------------------	

def allIsMinusOne(items):
    '''return True if all the elements in items are equal -1'''  
    return all(x == -1 for x in items)
    
###---------------------------------------------------------------------	

def allTheSame(items):
    '''return True if all the elements in items are similar'''  
    return all(x == items[0] for x in items)
    
###---------------------------------------------------------------------	

def getIndex(items, itm):    
    '''get indices of itm in items'''  
    return [nn for nn in xrange(len(items)) if items.find(itm, nn) == nn]

###---------------------------------------------------------------------	
    
def removeCaretDollar(items):
    '''remove $, ^ and character after ^ in items'''    
    itemstmp = re.subn('\^.', '', items)       ## remove ^ and char after that (this should be done before removing $ - special case error: depth of 1x with mapping quality of $, example: chr1	1219526	g	1	^$A	%	1)
    items = itemstmp[0]
    items = items.replace('$', '')             ## remove $
    return items
    
###---------------------------------------------------------------------	

def getMode(items):
    '''find the mode in items'''    
    return np.argmax(np.bincount(np.asarray(items)))

###---------------------------------------------------------------------	

def removeINDELnt(items):
    '''remove +[0-9]+[ACGTNacgtn] and -[0-9]+[ACGTNacgtn] in items'''			
    items_out = items
	
    ind_plus = getIndex(items, '+')                
    ind_minus =  getIndex(items, '-')      
	
    ind_p_int = map(lambda x:x+1, ind_plus)
    ind_m_int = map(lambda x:x+1, ind_minus)
	
    ind_p_int2 = map(lambda x:x+2, ind_plus)
    ind_m_int2 = map(lambda x:x+2, ind_minus)
	
    vec_p = list(range(len(ind_plus)))
    vec_m = list(range(len(ind_minus)))
	
    if ind_plus:
		for ii in range(len(ind_plus)):
			vec_p[ii] = itemgetter(ind_p_int[ii])(items)
			vec_p2 = itemgetter(ind_p_int2[ii])(items)
			if vec_p2.isdigit():
				vec_p[ii] = str(vec_p[ii]) + vec_p2
    else:
		vec_p = []		
			
    if ind_minus:
		for ii in range(len(ind_minus)):
			vec_m[ii] = itemgetter(ind_m_int[ii])(items)
			vec_m2 = itemgetter(ind_m_int2[ii])(items)
			if vec_m2.isdigit():
				vec_m[ii] = str(vec_m[ii]) + vec_m2
    else:
		vec_m = []		
									
    ind_rm_p = range(len(vec_p))
    ind_rm_m = range(len(vec_m))

    if vec_p:
		for ii in range(len(vec_p)):
			if int(vec_p[ii]) < 10:
				ind_rm_p[ii] = int(ind_p_int[ii]) + 1 + int(vec_p[ii])
				items_out = items_out.replace(str(items[int(ind_plus[ii]):ind_rm_p[ii]]), '') 
			elif (int(vec_p[ii]) >= 10 and int(vec_p[ii]) < 100):
				ind_rm_p[ii] = int(ind_p_int[ii]) + 2 + int(vec_p[ii])
				items_out = items_out.replace(str(items[int(ind_plus[ii]):ind_rm_p[ii]]), '') 
											
    if vec_m:
		for ii in range(len(vec_m)):
			if int(vec_m[ii]) < 10:
				ind_rm_m[ii] = int(ind_m_int[ii]) + 1 + int(vec_m[ii])
				items_out = items_out.replace(str(items[int(ind_minus[ii]):ind_rm_m[ii]]), '') 
			elif (int(vec_m[ii]) >= 10 and int(vec_m[ii]) < 100):
				ind_rm_m[ii] = int(ind_m_int[ii]) + 2 + int(vec_m[ii])	
				items_out = items_out.replace(str(items[int(ind_minus[ii]):ind_rm_m[ii]]), '')

    return items_out
    
###---------------------------------------------------------------------	

def getQualAndPos(ind_frw, ind_rev, ph, pos):
    '''given ind_frw and ind_rev, find the corresponding qualities and'''
    '''positions in ph and pos'''
    if ind_frw:
	qual_alt_frw = itemgetter(*ind_frw)(ph)
	qual_frw = [ord(i) - args.qual_offset for i in ''.join(qual_alt_frw)] 
	mean_qual_frw = np.mean(qual_frw)
	pos_f = itemgetter(*ind_frw)(pos)
	if (len(ind_frw) == 1):
	    pos_f = (pos_f,)	
	mode_pos_frw = getMode(map(int, list(pos_f)))
    else:
	mean_qual_frw = 0; mode_pos_frw = 0; pos_f = 0
    
    if ind_rev:
	qual_alt_rev = itemgetter(*ind_rev)(ph)
	qual_rev = [ord(i) - args.qual_offset for i in ''.join(qual_alt_rev)]
	mean_qual_rev = np.mean(qual_rev)
	pos_r = itemgetter(*ind_rev)(pos)
	if (len(ind_rev) == 1):
	    pos_r = (pos_r,)
	mode_pos_rev = getMode(map(int, list(pos_r)))
    else:
	mean_qual_rev = 0; mode_pos_rev = 0; pos_r = 0
	
    return mean_qual_frw, mean_qual_rev, mode_pos_frw, mode_pos_rev
    
###---------------------------------------------------------------------	

def markBases(items, del_frw, del_rev):
    '''mark bases before del_bases_frw and del_bases_rev with F for '''
    '''forward and R for reverse'''
    items_tmp = items
    ind_del_frw = getIndex(items, del_frw)
    ind_p_int = map(lambda x:x-1, ind_del_frw)
    
    if ind_p_int:
	for ii in range(len(ind_p_int)):
	    items_tmp = "".join((items_tmp[:int(ind_p_int[ii])], 'F' ,items_tmp[int(ind_p_int[ii])+1:]))
	    
    ind_del_rev = getIndex(items, del_rev)
    ind_m_int = map(lambda x:x-1, ind_del_rev)
	    
    if ind_m_int:
	for ii in range(len(ind_m_int)):
	    items_tmp = "".join((items_tmp[:int(ind_m_int[ii])], 'R' ,items_tmp[int(ind_m_int[ii])+1:]))
	    
    return items_tmp
    
###---------------------------------------------------------------------	

def getErrorRate(items):
    
    my_count = np.zeros([1,10]) 
    my_phat = np.zeros([1,9]) 
	    
    nzfrw = np.nonzero(items[:,2])
    nzrev = np.nonzero(items[:,3])
	
    my_count[0,0] = items[:,0].sum()
    my_count[0,1] = items[:,1].sum()
    my_count[0,2] = items[:,2].sum()
    my_count[0,3] = items[:,3].sum()
    my_count[0,4] = items[:,5].sum()
    my_count[0,5] = items[:,6].sum()
    
    if (np.size(nzfrw) != 0):
	my_count[0,6] = getMode((items[nzfrw[0],7]).astype(int))
    else:
	my_count[0,6] = 0
	
    if (np.size(nzrev) != 0):
	my_count[0,7] = getMode((items[nzrev[0],8]).astype(int))
    else:
	my_count[0,7] = 0
	
    my_count[0,8] = np.size(nzfrw)
    my_count[0,9] = np.size(nzrev)
    
    if (my_count[0,8] != 0):
    	my_count[0,4] = my_count[0,4]/my_count[0,8]
    if (my_count[0,9] != 0):
    	my_count[0,5] = my_count[0,5]/my_count[0,9]
    
    ### Error rate
    if (my_count[0,0] != 0):
    	my_phat[0,0] = my_count[0,2]/my_count[0,0]  ### frw
    if (my_count[0,1] != 0):	
    	my_phat[0,1] = my_count[0,3]/my_count[0,1]  ### rev
    my_phat[0,2] = (my_count[0,2] + my_count[0,3])/(my_count[0,0] + my_count[0,1]) ### frw + rev	
    my_phat[0,3:] = my_count[0,4:]	
    
    return my_count, my_phat

###---------------------------------------------------------------------	

def PileupExtractor(line, dic_mpileup):
    
	itr = 0;
	###start_time = time.time()
	cols = line.split("\t")
	var = ''.join(cols[4:5])
	ref = ''.join(cols[2:3])
	mpileup_line = dic_mpileup['-'.join((cols[0], cols[1]))]
	bases_org = mpileup_line[4]
	bases = removeCaretDollar(bases_org)    			### bases do not include Caret and Dollar 
	phreds = list(map(str, mpileup_line[5]))
	pos_in_read = (mpileup_line[6].rstrip()).split(",")
	###print 'Init completed in ',(time.time() - start_time), ' sec'
	indel_nt = '-'
	
	if itr % 100 == 0 and itr != 0:
		print "Adding quality ...", itr, "variants completed."
	
	itr += 1
	
	len_var = len(var)
	len_ref = len(ref)
	len_phreds = len(phreds)
	len_bases= len(bases)

	################ SNV with no indel at next position
	if (len_var == len_ref and len_phreds == len_bases): 

	    ind_alt_frw = getIndex(bases, var)
	    ind_alt_rev = getIndex(bases, var.lower())
	    
	    [mean_qual_frw, mean_qual_rev , mode_pos_frw, mode_pos_rev] = \
	    	getQualAndPos(ind_alt_frw, ind_alt_rev, phreds, pos_in_read)
	
	    ### Get total and allele depth on frw and rev strands 
	    d_var_frw = len(ind_alt_frw)
	    d_var_rev = len(ind_alt_rev)	
	    
	    [d_tot_frw, d_tot_rev] = getFrwRevTotalDepth(bases)
	
	################ SNV with indel at next position
	elif (len_var == len_ref and len_phreds != len_bases): 
	    
	    bases_tmp = removeINDELnt(bases)

	    ind_alt_frw = getIndex(bases_tmp, var)

	    ind_alt_rev = getIndex(bases_tmp, var.lower())

	    [mean_qual_frw, mean_qual_rev , mode_pos_frw, mode_pos_rev] = \
	    	getQualAndPos(ind_alt_frw, ind_alt_rev, phreds, pos_in_read)
	    	
	    ### Get total and allele depth on frw and rev strands 
	    d_var_frw = len(ind_alt_frw)
	    d_var_rev = len(ind_alt_rev)	
	    
	    [d_tot_frw, d_tot_rev] = getFrwRevTotalDepth(bases_tmp)
	    
	################ Deletion
	elif (len_var < len_ref): 
	    
	    bases_tmp  = bases
	    
	    ### Get deleted bases (del_nt)
	    altsuffix = getSuffix(ref, var)
	    del_nt = ref[1:len(ref) - len(altsuffix)]
	    indel_nt = del_nt
	    del_bases_frw = "-" + str(len(del_nt)) + del_nt
	    del_bases_rev = "-" + str(len(del_nt)) + del_nt.lower()
	    			
	    ### Get the position after deleted bases 
	    pos_n = int(''.join(cols[1:2])) + len(del_nt) + 1
	    next_pos_avail = pos_n in dic_mpileup.keys()
	    if next_pos_avail:
		mpileup_line_n = dic_mpileup[str(pos_n)]
		bases_org_n = mpileup_line_n[4]
		bases_n = removeCaretDollar(bases_org_n)    
		phreds_n = list(map(str, mpileup_line_n[5]))
		pos_in_read_n = (mpileup_line_n[6].rstrip()).split(",")	    
		
		### Inspect the positions with deleted bases (*)
		has_dollar = range(len(del_nt))
		
		dollar_count = 0
		caret_count_frw = 0
		caret_count_rev = 0
		dollar_count_frw = 0
		dollar_count_rev = 0
		
		for ilp in range(len(del_nt)):
		    mpileup_line_ilp = dic_mpileup[str(int(''.join(cols[1:2])) + ilp + 1)]
		    bases_ilp = mpileup_line_ilp[4]
		    has_dollar[ilp] = bases_ilp.find('$')
		    dollar_count += bases_ilp.count('$')
		    
		    if ilp == 0:
			ind_caret = getIndex(bases_ilp, '^')
			if ind_caret:
			    ind_first_base = map(lambda x:x+2, ind_caret)
			    first_bases =  itemgetter(*ind_first_base)(bases_ilp)
			    for i in range(len(frw_charsm)):
				caret_count_frw += first_bases.count(frw_charsm[i])
				caret_count_rev += first_bases.count(rev_charsm[i])
				
		    if ilp < len(del_nt):
			    ind_dollar = getIndex(bases_ilp, '$')
			    if ind_dollar:
				ind_last_base = map(lambda x:x-1, ind_dollar)
				last_bases =  itemgetter(*ind_last_base)(bases_ilp)
				for i in range(len(frw_charsm)):
				    dollar_count_frw += last_bases.count(frw_charsm[i])
				    dollar_count_rev += last_bases.count(rev_charsm[i])
	    else:
		
		dollar_count = 0
		caret_count_frw = 0
		caret_count_rev = 0
		dollar_count_frw = 0
		dollar_count_rev = 0
		has_dollar = range(len(del_nt))
		for ilp in range(len(del_nt)):
			has_dollar[ilp] = -1	    
				    
						
	    ### Mark bases before minus sign 
	    bases_marked = markBases(bases, del_bases_frw, del_bases_rev)
	    bases_marked_tmp = removeINDELnt(bases_marked)
	    
	    ### Get the index of bases before deleted bases
	    ind_del_frw_b = getIndex(bases_marked_tmp, 'F')
	    ind_del_rev_b = getIndex(bases_marked_tmp, 'R')
	
	    ### Get index of dollar signs in bases_org
	    ind_dollar = getIndex(bases_org, '$')
	    
	    ### Get the index of bases after deleted bases
	    ### Correct ind_del_frw_a and ind_del_rev_a for $ in bases
	    if ind_dollar:
			ind_del_frw_a = map(lambda x:x-len(ind_dollar), ind_del_frw_b)
			ind_del_rev_a = map(lambda x:x-len(ind_dollar), ind_del_rev_b)
	    else:
			ind_del_frw_a = ind_del_frw_b
			ind_del_rev_a = ind_del_rev_b
	    
	    ### Correct ind_del_frw_a and ind_del_rev_a for $ in positions with *
	    if not allIsMinusOne(has_dollar):
		ind_del_frw_a = map(lambda x:x-dollar_count, ind_del_frw_a)
		ind_del_rev_a = map(lambda x:x-dollar_count, ind_del_rev_a)
		
	    [mean_qual_frw_b, mean_qual_rev_b, mode_pos_frw_b, mode_pos_rev_b] = \
	    	getQualAndPos(ind_del_frw_b, ind_del_rev_b, phreds, pos_in_read)
		
	    if next_pos_avail:
		
		[mean_qual_frw_a, mean_qual_rev_a, mode_pos_frw_a, mode_pos_rev_a] = \
		    getQualAndPos(ind_del_frw_a, ind_del_rev_a, phreds_n, pos_in_read_n)
		    
		mean_qual_frw = (mean_qual_frw_b + mean_qual_frw_a)/2.0
		mean_qual_rev = (mean_qual_rev_b + mean_qual_rev_a)/2.0
		
		mode_pos_frw = (mode_pos_frw_b + mode_pos_frw_a)/2.0
		mode_pos_rev = (mode_pos_rev_b + mode_pos_rev_a)/2.0
		
	    else:
		
		mean_qual_frw = mean_qual_frw_b
		mean_qual_rev = mean_qual_rev_b 
		
		mode_pos_frw = mode_pos_frw_b
		mode_pos_rev = mode_pos_rev_b 	
	    
	    ### Get total and allele depth on frw and rev strands 
	    d_var_frw = len(ind_del_frw_b)
	    d_var_rev = len(ind_del_rev_b)	
	    
	    [d_tot_frw, d_tot_rev] = getFrwRevTotalDepth(bases_marked_tmp)
	    
	    ### Correct total depth
	    if ind_dollar:
		ind_last_base = map(lambda x:x-1, ind_dollar)
		last_bases =  itemgetter(*ind_last_base)(bases_org)
		for i in range(len(frw_charsm)):
			d_tot_frw -= last_bases.count(frw_charsm[i])
		 	d_tot_rev -= last_bases.count(rev_charsm[i])	
			
	    d_tot_frw -= dollar_count_frw
	    d_tot_rev -= dollar_count_rev		 

	    d_tot_frw += caret_count_frw
	    d_tot_rev += caret_count_rev
		
	################ Insertion
	elif (len_var > len_ref): 

	    bases_tmp  = bases
	    
	    ### Get deleted bases (del_nt)
	    altsuffix = getSuffix(var, ref)
	    ins_nt = var[1:len(var) - len(altsuffix)]
	    indel_nt = ins_nt
	    ins_bases_frw = "+" + str(len(ins_nt)) + ins_nt
	    ins_bases_rev = "+" + str(len(ins_nt)) + ins_nt.lower()
	    	    
	    ### Get the position after current one
	    pos_n = int(''.join(cols[1:2])) + 1
	    next_pos_avail = pos_n in dic_mpileup.keys()
	    if next_pos_avail:
		mpileup_line_n = dic_mpileup[str(pos_n)]
		bases_org_n = mpileup_line_n[4]
		bases_n = removeCaretDollar(bases_org_n)    
		phreds_n = list(map(str, mpileup_line_n[5]))
		pos_in_read_n = (mpileup_line_n[6].rstrip()).split(",")	
	    
	    ### Mark bases before plus sign 
	    bases_marked = markBases(bases, ins_bases_frw, ins_bases_rev)
	    bases_marked_tmp = removeINDELnt( bases_marked )
	    
	    ### Get the index of bases before inserted bases
	    ind_ins_frw_b = getIndex(bases_marked_tmp, 'F')
	    ind_ins_rev_b = getIndex(bases_marked_tmp, 'R')

	    ### Get index of dollar signs in bases_org
	    ind_dollar = getIndex(bases_org, '$')

	    ### Get the index of bases after inserted bases
	    ### Correct ind_ins_frw_a and ind_ins_rev_a for $ in bases
	    if ind_dollar:
			ind_ins_frw_a = map(lambda x:x-len(ind_dollar), ind_ins_frw_b)
			ind_ins_rev_a = map(lambda x:x-len(ind_dollar), ind_ins_rev_b)
	    else:
			ind_ins_frw_a = ind_ins_frw_b
			ind_ins_rev_a = ind_ins_rev_b
		
	    [mean_qual_frw_b, mean_qual_rev_b, mode_pos_frw_b, mode_pos_rev_b] = \
		getQualAndPos(ind_ins_frw_b, ind_ins_rev_b, phreds, pos_in_read)		
	    
	    if next_pos_avail:
		
		[mean_qual_frw_a, mean_qual_rev_a, mode_pos_frw_a, mode_pos_rev_a] = \
		    getQualAndPos(ind_del_frw_a, ind_del_rev_a, phreds_n, pos_in_read_n)
		    
		mean_qual_frw = (mean_qual_frw_b + mean_qual_frw_a)/2.0
		mean_qual_rev = (mean_qual_rev_b + mean_qual_rev_a)/2.0
		
		mode_pos_frw = (mode_pos_frw_b + mode_pos_frw_a)/2.0
		mode_pos_rev = (mode_pos_rev_b + mode_pos_rev_a)/2.0
		
	    else:
		
		mean_qual_frw = mean_qual_frw_b
		mean_qual_rev = mean_qual_rev_b 
		
		mode_pos_frw = mode_pos_frw_b
		mode_pos_rev = mode_pos_rev_b 	
	    
	    ### Get total and allele depth on frw and rev strands 
	    d_var_frw = len(ind_ins_frw_b)
	    d_var_rev = len(ind_ins_rev_b)	
	    
	    [d_tot_frw, d_tot_rev] = getFrwRevTotalDepth(bases_marked_tmp)
	    
	    ### Correct total depth	    
	    if ind_dollar:
		ind_last_base = map(lambda x:x-1, ind_dollar)
		last_bases =  itemgetter(*ind_last_base)(bases_org)
		for i in range(len(frw_charsm)):
			d_tot_frw -= last_bases.count(frw_charsm[i])
		 	d_tot_rev -= last_bases.count(rev_charsm[i])	

	######################
	newrow = cols[0:6] + (str(d_tot_frw)).split() + (str(d_tot_rev)).split() + (str(d_var_frw)).split() + (str(d_var_rev)).split() \
		+ (str((float(d_var_frw + d_var_rev)/float(d_tot_frw + d_tot_rev)*100))).split() + (str(mean_qual_frw)).split() + (str(mean_qual_rev)).split() \
		+ (str(mode_pos_frw)).split() + (str(mode_pos_rev)).split() +  (indel_nt).split() + cols[15:]
		
	
	return newrow

###---------------------------------------------------------------------	

def run_cmd(cmd, bool_verbose, bool_getstdout, stepi):
	"""Run a shell command"""

	# if verbose, print command
	if (bool_verbose):
	    	print(stepi)
		print("Step command: " + cmd)

	proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

	(stdout, stderr) =  proc.communicate()

	# if error, print it
	if (stderr and bool_verbose):
		print("[Step stdout] " + stderr),

	# return stdout
	if (bool_getstdout): 
		return stdout.rstrip()
	else:
		return "0" # note: this must return a str   
		
###---------------------------------------------------------------------	

class Steps(object):
	"""A parent step class from which step children inherit"""
	
	def __init__ (self, args, step_index):
		"""Initialize step object"""
		self.args = args 
		self.step_index = step_index 
		
	def set_descrip(self, mydescrip):
		"""Set the description"""
		self.description = mydescrip 
		
	def set_input(self, myfile):
		"""Set the input file for the step"""
		self.input = myfile

		### check if input file has nonzero size
		check_file_exists_and_nonzero(self.input)

	def set_output1(self, myfile):
		"""Set the output file 1 for the step"""
		self.output1 = myfile
		
	def set_output2(self, myfile):
		"""Set the output file 2 for the step"""
		self.output2 = myfile

	def set_output3(self, myfile):
		"""Set the output file 3 for the step"""
		self.output3 = myfile
		
	def run(self):
		"""The run method, meant to be overridden by the children"""
		print self.description		
		
###---------------------------------------------------------------------	
		
class Step1(Steps):
	"""Run samtools mpileup to detect variants and generate Pileup and BCF/VCF files"""

	def __init__ (self, args, step_index):
		"""Initialize step object"""

		### Parent's constructor
		Steps.__init__(self, args, step_index)
		
		### Define description as well as input and output attributes
		self.set_descrip("###### Step 1 - Run samtools mpileup to detect variants and generate Pileup and BCF/VCF. ######")
		self.set_input(args.BAM_dir + args.Name + ".sorted.bam")
		self.set_output1(self.args.RES_dir + "/" + self.args.Name + "/" + self.args.Name + ".mpileup")         ### Pileup output
		self.set_output2(self.args.RES_dir + "/" + self.args.Name + "/" +  self.args.Name + ".mp.bcf")         ### Multiallelic BCF output
		self.set_output3(self.args.RES_dir + "/" + self.args.Name + "/" +  self.args.Name + ".mp.vcf.gz")      ### Biallelic VCF output
		
		
	@mytimer	
	def run(self):	
		"""Run samtools mpileup"""

		### Run parent's run method
		super(Step1, self).run()
				
		### Set the flag for samtools mpileup to generate Pileup
		if (self.args.positions):   
			self.mpileupflag = "-q {} -Q {} -d {} -L {} -A -Bx -F {} --tandem-qual {} --ext-prob {} --open-prob {} --output-BP --output-tags DP,AD,ADF,ADR --positions {} -f {} {} > {}".format(self.args.min_MQ, self.args.min_BQ, self.args.max_depth, self.args.max_depth, self.args.gap_frac, self.args.tandem_qual, self.args.ext_prob, self.args.open_prob, self.args.positions, self.args.REF, self.input, self.output1)
					
		else:
			self.mpileupflag = "-q {} -Q {} -d {} -L {} -A -Bx -F {} --tandem-qual {} --ext-prob {} --open-prob {} --output-BP --output-tags DP,AD,ADF,ADR -f {} {} > {}".format(self.args.min_MQ, self.args.min_BQ, self.args.max_depth, self.args.max_depth, self.args.gap_frac, self.args.tandem_qual, self.args.ext_prob, self.args.open_prob, self.args.REF, self.input, self.output1)
			
		### Set the flag for samtools mpileup and bcftools to generate BCF/VCF
		self.mpileupflag_bcf = self.mpileupflag.replace('-f', '-gf')
		self.mpileupflag_bcf = self.mpileupflag_bcf.replace(self.output1, self.output2)
			
		### Set the flag for bcftools norm
		self.bcftoolsflag = "--do-not-normalize -m - -Oz -o {} {}".format(self.output3, self.output2)   	 #### ATT: this is correct flag, use it later
		###self.bcftoolsflag = "--do-not-normalize -m - -Oz -o {} > {}".format(self.output3, self.output3) 	 #### ATT: this is only for testing

	    	### Generating Pileup
		#time.sleep(1)
		self.cmd_mpileup = "samtools mpileup " + self.mpileupflag
		run_cmd(self.cmd_mpileup, self.args.verbose, 0, "### Step 1 - generating Pileup via samtools mpileup ###")
		
		### Generating BCF
		self.cmd_mpileup_bcf = "samtools mpileup " + self.mpileupflag_bcf
		run_cmd(self.cmd_mpileup_bcf, self.args.verbose, 0, "### Step 1 - generating BCF via samtools mpileup ###")
		
		### Split multiallelic sites in BCF file into biallelic records via bcftools norm
		self.cmd_bcftools = "bcftools norm " + self.bcftoolsflag
		run_cmd(self.cmd_bcftools, self.args.verbose, 0, "### Step 1 - converting multiallelic into biallelic records ###")		
				
		### Check if output files have nonzero size
		check_file_exists_and_nonzero(self.output1)  ### Pileup
		check_file_exists_and_nonzero(self.output2)  ### Multiallelic BCF
		check_file_exists_and_nonzero(self.output3)  ### Biallelic VCF
		
		### Delet multiallelic BCF
		os.remove(self.output2)		

###---------------------------------------------------------------------	

class Step2(Steps):
	"""Run snpEff and SnpSift for annotatation"""

	def __init__ (self, args, step_index):
		"""Initialize step object"""
		
		### Parent's constructor
		Steps.__init__(self, args, step_index)
		
		### Define description as well as input and output attributes
		self.set_descrip("###### Step 2 - Optinal annotation via snpEff and SnpSift. ######")
		self.set_input(self.args.RES_dir + "/" + self.args.Name + "/" + self.args.Name + ".mp.vcf.gz")             ### Biallelic VCF from step 1
		self.set_output1(self.args.RES_dir + "/" + self.args.Name + "/" + self.args.Name + ".ann.0.vcf.gz")        ### SnpEff output
		
	@mytimer
	def run(self):	
		"""Run snpEff and SnpSift"""
		
		### Run parent's run method
		super(Step2, self).run()	   
		 	
		### get paths
		self.effpath = spawn.find_executable('snpEff.jar')
		self.siftpath = spawn.find_executable('SnpSift.jar')           
		self.configpath = spawn.find_executable('snpEff.config')
				
		### set the flag for snpEff
		self.effopts = "-noLog -noHgvs -q -formatEff -noStats -lof -canon"
		
		### Run npEff 
		self.cmd_snpeff = "java -Xmx" + self.args.memory + "g -jar " + self.effpath + " ann hg19 " + self.effopts + " -c " + self.configpath + " " + self.input + " | gzip > " + self.output1
		run_cmd(self.cmd_snpeff, self.args.verbose, 1, "### Step 2 - annotation via snpEff ###")

		## Additional annotation with SnpSift
		if (self.args.ann):
			for j,i in enumerate(self.args.ann.split(",")):
				#### Set input and output VCF files
				self.in_vcf_snpsift = self.args.RES_dir + "/" + self.args.Name + "/" + self.args.Name + ".ann." + str(j) + ".vcf.gz"
				self.out_vcf_snpsift = self.args.RES_dir + "/" + self.args.Name + "/" + self.args.Name + ".ann." + str(j+1) + ".vcf.gz"
								
				#### Run SnpSift
				self.cmd_snpsift = "java -Xmx" + self.args.memory + "g -jar " + self.siftpath + " annotate " + i + " -noLog " +  self.in_vcf_snpsift + " | gzip > " + self.out_vcf_snpsift
				run_cmd(self.cmd_snpsift, self.args.verbose, 1, "### Step 2 - annotation via SnpSift: " + i + " ###")
				
				#### Delete intermediate VCF files
				os.remove(self.in_vcf_snpsift)
				
		### Check if output files have nonzero size
		if (self.args.ann):
		    	numann = len(self.args.ann.split(",")) 
			check_file_exists_and_nonzero(self.args.RES_dir + "/" + self.args.Name + "/" + self.args.Name + ".ann." + str(numann) + ".vcf.gz")   ### Annotated VCF via snpEff
		else:
		    	check_file_exists_and_nonzero(self.output1)	        ### SnpEff annotated VCF
				
###---------------------------------------------------------------------			
		
class Step3(Steps):
	"""Convert VCf to human readable file via bcftools query"""

	def __init__ (self, args, step_index):
		"""Initialize step object"""
		
		### Parent's constructor
		Steps.__init__(self, args, step_index)
		
		### Define description as well as input and output attributes
		self.set_descrip("###### Step 3 - Convert VCf to human readable file via bcftools query and generate .var ######")
		self.set_output1(self.args.RES_dir + "/" + self.args.Name + "/" + self.args.Name + ".var")                					### var output
				
		if ("2" in self.args.steps):
			if (self.args.ann):
					numann = len(self.args.ann.split(",")) 
					self.set_input(self.args.RES_dir + "/" + self.args.Name + "/" + self.args.Name + ".ann." + str(numann) + ".vcf.gz")             ### Annotated VCF via snpEff
			else:
				self.set_input(self.args.RES_dir + "/" + self.args.Name + "/" + self.args.Name + ".ann.0.vcf.gz")                               ### Annotated VCF via SnpSift			    
		else:
			self.set_input(self.args.RES_dir + "/" + self.args.Name + "/" +  self.args.Name + ".mp.vcf.gz")                                         ### Biallelic VCF (not annotated)
		
	@mytimer
	def run(self):	
		"""Run bcftools query"""
		
		### Run parent's run method
		super(Step3, self).run()
		
		### Write header for var		
		if ("2" in self.args.steps):
		    		if (self.args.ann):
				    	self.bcf_list = "'%CHROM %POS %REF %ALT %DP [%AD{0}] [%AD{1}] [%ADF{0}] [%ADR{0}] [%ADF{1}] [%ADR{1}] %QUAL %SNP %CDS %EFF\n'" 	
					self.header_var = ["#chr", "pos", "ref", "ref_ctx", "alt", "alt_ctx", "total_depth", "ref_depth", \
						"alt_depth_frw", "alt_depth_rev", "freq", "ave_alt_qual_frw", "ave_alt_qual_rev", \
						"pos_inread_alt_frw", "pos_inread_alt_rev", "snp", "CDS", "effect"]  
				else:	
				    	self.bcf_list = "'%CHROM %POS %REF %ALT %DP [%AD{0}] [%AD{1}] [%ADF{0}] [%ADR{0}] [%ADF{1}] [%ADR{1}] %QUAL %EFF %EFF %EFF\n'" 
					self.header_var = ["#chr", "pos", "ref", "ref_ctx", "alt", "alt_ctx", "total_depth", "ref_depth", \
						"alt_depth_frw", "alt_depth_rev", "freq", "ave_alt_qual_frw", "ave_alt_qual_rev", \
						"pos_inread_alt_frw", "pos_inread_alt_rev", "effect"]					
					
		else:
		    		self.bcf_list = "'%CHROM %POS %REF %ALT %DP [%AD{0}] [%AD{1}] [%ADF{0}] [%ADR{0}] [%ADF{1}] [%ADR{1}] %QUAL\n'"
				self.header_var = ["#chr", "pos", "ref", "ref_ctx", "alt", "alt_ctx", "total_depth", "ref_depth", \
					"alt_depth_frw", "alt_depth_rev", "freq", "ave_alt_qual_frw", "ave_alt_qual_rev", \
					"pos_inread_alt_frw", "pos_inread_alt_rev"]

		with open(self.output1, 'w') as fvar:
    			csv.writer(fvar, delimiter='\t').writerow(self.header_var)
							
		### Set the flag for bcftools query
		if ("2" in self.args.steps):
		    		if (self.args.ann):	
				    self.cmd_bcftools = "bcftools query -f {in1} {in2} | awk '{{if($4 !~ \"*\") \
					print $1   \"\t\"   $2   \"\t\"   $3   \"\t\"   \"-\"   \"\t\"   $4   \"\t\"   \"-\"   \"\t\" \
					$5   \"\t\"   $6   \"\t\"   $10   \"\t\"  $11   \"\t\"   \"-\"   \"\t\"   \"-\"   \"\t\"   \"-\"   \"\t\" \
					\"-\"  \"\t\"   \"-\"   \"\t\"   $13  \"\t\"  $14 \"\t\" $15}}' >> {in3}".format(in1=self.bcf_list, in2=self.input, in3=self.output1)  
				else:	
				    self.cmd_bcftools = "bcftools query -f {in1} {in2} | awk '{{if($4 !~ \"*\") \
					print $1   \"\t\"   $2   \"\t\"   $3   \"\t\"   \"-\"   \"\t\"   $4   \"\t\"   \"-\"   \"\t\" \
					$5   \"\t\"   $6   \"\t\"   $10   \"\t\"  $11   \"\t\"   \"-\"   \"\t\"   \"-\"   \"\t\"   \"-\"   \"\t\" \
					\"-\"  \"\t\"   \"-\"   \"\t\"  $15}}' >> {in3}".format(in1=self.bcf_list, in2=self.input, in3=self.output1)  					
					
		else:
				self.cmd_bcftools = "bcftools query -f {in1} {in2} | awk '{{if($4 !~ \"*\") \
					print $1   \"\t\"   $2   \"\t\"   $3   \"\t\"   \"-\"   \"\t\"   $4   \"\t\"   \"-\"   \"\t\" \
					$5   \"\t\"   $6   \"\t\"   $10   \"\t\"  $11   \"\t\"   \"-\"   \"\t\"   \"-\"   \"\t\"   \"-\"   \"\t\" \
					\"-\"  \"\t\"   \"-\" }}' >> {in3}".format(in1=self.bcf_list, in2=self.input, in3=self.output1)  
							
    
		### Run bcftools query			
		run_cmd(self.cmd_bcftools, self.args.verbose, 1, "### Step 3 - Converting VCF to human readable file via bcftools query and storing in .var ###")		
				
		### Check if output file have nonzero size
		check_file_exists_and_nonzero(self.output1)  	### var
		
###---------------------------------------------------------------------					
		
class Step4(Steps):
	"""Extract Phred scores and more from Pileup and include in var as var.qual"""

	def __init__ (self, args, step_index):
		"""Initialize step object"""
		
		### Parent's constructor
		Steps.__init__(self, args, step_index)
		
		### Define description as well as input and output attributes
		self.set_descrip("###### Step 4 - Extract Phred scores and more from Pileup and store in .var.qual ######")
		self.set_input(self.args.RES_dir + "/" + self.args.Name + "/" + self.args.Name + ".var")                                   		        ### var input
		self.set_output1(self.args.RES_dir + "/" + self.args.Name + "/" + self.args.Name + ".var.qual")                					### var.qual output
		
	@mytimer
	def run(self):	
		"""Run bcftools query"""
		### Run parent's run method
		super(Step4, self).run()
		
		### Read Pileup into dictionary 
		self.dict_pileup = {}	
		with open(self.args.RES_dir + "/" + self.args.Name + "/" + self.args.Name + ".mpileup", 'r') as self.fpileup:
			for self.line in self.fpileup:
				if ( self.line[0:1] != "#" ):
					self.cols_mp = self.line.split("\t")
					self.dict_pileup['-'.join((self.cols_mp[0], self.cols_mp[1]))] = self.cols_mp
					
		### Write header for var
		if ("2" in self.args.steps):
		    		if (self.args.ann):
					self.header_varqual = ["#chr", "pos", "ref", "ref_ctx", "alt", "alt_ctx", "total_depth_frw", "total_depth_rev", \
						"alt_depth_frw", "alt_depth_rev", "freq", "ave_alt_qual_frw", "ave_alt_qual_rev", \
						"pos_inread_alt_frw", "pos_inread_alt_rev", "indel_bases", "snp", "CDS", "effect"]
				else:	
					self.header_varqual = ["#chr", "pos", "ref", "ref_ctx", "alt", "alt_ctx", "total_depth_frw", "total_depth_rev", \
						"alt_depth_frw", "alt_depth_rev", "freq", "ave_alt_qual_frw", "ave_alt_qual_rev", \
						"pos_inread_alt_frw", "pos_inread_alt_rev", "indel_bases", "effect"]					
					
		else:
				self.header_varqual = ["#chr", "pos", "ref", "ref_ctx", "alt", "alt_ctx", "total_depth_frw", "total_depth_rev", \
					"alt_depth_frw", "alt_depth_rev", "freq", "ave_alt_qual_frw", "ave_alt_qual_rev", \
					"pos_inread_alt_frw", "pos_inread_alt_rev", "indel_bases"]		
						
		with open(self.output1, 'w') as self.fvarqual:
    			csv.writer(self.fvarqual, delimiter='\t').writerow(self.header_varqual)
			
		### Loop over variants in var
		with open(self.input, 'r') as fvar:	
			for self.line in fvar:
				if (self.line[0:1] != "#"):
					self.varqualrow = PileupExtractor(self.line, self.dict_pileup)
					with open(self.output1, 'a') as fvarqual:
						self.varqualrow[-1] = self.varqualrow[-1].rstrip()
						fvarqual.write("\t".join(self.varqualrow) + "\n")      

		### Check if output file have nonzero size
		check_file_exists_and_nonzero(self.output1)  	### var.qual
		
###---------------------------------------------------------------------					

class Step5(Steps):
	"""Adding context to variants"""

	def __init__ (self, args, step_index):
		"""Initialize step object"""
		
		### Parent's constructor
		Steps.__init__(self, args, step_index)
		
		### Define description as well as input and output attributes
		self.set_descrip("###### Step 5 - Adding context to variants and generating .var.ctx ######")
		self.set_input(self.args.RES_dir + "/" + self.args.Name + "/" + self.args.Name + ".var.qual")                	### var.qual input
		self.set_output1(self.args.RES_dir + "/" + self.args.Name + "/" + self.args.Name + ".var.ctx")                  ### var.ctx output
		self.dashlist = list('-----------------------------------------------------')	

		if ("2" in self.args.steps):
		    		if (self.args.ann):
					self.header_varctx = ["#chr", "pos", "ref", "ref_ctx", "alt", "alt_ctx", "total_depth_frw", "total_depth_rev", \
						"alt_depth_frw", "alt_depth_rev", "freq", "ave_alt_qual_frw", "ave_alt_qual_rev", \
						"pos_inread_alt_frw", "pos_inread_alt_rev", "indel_bases", "No. of hp repeats", "snp", "CDS", "effect"]
				else:	
					self.header_varctx = ["#chr", "pos", "ref", "ref_ctx", "alt", "alt_ctx", "total_depth_frw", "total_depth_rev", \
						"alt_depth_frw", "alt_depth_rev", "freq", "ave_alt_qual_frw", "ave_alt_qual_rev", \
						"pos_inread_alt_frw", "pos_inread_alt_rev", "indel_bases", "No. of hp repeats", "effect"]					
					
		else:
				self.header_varctx = ["#chr", "pos", "ref", "ref_ctx", "alt", "alt_ctx", "total_depth_frw", "total_depth_rev", \
					    "alt_depth_frw", "alt_depth_rev", "freq", "ave_alt_qual_frw", "ave_alt_qual_rev", \
					    "pos_inread_alt_frw", "pos_inread_alt_rev", "indel_bases", "No. of hp repeats"]
	
	@mytimer
	def run(self):	
		"""Get the context"""
		### Run parent's run method
		super(Step5, self).run()
		
		### Loop over variants in var
		with open(self.input, 'r') as fvarqual:	
		    	self.alt_chr_p = 0
			for self.line in fvarqual:
				if (self.line[0:1] == "#"):
					with open(self.output1, 'w') as self.fvarctx:
    						csv.writer(self.fvarctx, delimiter='\t').writerow(self.header_varctx)
				else:
							    
					self.cols = self.line.split("\t")

					self.alt_chr = self.cols[0][3:]
					self.alt_pos = int(self.cols[1])
					self.reff = self.cols[2]
					self.altt = self.cols[4]
					self.lref = len(self.reff)
					self.lalt = len(self.altt)
					self.ind_bases = self.cols[15].rstrip()

					if (self.alt_chr != self.alt_chr_p):

					    self.cmd_faidx = "samtools faidx {} chr{}".format(self.args.REF, self.alt_chr)
					    self.refseq = run_cmd(self.cmd_faidx, self.args.verbose, 1, "### Step 5 ###")
					    self.refseq = self.refseq.replace(">chr{}".format(self.alt_chr), '').replace('\n', '').replace('\r', '')

					### SNV
					if ( self.lref == 1 and self.lalt == 1):
					    self.ref_cntx = map(str.upper,self.refseq[self.alt_pos-2:self.alt_pos+1])
					    self.alt_cntx = self.ref_cntx
					    self.ref_cntx = ''.join(self.ref_cntx)
					    self.alt_cntx[1] = self.altt
					    self.alt_cntx = ''.join(self.alt_cntx)
					    
					    self.num_homo_nt = 0
					
					#############		
					### Deletion
					if (self.lref > self.lalt):
					    
					    self.ref_cntx = map(str.upper,self.refseq[self.alt_pos-1:self.alt_pos+self.lref])
					    self.ref_cntx = ''.join(self.ref_cntx)
					    self.alt_cntx = self.ref_cntx[0] + (self.ref_cntx[1:self.lref]).replace(self.ind_bases, ''.join(self.dashlist[0:len(self.ind_bases)]), 1) + self.ref_cntx[-1:]
					    self.homo_rmdr = (self.ref_cntx[1:self.lref]).replace(self.ind_bases, '')  
					    
					    if allTheSame(self.ind_bases):
						self.num_homo_nt = (self.alt_cntx[1:self.lalt]).count(self.ind_bases[0])		    
					    else:
						self.num_homo_nt = (self.alt_cntx[1:self.lalt]).count(self.ind_bases)
			    
					    if allTheSame(self.ind_bases):
						self.num_homo_nt = (self.ref_cntx[1:self.lref]).count(self.ind_bases[0])		    
					    else:
						self.num_homo_nt = (self.ref_cntx[1:self.lref]).count(self.ind_bases)
						
					    self.ind_bases = '-' + self.ind_bases
													    
					#############
					### Insertion 
					if (self.lref < self.lalt):
					    self.alt_cntx = map(str.upper,self.refseq[self.alt_pos-1:self.alt_pos+self.lref])
					    self.alt_cntx = ''.join(self.alt_cntx)
					    self.alt_cntx = self.alt_cntx[0] + self.ind_bases + self.alt_cntx[1:]
					    self.ref_cntx = self.alt_cntx[0] + (self.alt_cntx[1:self.lalt]).replace(self.ind_bases, ''.join(self.dashlist[0:len(self.ind_bases)]), 1) + self.alt_cntx[-1:]
					    self.homo_rmdr = (self.alt_cntx[1:self.lalt]).replace(self.ind_bases, '')  
			    
						
					    if allTheSame(self.ind_bases):
						self.num_homo_nt = (self.alt_cntx[1:self.lalt]).count(self.ind_bases[0])		    
					    else:
						self.num_homo_nt = (self.alt_cntx[1:self.lalt]).count(self.ind_bases)		    
			    
					    self.ind_bases = '+' + self.ind_bases
			    
					#############
					self.varctxrow = self.cols[0:3] + self.ref_cntx.split() + self.altt.split() + self.alt_cntx.split() + self.cols[6:15] +  self.ind_bases.split() + (str(self.num_homo_nt)).split() + self.cols[16:]
					    
					with open(self.output1, 'a') as self.fvarctx:	
						self.varctxrow[-1] = self.varctxrow[-1].rstrip()
						self.fvarctx.write("\t".join(self.varctxrow) + "\n")
					
					self.alt_chr_p = self.alt_chr		
					
		### Check if output file have nonzero size
		check_file_exists_and_nonzero(self.output1)  	### var.ctx
		     
###---------------------------------------------------------------------					

class Step6(Steps):
	"""Adding context to variants"""

	def __init__ (self, args, step_index):
		"""Initialize step object"""
		
		### Parent's constructor
		Steps.__init__(self, args, step_index)
		
		### Define description as well as input and output attributes
		self.set_descrip("###### Step 6 - Estimate the error rate. ######")
	    	self.set_input(self.args.RES_dir + "/" + self.args.Name + "/" + self.args.Name + ".var.ctx")                  ### var.ctx output
		
		self.dashstr = '-----------------------------------------------------'
		self.dotstr = '......................................................'

		self.ind1nt = ['T', 'G', 'C', 'A']
		self.ind2nt = [['G','T'], ['G','C'], ['G','A'], ['A','T'], ['A','G'], ['A','C'], ['C','A'], ['C','G'], ['C','T'], ['T','A'], ['T','C'], ['T','G']]
		
		
		self.header_snv_count = ["#ref_ctx", "alt_ctx", "total_depth_frw", "total_depth_rev", \
		"alt_depth_frw", "alt_depth_rev", "alt_qual_frw", "alt_qual_rev", \
		"alt_pos_in_read_frw", "alt_pos_in_read_rev", "occurrence_frw", \
		"occurrence_rev"]  
	
		self.header_snv_phat = ["#ref_ctx", "alt_ctx", "phat_PE_frw", "phat_PE_rev", \
				"phat_merged", "alt_qual_frw", "alt_qual_rev", \
				"alt_pos_in_read_frw", "alt_pos_in_read_rev", "occurrence_frw", \
				"occurrence_rev"] 
				
		self.header_ins_count = ["#inserted_bases", "total_depth_frw", "total_depth_rev", \
				"alt_depth_frw", "alt_depth_rev", "alt_qual_frw", "alt_qual_rev", \
				"alt_pos_in_read_frw", "alt_pos_in_read_rev", "occurrence_frw", \
				"occurrence_rev"]  
			
		self.header_ins_phat = ["#inserted_bases", "phat_PE_frw", "phat_PE_rev", \
				"phat_merged", "alt_qual_frw", "alt_qual_rev", \
				"alt_pos_in_read_frw", "alt_pos_in_read_rev", "occurrence_frw", \
				"occurrence_rev"] 
				
		self.header_del_count = ["#deleted_bases", "total_depth_frw", "total_depth_rev", \
				"alt_depth_frw", "alt_depth_rev", "alt_qual_frw", "alt_qual_rev", \
				"alt_pos_in_read_frw", "alt_pos_in_read_rev", "occurrence_frw", \
				"occurrence_rev"]  
			
		self.header_del_phat = ["#deleted_bases", "phat_PE_frw", "phat_PE_rev", \
				"phat_merged", "alt_qual_frw", "alt_qual_rev", \
				"alt_pos_in_read_frw", "alt_pos_in_read_rev", "occurrence_frw", \
				"occurrence_rev"] 
				
		### Create reuired directories 		
		self.err_count_res_dir = self.args.RES_dir + self.args.Name + "/err_count/"    ### directory for counting errors
		if not os.path.exists(self.err_count_res_dir):
			os.makedirs(self.err_count_res_dir)
			
		self.err_rate_res_dir = self.args.RES_dir + self.args.Name + "/err_rate/"      ### directory for error rates
		if not os.path.exists(self.err_rate_res_dir):
			os.makedirs(self.err_rate_res_dir)
		
		### Remove DEL_1nt_hp_err_count and INS_1nt_hp_err_count if they already exist
		self.countfiles = ['DEL_1nt_hp_err_count', 'INS_1nt_hp_err_count']
		for i_cf in range(len(self.countfiles)):
		    self.outputcsv = self.err_count_res_dir + self.countfiles[i_cf] + ".csv"
		    if os.path.exists(self.outputcsv):
				self.cmd_rm = "rm {}".format(self.outputcsv)
				os.system(self.cmd_rm)
			    
		self.ratefiles = ['DEL_1nt_hp_err_rate', 'INS_1nt_hp_err_rate']
		for i_rf in range(len(self.ratefiles)):
		    self.outputcsv = self.err_rate_res_dir + self.ratefiles[i_rf] + ".csv"
		    if os.path.exists(self.outputcsv):
				self.cmd_rm = "rm {}".format(self.outputcsv)
				os.system(self.cmd_rm)
		
		### Remove DEL_2nt_hp_err_count and INS_2nt_hp_err_count if they already exist	
		if self.args.two_nt_indel:
			self.countfiles = ['DEL_2nt_hp_err_count', 'INS_2nt_hp_err_count']
			for i_cf in range(len(self.countfiles)):	
			    self.outputcsv = self.err_count_res_dir + self.countfiles[i_cf] + ".csv"
			    if os.path.exists(self.outputcsv):
					self.cmd_rm = "rm {}".format(self.outputcsv)
					os.system(self.cmd_rm)
			    
			self.ratefiles = ['DEL_2nt_hp_err_rate', 'INS_2nt_hp_err_rate']
			for i_rf in range(len(self.ratefiles)):
			    self.outputcsv = self.err_rate_res_dir + self.ratefiles[i_rf] + ".csv"
			    if os.path.exists(self.outputcsv):
					self.cmd_rm = "rm {}".format(self.outputcsv)
					os.system(self.cmd_rm)
		    
		###  Generate 96 context-specific SNVs
		self.snv_mut_list = np.chararray([192,2], itemsize=16) 
		
		self.itr = 0
		for k in range(len(self.ind2nt)):
		    for j in range(len(self.ind1nt)):
				for i in range(len(self.ind1nt)):
				    if self.itr < 96:
						self.ref = self.ind1nt[i] + self.ind2nt[k][0] + self.ind1nt[j]
						self.var = self.ind1nt[i] + self.ind2nt[k][1] + self.ind1nt[j]
				    else:
						self.ref = self.ind1nt[3-j] + self.ind2nt[k][0] + self.ind1nt[3-i]
						self.var = self.ind1nt[3-j] + self.ind2nt[k][1] + self.ind1nt[3-i]	
					
				    self.snv_mut_list[self.itr,0] = self.ref
				    self.snv_mut_list[self.itr,1] = self.var
				    self.itr += 1
				    
		###  Generate 4 INDELs - 1 nt
		self.ind1_mut_list = np.chararray([4,1], itemsize=16) 
		for i in range(len(self.ind1nt)):
		    self.ind1_mut_list[i,0] = self.ind1nt[i]
		    
		###  Generate 16 INDELs - 2 nt		    
		self.ind2_mut_list = np.chararray([16,1], itemsize=16) 
		self.itr = 0
		for j in range(len(self.ind1nt)):
			for i in range(len(self.ind1nt)):
			    self.indel_bases = self.ind1nt[i] + self.ind1nt[j]
			    self.ind2_mut_list[self.itr, 0] = self.indel_bases
			    self.itr += 1
		
		### Read var.ctx into dictionary 
		self.dctx_snv = {}
		self.dctx_ins = {}
		self.dctx_del = {}
	    	self.itr = 1;
		with open(self.input, 'r') as self.fvarctx:
			for self.line in self.fvarctx:
				if (self.line[0:1] != "#"):
				    self.line = self.line.rstrip()
				    self.cols = self.line.split("\t")
			
				    if (len(self.cols[2]) == len(self.cols[4]) and float(self.cols[10]) < self.args.cut_off_freq and (int(self.cols[6]) + int(self.cols[7])) > self.args.min_depth):
					self.mykey = str(self.itr) + ',' + self.cols[3] + ',' + self.cols[5]
					self.dctx_snv[self.mykey] = self.cols[0:]
						    
				    if (len(self.cols[2]) < len(self.cols[4]) and float(self.cols[10]) < self.args.cut_off_freq and (int(self.cols[6]) + int(self.cols[7])) > self.args.min_depth):
					self.mykey = str(self.itr) + ',' + str(len(self.cols[15])-1) + ',' + self.cols[15] + ',HP:' + self.cols[16] + '+'
					self.dctx_ins[self.mykey] = self.cols[0:]
			
				    if (len(self.cols[2]) > len(self.cols[4]) and float(self.cols[10]) < self.args.cut_off_freq and (int(self.cols[6]) + int(self.cols[7])) > self.args.min_depth):
					self.mykey = str(self.itr) + ',' + str(len(self.cols[15])-1) + ',' + self.cols[15] + ',HP:' + self.cols[16] + '-'
					self.dctx_del[self.mykey] = self.cols[0:]
					
				self.itr += 1		    

		
	@mytimer
	def run(self):	
		"""Get the error rate"""
		### Run parent's run method
		super(Step6, self).run()		
		
		#############
		### 96 context-specific SNVs
		self.snv_count = np.zeros([96,10]) 
		self.snv_phat = np.zeros([96,9]) 
		
		for i_snv in range(96):
		
			self.ref_cntx = self.snv_mut_list[i_snv,0]
			self.alt_cntx = self.snv_mut_list[i_snv,1] 
			self.ref_cntx_c = self.snv_mut_list[i_snv + 96,0]
			self.alt_cntx_c = self.snv_mut_list[i_snv + 96,1] 
			
			self.mut_pos_tmp1 = list(map(float, self.myval[6:15]) for self.mkkey, self.myval in self.dctx_snv.iteritems() if (self.ref_cntx + ',' + self.alt_cntx) in self.mkkey) 
			self.mut_pos_tmp2 = list(map(float, self.myval[6:15]) for self.mkkey, self.myval in self.dctx_snv.iteritems() if (self.ref_cntx_c + ',' + self.alt_cntx_c) in self.mkkey) 
			
			self.mut_pos_tmp = self.mut_pos_tmp1 + self.mut_pos_tmp2
			self.mut_pos = np.array(self.mut_pos_tmp)
		    
			if self.mut_pos_tmp:
				[self.snv_count[i_snv,:], self.snv_phat[i_snv,:]] = getErrorRate(self.mut_pos)
				
		### write snv_count to file
		self.snv_count = np.concatenate((self.snv_mut_list[96:,:], self.snv_count), axis=1)	
		with open(self.err_count_res_dir + "SNV_err_count.csv", 'w') as self.f100:
			csv.writer(self.f100).writerow(self.header_snv_count)
			csv.writer(self.f100).writerows(self.snv_count)
			
		### write snv_phat to file
		self.snv_phat = np.concatenate((self.snv_mut_list[96:,:], self.snv_phat), axis=1)	
		with open(self.err_rate_res_dir + "SNV_err_rate.csv", 'w') as self.f101:
			csv.writer(self.f101).writerow(self.header_snv_phat)
			csv.writer(self.f101).writerows(self.snv_phat)

		#############
		### Single nt indel - no homopolymeric 
		self.del_count = np.zeros([4,10]) 
		self.ins_count = np.zeros([4,10]) 
		
		self.del_phat = np.zeros([4,9]) 
		self.ins_phat = np.zeros([4,9]) 
		
		for i_ind1 in range(len(self.ind1_mut_list)):
		    
			self.indel_nt = self.ind1_mut_list[i_ind1,0]
		    
			### INS
			self.mut_pos_ins1 = list(map(float, self.myval[6:15]) for self.mkkey, self.myval in self.dctx_ins.iteritems() if (',+' + self.indel_nt + ',HP:1+') in self.mkkey) 
			
			self.mut_pos_tmp_ins = self.mut_pos_ins1 
			self.mut_pos_ins = np.array(self.mut_pos_tmp_ins)
			    
			if self.mut_pos_tmp_ins:
				[self.ins_count[i_ind1,:], self.ins_phat[i_ind1,:]] = getErrorRate(self.mut_pos_ins)
			
			### DEL
			self.mut_pos_del1 = list(map(float, self.myval[6:15]) for self.mkkey, self.myval in self.dctx_del.iteritems() if (',-' + self.indel_nt + ',HP:1-') in self.mkkey) 
		    
			self.mut_pos_tmp_del = self.mut_pos_del1 
			self.mut_pos_del = np.array(self.mut_pos_tmp_del)
			
			if self.mut_pos_tmp_del:
				[self.del_count[i_ind1,:], self.del_phat[i_ind1,:]] = getErrorRate(self.mut_pos_del)
		
		### write del_count and ins_count to file
		self.ins_count = np.concatenate((self.ind1_mut_list[:,:], self.ins_count), axis=1)	
		with open(self.err_count_res_dir + "INS_1nt_err_count.csv", 'w') as self.f200:
			csv.writer(self.f200).writerow(self.header_ins_count)
			csv.writer(self.f200).writerows(self.ins_count)
		    
		self.del_count = np.concatenate((self.ind1_mut_list[:,:], self.del_count), axis=1)	    
		with open(self.err_count_res_dir + "DEL_1nt_err_count.csv", 'w') as self.f201:
			csv.writer(self.f201).writerow(self.header_del_count)
			csv.writer(self.f201).writerows(self.del_count)	
		    
		### write del_phat and ins_phat to file
		self.ins_phat = np.concatenate((self.ind1_mut_list[:,:], self.ins_phat), axis=1)	
		with open(self.err_rate_res_dir + "INS_1nt_err_rate.csv", 'w') as self.f202:
			csv.writer(self.f202).writerow(self.header_ins_phat)
			csv.writer(self.f202).writerows(self.ins_phat)
		    
		self.del_phat = np.concatenate((self.ind1_mut_list[:,:], self.del_phat), axis=1)	
		with open(self.err_rate_res_dir + "DEL_1nt_err_rate.csv", 'w') as self.f203:
			csv.writer(self.f203).writerow(self.header_del_phat)
			csv.writer(self.f203).writerows(self.del_phat)
    
		#############
		### Single nt indel - homopolymeric 
								
		for kl in range(2,self.args.max_hp_no+1):
			    
			self.del_count = np.zeros([4,10]) 
			self.ins_count = np.zeros([4,10]) 
			
			self.del_phat = np.zeros([4,9]) 
			self.ins_phat = np.zeros([4,9]) 
			
			for ij in range(len(self.ind1_mut_list)):
			    
			    self.indel_nt = self.ind1_mut_list[ij,0]
			
			    ### INS
			    self.mut_pos_ins1 = list(map(float, self.myval[6:15]) for self.mkkey, self.myval in self.dctx_ins.iteritems() if (',+' + self.indel_nt + ',HP:' + str(kl) + '+') in self.mkkey) 
			    
			    self.mut_pos_tmp_ins = self.mut_pos_ins1 
			    self.mut_pos_ins = np.array(self.mut_pos_tmp_ins)
				
			    if self.mut_pos_tmp_ins:
				    [self.ins_count[ij,:], self.ins_phat[ij,:]] = getErrorRate(self.mut_pos_ins)
			    
			    ### DEL
			    self.mut_pos_del1 = list(map(float, self.myval[6:15]) for self.mkkey, self.myval in self.dctx_del.iteritems() if (',-' + self.indel_nt + ',HP:' + str(kl) + '-') in self.mkkey) 
			
			    self.mut_pos_tmp_del = self.mut_pos_del1 
			    self.mut_pos_del = np.array(self.mut_pos_tmp_del)
			    
			    if self.mut_pos_tmp_del:
				    [self.del_count[ij,:], self.del_phat[ij,:]] = getErrorRate(self.mut_pos_del)
			
			#### write del_count and ins_count to file
			self.ins_count = np.concatenate((self.ind1_mut_list[:,:], self.ins_count), axis=1)	
			with open(self.err_count_res_dir + "INS_1nt_hp_err_count.csv", 'a') as self.f300:
			    csv.writer(self.f300).writerow(['##Number of homopolymer repeats:', str(kl)])
			    csv.writer(self.f300).writerow(self.header_ins_count)
			    csv.writer(self.f300).writerows(self.ins_count)
			    
			self.del_count = np.concatenate((self.ind1_mut_list[:,:], self.del_count), axis=1)	    
			with open(self.err_count_res_dir + "DEL_1nt_hp_err_count.csv", 'a') as self.f301:
			    csv.writer(self.f301).writerow(['##Number of homopolymer repeats:', str(kl)])
			    csv.writer(self.f301).writerow(self.header_del_count)
			    csv.writer(self.f301).writerows(self.del_count)	
			    
			#### write del_phat and ins_phat to file
			self.ins_phat = np.concatenate((self.ind1_mut_list[:,:], self.ins_phat), axis=1)	
			with open(self.err_rate_res_dir + "INS_1nt_hp_err_rate.csv", 'a') as self.f302:
			    csv.writer(self.f302).writerow(['##Number of homopolymer repeats:', str(kl)])
			    csv.writer(self.f302).writerow(self.header_ins_phat)
			    csv.writer(self.f302).writerows(self.ins_phat)
			    
			self.del_phat = np.concatenate((self.ind1_mut_list[:,:], self.del_phat), axis=1)	
			with open(self.err_rate_res_dir + "DEL_1nt_hp_err_rate.csv", 'a') as self.f303:
			    csv.writer(self.f303).writerow(['##Number of homopolymer repeats:', str(kl)])
			    csv.writer(self.f303).writerow(self.header_del_phat)
			    csv.writer(self.f303).writerows(self.del_phat)
		    
		
		#############
		### Double nt indel - no homopolymeric 						
		if self.args.two_nt_indel:
			self.del_2_count = np.zeros([16,10]) 
			self.ins_2_count = np.zeros([16,10]) 
			
			self.del_2_phat = np.zeros([16,9]) 
			self.ins_2_phat = np.zeros([16,9]) 
			
			for ij in range(len(self.ind2_mut_list)):
			
			    self.indel_nt = self.ind2_mut_list[ij,0]
			
			    ### INS
			    self.mut_pos_ins1 = list(map(float, self.myval[6:15]) for self.mkkey, self.myval in self.dctx_ins.iteritems() if (',+' + self.indel_nt + ',HP:1+') in self.mkkey) 
			    
			    self.mut_pos_tmp_ins = self.mut_pos_ins1 
			    self.mut_pos_ins = np.array(self.mut_pos_tmp_ins)
				
			    if self.mut_pos_tmp_ins:
				    [self.ins_2_count[ij,:], self.ins_2_phat[ij,:]] = getErrorRate(self.mut_pos_ins)
			    
			    ### DEL
			    self.mut_pos_del1 = list(map(float, self.myval[6:15]) for self.mkkey, self.myval in self.dctx_del.iteritems() if (',-' + self.indel_nt + ',HP:1-') in self.mkkey) 
			
			    self.mut_pos_tmp_del = self.mut_pos_del1 
			    self.mut_pos_del = np.array(self.mut_pos_tmp_del)
			    
			    if self.mut_pos_tmp_del:
				    [self.del_2_count[ij,:], self.del_2_phat[ij,:]] = getErrorRate(self.mut_pos_del)
			
			#### write del_2_count and ins_2_count to file
			self.ins_2_count = np.concatenate((self.ind2_mut_list[:,:], self.ins_2_count), axis=1)	
			with open(self.err_count_res_dir + "INS_2nt_err_count.csv", 'w') as self.f400:
			    csv.writer(self.f400).writerow(self.header_ins_count)
			    csv.writer(self.f400).writerows(self.ins_2_count)
			    
			self.del_2_count = np.concatenate((self.ind2_mut_list[:,:], self.del_2_count), axis=1)	    
			with open(self.err_count_res_dir + "DEL_2nt_err_count.csv", 'w') as self.f401:
			    csv.writer(self.f401).writerow(self.header_del_count)
			    csv.writer(self.f401).writerows(self.del_2_count)	
			    
			#### write del_2_phat and ins_2_phat to file
			self.ins_2_phat = np.concatenate((self.ind2_mut_list[:,:], self.ins_2_phat), axis=1)	
			with open(self.err_rate_res_dir + "INS_2nt_err_rate.csv", 'w') as self.f402:
			    csv.writer(self.f402).writerow(self.header_ins_phat)
			    csv.writer(self.f402).writerows(self.ins_2_phat)
			    
			self.del_2_phat = np.concatenate((self.ind2_mut_list[:,:], self.del_2_phat), axis=1)	
			with open(self.err_rate_res_dir + "DEL_2nt_err_rate.csv", 'w') as self.f403:
			    csv.writer(self.f403).writerow(self.header_del_phat)
			    csv.writer(self.f403).writerows(self.del_2_phat)			
			
			
		#############
		### Double nt indel - homopolymeric 						
		if self.args.two_nt_indel:			
			for kl in range(2,self.args.max_hp_no+1):
				self.del_2_count = np.zeros([16,10]) 
				self.ins_2_count = np.zeros([16,10]) 
				
				self.del_2_phat = np.zeros([16,9]) 
				self.ins_2_phat = np.zeros([16,9]) 
				
				for ij in range(len(self.ind2_mut_list)):
				
				    self.indel_nt = self.ind2_mut_list[ij,0]
				
				    ### INS
				    self.mut_pos_ins1 = list(map(float, self.myval[6:15]) for self.mkkey, self.myval in self.dctx_ins.iteritems() if (',+' + self.indel_nt + ',HP:' + str(kl) + '+') in self.mkkey) 
				    
				    self.mut_pos_tmp_ins = self.mut_pos_ins1 
				    self.mut_pos_ins = np.array(self.mut_pos_tmp_ins)
					
				    if self.mut_pos_tmp_ins:
					    [self.ins_2_count[ij,:], self.ins_2_phat[ij,:]] = getErrorRate(self.mut_pos_ins)
				    
				    ### DEL
				    self.mut_pos_del1 = list(map(float, self.myval[6:15]) for self.mkkey, self.myval in self.dctx_del.iteritems() if (',-' + self.indel_nt + ',HP:' + str(kl) + '-') in self.mkkey) 
				
				    self.mut_pos_tmp_del = self.mut_pos_del1 
				    self.mut_pos_del = np.array(self.mut_pos_tmp_del)
				    
				    if self.mut_pos_tmp_del:
					    [self.del_2_count[ij,:], self.del_2_phat[ij,:]] = getErrorRate(self.mut_pos_del)
				
				#### write del_2_count and ins_2_count to file
				self.ins_2_count = np.concatenate((self.ind2_mut_list[:,:], self.ins_2_count), axis=1)	
				with open(self.err_count_res_dir + "INS_2nt_hp_err_count.csv", 'a') as self.f:
				    csv.writer(self.f).writerow(['##Number of repeats:', str(kl)])
				    csv.writer(self.f).writerow(self.header_ins_count)
				    csv.writer(self.f).writerows(self.ins_2_count)
				    
				self.del_2_count = np.concatenate((self.ind2_mut_list[:,:], self.del_2_count), axis=1)	    
				with open(self.err_count_res_dir + "DEL_2nt_hp_err_count.csv", 'a') as self.f:
				    csv.writer(self.f).writerow(['##Number of repeats:', str(kl)])
				    csv.writer(self.f).writerow(self.header_del_count)
				    csv.writer(self.f).writerows(self.del_2_count)	
				    
				#### write del_2_phat and ins_2_phat to file
				self.ins_2_phat = np.concatenate((self.ind2_mut_list[:,:], self.ins_2_phat), axis=1)	
				with open(self.err_rate_res_dir + "INS_2nt_hp_err_rate.csv", 'a') as self.f:
				    csv.writer(self.f).writerow(['##Number of repeats:', str(kl)])
				    csv.writer(self.f).writerow(self.header_ins_phat)
				    csv.writer(self.f).writerows(self.ins_2_phat)
				    
				self.del_2_phat = np.concatenate((self.ind2_mut_list[:,:], self.del_2_phat), axis=1)	
				with open(self.err_rate_res_dir + "DEL_2nt_hp_err_rate.csv", 'a') as self.f:
				    csv.writer(self.f).writerow(['##Number of repeats:', str(kl)])
				    csv.writer(self.f).writerow(self.header_del_phat)
				    csv.writer(self.f).writerows(self.del_2_phat)				    
    		
args = get_arg()
main()
