import sys
import json
import numpy as np
import random
import bisect
import matplotlib.pyplot as plt
from collections import defaultdict

# Copyright 2018.
#    Yi-Pin Lai and Thomas R. Ioerger.
#
#    ptACR is a free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License.
#
#    ptACR is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with the ptACR. If not, see <http://www.gnu.org/licenses/>.

#__author__ = "Yi-Pin Lai and Thomas R. Ioerger, https://github.com/yplai/ptACR"
#__version__ = "1.0.0"
#__license__ = "GNU"


# input format: similar to the fasta format, the number of taxa and the length of the alignment are shown in the first line, followed by lines where one line for each strain/isolate with its ID and DNA sequence
def parsing_input_SNP(infile):
	strain_id = []
        strain_snp = []
        strain_snp_dict = {}
        for ii, line in enumerate(open(infile)):
                if '#' in line: continue
		line = line.rstrip()
                line = line.split()
		if ii is 0:
                        num_strain = line[0]
                        num_snp = line[1]
                else:
                        if '>' in line[0]:
                                strain_id.append(line[0][1:])
                                strain_snp_dict[line[0][1:]] = line[1]
                        else:
                                strain_id.append(line[0])
                                strain_snp_dict[line[0]] = line[1]
                        strain_snp.append(line[1])
        return strain_id, strain_snp


# determine the compatibility for a pair of binary-state characters
def compatible_binary(C, D):
        m = len(C)
        pairs = []
        allSets = True
        for k in range(m):
                pair = C[k] + D[k]
                if pair not in pairs:
                        pairs.append(pair)
        if len(pairs) == 4:
                allSets = False
        return allSets


# dictionary of character patterens 
def site_char(strain_snp):
        site_dict = {}
        site_dict_set = {}
        for ii in range(len(strain_snp[0])):
                char = [cc[ii] for cc in strain_snp]
                char1 = np.unique(char)
		str_char = ''
		char2 = list(char1)
		for hh in char:
			aa = char2.index(hh)
			str_char += str(aa)
		site_dict[ii] = str_char
                idx1 = []
                for uu in char1:
                        idx1.append([pp for pp,cc in enumerate(char) if cc == uu])
		site_dict_set[str_char] = idx1
	return site_dict, site_dict_set


# determine the compatibility for a pair of multi-state characters
def compatible_multi(P, Q):
        nodeP = ['a', 'b', 'c', 'd']
        nodeQ = ['m', 'n', 'o', 'p']
        P = [set(w) for w in P]
        Q = [set(w) for w in Q]
        graph = defaultdict(list)
        visited = {}
        for i1, ss1 in enumerate(P):
                for i2, ss2 in enumerate(Q):
                        if len(ss1.intersection(ss2)) > 0:  # ss1 & ss2   
                                graph[nodeP[i1]].append(nodeQ[i2])
                                graph[nodeQ[i2]].append(nodeP[i1])
                        visited[nodeQ[i2]] = False
                visited[nodeP[i1]] = False
        for v in visited.keys():
                if visited[v] == False:
                        if detect_cycle(v,visited,v,graph) == True:
                                return False
        return True


# supporting function for "compatible_multi"
def detect_cycle(u, visited, parent, graph):
        visited[u] = True
        iscycle = False
        for v in graph[u]:
                if visited[v] == False:
                        if detect_cycle(v,visited,u,graph) == True:
                                iscycle = True
                                return iscycle
                elif parent != v:
                        iscycle = True
                        return iscycle
        return iscycle


# dictionay of compatibility for pairs of unique patterns
def pairwise_compat_dict(site_snp_dict, site_dict_set):
        all_sites = [site_snp_dict[pp] for pp in range(len(site_snp_dict))]
        unique_patterns = list(set(all_sites))
        pw_compat = {}
        for tt in unique_patterns:
                for tt2 in unique_patterns:
			uu = list(site_dict_set[tt])
			uu2 = list(site_dict_set[tt2])
                        if compatible_multi(uu,uu2) == True:
                                pw_compat[(tt,tt2)] = True
                                pw_compat[(tt2,tt)] = True
                        else:
                                pw_compat[(tt,tt2)] = False
                                pw_compat[(tt2,tt)] = False
        return pw_compat


# calculate compatibility score of all paired patterns of one from upstream regions and the other from downstream regions for a given site with a given window size
def compat_pattern_between(site_index_list_left,site_index_list_right,site_snp_dict,pw_compat_dict):
        patterns_left, patterns_right = {},{}
        pat_list_left, pat_list_right = [],[]
        for aa in site_index_list_left:
                snp = site_snp_dict[aa]
                if snp not in patterns_left:
                        patterns_left[snp] = 1
                        pat_list_left.append(snp)
                else:
                        patterns_left[snp] += 1
        for bb in site_index_list_right:
                snp = site_snp_dict[bb]
                if snp not in patterns_right:
                        patterns_right[snp] = 1
                        pat_list_right.append(snp)
                else:
                        patterns_right[snp] += 1
        scores = 0
        for cc in pat_list_left:
                for dd in pat_list_right:
                        if pw_compat_dict[(cc,dd)] == True:
                                scores += 1 * patterns_left[cc] * patterns_right[dd]
        return scores


# for a given local minimum, calculate p-value by performing permutaion
def between_permutation_dict(site_snp_dict, site_index, win_size, times, pw_compat_dict):
        upstream = range(int(site_index) - int(win_size)/2, int(site_index))
        downstream = range(int(site_index) + 1, int(site_index) + int(win_size)/2 + 1)
        compatible_score = compat_pattern_between(upstream, downstream, site_snp_dict, pw_compat_dict)
        region_sites = upstream + downstream
        between_regions= []
        for tt in range(times):
                shuffled_index = random.sample(region_sites, len(region_sites))
                upstream = shuffled_index[:len(region_sites)/2]
                downstream = shuffled_index[len(region_sites)/2:]
                ss = compat_pattern_between(upstream, downstream, site_snp_dict, pw_compat_dict)
                between_regions += [ss]

        between_regions = sorted(between_regions)
        if compatible_score < between_regions[0]:
                p2 = 1/float(times)
        else:
                p2 = bisect.bisect(between_regions, compatible_score)/float(times)
        return p2


# estimate overall compatibility ratio within a region
def compat_pattern_within(site_index_list, site_snp_dict, pw_compat_dict):
        patterns = {}
        pat_list = []
        for ss in site_index_list:
                snp = site_snp_dict[ss]
                if snp not in patterns:
                        patterns[snp] = 1
                        pat_list.append(snp)
                else:
                        patterns[snp] += 1
        scores = 0
        for ii, aa in enumerate(pat_list):
                for bb in pat_list[ii:]:
                        if aa == bb:
                                scores += 1 * patterns[aa] * (patterns[aa]-1)/2
                        else:
                                if pw_compat_dict[(aa,bb)] == True:
                                        scores += 1 * patterns[aa] * patterns[bb]
        nn = len(site_index_list) * (len(site_index_list)-1)/2
        ratio = scores / float(nn)
        return ratio


# apply Bonferroni correction for multiple tests
def bonferroni(pval):
	ntests = float(len(pval))
	corrected_pval = [pp*ntests for pp in pval]
	return corrected_pval


# detect valleys in a distribution (1-D array) that are at least seperated by a given interval
def detect_troughs(arr, min_itv):
        diff = [arr[ii]-arr[ii+1] for ii in range(len(arr)-2)]
        min_itv = int(min_itv)
        trough_map = {}
        for jj, dd in enumerate(diff):
                if dd <= 0:
                        trough_map[jj] = arr[jj]
        trough_idx = []
        for w in sorted(trough_map, key=trough_map.get, reverse=False):
                left = arr[w-min_itv/2:w]
                right = arr[w+1:w+min_itv/2+1]
                if trough_map[w] <= np.mean(left) and trough_map[w] <= np.mean(right):
                        if len(trough_idx)==0:
                                trough_idx.append(w)
                        else:
                                if check_interval(w, trough_idx, min_itv)==True:
                                        trough_idx.append(w)
        return sorted(trough_idx)


# supporting function for "detect_troughs" 
def check_interval(idx,trough_idx, min_itv):
        for dd in trough_idx:
                if idx >= dd-min_itv and idx <= dd+min_itv:
                        return False
        return True


# find the potential breakpoints and perform permutation test on them, and obtain breakpoints with adjusted p-value
def find_bk(start, end, win_size, site_snp_dict, pw_compat_dict, times=10000, empval=0.05, correction=True):	
        win_size = int(win_size)
	all_comp_ratio = []
	x_all = range(start,end)
        for xx in x_all:
		if xx-win_size/2 >= 0 and xx+win_size/2 <= len(site_snp_dict):
                	site_index_list = range(xx-win_size/2, xx+win_size/2)
                	rr = compat_pattern_within(site_index_list, site_snp_dict, pw_compat_dict)
                	all_comp_ratio.append(rr)

        with open('all_comp_ratio_multi.json','w') as ff:
                json.dump(all_comp_ratio,ff)

	localmin = detect_troughs(all_comp_ratio, float(win_size))
	#print len(localmin), 'local minima are found.'
	#print localmin

	if start < win_size/2:
		start2 = win_size/2
        	valley_idx = [ll+win_size/2 for ll in localmin]
        else:
		start2 = start
		valley_idx = [ll+start2 for ll in localmin]
	valley_val = [all_comp_ratio[ii] for ii in localmin]

        valley_idx2, valley_val2 = [],[]

        acr_ratio_cutoff = np.mean(all_comp_ratio)-np.std(all_comp_ratio)
        for ii, vv in enumerate(valley_val):
                if vv < acr_ratio_cutoff:
                        valley_idx2.append(valley_idx[ii])
                        valley_val2.append(vv)
	if len(valley_idx2) > 0:
		print
		print len(valley_idx2), 'potential breakpoints are found.'
        	print valley_idx2
	else:
		print
		print 'No potential breakpoint is found.'
		print 'Compatibility ratio:', compat_pattern_within(range(start,end), site_snp_dict, pw_compat_dict)
		return 
	
	#run permutation
        valley_idx3, valley_val3 = [],[]
        between_p = []
        for dd, site_index in enumerate(valley_idx2):
                pval = between_permutation_dict(site_snp_dict, site_index, win_size, times, pw_compat_dict)
                between_p.append(pval)
		if correction == False:
                	if pval < empval:
                        	valley_idx3.append(site_index)
        			valley_val3.append(valley_val2[dd]-0.005)

	if correction == True:
		empval_bonferroni = empval/float(len(between_p))
		#corrected_p = bonferroni(between_p)
                print 
		print 'Adjusted p-value cutoff:', empval_bonferroni
                cnt_bks = 0
		for dd, site_index in enumerate(valley_idx2):
			if between_p[dd] < empval_bonferroni:
                                cnt_bks += 1
				if cnt_bks == 1:
					print 'order', 'position', 'p-value'
				valley_idx3.append(site_index)
                        	valley_val3.append(valley_val2[dd]-0.005)
				print dd+1, site_index, between_p[dd]
	
	if len(valley_idx3) > 0:
		print
		print len(valley_idx3), 'recombination breakpoints are found.'	
		print valley_idx3
	else:
		print
		print 'No recombination breakpoint is found.'
		print 'Compatibility ratio:', compat_pattern_within(range(start,end), site_snp_dict, pw_compat_dict)
		return

	#overall compatibility and regional compatibility
	if start-win_size/2 <= 0:
                start = 0
        if end+win_size/2 >= len(site_snp_dict)-1:
                end = len(site_snp_dict)
	print
        print 'Overall compatibility', '[', start, ',', end, ']:', compat_pattern_within(range(start,end), site_snp_dict, pw_compat_dict)
        valley_idx_all = [start] + valley_idx3 + [end]
	for dd, xx in enumerate(valley_idx_all[:-1]):
		site_index_list = range(xx, valley_idx_all[dd+1])
                rr = compat_pattern_within(site_index_list, site_snp_dict, pw_compat_dict)
		print dd+1, '[', xx, ',', valley_idx_all[dd+1]-1, ']:', rr

	x_all2 = range(start2, start2+len(all_comp_ratio))
	return x_all2, all_comp_ratio, valley_idx2, valley_val2, valley_idx3, valley_val3


def main():
        #convert input file to dictionary
        infile = sys.argv[1]
        strain_id, strain_snp = parsing_input_SNP(infile)
        site_snp_dict, site_dict_set = site_char(strain_snp)
	print 'Number of strains:', len(strain_id)
	print 'Number of sites:', len(site_snp_dict)

        pw_compat_dict = pairwise_compat_dict(site_snp_dict,site_dict_set)
	print 'Number of unique paired patterns:', len(pw_compat_dict)

        start = int(sys.argv[2])
        end = int(sys.argv[3])
        print 'Starting position', start
	print 'Ending position:', end
	
	# up to five sliding window size can be applied to the ptACR at a time
        color = ['b','m','c','k','y']
        _, ax = plt.subplots(1, 1, figsize=(8, 6))
        for ii,ww in enumerate(sys.argv[4:]):
                win_size = int(ww)
		print 
		print 'Sliding window size:', win_size
                x_all, all_comp_ratio, valley_idx2, valley_val2, valley_idx3, valley_val3 = find_bk(start, end, win_size, site_snp_dict, pw_compat_dict, correction=True)
		ax.plot(x_all, all_comp_ratio , color[ii], lw=1, label=str(win_size))
                ax.plot(valley_idx2, valley_val2, '*', mfc=None, mec='r', mew=2, ms=8)
                ax.plot(valley_idx3, valley_val3, '+', mfc=None, mec='g', mew=2, ms=8)
        
        ax.set_ylabel('Average compatibility ratio', fontsize=14)
        ax.set_xlabel('Position of site', fontsize=14)
        ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
        plt.show()
        
if __name__ == '__main__':
        main()
