# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 15:04:19 2015

Based on the hyper.enrichment R script excerpt

@author: dbray
"""
import scipy.stats as stats
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

def read_gene_sigs(gene_sigs_file):
    """ returns the gene signatures in the input file """
    
    # open the gene signature file and read the header line
    f = open(gene_sigs_file, 'r')
    line = f.readline()
    
    # initialize a nested list to hold the gene signatures and an int to store cluster number
    gene_sigs = []
    cluster_number = len(line.split())
    for c in range(cluster_number):
        gene_sigs.append([])
    
    # iterate through the rest of the file to read the entire gene signature
    while True:
        line = f.readline()
        if not line:
            break # break out of the loop at the end of the file
        gene_parts = line.split()
        
        # sort each gene to its proper cluster signature
        for c in range(cluster_number):
            gene_sigs[c].append(gene_parts[c])
    
    # close the input file
    f.close()
    
    # return the gene signatures
    return gene_sigs


def hypergeo_test(train_sig, test_sig, total_genes):
    """ Conducts a hypergeometric test to obtain p-values for cluster enrichment """
    
    # compute the intersection between the train set signature and the test set signature
    n_hits = len([gene for gene in train_sig if gene in test_sig])
    
    # create variables for the rest of the parameters of the hypergeometric test
    n_drawn = len(train_sig)
    n_cats = len(test_sig)
    
    # compute the p-value Pr(X >= n_hits)
    p = stats.hypergeom.sf(n_hits, total_genes, n_cats, n_drawn)
    
    # return the p-value computed
    return p


def correct_p_values(pvalues):
    """ apply Benjamini-Hochberg FDR correction to p-values """
    
    # import the stats package    
    stats = importr('stats')       
    
    # collapse the pvalues nested list into a single list
    temp_p_values = []
    for i in range(len(pvalues)):
        for j in range(len(pvalues[0])):
            temp_p_values.append(pvalues[i][j])
    
    # adjust the p values from the collapsed list using the Benjamini-Hochberg correction
    p = stats.p_adjust(FloatVector(temp_p_values), method = 'BH')
    
    # rebuild the original list
    iter_p = iter(p)
    for i in range(len(pvalues)):
        for j in range(len(pvalues[0])):
            pvalues[i][j] = next(iter_p)
    
    # return the corrected p values
    return pvalues


def write_p_values(all_p_values):
    """ writes the corrected p values to an output file """
    
    # open an output file
    output_file = "C:/Users/dbray/Dropbox/Boston University/Fall 2015/Challenge Project/TrainingSet/hypergeo_enrich_test.txt"
    o = open(output_file, 'w')
    header = ""   
    for i in range(len(all_p_values)):
        header += str(i+1) + '\t'
    o.write(header.rstrip() + '\n')
    for n in range(len(all_p_values[0])):
        line = ""
        for c in range(len(all_p_values)):
            line += (str(all_p_values[c][n]) + '\t')
        o.write(line.rstrip() + '\n')
    o.close()
    
    # return back
    return 1


def main():
    # assign variables to the file names of the gene signature files
    train_sigs_file = "C:/Users/dbray/Dropbox/Boston University/Fall 2015/Challenge Project/TrainingSet/sample_cluster_gene_sigs_TEST.txt"
    test_sigs_file = "C:/Users/dbray/Dropbox/Boston University/Fall 2015/Challenge Project/TrainingSet/sample_cluster_gene_sigs_TEST.txt"
    
    # read the gene signatures file and assign a variable name to each
    train_sigs = read_gene_sigs(train_sigs_file)
    test_sigs = read_gene_sigs(test_sigs_file)
    
    # store the number of genes used to determine the cluster
    gene_number_file = "C:/Users/dbray/Dropbox/Boston University/Fall 2015/Challenge Project/TrainingSet/genes4clustering.txt"
    f = open(gene_number_file, 'r')
    lines = f.readlines()
    total_genes = len(lines) - 1
    f.close()
    
    # pass the signature lists to the hypergeometric test function to obtain p-values
    all_p_values = []    
    for train_sig in train_sigs:
        p_values = []
        for test_sig in test_sigs:
            p = hypergeo_test(train_sig, test_sig, total_genes)
            p_values.append(p)
        all_p_values.append(p_values)
    
    # adjust p-values for FDR
    all_p_values = correct_p_values(all_p_values)
  
    # print the p-values to a table
    write_p_values(all_p_values)
    
    
if __name__ == "__main__":
    main()
