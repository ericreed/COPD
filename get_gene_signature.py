# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 20:31:41 2015

@author: dbray
"""
import itertools

def get_gene_sig(sample_cluster_file, gene_file, norm_expr_file, g):
    """ gets gene signature from a sample cluster """
    
    # initialize variables needed from the sample cluster file
    samples = []            # total of i samples
    sample_clusters = []    # 1 cluster per sample (total = i)
    
    # open the sample cluster input file
    with open(sample_cluster_file, 'r') as f:
        line = f.readline() # read the header
        while True:
            # read the current line
            line = f.readline()
            if not line:
                break # break at the end of the file
           
            # split the line into its parts
            line_parts = line.split()
           
            # store the relevant parts
            samples.append(line_parts[0])
            sample_clusters.append(line_parts[4])
            
    # initialize a variable to conveniently hold the number of clusters
    cluster_number = len(set(sample_clusters))
            
    # initialize variables needed from the normalized expression matrix
    expr_samples = []   # total of m samples
    genes = []          # total of n genes
    expr_matrix = []    # n genes x i (subset of m) matrix
    
    # open the DGE list of genes used to cluster
    with open(gene_file, 'r') as f:
       line = f.readline() # read the header
       while True:
           # read the current line
           line = f.readline()
           if not line:
               break # break at the end of the file
           
           # split the line into its parts
           line_parts = line.split()
           
           # store the gene on the line
           genes.append(line_parts[1])
    
    # open the normalized expression input file
    with open(norm_expr_file, 'r') as f:
        sample_line = f.readline().replace('"', '') 
        expr_samples = sample_line.split()  # store samples from first line
        
        # determine which samples to store gene expression values from
        relevant_indeces = []
        for expr_sample in expr_samples:
            for sample in samples:
                if str(expr_sample) == str(sample):
                    relevant_indeces.append(expr_samples.index(expr_sample))
                    continue
                
        # begin iterating through rest of the file
        while True:
            line = f.readline()
            if not line:
                break # break at the end of the file
            
            # split the line into its parts
            line_parts = line.split()
            
            # store the gene name from the first part
            gene_name = line_parts[0].replace('"', '')
            
            # SKIP GENE IF IT IS NOT USED IN THE CLUSTERING
            if gene_name not in genes:
                continue
            
            # if the gene is used in the clustering
            line_parts = line_parts[1:] # skip the gene name at the beginning
            
            # list comprehension to get expression values from relevant samples
            expr_values = [line_parts[i] for i in relevant_indeces]
            
            # add the list of values to the expression matrix
            expr_matrix.append(expr_values)
        
        # print dimensions of the expression matrix as check
        print("number of samples: " + str(len(expr_matrix[0])))
        print("number of genes: " + str(len(expr_matrix)))

    # initialize gene signature variables needed
    gene_freq_per_cluster = [{} for c in range(cluster_number)] # defines dictionaries (indexed by cluster number - 1)
    cluster_dict = {samples[i]: int(sample_clusters[i]) for i in range(len(samples))}
        
    # crawl through the expression matrix one sample (column) at a time
    for i in range(len(samples)):
        gene_expr_per_sample = {} # new one needed each time sample is incremented
        for n in range(len(genes)):
            gene_expr_per_sample[genes[n]] = expr_matrix[n][i]
        c = cluster_dict[samples[i]] - 1 # get cluster assignment of current sample
        top_genes = ((gene, gene_expr_per_sample) for gene in sorted(gene_expr_per_sample, key=gene_expr_per_sample.get, reverse = True))
        for gene, expr in itertools.islice(top_genes, g): # get how many times a gene appears in the top 200
            gene_freq_per_cluster[c][gene] = gene_freq_per_cluster[c].get(gene, 0) + 1
    
    # sort the most frequence genes in each cluster
    top_cluster_genes = []  # initialize a list of tuples to hold the top genes in a cluster
    for c in range(cluster_number):
        top_genes = [(gene, gene_freq_per_cluster[c][gene]) for gene in sorted(gene_freq_per_cluster[c], key=gene_freq_per_cluster[c].get, reverse = True)]
        top_cluster_genes.append([gene_tuple[0] for gene_tuple in top_genes])
        
    # determine a set of banned genes (patent pending) that appear in more than 1 cluster consensus top 200 (via intersection)
    banned_gene_set = set()
    for c in range(cluster_number):
        for q in range(cluster_number):
            if c != q: # clusters are not being compared to themselves
                for x in range(len(top_cluster_genes[c])):
                    for y in range(len(top_cluster_genes[q])):
                        if top_cluster_genes[c][x] == top_cluster_genes[q][y]:
                            if abs(x-y) < 10: # makes sure each gene is at least ranked 10 positions away in different clusters
                                banned_gene_set.add(top_cluster_genes[c][x]) # ban that sick filth
                                
    # determine the unique gene signatures based on the top cluster genes previously determined
    gene_sigs = []
    for c in range(cluster_number):
        gene_sigs.append([])
    for x in range(g):
        for c in range(cluster_number):
            if top_cluster_genes[c][x] not in banned_gene_set:
                gene_sigs[c].append(top_cluster_genes[c][x])
                
    # determine the minimum number of genes in a gene signature
    min_number = 10000 # use arbitrary high number that will be replaced
    for c in range(cluster_number):
        if len(gene_sigs[c]) < min_number:
            min_number = len(gene_sigs[c])
                
    # return the gene signatures (1 for each cluster)
    for c in range(cluster_number):
        gene_sigs[c] = gene_sigs[c][:min_number] # keeps number symmetric
    return gene_sigs


def print_sigs_to_file(gene_sigs, dir_path, g):
    o = open(dir_path + "sample_cluster_gene_sigs_v2.txt", 'w')
    header = ""   
    for i in range(len(gene_sigs)):
        header += str(i+1) + '\t'
    o.write(header.rstrip() + '\n')
    for n in range(len(gene_sigs[0])):
        line = ""
        for c in range(len(gene_sigs)):
            line += (gene_sigs[c][n] + '\t')
        o.write(line.rstrip() + '\n')
    o.close()

    
def main():
    
    dir_path = "C:/Users/dbray/Dropbox/Boston University/Fall 2015/Challenge Project/TrainingSet/"
    sample_cluster_file = dir_path + "sampK3_clusters.txt"
    norm_expr_file = dir_path + "NormalizedExpMat.txt"
    gene_file = dir_path + "genes4clustering.txt"
    g = 200 # desired size of gene signatures
    
    gene_sigs = get_gene_sig(sample_cluster_file, gene_file, norm_expr_file, g)
    print_sigs_to_file(gene_sigs, dir_path, g)

if __name__ == "__main__":
    main()
