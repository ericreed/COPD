# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 17:43:27 2015

@author: dbray
"""

def reformat_entrez(in_file, out_file):
    """ reformats an entrez gene ID to be compatible with GSEA """
    
    # open the input file
    f = open(in_file, 'r')
    o = open(out_file, 'w')
    
    # iterate through the lines in the input file
    while True:
        line = f.readline() # read a line
        if not line:
            break # break out of the file at the end
        line = line.replace('_at', '') # strip the '_at' suffix
        o.write(line) # write the line to file
        
    # close the files
    f.close()
    o.close()
    
    return 1
    
def main():
    
    # set a current directory
    curr_dir = "C:/Users/dbray/Google Drive/COPD_project/FormattedData/GSEA_inputs/"
    
    # define the input file
    in_file = curr_dir + "pooledExpMat_GSEA.txt"
    
    # define an output file
    out_file = curr_dir + "pooledExpMat_GSEA_entrez.txt"
    
    # pass the directory, files to the script
    reformat_entrez(in_file, out_file)
    
    
if __name__ == "__main__":
    main()
