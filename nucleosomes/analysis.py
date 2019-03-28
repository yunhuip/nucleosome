# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os
import numpy as np
import pandas as pd 

def convert(list): 
      
    # Converting integer list to string list 
    s = [str(i) for i in list] 
      
    # Join list items using join() 
    res = str("|".join(s)) 
      
    return(res) 
# load files 
def generate_contact(inputfile, outputfile,histone_type, domain = "False"):
    input = pd.read_table(inputfile)
    uniprot = input.iloc[:,0:1]
    uniprot = np.unique(uniprot)
    f = open(outputfile,'a')
  
    for ID in uniprot:
        his =[]
        site = []
        organism = []
        for i in range(0,len(input)):
            if domain == "True":
                if (input.iloc[i,0] == ID) and (input.iloc[i,8] == histone_type) and (input.isnull().iloc[i,6] != True):
                    his.append(input.iloc[i,1])
                    organism.append(input.iloc[i,2])
                    site.append(input.iloc[i,9])                
            else:
                if (input.iloc[i,0] == ID) and (input.iloc[i,8] == histone_type):
                    his.append(input.iloc[i,1])
                    organism.append(input.iloc[i,2])
                    site.append(input.iloc[i,9])
        his_counts = np.unique(his,return_counts=True)
        organism_counts = np.unique(organism,return_counts=True)
        site_counts = np.unique(site,return_counts=True)
        if convert(site_counts[0]) != "":
            f.write(ID + "\t" + convert(his_counts[0]) + "\t"+ convert(organism_counts[0]) + "\t" + convert(site_counts[0]) + "\t" + convert(site_counts[1]) + "\n")
    f.close

generate_contact("all_nucleosomes_results.tsv", "H2A_contact_combine.txt", "H2A")
generate_contact("all_nucleosomes_results.tsv", "H2A_contact_domain.txt", "H2A", "True")

generate_contact("all_nucleosomes_results.tsv", "H2B_contact_combine.txt", "H2B")
generate_contact("all_nucleosomes_results.tsv", "H2B_contact_domain.txt", "H2B", "True")

generate_contact("all_nucleosomes_results.tsv", "H3_contact_combine.txt", "H3")
generate_contact("all_nucleosomes_results.tsv", "H3_contact_domain.txt", "H3", "True")

generate_contact("all_nucleosomes_results.tsv", "H4_contact_combine.txt", "H4")
generate_contact("all_nucleosomes_results.tsv", "H4_contact_domain.txt", "H4", "True")


def generate_annotation(inputfile, outputfile):
    
    input = pd.read_table(inputfile, names = ["uniprot_ID", "Histone_name","organism", "binding_site", "counts"])  
    f = open(outputfile,"a")
    f.write("JALVIEW_ANNOTATION\n")  
    f.write("\n")     
    for i in range(0,len(input)):
        list_null = [0]* 400 
        position = input.iloc[i,3].split("|")
        count = input.iloc[i,4].split("|")
        for j in range(0,len(position)):
            list_null[int(position[j])-1] = count[j]
        contacts_list = convert(list_null)
        f.write("SEQUENCE_REF" + "\t" + input.iloc[i,0].upper() + "\t" + "0\n")
        f.write("BAR_GRAPH" + "\t" + "Contacts" + "\t" + "Contacts in " + input.iloc[i,0].upper() + "_"+ input.iloc[i,1] \
                + "\t" + contacts_list + "\n")
        f.write("\n")    
    f.close
        
generate_annotation("H2A_contact_combine.txt", "H2A_annotation_combine.txt") 
generate_annotation("H2A_contact_domain.txt", "H2A_annotation_domain.txt") 

generate_annotation("H2B_contact_combine.txt", "H2B_annotation_combine.txt") 
generate_annotation("H2B_contact_domain.txt", "H2B_annotation_domain.txt") 

generate_annotation("H3_contact_combine.txt", "H3_annotation_combine.txt") 
generate_annotation("H3_contact_domain.txt", "H3_annotation_domain.txt") 

generate_annotation("H4_contact_combine.txt", "H4_annotation_combine.txt") 
generate_annotation("H4_contact_domain.txt", "H4_annotation_domain.txt")           
        
          
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
 
            
