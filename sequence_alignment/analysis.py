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
def generate_contact(inputfile, outputfile):
    input = pd.read_table(inputfile, names = ["uniprot_ID", "Histone_name", "binding_site"])
    uniprot = input.iloc[:,0:1]
    uniprot = np.unique(uniprot)
    f = open(outputfile,'a')
  
    print(convert(list)) 
    for ID in uniprot:
        his =[]
        site = []
        for i in range(0,len(input)):
            if (input.iloc[i,0] == ID):
                his.append(input.iloc[i,1])
                site.append(input.iloc[i,2])
        his_counts = np.unique(his,return_counts=True)
        site_counts = np.unique(site,return_counts=True)
        f.write(ID + "\t" + convert(his_counts[0]) + "\t" + convert(site_counts[0]) + "\t" + convert(site_counts[1]) + "\n")
    f.close

#generate_contact("H1.txt", "H1_contact.txt")
#generate_contact("H2.txt", "H2_contact.txt")
#generate_contact("H3.txt", "H3_contact.txt")
#generate_contact("H4.txt", "H4_contact.txt")

def generate_annotation(inputfile, outputfile):
    
    input = pd.read_table(inputfile, names = ["uniprot_ID", "Histone_name", "binding_site", "counts"])  
    f = open(outputfile,"a")
    f.write("JALVIEW_ANNOTATION\n")  
    f.write("\n")     
    for i in range(0,len(input)):
        list_null = [0]* 400 
        position = input.iloc[i,2].split("|")
        count = input.iloc[i,3].split("|")
        for j in range(0,len(position)):
            index = position[j].split("_")[1]
            list_null[int(index)-1] = count[j]
        contacts_list = convert(list_null)
        f.write("SEQUENCE_REF" + "\t" + input.iloc[i,0] + "\n")
        f.write("BAR_GRAPH" + "\t" + "Contacts" + "\t" + "Contacts in " + input.iloc[i,0] \
                + "\t" + contacts_list + "\n")
        f.write("\n")    
    f.close
        
generate_annotation("H1_contact.txt", "H1_annotation.txt") 
generate_annotation("H2_contact.txt", "H2_annotation.txt")           
generate_annotation("H3_contact.txt", "H3_annotation.txt")           
generate_annotation("H4_contact.txt", "H4_annotation.txt")           
          
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
 
            
