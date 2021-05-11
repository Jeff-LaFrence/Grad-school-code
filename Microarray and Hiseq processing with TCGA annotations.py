# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 22:12:26 2021

@author: golde
"""

##initialize open structure for data
HiSeqV2 = []
microarray = []
ov_tcga_clinical_data = []

##read in HiSeq data
HiSeqV2file = open('HiSeqV2','r')


for line in HiSeqV2file:
    line = line.rstrip()
    HiSeqV2.append(line)
    
HiSeqV2file.close()
    
##read in microarray data
microarrayfile = open('HT_HG-U133A','r')

for line in microarrayfile:
    line = line.rstrip()
    microarray.append(line)  

microarrayfile.close()

## read in TCGA data
ov_tcga_clinical_data_file = open('ov_tcga_clinical_data.tsv','r')

#format data, remocving special characters to make coding easier
for (row, line) in enumerate(ov_tcga_clinical_data_file):
    line = line.rstrip()
    if row == 0:
        line = line.replace(" ","_")
        line = line.replace("(","")
        line = line.replace(")","")
    ov_tcga_clinical_data.append(line)

ov_tcga_clinical_data_file.close()


#pull column from ov_tcga_clinical_data_file containing patient ID's
for (rownum, rowname) in enumerate(ov_tcga_clinical_data):
    #take rowname[1] (header), split to list
    if (rownum == 0):
        clindataheader = rowname.split("\t")
        for i in clindataheader:
            listname = i
            vars()[listname]  = []
    #pull non-header rows        
    if (rownum > 0):
        #split each row
        current_row = rowname.split("\t")
        for (current_num, current_dat) in enumerate(current_row):
            listname = clindataheader[current_num]
            vars()[listname].append(current_dat)
            
#generate group statuses, and report whicg loop was used to make call
groupID = []

for (stat_integer,this_stat) in enumerate(Overall_Survival_Status):
    if this_stat == "NA":
        groupID.append("NA")
        print("loop 1")
    elif this_stat == "1:DECEASED":
        print("loop 2")
        if Overall_Survival_Months[stat_integer] == "NA":
           groupID.append("NA")
           print("loop 2.1")
        elif int(float(Overall_Survival_Months[stat_integer])) <= 36:
           groupID.append("1")
           print("loop 2.2")
        else:
           groupID.append("2")
           print("loop 2.3")
    else:
        print("loop 3")
        if Overall_Survival_Months[stat_integer] == "NA":
           groupID.append("NA")
           print("loop 3.1")
        elif int(float(Overall_Survival_Months[stat_integer])) <= 36:
            groupID.append("NA")
            print("loop 3.2")
        else:
           groupID.append("2")
           print("loop 2.3")
           

#create dictionary to add group Ids to
groupIDdict = {}

for (sample, values) in enumerate(Sample_ID):
    groupIDdict[values] = groupID[sample]

##format Hiseq data into dictionary of lists
Hiseqdict = {}

for (rownum, rowname) in enumerate(HiSeqV2):
    #take rowname[1] (header), split to list
    if (rownum == 0):  
        hiseqheader = rowname.split("\t")
        for (index_num, index_val) in enumerate(hiseqheader):
            Hiseqdict[index_num] = [index_val]
    #make dictionary of lists for the rest of data table
    if (rownum > 0):
        hiseqheader = rowname.split("\t")
        for (index_num, index_val) in enumerate(hiseqheader):
            Hiseqdict[index_num].append(index_val)
            
#go through dictionary, reformat to matrix by grabbing each list from dictionary
final_file_hiseq = []
currenthiseqlist = []
for key in Hiseqdict:
    currenthiseqlist = Hiseqdict[key]
    #make first list into header, add line for group ID
    if key == 0:
        currenthiseqlist.insert(1,"group_ID")
        final_file_hiseq.append(currenthiseqlist)
    #add remaining lines to matrix, add groupId notation to each
    if currenthiseqlist[0] in Sample_ID:
        #pull group ID for sample_ID
        currenthiseqlist.insert(1,groupIDdict[currenthiseqlist[0]])
        final_file_hiseq.append(currenthiseqlist)

#write HiSeq with Id annotations data to file       
SeqData_file = open("SeqData.txt", "w")
for line in final_file_hiseq:
    print(*line, sep = "\t" , end="\n", file = SeqData_file)
SeqData_file.close()
        

microarraydict = {}

for (rownum, rowname) in enumerate(microarray):
    #take rowname[1] (header), split to list
    if (rownum == 0):  
        microarrayheader = rowname.split("\t")
        for (index_num, index_val) in enumerate(microarrayheader):
            microarraydict[index_num] = [index_val]
    #make dictionary of lists for the rest of data table
    if (rownum > 0):
        microarrayheader = rowname.split("\t")
        for (index_num, index_val) in enumerate(microarrayheader):
            microarraydict[index_num].append(index_val)

#go through dictionary, reformat to matrix by grabbing each list from dictionary            
final_file_array = []
currentarraylist = []
for key in microarraydict:
    currentarraylist = microarraydict[key]
    #make first list into header, add line for group ID
    if key == 0:
        currentarraylist.insert(1,"group_ID")
        final_file_array.append(currentarraylist)
    #add remaining lines to matrix, add groupId notation to each
    if currentarraylist[0] in Sample_ID:
        #pull group ID for sample_ID
        currentarraylist.insert(1,groupIDdict[currentarraylist[0]])
        final_file_array.append(currentarraylist)

#write array with Id annotations data to file        
ArrayData_file = open("ArrayData.txt", "w")
for line in final_file_array:
    print(*line, end="\n", sep = "\t", file = ArrayData_file)
ArrayData_file.close()





