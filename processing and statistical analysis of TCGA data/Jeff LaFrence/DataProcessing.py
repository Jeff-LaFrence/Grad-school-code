# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 22:12:26 2021

@author: golde
"""
HiSeqV2 = []
microarray = []
ov_tcga_clinical_data = []

HiSeqV2file = open('HiSeqV2','r')


for line in HiSeqV2file:
    line = line.rstrip()
    HiSeqV2.append(line)
    
HiSeqV2file.close()
    

microarrayfile = open('HT_HG-U133A','r')

for line in microarrayfile:
    line = line.rstrip()
    microarray.append(line)  

microarrayfile.close()


ov_tcga_clinical_data_file = open('ov_tcga_clinical_data.tsv','r')

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
            
#generate group statuses
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
           

groupIDdict = {}

for (sample, values) in enumerate(Sample_ID):
    groupIDdict[values] = groupID[sample]

##First Pathologic Diagnosis Biospecimen Acquisition Other Method Type occurs twice in ov_TCGA##

#pull 
# HiSeqV2data = []

# for (rownum, rowname) in enumerate(HiSeqV2):
#     #take rowname[1] (header), split to list
#     if (rownum == 0):  
#         HiSeqV2header = rowname.split("\t")
#     if (rownum > 0): 
#         HiSeqV2data.append(rowname)
        
Hiseqdict = {}

for (rownum, rowname) in enumerate(HiSeqV2):
    #take rowname[1] (header), split to list
    if (rownum == 0):  
        hiseqheader = rowname.split("\t")
        for (index_num, index_val) in enumerate(hiseqheader):
            Hiseqdict[index_num] = [index_val]
    if (rownum > 0):
        hiseqheader = rowname.split("\t")
        for (index_num, index_val) in enumerate(hiseqheader):
            Hiseqdict[index_num].append(index_val)
            
            
final_file_hiseq = []
currenthiseqlist = []
for key in Hiseqdict:
    currenthiseqlist = Hiseqdict[key]
    if key == 0:
        currenthiseqlist.insert(1,"group_ID")
        final_file_hiseq.append(currenthiseqlist)
    if currenthiseqlist[0] in Sample_ID:
        #pull group ID for sample_ID
        currenthiseqlist.insert(1,groupIDdict[currenthiseqlist[0]])
        final_file_hiseq.append(currenthiseqlist)
        
SeqData_file = open("SeqData.txt", "w")
for line in final_file_hiseq:
    print(*line, sep = "\t" , end="\n", file = SeqData_file)
SeqData_file.close()
        
# ArrayData_file = open("ArrayData.txt", "w")
# print(*microarrayheadercut, end="\n", sep = "\t", file = ArrayData_file)
# print(*group_ID_micro, end="\n", sep = "\t", file = ArrayData_file)
# print(*microarraydata, end="\n", sep = "\n", file = ArrayData_file)
# ArrayData_file.close() 



microarraydict = {}

for (rownum, rowname) in enumerate(microarray):
    #take rowname[1] (header), split to list
    if (rownum == 0):  
        microarrayheader = rowname.split("\t")
        for (index_num, index_val) in enumerate(microarrayheader):
            microarraydict[index_num] = [index_val]
    if (rownum > 0):
        microarrayheader = rowname.split("\t")
        for (index_num, index_val) in enumerate(microarrayheader):
            microarraydict[index_num].append(index_val)
            
final_file_array = []
currentarraylist = []
for key in microarraydict:
    currentarraylist = microarraydict[key]
    if key == 0:
        currentarraylist.insert(1,"group_ID")
        final_file_array.append(currentarraylist)
    if currentarraylist[0] in Sample_ID:
        #pull group ID for sample_ID
        currentarraylist.insert(1,groupIDdict[currentarraylist[0]])
        final_file_array.append(currentarraylist)
        
ArrayData_file = open("ArrayData.txt", "w")
for line in final_file_array:
    print(*line, end="\n", sep = "\t", file = ArrayData_file)
ArrayData_file.close()


# for (rownum, rowname) in enumerate(microarray):
#     #take rowname[1] (header), split to list
#     if (rownum == 0):  
#         microarrayheader = rowname.split("\t")
#     if (rownum > 0):
#         #make into list of lists
#         microdatanames = rowname.split("\t")
#         microarraydata.append(microdatanames)

# for (rownum, rowname) in enumerate(microarraydata):
#         listname = str(rownum)
#         vars()[listname]  = []
#         vars()[listname].append(rowname)    

    
    
        
###
# #pull column from ov_tcga_clinical_data_file containing patient ID's
# for (rownum, rowname) in enumerate(microarray):
#     #take rowname[1] (header), split to list
#     if (rownum == 0):
#         microarrayheader = rowname.split("\t")
#         for i in microarrayheader:
#             listname = i
#             vars()[listname]  = []
#     #pull non-header rows        
#     if (rownum > 0):
#         #split each row
#         current_row = rowname.split("\t")
#         for (current_num, current_dat) in enumerate(current_row):
#             listname = current_num
#             vars()[listname].append(current_dat)
###       
 #need to leave out if not in ov_tcga
#this means processing in and out entire data set
###split, grab rows of interest, move to next one###


# for (rownum, rowname) in enumerate(microarray):
#     #take rowname[1] (header), split to list
#     if (rownum == 0):
#         microarrayheader = rowname.split("\t")
#     #pull non-header rows        
#     if (rownum > 0):
#         #split each row
#         current_row = rowname.split("\t")
#         listname = str(rownum)
#         vars()[listname].append(current_row)    

        
####       
        
        
   
        
        
#make row variable groupID for HVseq
# group_ID_seq = []
# group_IDnum_seq = []
# group_IDDFM_seq = []
# group_IDDFS_seq = []
# clindat_sample_seq = []
# hiseq_sample_seq = []

# for (sampleIDnum, sampleIDname) in enumerate(HiSeqV2header):
#     if (sampleIDname == "sample"):
#         group_ID_seq.append("group_ID")
#     else:
#         for (clinnum, clinname) in enumerate(Sample_ID):
#             if sampleIDname == clinname:
#                 group_ID_seq.append(groupID[clinnum])
#                 #analytics
#                 group_IDnum_seq.append(clinnum)
#                 group_IDDFM_seq.append(Disease_Free_Months[clinnum])
#                 group_IDDFS_seq.append(Disease_Free_Status[clinnum])
#                 clindat_sample_seq.append(clinname)
#                 hiseq_sample_seq.append(sampleIDname)
                
# SeqData_file = open("SeqData.txt", "w")
# print(*HiSeqV2header, end="\n", sep = "\t", file = SeqData_file)
# print(*group_ID_seq, end="\n", sep = "\t", file = SeqData_file)
# print(*HiSeqV2data, end="\n", sep = "\n", file = SeqData_file)
# SeqData_file.close() 



# #make row variable groupID for clindataseq

# group_ID_micro = []
# group_IDnum_micro = []
# group_IDDFM_micro = []
# group_IDDFS_micro = []
# clindat_sample_micro = []
# hiseq_sample_micro = []
# sample_ID_micro = []
# microarrayheaderprint = []
# microarraydataprint = []
# temparraydata = []
# microarrayheadercut = []


# for (sampleIDnum, sampleIDname) in enumerate(microarrayheader):
#     #set rowname for first column   
#     if (sampleIDname == "sample"):
#         group_ID_micro.append("group_ID")
#     else:
#         for (clinnum, clinname) in enumerate(Sample_ID):
#             if sampleIDname == clinname:
#                  #microarraydata   

# #same as hiseqv2
# for (sampleIDnum, sampleIDname) in enumerate(microarrayheader):
#     #set rowname for first column   
#     if (sampleIDname == "sample"):
#         group_ID_micro.append("group_ID")
#     else:
#         for (clinnum, clinname) in enumerate(Sample_ID):
#             if sampleIDname == clinname:
#                   #grab relevant
#                 group_ID_micro.append(groupID[clinnum])
#                     # #grab specific items from list for new list
#                     # arrayrowinteration.append(i[sampleIDnum])
#                     # group_ID_micro.append(groupID[clinnum])
#                     #analytics
#                 group_IDnum_micro.append(clinnum)
#                 group_IDDFM_micro.append(Disease_Free_Months[clinnum])
#                 group_IDDFS_micro.append(Disease_Free_Status[clinnum])
#                 clindat_sample_micro.append(clinname)
#                 hiseq_sample_micro.append(sampleIDname)


                
           
                
# ArrayData_file = open("ArrayData.txt", "w")
# print(*microarrayheadercut, end="\n", sep = "\t", file = ArrayData_file)
# print(*group_ID_micro, end="\n", sep = "\t", file = ArrayData_file)
# print(*microarraydata, end="\n", sep = "\n", file = ArrayData_file)
# ArrayData_file.close() 

# group_IDnum_micro.sort()
        


