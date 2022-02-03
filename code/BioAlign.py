#!/usr/bin/env python
# coding: utf-8

import copy
import numpy as np
import os
import random
import time
import sys

# ## Stage-1

def get_alignment_stage_1(specie1, specie2):
    
    alignment_list = []
    aligned_nodes = []
    
    f = open("../3D-structure-similarity/"+specie1+"-"+specie2+".tmscore")
    structure = f.readlines()
    f.close()
    
    f = open("../global-sequence-similarity/"+specie1+"-"+specie2+".blast")
    sequence = f.readlines()
    f.close()
    
    f = open("../local-sequence-similarity/"+specie1+"-"+specie2+".local")
    local = f.readlines()
    f.close()
    
    # Top-Quality Nodes (TOP)
    for i in range(len(structure)):
        tokens = (structure[i].strip("\n")).split("\t")
        if float(tokens[-1]) >= 0.8:
            if tokens[0] not in aligned_nodes and tokens[1] not in aligned_nodes:
                aligned_nodes.append(tokens[0])
                aligned_nodes.append(tokens[1])
                alignment_list.append(tokens[0]+"\t"+tokens[1])
                
    for i in range(len(sequence)):
        tokens = (sequence[i].strip("\n")).split("\t")
        if float(tokens[-1]) >= 200:
            if tokens[0] not in aligned_nodes and tokens[1] not in aligned_nodes:
                aligned_nodes.append(tokens[0])
                aligned_nodes.append(tokens[1])
                alignment_list.append(tokens[0]+"\t"+tokens[1])
    
    for i in range(len(local)):
        tokens = (local[i].strip("\n")).split("\t")
        if float(tokens[-1]) >= 4:
            if tokens[0] not in aligned_nodes and tokens[1] not in aligned_nodes:
                aligned_nodes.append(tokens[0])
                aligned_nodes.append(tokens[1])
                alignment_list.append(tokens[0]+"\t"+tokens[1])
                
    top_nodes = copy.deepcopy(alignment_list)

    # Stage-1
    for i in range(len(structure)):
        tokens = (structure[i].strip("\n")).split("\t")
        if float(tokens[-1]) >= 0.5:
            if tokens[0] not in aligned_nodes and tokens[1] not in aligned_nodes:
                aligned_nodes.append(tokens[0])
                aligned_nodes.append(tokens[1])
                alignment_list.append(tokens[0]+"\t"+tokens[1])
                
    for i in range(len(sequence)):
        tokens = (sequence[i].strip("\n")).split("\t")
        if float(tokens[-1]) >= 50:
            if tokens[0] not in aligned_nodes and tokens[1] not in aligned_nodes:
                aligned_nodes.append(tokens[0])
                aligned_nodes.append(tokens[1])
                alignment_list.append(tokens[0]+"\t"+tokens[1])
    
    for i in range(len(local)):
        tokens = (local[i].strip("\n")).split("\t")
        if float(tokens[-1]) >= 2:
            if tokens[0] not in aligned_nodes and tokens[1] not in aligned_nodes:
                aligned_nodes.append(tokens[0])
                aligned_nodes.append(tokens[1])
                alignment_list.append(tokens[0]+"\t"+tokens[1])
                
    return top_nodes, alignment_list

# ## Stage-2 Topology

def get_unique_nodes(specie):

    f = open("../networks/"+specie+".interaction")
    specie_lines = f.readlines()
    f.close()
    specie_nodes=[]

    for i in range(len(specie_lines)):
        tokens = specie_lines[i].split("\t")
        tokens[1] = tokens[1].strip("\n")
        if tokens[0] not in specie_nodes:
            specie_nodes.append(tokens[0])
        if tokens[1] not in specie_nodes:
            specie_nodes.append(tokens[1])
            
    return specie_nodes

def get_remaining_nodes_top(align_nodes,net1,net2):

    for i in align_nodes:
        tokens = i.split("\t")
        if tokens[0] in net1:
            net1.remove(tokens[0])
        if tokens[1] in net2:
            net2.remove(tokens[1])
        
    return net1, net2

def get_neighbors(node, network):
    
    neighbor = []
    for i in network:
        tokens = i.strip("\n").split("\t")
        if node == tokens[0]:
            neighbor.append(tokens[1])
        if node == tokens[1]:
            neighbor.append(tokens[0])
            
    return neighbor

def get_alignment_top(aligned,remaining1,remaining2,specie1,specie2):
    
    f = open("../networks/"+specie1+".interaction")
    network1 = f.readlines()
    f.close()
    
    f = open("../networks/"+specie2+".interaction")
    network2 = f.readlines()
    f.close()
    
    net_a = []
    net_b = []
    for i in aligned:
        tokens = i.strip("\n").split("\t")
        net_a.append(tokens[0])
        net_b.append(tokens[1])
        
    neighbors1 = [] # list of neighbors of all unaligned nodes
    aligned_neighbors1 = [] # list of aligned-neighbors of all unaligned nodes
    
    for i in range(len(remaining1)):
        n = get_neighbors(remaining1[i],network1) # get neighbors of un-aligned nodes
        neighbors1.append(n)
        
        node_aligned_neighbor = []
        for j in n:
            if j in net_a: # neighbors that are aligned
                node_aligned_neighbor.append(j)
        aligned_neighbors1.append(node_aligned_neighbor)
            
    neighbors2 = [] # list of neighbors of all unaligned nodes
    aligned_neighbors2 = [] # list of aligned-neighbors of all unaligned nodes
    
    for i in range(len(remaining2)):
        n = get_neighbors(remaining2[i],network2) # get neighbors of un-aligned nodes
        neighbors2.append(n)
        
        node_aligned_neighbor = []
        for j in n:
            if j in net_b: # neighbors that are aligned
                node_aligned_neighbor.append(j)
        aligned_neighbors2.append(node_aligned_neighbor)
    
    # here we check the count of aligned_neighbors of u 
      #  that are aligned with aligned_neighbors of v
    scores=[]
    candidates_a = []
    candidates_b = []
    for i in range(len(aligned_neighbors1)):
        
        if len(aligned_neighbors1[i]) > 0:
            for j in range(len(aligned_neighbors2)):
                if len(aligned_neighbors2[j]) > 0: 
                    count=0
                    for a in range(len(aligned_neighbors1[i])):
                        for b in range(len(aligned_neighbors2[j])):
                            align_pair = aligned_neighbors1[i][a]+"\t"+aligned_neighbors2[j][b]
                            if align_pair in aligned:
                                count+=1
                    if count>=1: ### threhold on neighbors (r)
                        candidates_a.append(remaining1[i])
                        candidates_b.append(remaining2[j])
                        scores.append(count)

    scores, candidates_a, candidates_b = (list(t) for t in zip(*sorted(zip(scores, candidates_a, candidates_b),reverse=True)))

    list_ = []
    
    for i in range(len(scores)):
        if candidates_a[i] not in list_ and candidates_b[i] not in list_:
            list_.append(candidates_a[i])
            list_.append(candidates_b[i])
            aligned.append(candidates_a[i]+"\t"+candidates_b[i])
            
    return aligned


# ## Stage-2-i Remote Homology

def get_remaining_nodes(net1,net2,align_nodes):

    net_a = []
    net_b = []

    for i in align_nodes:
        tokens = i.split("\t")
        net_a.append(tokens[0])
        net_b.append(tokens[1].strip('\n'))


    f = open("../networks/"+net1+".interaction") # first network
    net1_lines = f.readlines()
    f.close()

    f = open("../networks/"+net2+".interaction") # second network
    net2_lines = f.readlines()
    f.close()

    # getting the remaining nodes
    net1_remaining = []
    for i in net1_lines:
        tokens = i.split()
        if tokens[0] not in net_a and tokens[0] not in net1_remaining:
            net1_remaining.append(tokens[0])
        if tokens[1] not in net_a and tokens[1] not in net1_remaining:
            net1_remaining.append(tokens[1])

    net2_remaining = []
    for i in net2_lines:
        tokens = i.split()
        if tokens[0] not in net_b and tokens[0] not in net2_remaining:
            net2_remaining.append(tokens[0])
        if tokens[1] not in net_b and tokens[1] not in net2_remaining:
            net2_remaining.append(tokens[1])

    return net1_remaining, net2_remaining

def get_values(specie,net1):
    
    net1_proteins = []

    for i in range(len(net1)):
        
        protiens = []
        # get the file
        f = open("../homologs/"+specie+"/"+net1[i]+".txt")
        lines = f.readlines()
        f.close()

        for j in range(len(lines)):
            protiens.append(lines[j].strip("\n"))
        net1_proteins.append(protiens)
    return net1_proteins

def alignment_remote(net1_proteins,net2_proteins,aligned, net1, net2):

    mylist = []
    count=0
    for i in range(len(net1_proteins)):
        for j in range(len(net2_proteins)):
            common = len(set(net1_proteins[i]).intersection(net2_proteins[j]))
            if common > 0:
                count+=1
            mylist.append(common)

    mylist = np.asarray(mylist).reshape((len(net1),len(net2)))

    list = []
    for i in range(mylist.shape[0]):
        index = np.argmax(mylist[i])
        if index not in list: # one index can align only once --> global alignment
            list.append(index)
            my_str = net1[i] + "\t" + net2[index]
            aligned.append(my_str)
            
    return aligned

# ## Stage-2-ii Secondary Structure

def motif_count(string):

    i=-1
    turn_count=0
    loop_count=0

    while (i < len(string)-2):

        start_h = 0
        end_h = 0
        end_loop = 0
        end_turn = 0
        end_h2 = 0
        
        flag_t = 0
        flag_l = 0

        i+=1
        if i+2 > len(string):
            break;

        if string[i] == 'H':
            # might be a start of motif -- check length of Hs
            for j in range(i+1,len(string)):
                if string[j] != 'H':
                    len_h = j-i
                    break

            if len_h > 3: # Hs > 3 can be considerd as Helix
                start_h = i
                end_h = j

                if end_h+2 > len(string):
                    break;
                # here we get H , now find Loop / Turn
                len_loop = 0
                if string[end_h+1] == 'L':# or (string[end_h+1] == 'B' and string[end_h+2] == 'L'):
                    for j in range(end_h+1, len(string)):
                        if string[j] != 'L':
                            len_loop = j - end_h # L start from enh_h and end at j
                            break;

                    # between 2-5 -- > turn
                    if len_loop > 0 and len_loop < 6:
                        end_turn = j
                        if end_turn+2 > len(string):
                            break;
                        if string[end_turn+1] == 'H' and end_turn+2 < len(string):
                            for j in range(end_turn+1, len(string)):
                                if string[j] != 'H':
                                    len_h2 = j - end_turn # L start from enh_h and end at j
                                    break;
                            if len_h2 > 3:
                                end_h2 = j
                                i = end_h2
                                flag_t = 1

                    elif len_loop > 5: # its loop
                        end_loop = j
                        if end_loop+2 > len(string):
                            break;
                        if string[end_loop+1] == 'H':
                            for j in range(end_loop+1, len(string)):
                                if string[j] != 'H':
                                    len_h2 = j - end_loop # L start from enh_h and end at j
                                    break;
                            if len_h2 > 3:
                                end_h2 = j
                                i = end_h2
                                flag_l = 1

            if flag_t != 0:
                turn_count+=1
            elif flag_l !=0:
                loop_count+=1
                
    return turn_count, loop_count

def alignment_sec(net1,net2,specie1,specie2,aligned):

    network1=[]
    network2=[]

    sum_1=[]
    sum_2=[]

    proteins1=[]
    proteins2=[]


    for j in range(len(net1)):

        f = open("../sec/"+specie1+"/"+net1[j]+"_sec.txt")
        linei = str(f.readlines())
        f.close()
        line = []

        for i in range(len(linei)):
            if linei[i] == 'L' or linei[i] == 'T' or linei[i] == 'S':
                line.append('L')
            elif linei[i] == 'H' or linei[i] == 'G' or linei[i] == 'I':
                line.append('H')
            elif linei[i] == 'B' or linei[i] == 'E':
                line.append('B')

        turni, loopi = motif_count(line)
        if turni > 0 or loopi > 0: # atleast one motif exist 
            network1.append((turni,loopi))
            sum_1.append(turni+loopi)
            proteins1.append(net1[j])

    for j in range(len(net2)):
        try:
            f = open("../sec/"+specie2+"/"+net2[j]+"_sec.txt")
            linei = str(f.readlines())
            f.close()
            line = []
            for i in range(len(linei)):
                if linei[i] == 'L' or linei[i] == 'T' or linei[i] == 'S':
                    line.append('L')
                elif linei[i] == 'H' or linei[i] == 'G' or linei[i] == 'I':
                    line.append('H')
                elif linei[i] == 'B' or linei[i] == 'E':
                    line.append('B')
            turni, loopi = motif_count(line)
            if turni > 0 or loopi > 0: # atleast one motif exist 
                network2.append((turni,loopi))
                sum_2.append(turni+loopi)
                proteins2.append(net2[j])

        except:
            continue
    
    # sort on the basis of sum of turn and helix
    sum_1, proteins1, network1 = (list(t) for t in zip(*sorted(zip(sum_1, proteins1, network1),reverse=True)))
    dist = []
    aligni = []
    alignj = []
    for i in range(len(proteins1)):
        dist = []
        for j in range(len(proteins2)):
            dist.append(sum(abs (np.asarray(network1[i]) - np.asarray(network2[j]) )))
        index = np.argmin(dist)

        while(True):
            index = np.argmin(dist)
            if i not in aligni and index not in alignj:
                aligni.append(i)
                alignj.append(index)
                break;
            else:
                dist.remove(dist[index])
                
    for i in range(len(proteins1)):
        aligned.append(proteins1[i]+"\t"+proteins2[i])
        
    return aligned

def main(specie1, specie2, choice="bt"):
    
    print ("***************************************** BioAlign *****************************************")
    print ("You are Aligning", specie1, "Network with", specie2, "Network")
    print ("BioAlign will Use 3D-Strcuture, Global-Sequence, and Local-Sequence Similarities in Stage-1")

    if choice == 'T' or choice == 't':
        print ("Method: Stage-1 + Topology -- Neighbourhood Expansion")

    elif choice == 'B' or choice == 'b':
        print ("Method: Stage-1 + Remote Homology + Secondary Strcuture")

    else:
        print ("Method: Stage-1 + Remote Homology + Secondary Strcuture + Topology")
        
    print ("********************************************************************************************")

    specie1_nodes = get_unique_nodes(specie1)
    specie2_nodes = get_unique_nodes(specie2)
    
    if len(specie1_nodes) <= len(specie2):
        print ("Nodes in Network 1 are greater than Network 2, Swap the Networks")
    
    print ("\nStage-1")
    top_nodes, aligned_stage1 = get_alignment_stage_1(specie1, specie2)
    print ("Top Aligned Nodes", len(top_nodes))
    print ("Nodes Aligned After Stage-1", len(aligned_stage1))
    
    max_length = min(len(specie1_nodes), len(specie2_nodes))
    
    if len(aligned_stage1) < max_length*0.99 and (choice == "t" or choice == "T"): #t
        
        print ("\nStage-2 [Topology -- Neighbourhood Expansion]")
        remaining1, remaining2 = get_remaining_nodes_top(aligned_stage1,specie1_nodes,specie2_nodes)
        aligned_stage2 = get_alignment_top(aligned_stage1,remaining1,remaining2,specie1,specie2)
        print ("Nodes Aligned After Stage-2", len(aligned_stage2))
    
    elif len(aligned_stage1) < max_length*0.99 and (choice == "b" or choice == "B"): #b
        
        print ("\nStage-2-i [Remote Homology]")
        net1,net2 = get_remaining_nodes(specie1,specie2, aligned_stage1)
        net1_proteins = get_values(specie1,net1)
        net2_proteins = get_values(specie2,net2)
        aligned_stage2 = alignment_remote(net1_proteins,net2_proteins,aligned_stage1,net1,net2)
        print ("Nodes Aligned After Stage-2-i: ",len(aligned_stage2))
        
        print ("\nStage-2-ii [Secondary Structure]")
        net1,net2 = get_remaining_nodes(specie1,specie2, aligned_stage2)
        aligned_stage2 = alignment_sec(net1,net2,specie1,specie2,aligned_stage2)
        print ("Nodes Aligned After Stage-2-ii: ",len(aligned_stage2))

    elif len(aligned_stage1) < max_length*0.99: #bt

        print ("\nStage-2-i [Remote Homology]")
        net1,net2 = get_remaining_nodes(specie1,specie2, aligned_stage1)
        net1_proteins = get_values(specie1,net1)
        net2_proteins = get_values(specie2,net2)
        aligned_stage2 = alignment_remote(net1_proteins,net2_proteins,aligned_stage1,net1,net2)
        print ("Nodes Aligned After Stage-2-i: ",len(aligned_stage2))
        
        print ("\nStage-2-ii [Secondary Structure]")
        net1,net2 = get_remaining_nodes(specie1,specie2, aligned_stage2)
        aligned_stage2 = alignment_sec(net1,net2,specie1,specie2,aligned_stage2)
        print ("Nodes Aligned After Stage-2-ii: ",len(aligned_stage2))


        print ("\nStage-2-iii [Topology -- Neighbourhood Expansion]")
        remaining1, remaining2 = get_remaining_nodes_top(aligned_stage2,specie1_nodes,specie2_nodes)
        aligned_stage2 = get_alignment_top(aligned_stage2,remaining1,remaining2,specie1,specie2)
        print ("Nodes Aligned After Stage-2", len(aligned_stage2))

    else:

        aligned_stage2 = copy.deepcopy(aligned_stage1)
        print ("\nAlignment is Completed in Stage-1")
        
    f = open("../alignments/"+specie1+"-"+specie2+"-"+choice+".alignment",'w')
    for i in range(len(aligned_stage2)):
        f.write(aligned_stage2[i]+"\n")
    f.close()

specie1 = sys.argv[1]
specie2 = sys.argv[2]
choice = sys.argv[3]

a = main(specie1,specie2,choice)

