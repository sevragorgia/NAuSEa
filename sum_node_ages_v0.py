#!/usr/bin/env python

#Copyright 2018 Sergio Vargas
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys

import numpy as np

from collections import defaultdict

import dendropy
from dendropy.calculate import treecompare




age_outfile = open("/mnt/localspace/s.vargas/Data/Synechococcus_Phylogenomics/Alignments/Geneious/Reduced_dataset/BEAST2/fixed_topology/Correct_node_constraints/20170531_All_Cyanobacteria_Phylogenomics_04_manually_reduced_AllRuns_Node_Ages.csv", "w")
edge_outfile = open("/mnt/localspace/s.vargas/Data/Synechococcus_Phylogenomics/Alignments/Geneious/Reduced_dataset/BEAST2/fixed_topology/Correct_node_constraints/20170531_All_Cyanobacteria_Phylogenomics_04_manually_reduced_AllRuns_Edge_Lengths.csv", "w")

taxa = dendropy.TaxonNamespace()

guide_tree = dendropy.Tree.get_from_path("/mnt/localspace/s.vargas/Data/Synechococcus_Phylogenomics/Alignments/Geneious/Reduced_dataset/BEAST2/fixed_topology/Correct_node_constraints/target_Cyano_phylo.tre", "nexus", taxon_namespace=taxa, rooting="default-rooted")
guide_tree.encode_bipartitions()

guide_tree.calc_node_ages()#this will add an age attribute to the guide tree, we will turn this attribute into a list

node_ages = defaultdict(list)
node_subtending_edge_lengths = defaultdict(list)

source_files = [open("/mnt/localspace/s.vargas/Data/Synechococcus_Phylogenomics/Alignments/Geneious/Reduced_dataset/BEAST2/fixed_topology/Correct_node_constraints/All_Cyanobacteria_Phylogenomics_04_manually_reduced_run1.trees", "r")] # Note: for 'Tree.yield_from_files',

tree_yielder = dendropy.Tree.yield_from_files(files=source_files,schema='nexus',taxon_namespace=taxa, rooting="default-rooted")

tree_index = 0

for mcmc_tree in tree_yielder:
  
  mcmc_tree.encode_bipartitions()
  
  if treecompare.symmetric_difference(guide_tree, mcmc_tree) == 0:#if tree is equal to the guide tree, which should be the case in our case
    
    mcmc_tree.calc_node_ages() #adds an atribute age to each node in the tree
    mcmc_node_number = 0
    
    for mcmc_node in mcmc_tree.postorder_node_iter():#not sure if this for is required
      #print the node age in postorder to the outfile
      age_outline = str(mcmc_node.age) + "\t"
      age_outfile.write(age_outline)
      
      edge_outline = str(mcmc_node.edge_length) + "\t"
      edge_outfile.write(edge_outline)
      
      #now store the value in a list in memory
      node_ages[mcmc_node_number].append(mcmc_node.age)
      
      #need to summarize branch length as well to get the mean/media branch length based on the mcmc chain results...
      node_subtending_edge_lengths[mcmc_node_number].append(mcmc_node.edge_length)
      
      mcmc_node_number += 1
    
    age_outfile.write("\n")
    edge_outfile.write("\n")
    tree_index += 1
  
  print(tree_index)
  

####################################################################
####################################################################


node_number = 0

for guide_node in guide_tree.postorder_node_iter():
  
  age_outfile.write(str(node_number) + "\t")
  edge_outfile.write(str(node_number) + "\t")
  
  guide_node.annotations.add_new(name="&node_number", value=str(node_number))
  
  #calculate mean node age
  mean_node_age = np.mean(node_ages[node_number])
  guide_node.annotations.add_new(name="&mean_node_age", value=str(mean_node_age))
  
  #calculate stdev of node age
  stdev_node_age = np.std(node_ages[node_number])
  stdev_str = "{" + str(mean_node_age+stdev_node_age) + "," + str(mean_node_age-stdev_node_age) + "}"
  guide_node.annotations.add_new(name="&stdev_node_age", value=stdev_str)
  
  #calculate median of node age
  median_node_age = np.median(node_ages[node_number])
  guide_node.annotations.add_new(name="&median_node_age", value=str(median_node_age))
  
  #calculate mean edge lengths
  mean_subtending_edge_length = np.mean(node_subtending_edge_lengths[node_number])
  guide_node.edge.length = mean_subtending_edge_length
  
  node_number += 1
  

age_outfile.write("\n")
edge_outfile.write("\n")
guide_tree.write(path="/mnt/localspace/s.vargas/Data/Synechococcus_Phylogenomics/Alignments/Geneious/Reduced_dataset/BEAST2/fixed_topology/Correct_node_constraints/Annotated_summary_tree.tre", schema='nexus')


####################################################################
####################################################################

age_outfile.close()
edge_outfile.close()

sys.exit(0)
