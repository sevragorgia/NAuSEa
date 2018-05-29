#!/usr/bin/env python3

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

import argparse

import dendropy
#from dendropy.calculate import treecompare

###################
#
# set the argument list.
#
###################################


parser = argparse.ArgumentParser()

parser.add_argument("tree_file", help="file containing trees sampled during the MCMC")
parser.add_argument("-b","--beast", help="beast MCMC", action="store_true")

args = parser.parse_args()

######################
# 
# Read the trees
# 
############################################

mcmc_trees = dendropy.TreeList()

if args.beast:
  
  mcmc_trees.read_from_path(args.tree_file, "nexus")
  
else:
  
  mcmc_trees.read_from_path(args.tree_file, "newick")
  

######
#
# loop through trees
#
########################

tree_index = 0

for tree in mcmc_trees:
  
  print(tree_index)

  mcmc_node_number = 0

  for mcmc_node in tree.postorder_node_iter():#not sure if this for is required
    #print the node age in postorder to the outfile
    age_outline = str(mcmc_node.age) + "\t"
    print(age_outline)
#    age_outfile.write(age_outline)
    edge_outline = str(mcmc_node.edge_length) + "\t"
    print(edge_outline)
#    edge_outfile.write(edge_outline)
      
    #now store the value in a list in memory
#    node_ages[mcmc_node_number].append(mcmc_node.age)
      
    #need to summarize branch length as well to get the mean/media branch length based on the mcmc chain results...
#    node_subtending_edge_lengths[mcmc_node_number].append(mcmc_node.edge_length)
      
    mcmc_node_number += 1
    
#    age_outfile.write("\n")
#    edge_outfile.write("\n")
  
  tree_index += 1
  


exit(0)
