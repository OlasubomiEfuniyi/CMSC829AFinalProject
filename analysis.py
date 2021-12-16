# This script can be used to find the percentage of edges in a tree
# with bootstrap support values >= a provided minimum.

# It can also be used to compare the bipartitions of two trees.

# Only newick tree files are supported.

import dendropy
import sys

# Returns the percentage of edges with bootstrap support >= minimum
def bootstrapPercentageAnalysis(path_to_file="", minimum=0):
    tree = dendropy.Tree.get(path=path_to_file,schema="newick") 
    count = 0
    count_internal = 0

    # Bootstrap support values are stored as internal node labels
    for node in tree.nodes():
        if node.level() != 0 and node.is_internal(): # The imaginary root is considere an internal node so we avoid it.
            count_internal += 1
            assert(node.label != None) # Every internal node except from the imaginary root should have a bootstrap branch support value labelling it
            if int(node.label) >= minimum:
                count += 1
            


    if count_internal > 0:
        return (count/count_internal) * 100
    else:
        return 0

# Compares two trees to see how many non trivial bipartitions they have in common.
# Returns a tuple whose first item is the count of non trivial bipartitions in common,
# the second item is the number of non trivial bipartitions in the first tree,
# and the third item is the number of non trivial bipartitions in the second tree.
def treeComparison(path_to_file_1="", path_to_file_2=""):
    # Use the same taxon namespace for both trees so that the interpretations of bipartition bitmasks can be consistent
    tree1 = dendropy.Tree.get(path=path_to_file_1,schema="newick")
    tree2 = dendropy.Tree.get(path=path_to_file_2,schema="newick", taxon_namespace=tree1.taxon_namespace) 

    tree1_bipartitions = list(filter(lambda bip: not bip.is_trivial(), tree1.encode_bipartitions()))
    tree2_bipartitions = list(filter(lambda bip: not bip.is_trivial(), tree2.encode_bipartitions()))

    count = 0

    for bipartition1 in tree1_bipartitions: # This collapses unrooted basal bifurcations by default, thereby getting rid fo the imaginary root
            bitstring1 = bipartition1.leafset_as_bitstring()
            for bipartition2 in tree2_bipartitions:
                bitstring2 = bipartition2.leafset_as_bitstring()
                if bitstring1 == bitstring2:
                    count += 1
    
    return (count, len(tree1_bipartitions), len(tree2_bipartitions))
            




if len(sys.argv) != 3:
    print("Usage: python3 analysis.py <path_to_file> <minimum_bootstrap_value>")
    print("Usage: python3 analysis.py <path_to_file> <path_to_file>")
    exit()
else:
    if sys.argv[2].isnumeric(): # Then this could only work for bootstrap analysis per usage
        minimum = int(sys.argv[2])
        result = bootstrapPercentageAnalysis(path_to_file = sys.argv[1], minimum = minimum)
        print(str(result) + "% of internal nodes have a bootstrap branch support of " + str(minimum) + " or more.")
    else:
        (count_in_common, num_non_trivial_bip_tree1, num_non_trivial_bip_tree2) = treeComparison(sys.argv[1], sys.argv[2])
        print("Both trees have " + str(count_in_common) + " non-trivial bipartitions in common.")
        print("First tree has " + str(num_non_trivial_bip_tree1) + " non-trivial biparitions.")
        print("Second tree has " + str(num_non_trivial_bip_tree2) + " non-trivial bipartitions.")

