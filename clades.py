# This script operates in multiple modes
# Usage 1: python3 clades.py -b <path_to_tree_file> outputs the percentage of bootsrap values that are 75 or more assuming a valid Newick tree file.
# Usage 2: python3 clades.py -c <path_to_tree_file_1> <path_to_tree_file_2> outputs the list of the clades both tree files have in common, with
#          their bootstrap support values. Assumes both tree files were built on the same taxas.


import sys

def isUpperCaseAlphabet(val):
    if (val == "A" or val == "B" or val == "C" or val == "D" or val == "E" or val == "F" or val == "G" or val == "H" or val == "I" or val == "J" or
        val == "K" or val == "L" or val == "M" or val == "N" or val == "O" or val == "P" or val == "Q" or val == "R" or val == "S" or val == "T" or
        val == "U" or val == "V" or val == "W" or val == "X" or val == "Y" or val == "Z"):
        return True
    else:
        return False

def isDigit(d):
    if d == "0" or d == "1" or d == "2" or d == "3" or d == "4" or d == "5" or d == "6" or d == "7" or d == "8" or d == "9":
        return True
    else:
        return False

# Given a map from clade ids to their bootstrap values, return the percentage of the clades that have a bootsrap support value >= 75
def percentAbove(bootstrapMap):
    count = 0
    for k in bootstrapMap:
        if (bootstrapMap[k] >= 75):
            count = count + 1
    return (count/len(bootstrapMap) * 100)

# Given a map from clade ids to their bootstrap values and a list of clade ids to consider, return the percentage of the clades that have a bootstrap 
# support value >= 75
def percentageAboveWithCladeIds(bootstrapMap, cladeIds):
    count = 0
    cladeIds = set(cladeIds)

    for k in bootstrapMap:
        if k in cladeIds and bootstrapMap[k] >= 75:
            count = count + 1
    return (count/len(cladeIds) * 100)

def jointPercentageAbove(bootstrapMap1, bootstrapMap2, cladeIds1, cladeIds2):
    assert(len(cladeIds1) == len(cladeIds2))
    count = 0
    for i in range(0,len(cladeIds1)):
        if bootstrapMap1[cladeIds1[i]] >= 75 and bootstrapMap2[cladeIds2[i]] >= 75:
            count = count + 1
    return (count/len(cladeIds1) * 100)

def usageError():
    print("Usage 1: python3 clades.py -b <path_to_tree_file> outputs the percentage of bootsrap values that are 75 or more assuming a valid Newick tree file.")
    print("""Usage 2: python3 clades.py -c <path_to_tree_file_1> <path_to_tree_file_2> outputs the list of the clades both tree files have in common, with 
    their bootstrap support values. Assumes both tree files were built on the same taxas.""")
    exit()

# Returns true if list l1 is equal to list l2. Returns false otherwise
def listsEqual(l1, l2):
    l1set = set(l1)
    l2set = set(l2)

    if len(l1set) == len(l2set):
        for v in l1set:
            if not v in l2set:
                return False
        
        return True
    else:
        return False

    
# Given a file object, extract the clades and return a tuple whose
# first item is a map from clade ids to a list of taxa in the clade,
# whose second item is a map from clade ids to the bootstrap support value for the clade
# and whose third item is the cladeId corresponding to the entire tree
def extractClades(f):
    # Read the content of the file
    content = f.readline() # string to parse
    stack = [] 
    i = 0
    accum = ""
    digitAccum = ""
    considerDigit = False
    cladeMap = dict() # a map from clade ids to list of taxa in clade/subclades
    bootstrapMap = dict() # a map from clade id to the bootstrap branch value for the clade
    cladeId = 1


    while i < len(content) - 1:
        if content[i] == "(":
            considerDigit = False
            stack.append(content[i])
        elif isUpperCaseAlphabet(content[i]): # Build up taxa name
            accum += content[i]
        elif content[i] == ":" and accum != "": # Since every taxa name is followed by branch length, stop when you see : and a taxa is being built
            stack.append(accum) # so that it will be part of the next clade
            accum = "" # reset for next time
        elif content[i] == ")": # a clade is being closed
            # stack has all the members of this clade (some may be sub-clade ids) until we hit its corresponding "("
            members = []
            val = stack.pop()
            while val != "(":
                if isinstance(val, int): # This is a clade id so merge its contents with the contents of the new clade
                    members = members + cladeMap[val]
                else:
                    members.append(val)
                val = stack.pop()

            cladeMap[cladeId] = members # Map clade id to the members of the clade
            stack.append(cladeId) # It will be a subclade in any outer clades
            considerDigit = True # From now until next non ":", build up the bootstrap support value for this clade
        elif (isDigit(content[i]) or content[i] == ".") and considerDigit == True:
            digitAccum += content[i]
        elif content[i] == ":" and considerDigit == True: # The end of a the bootsrap support value for the current clade
            considerDigit = False
            bootstrapMap[cladeId] = int(digitAccum)
            digitAccum = ""  # Reset for next bootstrap support value
            cladeId = cladeId + 1 # I am done mapping this clade to anything. Move on to new clade id
        

        i = i + 1



    return (cladeMap, bootstrapMap, cladeId)

# Check for the different usages and handle them accordingly
if len(sys.argv) <= 2:
    usageError()
else:
    if sys.argv[1] == "-b" and len(sys.argv) == 3: # percentage of bootstrap values above 75 mode
        filename = sys.argv[2]
        f = open(filename)
        (_, bootstrapMap, _) = extractClades(f)
        print(str(percentAbove(bootstrapMap)) + "% of clades have branch support value >= 75")
    elif sys.argv[1] == "-c" and len(sys.argv) == 4: # Clade comparison mode
        filename1 = sys.argv[2]
        filename2 = sys.argv[3]
        (cladeMap1, bootstrapMap1, cladeId1) = extractClades(open(filename1))
        (cladeMap2, bootstrapMap2, cladeId2) = extractClades(open(filename2))
        ids_in_common_1 = []
        ids_in_common_2 = []

        # assume they have no clades in common until proven otherwise
        ids_not_in_common_1 = list(cladeMap1.keys()) 
        ids_not_in_common_2 = list(cladeMap2.keys())

        # they have the trivial clade in common so remove it from both lists
        ids_not_in_common_1.remove(cladeId1)
        ids_not_in_common_2.remove(cladeId2)

        # Find the clades they have in common and not in common
        print("file_1 has " + str(len(cladeMap1)) + " clades")
        print("file_2 has " + str(len(cladeMap2)) + " clades")
        print()
        inCommon = 0
        for k1 in cladeMap1:
            if k1 != cladeId1: # Avoid the clade for the entire tree
                for k2 in cladeMap2:
                    if k2 != cladeId2 and listsEqual(cladeMap1[k1], cladeMap2[k2]):
                        inCommon += 1
                        # Record the ids of the clades that they have in common. Although the clades are the same, their ids may not be the same
                        ids_in_common_1.append(k1)
                        ids_in_common_2.append(k2)

                        ids_not_in_common_1.remove(k1)
                        ids_not_in_common_2.remove(k2)
                        break # no duplicate clades so no need to go on
            
        print("Both files have " + str(inCommon)  +  " non trivial clades in common")
        print(str(percentageAboveWithCladeIds(bootstrapMap1, ids_in_common_1)) + "% of them have branch support value >= 75 in file_1")
        print(str(percentageAboveWithCladeIds(bootstrapMap2, ids_in_common_2)) + "% of them have branch support value >= 75 in file_2")
        print(str(jointPercentageAbove(bootstrapMap1, bootstrapMap2, ids_in_common_1, ids_in_common_2)) + 
        """% of them have joint branch support value >= 75 in file_1 and file_2""")
        print()

        if len(ids_not_in_common_1) != 0:
            print("file_1 has " + str(len(ids_not_in_common_1)) + " clades not in file_2")
            print(str(percentageAboveWithCladeIds(bootstrapMap1, ids_not_in_common_1)) + "% of them have branch support value >= 75")
        if len(ids_not_in_common_2) != 0:
            print()
            print("file_2 has " + str(len(ids_not_in_common_2)) + " clades not in file_1")
            print(str(percentageAboveWithCladeIds(bootstrapMap2, ids_not_in_common_2)) + "% of them have branch support value >= 75")         
    else:
        usageError()




