# Chou and Fasman parameters for Alpha Helix detection
palpha = {
    'E' : [1.53, 1], 
    'A' : [1.45, 1],
    'L' : [1.34, 1],
    'H' : [1.24, 1],
    'M' : [1.20, 1],
    'Q' : [1.17, 1],
    'W' : [1.14, 1],
    'V' : [1.14, 1],
    'F' : [1.12, 1],
    'K' : [1.07, 0.5],
    'I' : [1.00, 0.5],
    'D' : [0.98, 0],
    'T' : [0.82, 0],
    'S' : [0.79, 0],
    'R' : [0.79, 0],
    'C' : [0.77, 0],
    'N' : [0.73, -1],
    'Y' : [0.61, -1],
    'P' : [0.59, -1],
    'G' : [0.53, -1],
}

# Chou and Fasman parameters for Beta strand detection
pbeta = {
    'M' : [1.67, 1],
    'V' : [1.65, 1],
    'I' : [1.60, 1],
    'C' : [1.30, 1],
    'Y' : [1.29, 1],
    'F' : [1.28, 1],
    'Q' : [1.23, 1],
    'L' : [1.22, 1],
    'T' : [1.20, 1],
    'W' : [1.19, 1],
    'A' : [0.97, 0],
    'R' : [0.90, 0],
    'G' : [0.81, -1],
    'D' : [0.80, -1],
    'K' : [0.74, -1],
    'S' : [0.72, -1],
    'H' : [0.71, 0],
    'N' : [0.65, 0],
    'P' : [0.62, -1],
    'E' : [0.26, -1],
}

# output of STRIDE webserver
stride = 'TTTT     HHHHHH EEEEEETTEEEEEEEETTEEEEEGGGG  HHHHH   HHHHHHH  GGG EEEETTEEE EEEEEEETTEEEEEE   TTTT        TTTEEEEEEEEETTEEEEEEEEEETTTT B    TTTTTTTEE '
# given input sequence
seq = 'SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDTVYCPRHVICTAEDMLNPNYEDLLIRKSNHSFLVQA\
GNVQLRVIGHSMQNCLLRLKVDTSNPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNHTIKGSFLNGSCGSVGF'
#converting string to list
seq = list(seq) 
stride = list(stride)
# initialise output strings to all turns
str1 = []
for x in range(len(seq)):
    str1.append('T')
str2 = []
for x in range(len(seq)):
    str2.append('T')
final = []
for x in range(len(seq)):
    final.append('T')

# identifying alpha helices
# traverse through sequence forming windows of six residues and extending in both directions
for x in range(len(seq)-5):
    ctr1 = 0
    site = []
    #checking windows of six residues for nucleation sites
    for i in range(x, x+6):
        site.append(seq[i])
        if palpha[seq[i]][0] >= 1:
            ctr1 +=1

    if ctr1 >= 4:
        
        #assigning alpha helix structure to residues in identified nucleation sites
        for j in range(x, x+6):
            str1[j] = 'H'
        
        #extending in right direction from nucleation sites
        ctr2 = 0
        for y in range(x+3, len(seq)-3):
            ctr2 = 0
            for i in range(y, y+4):
                ctr2 += palpha[seq[i]][0]
            if ctr2 >= 4:
                for j in range(y, y+4):
                    str1[j] = 'H'
            else:
                break
        
        #extending in left direction from nucleation sites
        ctr2 = 0
        for y in range(x+2, -1, -1):
            ctr2 = 0
            for i in range(y, y-4, -1):
                ctr2 += palpha[seq[i]][0]
            if ctr2 >= 4:
                 for j in range(y, y-4, -1):
                    str1[j] = 'H'
            else:
                break
                
# identifying beta strands
# traverse through sequence forming windows of fuve residues and extending in both directions
for x in range(len(seq)-4):
    ctr1 = 0
    site = []
    
    #checking windows of five residues for nucleation sites
    for i in range(x, x+5):
        site.append(seq[i])
        if pbeta[seq[i]][0] > 1:
            ctr1 +=1

    if ctr1 >= 3:
        
        #assigning beta strand structure to residues in identified nucleation sites
        for j in range(x, x+5):
            str2[j] = 'S'
        
        #extending in right direction from nucleation sites
        ctr2 = 0
        for y in range(x+3, len(seq)-3):
            ctr2 = 0
            for i in range(y, y+4):
                ctr2 += pbeta[seq[i]][0]
            if ctr2 > 4:
                for j in range(y, y+4):
                    str2[j] = 'S'
            else:
                break
        
        #extending in left direction from nucleation sites
        ctr2 = 0
        for y in range(x+2, -1, -1):
            ctr2 = 0
            for i in range(y, y-4, -1):
                ctr2 += pbeta[seq[i]][0]
            if ctr2 > 4:
                 for j in range(y, y-4, -1):
                    str2[j] = 'S'
            else:
                break        
# conflict resolution by comparing average P(H) and P(S) in conflicting regions
totalpha = 0
totbeta = 0
flag = 0
pos = 0
for i in range(len(seq)):
    
#     adding values of alpha helix and beta strand propensities at conflicting region 
    if str1[i] == 'H' and str2[i] == 'S':
        totalpha += palpha[seq[i]][0]
        totbeta += pbeta[seq[i]][0]
        
#         keeping track of position where conflicting region starts
        if flag == 0:
            flag = 1
            pos = i
    else:
        if str1[i] == 'H' and str2[i] == 'T':
            final[i] = 'H'
        elif str1[i] == 'T' and str2[i] == 'S':
            final[i] = 'S'
            
#         assigning secondary structure based on greater average propensity value after conflicting region ends
        if flag == 1:
            flag = 0
            if totalpha > totbeta:
                for j in range(pos, i):
                    final[j] = 'H'
            else:
                for j in range(pos, i):
                    final[j] = 'S'
            totalpha = 0
            totbeta = 0

if flag == 1:
    flag = 0
    if totalpha > totbeta:
        for j in range(pos, i+1):
            final[j] = 'H'
    else:
        for j in range(pos, i+1):
            final[j] = 'S'
            
# printing outputs
# Output is displayed using notation: H: Helix, S: Beta strand, and T: turn
for ctr in range(len(seq)): # Input sequence
    print(seq[ctr], end = '')
print()
for ctr in range(len(seq)): # Secondary structure prediction output from implemented Chou-Fasman algorithm
    print(final[ctr], end = '')
print()
for ctr in range(len(seq)): # Secondary structure prediction output from STRIDE web server
    print(stride[ctr], end = '')
print()

diff = 0
for ctr in range(len(seq)): # Comparing outputs from implemented Chou-Fasman algorithm and STRIDE web server where * denotes differing region
#     comparing all types of assigned helix structures
    if (stride[ctr] == 'G' or stride[ctr] == 'H') and final[ctr] != 'H': # H(alpha helix) and G(310helix) on STRIDE are both assigned helix structure(H) by our algorithm.
            print('*', end = '')
            diff+=1
#     comparing all types of beta structures
    elif (stride[ctr] == 'E' or stride[ctr] == 'B') and final[ctr] != 'S':  # E(beta strand) and B(beta bridge) on STRIDE are both assigned beta strand structure(S) by our algorithm.
            print('*', end = '')
            diff+=1
#     comparing assigned turn and coil structures
    elif (stride[ctr] == 'T' or stride[ctr] == ' ') and final[ctr] != 'T':  # T(turn) and 'blank space'(coil) on STRIDE are both assigned turn structure(T) by our algorithm.
            print('*', end = '')
            diff+=1
    else:
        print(' ', end = '')
print()
# printing percentage match between outputs from implemented Chou-Fasman algorithm and STRIDE web server
print(100-(diff/len(seq)*100))            
