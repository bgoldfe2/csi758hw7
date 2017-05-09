import genbank as gb
import math
import numpy as np

def readGenesIn(fname):
    # Read in the bacteria file and parse the dna
    data = gb. ReadFile ( fname )
    dna = gb. ParseDNA ( data )
    return data,dna

def getGeneData():    
    
    # Look for non-complements second arg in gb.GenLocs is False
    path = 'bacteriaGB/Genbank/'

    # lists for starts and nonStarts
    startList = []
    nonStartList = []
    for filename in os.listdir(path):
        #print(path + filename)
        data,dna = readGenesIn(path+filename)
        
        # Find the keywoed locations and gene locations
        klocs = gb.FindKeywordLocs ( data )
        glocs = gb.GeneLocs ( data , klocs )
    
        for start in glocs:
            #print("=====================================")
            # check for valid gene
            goodness,stNum,endNum = preCheck(start,dna)
            if (goodness):
                # Check that the complement is False aka a Non-Complement
                isComp = start[1]
                codon = dna[stNum:stNum+3]
                #print("codon",codon)
                if (not isComp and codon == 'atg'): # these are the START 
                    #print("START Gene",isComp)
                    startList.append(dna[stNum:stNum+53])
                else:  # This is NOT a START
                    #print("NOT start Gene",isComp)
                    nonStartList.append(dna[stNum:stNum+53])

    #print("start size",len(startList),"nst size",len(nonStartList))
    return startList,nonStartList

def preCheck(start,dna):
    # check for chars 30 in front 20 after and valid A, T, G, C
    # check gene letters
    # remember to -1 in the start  position
    goodness = False
    stNum = (start[0][0][0] + 0)-1
    endNum = start[0][0][1]
    #print("start",start)
    if (stNum > 30 and len(dna) - endNum > 23):
        #print("good lenght",str(stNum),str(endNum))
        goodness = True
    #else:
        #print("bad length",str(stNum),str(endNum))

    #print ("Gene start",str(stNum),"gene end",str(endNum))

    # Check that all the characters are valid
    if (stNum+53 < len(dna)):
        okchars = ['a','c','g','t']
        gene = dna[stNum:stNum+53]
        #print ("gene chars ",gene)
        for i in range(0, len(gene)):
            if (gene[i] not in okchars):
                #print("BAD CHARACTER",i)
                goodness = False

    return goodness,stNum,endNum

def splitTrain(start):
    total = len(start)
    trainSize = math.floor(total * .2)
    stTrain = start[:trainSize]
    stOp = start[trainSize:]
    return stTrain,stOp

def getCounts(train):
    M = []
    # allocate 52 arrays
    for j in range(52):
        M.append(np.zeros((4,4)))

#TODO need to convert letters to index numbers

    for k in range(52):    
        for h in range(len(train[:5])):
            print(train[h][k],train[h][k+1])
                
def Driver():
    start, nonStart = getGeneData()
    stTrain,stOp = splitTrain(start)
    getCounts(stTrain)

    
    return stTrain
                

                
>>>>>>> f228cc02ec64aabdb5d53ca903f9cac9bc71fcd3
