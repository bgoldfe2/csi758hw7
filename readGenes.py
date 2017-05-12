import genbank as gb
import math
import numpy as np
import csv
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from operator import itemgetter, attrgetter

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
        
        # Find the keyword locations and gene locations
        klocs = gb.FindKeywordLocs ( data )
        glocs = gb.GeneLocs ( data , klocs )
    
        for start in glocs:
            # check for valid gene
            goodness,stNum,endNum = preCheck(start,dna)
            if (goodness):
                # Check that the complement is False aka a Non-Complement
                isComp = start[1]
                codon = dna[stNum:stNum+3]
                
                if (not isComp and codon == 'atg'): # these are the START 
                    startList.append(dna[stNum-30:stNum+23])
                else:  # This is NOT a START
                    nonStartList.append(dna[stNum-30:stNum+23])

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
    # Minus the last one due to k+1 evaluation
    seqCnt = len(train[0]) - 1
    # allocate 52 arrays
    for i in range(seqCnt):
        M.append(np.zeros((4,4)))

    #print(len(M),M[0].shape)
    #index = 4*k + 
    for h in train:
        clist = list(h)
        #print(clist[0],clist[1],clist[2])
        for j in range(len(clist)-1):
            k = sc2n(clist[j])
            kp1 = sc2n(clist[j+1])
            #print('j',j,'k',k,'kp1',kp1)
            M[j][k][kp1] += 1

    return M
                       
def sc2n(c):  # Stands for simple char to number
    # sends back the int for the char a=0,c=1,g=2,t=3
    ic = -1
    if (c == 'a'):
        ic = 0
    elif(c == 'c'):
        ic = 1
    elif(c == 'g'):
        ic = 2
    elif(c == 't'):
        ic = 3

    return ic

def convProb(M):

    M2 = []
    
    for MM in M:
        #print(MM.shape)
        rs = MM.sum(axis=1).astype(float)
        #print('rowsum',rs)
        #rowSum.append(rs)
        #AAAAHhhhhh this math is wrong as discovered on final!!!!
        #MM = MM/rs
        MM = (MM.T / rs).T
        #print(MM)
        M2.append(MM)

    return M2

def convOdds(M2):
    M3 = []
    odds = 0.25
    
    for MM in M2:
        MM = MM/odds
        #print(MM)
        M3.append(MM)

    return M3

def logOdds(M3):
    M4 = []
    replace = np.zeros((4,4))+0.0
    
    for MM in M3:
        MM = np.log(MM)
        
        ind = np.isnan(MM)
        #print("ind",ind)
        MM[ind] = replace[ind]
        ind = np.isinf(MM)
        #print("ind",ind)
        MM[ind] = replace[ind]
        #print(MM)
        M4.append(MM)

    return M4

def scoreString(M4,stringList):
    # Question 4 score the distributions
    sttScore = []
    score = 0.0
    geneNum = 0
    # loop through the strings
    for h in stringList:
        
        # convert string to list of chars
        clist = list(h)
        #print(clist[0],clist[1],clist[2])
        for j in range(len(clist)-1):
            k = sc2n(clist[j])
            kp1 = sc2n(clist[j+1])
            #print('gn',geneNum,'j',j,'k',k,'kp1',kp1)
            #read and sum the values
            val = float(M4[j][k][kp1])
            score += val
            
            #print('gn',geneNum,'score',j,score)
        
        # Append the total score to list
        #print('gn',geneNum,'final score',score)
        sttScore.append(score)
        score = 0.0
        geneNum += 1

    return sttScore

def findMaxScores(strList):
    #My creative solution
    out = []
    ind = 0

    for a in strList:
        mxInd = np.unravel_index(a.argmax(), a.shape)
        mxVal = a.max()
        out.append((mxVal,mxInd,ind))
        ind += 1
        
    sout = sorted(out, key=lambda x: x[0],reverse=True)

    return sout

def makeImportant(M4,hiVals):
    # make all but the important zeros
    M5 = []
    replace = np.zeros((4,4))+0.0
    ind = 0
    for MM in M4:
        if (ind in hiVals):
            MM = M4[ind]
            M5.append(MM)
        else:
            MM = replace
            M5.append(MM)
            
        ind += 1

    return M5    

def list2File(fname,data):
    with open(fname, 'w') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        wr.writerow(data)
        


def plotG(mu,sd,labs):
    # plot the mu,sd
    x = np.linspace(-10, 20, 100)
    fig = plt.figure()
    ax = plt.subplot(111)
    for i in range(3):
        ax.plot(x,mlab.normpdf(x, mu[i], sd[i]),label=labs[i])

    ax.legend()
    plt.show()

                
def Driver():
    start, nonStart = getGeneData()
    stTrain,stOp = splitTrain(start)

    # Get the initial counts matrix    
    M1 = getCounts(stTrain)

    # convert to a probability matrix
    M2 = convProb(M1)

    # convert to odds by dividing by 0.25
    M3 = convOdds(M2)

    # get log of odds matrix
    M4 = logOdds(M3)

    # score all of the strings
    sttScore = np.array(scoreString(M4,stTrain))
    stOpScore = np.array(scoreString(M4,stOp))
    nonStScore = np.array(scoreString(M4,nonStart))
    
    mu = [np.mean(sttScore),np.mean(stOpScore),np.mean(nonStScore)]
    sd = [np.std(sttScore),np.std(stOpScore),np.std(nonStScore)]
    print(mu)
    print(sd)
    labs = ['Training Starts', 'Non Training Starts','Non Starts']
    plotG(mu,sd,labs)

    hiValCodon = findMaxScores(M4)
    print(hiValCodon[:5])

    # this resulted in (29,33,45,36,22)
    hiVals = [x[2] for x in hiValCodon[:5]]
    M5 = makeImportant(M4,hiVals)

    sttSc2 = np.array(scoreString(M5,stTrain))
    stOpSc2 = np.array(scoreString(M5,stOp))
    nonStSc2 = np.array(scoreString(M5,nonStart))
    
    mu2 = [np.mean(sttSc2),np.mean(stOpSc2),np.mean(nonStSc2)]
    sd2 = [np.std(sttSc2),np.std(stOpSc2),np.std(nonStSc2)]
    print(mu2)
    print(sd2)
    labs = ['Training Starts', 'Non Training Starts','Non Starts']
    plotG(mu2,sd2,labs)

                

