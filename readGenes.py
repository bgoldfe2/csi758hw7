import genbank as gb

def readGenesIn(fname):
    # Read in the bacteria file and parse the dna
    data = gb. ReadFile ( fname )
    dna = gb. ParseDNA ( data )

    return data,dna

def getGeneData():    
    
    # Look for non-complements second arg in gb.GenLocs is False
    path = 'bacteriaGB/Genbank/'
    for filename in os.listdir(path):
        print(path + filename)
        data,dna = readGenesIn(path+filename)
        
        # Find the keywoed locations and gene locations
        klocs = gb.FindKeywordLocs ( data )
        glocs = gb.GeneLocs ( data , klocs )

        for start in glocs[:1]:
            # Check that the complement is False aka a Non-Complement
            if (not start[1]):
                # check gene letters
                # remember to -1 in the start  position
                stNum = start[0][0][0] + 0
                endNum = start[0][0][1]
                print ("Gene start",str(stNum),"gene end",str(endNum))
                print ("   ",dna[stNum-1:stNum+20])


                
