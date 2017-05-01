import genbank as gb

def readGenesIn(fname):
    #fname = 'bacteriaGB/Genbank/AE000512.gb.txt'
    data = gb. ReadFile ( fname )
    dna = gb. ParseDNA ( data )
    #len( dna )
    klocs = gb. FindKeywordLocs ( data )
    glocs = gb. GeneLocs ( data , klocs )
    #glocs [0]
    return glocs

def GetData( indir ):
    a = os.listdir( indir )
    bglocs = []
    for i in a:
        fbact = indir + i
        print(fbact) 
        bglocs.append(readGenesIn(fbact))
    return bglocs

        
