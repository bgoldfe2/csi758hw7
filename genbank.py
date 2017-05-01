# genbank.py
# Python 3.4
# Python scripts are provided as an educational tool. They are offered without guarantee of effectiveness or accuracy. Python scripts composed by the author may not be used for commercial uses without the author's explicit written permission.
# (c) 2016 Jason M. Kinser

import string

def ReadFile( filename ):
    fp = open( filename )
    data = fp.read()
    fp.close()
    return data

def ParseDNA( data ):
     # find ORIGIN
    orgn = data.find( 'ORIGIN' )
    # find the first row after that
    start = data.find( '1', orgn )
    # find the ending slashes
    end = data.find( '//', orgn )
    # split the substring into lines
    a = data[start:end].split('\n')
    dna = []   # blank string
    for i in a:
        spl = i.split()
        dna.append( ''.join(spl[1:]) )
    answ =  ''.join(dna)
    return answ

def Complement( st ):
    table = st.maketrans( 'acgt', 'tgca' )
    st1 = st.translate(table)
    st1 = st1[::-1]
    return st1

def FindKeywordLocs( data, keyword = '   CDS   ' ):
    # number of occurences
    ngenes = data.count( keyword )
    # find the occurences
    keylocs = []
    k =0
    for i in range( ngenes ):
        t = data.find( keyword,k)
        keylocs.append( t )
        k = t + 10
    return keylocs

# get the start and end location for a normal or complement-only
def EasyStartEnd( data, loc, cflag=0 ):
    # find the ..
    dots = data.find( '..', loc )
    # targets
    if cflag==0:
        t1, t2 = ' ','\n'
    if cflag==1:
        t1, t2 = '(',')'
    # find the first preceding blank
    b1 = data.rfind( t1,loc, dots )
    # find the first following newline
    b2 = data.find( t2, dots )
    # extract the numbers
    temp = data[b1+1:dots]
    temp = temp.replace( '>', '' )
    temp = temp.replace( '<', '' )
    start = int( temp )
    temp = data[dots+2:b2]
    temp = temp.replace( '>', '' )
    temp = temp.replace( '<', '' )
    end = int( temp )
    return start, end

# possible that there is a complement INSIDe of a JOIN.
# following will find the limits and remove any complements
def ShortString( indata, loc ):
    data = indata[loc:loc+1000]
    data = data.replace( '\n','' )
    data = data.replace( '\r','' )
    data = data.replace( ' ','')
    p1 = data.find( '(' )+1
    p2 = data.find( ')', p1 )    # this may be incorrect place
    shst = data[p1:p2]  # short string answer
    if 'complement' in shst and 'join' in shst:
        # there is a complement.  End of string will be ))
        p2 = data.find( '))', p1 )
        # remove complements
        shst = data[p1:p2]
        N = shst.count ('complement')
        for i in range( N ):
            a1 = shst.find( 'complement')
            a2 = shst.find(')', a1 )
            shst = shst[:a2] + shst[a2+1:]
            shst = shst[:a1] + shst[a1+11:]
    return shst, p1, p2

# extract splice locations
def Splices( data, loc ):
    # loc is the location of '
    join = data.find( 'join', loc )
    # find the parentheses
    shst, p1, p2 = ShortString( data, join )
    # count the number of dots
    ndots = shst.count( '..' )
    # extract the numbers
    numbs = []
    # the first has a ( .. ,
    dots = shst.find( '..' )
    t = shst[:dots]
    t = t.replace( '<','')
    st = int( t )
    comma = shst.find( ',' )
    en = int( shst[dots+2:comma] )
    numbs.append( (st,en) )
    # consider the rest except first and last
    # code is , .. ,
    for i in range( ndots - 2):
        dots = shst.find( '..', comma )
        comma2 = shst.find( ',', dots )
        st = int( shst[comma+1:dots] )
        en = int( shst[dots+2:comma2] )
        numbs.append( (st,en) )
        comma = comma2
    # last one has code , .. )
    shst = shst.replace('>','')
    dots = shst.find( '..', comma )
    st = int( shst[comma+1:dots] )
    en = int( shst[dots+2: ] )
    numbs.append( (st,en) )
    return numbs

# just extract the locations of the genes
def GeneLocs( data, keylocs  ):
    genes = []
    for i in range( len( keylocs  )):
        # get this line and look for join or complement
        end = data.find( '\n', keylocs [i] )
        temp = data[keylocs [i]: end ]
        joinf = 'join' in temp
        compf = 'complement' in temp
        # get the extracts
        c = 0
        if compf == True: c = 1
        if joinf:
            numbs = Splices( data, keylocs[i] )
            genes.append( (numbs, compf ))
        else:
            sten = EasyStartEnd( data, keylocs[i], c )
            genes.append( ([sten],compf) )
    return genes

# using a gene location from GeneLocs, extract the DNA
def GetCodingDNA( dna, genesi ):
    # dna from ParseDNA
    # genesi is genes[i] from GeneLocs
    codedna = ''
    N = len( genesi[0] ) # number of splices
    for i in range( N ):
        st, en = genesi[0][i]
        codedna += dna[st-1:en]
    # complment flag
    if genesi[1]:
        codedna = Complement( codedna )
    return codedna
    

# codon dictionary
# stop = 'p'
# start = 'M'
def Codons( ):
    codons = {}
    CODONS = ( ('TTT','F'), ('TTC','F'), ('TTA','L'),('TTG','L'),\
               ('TCT','S'), ('TCC','S'), ('TCA','S'),('TCG','S'),\
               ('TAT','Y'), ('TAC','Y'), ('TAA','p'),('TAG','p'),\
               ('TGT','C'), ('TGC','C'), ('TGA','p'),('TGG','W'),\
               ('CTT','L'), ('CTC','L'), ('CTA','L'),('CTG','L'),\
               ('CCT','P'), ('CCC','P'), ('CCA','P'),('CCG','P'),\
               ('CAT','H'), ('CAC','H'), ('CAA','Q'),('CAG','Q'),\
               ('CGT','R'), ('CGC','R'), ('CGA','R'),('CGG','R'),\
               ('ATT','I'), ('ATC','I'), ('ATA','I'),('ATG','M'),\
               ('ACT','T'), ('ACC','T'), ('ACA','T'),('ACG','T'),\
               ('AAT','N'), ('AAC','N'), ('AAA','K'),('AAG','K'),\
               ('AGT','S'), ('AGC','S'), ('AGA','R'),('AGG','R'),\
               ('GTT','V'), ('GTC','V'), ('GTA','V'),('GTG','V'),\
               ('GCT','A'), ('GCC','A'), ('GCA','A'),('GCG','A'),\
               ('GAT','D'), ('GAC','D'), ('GAA','E'),('GAG','E'),\
               ('GGT','G'), ('GGC','G'), ('GGA','G'),('GGG','G'))
    for i in CODONS:
        codons[i[0].lower()] = i[1]
    return codons
    
def TranslateCDNA( cdna, codontable ):
    ## codontable from function Codons()
    # answer empty string
    answ =[]
    # for each three letters
    for i in range( 0,len(cdna),3):
        # get from dict 
        answ.append( codontable[cdna[i:i+3]] )
    # return 
    answ = ''.join( answ)
    return answ

# get the protein translation data from the file
def Translation( data, loc ):
    # get the amino acids
    trans = data.find( '/translation', loc )
    # find the second "
    quot = data.find( '"', trans + 15 )
    # extract
    prot = data[trans+14:quot]
    # remove newlines and blanks
    prot = prot.replace( '\n','' )
    prot = prot.replace( '\r','' )
    prot = prot.replace( ' ','')
    return prot
