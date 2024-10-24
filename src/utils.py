import sys
sys.path.append('/home/sfromm/git/Mosaist/lib')
import mstpython as mst
from anarci import anarci
from Bio import Align





globalXXaligner = Align.PairwiseAligner()
localMSaligner = Align.PairwiseAligner()
localMSaligner.mode = "local"
localMSaligner.match = 2
localMSaligner.mismatch_score = -1
localMSaligner.open_gap_score = -0.5
localMSaligner.extend_gap_score = -0.1


def seqPer(seqA, seqB, alignType):
    try:
        if alignType == 'local':
            alignmentList = localMSaligner.align(seqA, seqB)
            alignment = alignmentList[0]
            string = alignment.format("psl")
            lst = string.split('\t')
            matchNum = int(lst[0])
            gaps1num = int(lst[5])
            gaps2num = int(lst[7])
            denom = matchNum + gaps1num + gaps2num
            if denom < 30:
                return 0
                
        elif alignType == 'global':
            matchNum = globalXXaligner.score(seqA, seqB)
            denom = min(len(seqA),len(seqB))
        
        else:
            sys.exit("can only provide global or local as the alignment type...")

        per = float(matchNum/denom)
        return per
    
    except IndexError:
        # accounts for cases like a: 'RSISSINIHTRDGSTTTLTGFPRIRS', aa: 'AAAAAAAAAAAAAAAAAAAAAAAAAAY' wherein there's no decent local alignment at all... in which case you should just skip checking the sequence ID % and return 0
        print(f"no alignment made...\nseqA was:\n{seqA}\nseqB was:\n{seqB}")
        return 0.0
    
    except Exception as error:
        sys.exit(f"An unexpected error occured:\n{error}")



def isScFv(hChain,lChain):
    hChainCapital = hChain.upper()
    lChainCapital = lChain.upper()
    if hChain == lChainCapital:
        return True
    elif lChain == hChainCapital:
        return True
    else:
        return False



def tripleToSingle(aa):
    if aa == "HIS":
        return "H"
    elif aa == "ALA":
        return "A"
    elif aa == "GLY":
        return "G"
    elif aa == "ILE":
        return "I"
    elif aa == "LEU":
        return "L"
    elif aa == "PRO":
        return "P"
    elif aa == "VAL":
        return "V"
    elif aa == "PHE":
        return "F"
    elif aa == "TRP":
        return "W"
    elif aa == "TYR":
        return "Y"
    elif aa == "ASP":
        return "D"
    elif aa == "GLU":
        return "E"
    elif aa == "ARG":
        return "R"
    elif aa == "LYS":
        return "K"
    elif aa == "SER":
        return "S"
    elif aa == "THR":
        return "T"
    elif aa == "CYS":
        return "C"
    elif aa == "MET":
        return "M"
    elif aa == "ASN":
        return "N"
    elif aa == "GLN":
        return "Q"
    else:
        #print("unknown amino acid " + aa)
        return "X"

# parses sequence from SEQRES section of a PDB file
def seqResParser(filePath):
    seqDict = {}
    pdbFile = open(filePath, 'r')
    fileList = pdbFile.readlines()
    for fileLine in fileList:
        lineSplit = fileLine.split()
        if lineSplit[0] == 'SEQRES':
            lineChain = fileLine[11]
            lineSeq = ''
            lineLen = len(lineSplit)
            a = 4
            while a < lineLen:
                lineRes = lineSplit[a]
                shortRes =tripleToSingle(lineRes)
                lineSeq += shortRes
                a += 1
            if lineChain in seqDict:
                seqDict[lineChain] += lineSeq
            else:
                seqDict[lineChain] = lineSeq
    return seqDict



# anarci parser
def anarciParser(oriSeqRes,anarciOut):

    # go along the tuples of the original sequence, and check the next residue in the anarci output; if it matches, the IMGT numbering to the tuple, and incriment to move on to the next in the original sequence; if it's a -, continue, and if it's neither... report an error cuz something's gone wrong

    # oriI & oriEnd will be values like 0 & 115, where oriI is "The index in the sequence of the first numbered residue", so essentially - when you put in the sequence AMG... if the numbering starts at M, oriI = 1, and oriEnd is "The index in the sequence of the last numbered residue", so if your sequence is 120 residues long but the last thing numbered by anarci is at index 114, oriEnd will be 114.

    # note: as anNumbering (the numbering according to anarci) will sometimes have lots of dashes in it (that is, instead of matching a residue you put in to a number, it has a gap / dash there that matches to an number) it will not necessarily be the size of the span of indexes (like 0-115 = 116), it will be at least that length or higher

    try:
        anNumbering = anarciOut[0]
        oriI = anarciOut[1]
        oriEnd = anarciOut[2]
    except:
        print("oriSeqRes: ")
        print(oriSeqRes)
        print("anNumbering: ")
        print(anNumbering)

    oriSeqIMGTtuples = []

    for anRes in anNumbering:
        anNumber = anRes[0]
        anName = anRes[1]
        if anName == '-':
            continue
        elif anName == oriSeqRes[oriI]: #sanity check to make sure it's the right residue
            tupleToAdd = (anNumber,anName)
            oriSeqIMGTtuples.append(tupleToAdd)
            oriI += 1
        else:
            print("error: anarci numbering & original pdb numbering off somehow")
            print("oriI:")
            print(oriI)
            print("oriEnd")
            print(oriEnd)
            print("anNumbering:")
            print(anNumbering)
            print("anarciOut:")
            print(anarciOut)
            print("oriSeqRes")
            print(oriSeqRes)
            print()
            quit()

    return oriSeqIMGTtuples



# structs to seqs function
def abAgStructs2Seqs(listArg,inputDb,strip_res_x=False):

    infoAndSeqList = []

    for i in listArg:

        currFilePath = inputDb + "/" + str(i[1]) + ".pdb"

        #currS = str(i[0])
        currPDB = str(i[1])
        currAg = str(i[4])
        currAbH = str(i[2])
        currAbL = str(i[3])
        currComplexType = str(i[5])
        currRes = i[6]

        if isScFv(currAbH,currAbL):
            abScFv = True
            currAbH = currAbH.upper()
            currAbL = currAbL.upper()
        else:
            abScFv = False

        # load the struct

        try:
            structAbAg = mst.Structure(currFilePath, "QUIET")
        except:
            print("testing failed file path to get structure from")
            print(currFilePath)
            continue

        seqResDict = seqResParser(currFilePath)

        # 1. get the sequences for the antibody and antigen [antigen seq,HloopSeqs,LloopSeqs] format

        if currAbH == "NA" and currAbL == "NA":
            print("WARNING: skipping entry" + currPDB + " because it apparently has neither heavy or light chains")
            continue

        seqFail = False

        agList = []
        aSeqReses = []
        aSeqs = []
        fractionAtomSeqres_a = 1

        if currAg != "NA":
            if " | " in currAg:
                agList = currAg.split(" | ")
            else:
                agList = [currAg]
            for ag in agList:
                aChain = structAbAg.getChainByID(ag)
                aSeq = ''
                try:
                    resA = aChain.getResidues()
                except:
                    seqFail = True
                    print("WARNING: failed at getting antigen residues from entry " + currPDB + "with chain ID(s): ")
                    print(agList)
                    continue

                try:
                    aSeqRes = seqResDict[ag]
                except:
                    seqFail = True
                    print("WARNING: failed at getting antigen sequence from entry " + currPDB + "with ag chain: " + ag)
                    print("dictionary of all sequences is as follows: ")
                    print(seqResDict)
                    continue

                #aNames = []
                for ii in range(len(resA)):

                    aRes = resA[ii]
                    aName = aRes.name
                    aNameS = tripleToSingle(aName)
                    aSeq += aNameS

                aSeqReses.append(aSeqRes)
                aSeqs.append(aSeq)

            a = 0
            fractionAtomSeqres_as = []
            while a < len(aSeqReses):
                fractionAtomSeqres_a = len(aSeqs[a])/len(aSeqReses[a])
                fractionAtomSeqres_as.append(fractionAtomSeqres_a)
                a +=1
            fractionAtomSeqres_a = min(fractionAtomSeqres_as)
        else:
            print("WARNING: no antigen found for entry " + currPDB)

        if seqFail == True:
            continue

        if currAbH != "NA":
            chainAbH = structAbAg.getChainByID(currAbH)
            try:
                resH = chainAbH.getResidues()
            except:
                print("WARNING: getting heavy chain residues failed for entry " + currPDB + " with chain ID: " + currAbH)
                continue
            hNames = []
            hSeq = ''
            if currAbH not in seqResDict:
                #print("chains off in file; skipping")
                continue
            hSeqRes = seqResDict[currAbH]

            for ii in range(len(resH)):
                hRes = resH[ii]
                hName = hRes.name
                hNum = hRes.num
                hCode = hRes.iCode
                hIndex = hRes.getResidueIndex()
                hNameS =tripleToSingle(hName)
                hNames.append([hIndex,(hNum,hCode),hNameS])
                hSeq += hNameS

            fractionAtomSeqres_h = len(hSeq)/len(hSeqRes)
        else:
            hSeqRes = ''
            fractionAtomSeqres_h = 1

        if abScFv:
            lSeqRes = hSeqRes
            fractionAtomSeqres_l = 1
        elif currAbL != "NA":
            chainAbL = structAbAg.getChainByID(currAbL)
            try:
                resL = chainAbL.getResidues()
            except:
                print("WARNING: getting light chain residues failed for entry " + currPDB + " with chain ID: " + currAbL)
                continue
            lNames = []
            lSeq = ''

            if currAbL not in seqResDict:
                continue

            lSeqRes = seqResDict[currAbL]

            for ii in range(len(resL)):
                lRes = resL[ii]
                lName = lRes.name
                lNum = lRes.num
                lCode = lRes.iCode
                lIndex = lRes.getResidueIndex()
                lNameS =tripleToSingle(lName)
                lNames.append([lIndex,(lNum,lCode),lNameS])
                lSeq += lNameS

            fractionAtomSeqres_l = len(lSeq)/len(lSeqRes)
        else:
            lSeqRes = ''
            fractionAtomSeqres_l = 1

        # 2. run anarci to get alignment, and 3. correspond the old numbering & indexes to Anarci numbering so loops can be identified later

        hNumbering = "NA"
        if currAbH != "NA":
            anarciH, anarciHextra, _ = anarci([("placeholderName",hSeqRes)], scheme="imgt", output=False)

            anarciHdicts = anarciHextra[0]
            if anarciHdicts is None:
                #print('anarciHdicts empty; skipping')
                continue

            d = 0
            anarciHl = 'NA'
            while d < len(anarciHdicts):
                anarciDict = anarciHdicts[d]
                if "chain_type" in anarciDict and anarciDict['chain_type'] == 'H':
                    anarciHl = anarciH[0][d]
                    break
                else:
                    d += 1

            if anarciHl == 'NA':
                print("H chain in " + currPDB +"'s anarci output not found, this structure will be skipped'")
                print("anarciHdicts were:")
                print(anarciHdicts)
                print("from h chain sequence:")
                print(hSeqRes)
                continue

            hNumbering = anarciParser(hSeqRes,anarciHl) # removed sorted() around anarciParser output

        lNumbering = "NA"
        if currAbL != "NA":
            anarciL, anarciLextra, _ = anarci([("placeholderName",lSeqRes)], scheme="imgt", output=False)
            anarciLdicts = anarciLextra[0]

            if anarciLdicts is None:
                continue

            d = 0
            anarciLl = 'NA'
            while d < len(anarciLdicts):
                anarciDict = anarciLdicts[d]
                if "chain_type" in anarciDict and ((anarciDict['chain_type'] == 'L') or (anarciDict['chain_type'] == 'K')):
                    anarciLl = anarciL[0][d]
                    break
                else:
                    d += 1

            if anarciLl == 'NA':
                print("L chain in " + currPDB +"'s anarci output not found, this structure will be skipped'")
                print("anarciLdicts were:")
                print(anarciLdicts)
                print("from l chain sequence:")
                print(lSeqRes)
                continue

            lNumbering = anarciParser(lSeqRes,anarciLl) # removed sorted() around anarciParser output

        # parse down the sequences to only be IMGT numbered ones

        hSeqShort = ''
        lSeqShort = ''
        if hNumbering != 'NA':
            for h in hNumbering:
                if len(h) == 1:
                    continue
                else:
                    hSeqShort += h[1]

        if lNumbering != 'NA':
            for l in lNumbering:
                if len(l) == 1:
                    continue
                else:
                    lSeqShort += l[1]

        # get all the loops by IMGT; CDR1 = 27-38, CDR2 = 56-65, CDR3 = 105-117

        currCDRH1 = ''
        currCDRH2 = ''
        currCDRH3 = ''
        currCDRL1 = ''
        currCDRL2 = ''
        currCDRL3 = ''
        if currAbH != "NA":
            ni = 27
            while ni < 118:
                for hTuple in hNumbering:
                    if len(hTuple) == 2:
                        hIMGTnum = hTuple[0][0]
                        if float(hIMGTnum) == ni:
                            if 27 <= float(hIMGTnum) and float(hIMGTnum) <= 38:
                                currCDRH1 += hTuple[1]
                            if 56 <= float(hIMGTnum) and float(hIMGTnum) <= 65:
                                currCDRH2 += hTuple[1]
                            if 105 <= float(hIMGTnum) and float(hIMGTnum) <= 117:
                                currCDRH3 += hTuple[1]
                ni += 1
            hLoopSeqs = currCDRH1 + currCDRH2 + currCDRH3
            hLoopSeqs = hLoopSeqs.strip("X")
        else:
            hLoopSeqs = ''

        if currAbL != "NA":
            ni = 27
            while ni < 118:
                for lTuple in lNumbering:
                    if len(lTuple) == 2:
                        lIMGTnum = ''
                        lIMGTnum = lTuple[0][0]
                        if float(lIMGTnum) == i:
                            if 27 <= float(lIMGTnum) and float(lIMGTnum) <= 38:
                                currCDRL1 += lTuple[1]
                            if 56 <= float(lIMGTnum) and float(lIMGTnum) <= 65:
                                currCDRL2 += lTuple[1]
                            if 105 <= float(lIMGTnum) and float(lIMGTnum) <= 117:
                                currCDRL3 += lTuple[1]
                ni += 1
            lLoopSeqs = currCDRL1 + currCDRL2 + currCDRL3
            lLoopSeqs = lLoopSeqs.strip("X")
        else:
            lLoopSeqs = ''

        minAtomSeqresFraction = min(fractionAtomSeqres_a, fractionAtomSeqres_h, fractionAtomSeqres_l)

        if strip_res_x:
            newaSeqs = []
            for aSeq in aSeqReses:
                newSeq = aSeq.strip("X")
                if newSeq != '':
                    newaSeqs.append(newSeq)
            hSeqShort = hSeqShort.strip("X")
            lSeqShort = lSeqShort.strip("X")

        infoAndSeqList.append({'currPDB':currPDB, 'aSeq':newaSeqs, 'hSeq': hSeqShort, 'lSeq': lSeqShort, 'hLoopSeqs':hLoopSeqs, 'lLoopSeqs':lLoopSeqs, 'complexType':currComplexType, 'currAg':currAg, 'currAbH':currAbH, 'currAbL':currAbL, 'minAtomSeqresFraction':minAtomSeqresFraction, 'currRes':currRes})

    return infoAndSeqList
