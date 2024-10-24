import datetime
from Bio.Seq import Seq
import pandas as pd

from .utils import seqPer


def prefilterdataset(sumFileName,methodsList,resCut,splitdate,cutoffDate=None):

    def convert_resolution(string):
        try:
            r = float(string)
            return r
        except ValueError:
            return None

    onOrAfterDateList = []
    beforeDateProtPepComplexesList = []

    cutoffFilterCount = 0
    methodFilterCount = 0
    qualFilterCount = 0
    
    summary_df = pd.read_csv(sumFileName, sep='\t', na_filter=False)
    for row in summary_df.itertuples(index=True, name='summaryrow'):

        currPdb = row.pdb
        currAbH = row.Hchain
        currAbL = row.Lchain
        currAg = row.antigen_chain
        currComplexType = row.antigen_type
        currDate = row.date
        currResolution = convert_resolution(row.resolution)
        currMethod = row.method

        dateSplit = currDate.split('/')
        dateYearEnd = int(dateSplit[2])
        if dateYearEnd > 50:
            dateYear = dateYearEnd + 1900
        else:
            dateYear = dateYearEnd + 2000
        try:
            d1 = datetime.datetime(dateYear, int(dateSplit[0]), int(dateSplit[1]))
        except:
            #print("date out of range so skipping line:")
            #print(line)
            continue


        # make a list A of all the PDB IDs on or after Date

        currList = [currDate,currPdb,currAbH,currAbL,currAg,currComplexType,currResolution]

        if cutoffDate:
            if d1 >= cutoffDate:
                continue

        if d1 >= splitdate:
            
            if currMethod not in methodsList:
                methodFilterCount += 1
                continue

            if currResolution and currResolution > float(resCut):
                #print("skipping " + currPdb + " because resolution too poor")
                qualFilterCount += 1
                continue
            if 'protein' in currComplexType:
                onOrAfterDateList.append(currList)

        # make a list B of all the PDB IDs before Date

        else:
            if 'protein' in currComplexType:
                cutoffFilterCount += 1
            if 'protein' in currComplexType or 'peptide' in currComplexType: # TODO allow all?
                beforeDateProtPepComplexesList.append(currList)

    return onOrAfterDateList, beforeDateProtPepComplexesList, cutoffFilterCount, methodFilterCount, qualFilterCount


def filterbysimilarity(onOrAfterDateInfoAndSeqsList, percentWinDatasetCut, pAbSeqCut, beforeDatePpInfoAndSeqsList=[],disallow_res_x=False, use_globalSeqID=False, whiteList=[],use_cutoffStrict=False, minAtomSeqresFraction=0):

    acceptedAbProtStructs = []
    prevSimilarityCount = 0
    knockoutWithinCount = 0

    for i in onOrAfterDateInfoAndSeqsList:
        toRemoveIfAdding = []
        keepI = True
        acceptedStructInfo = {}

        if disallow_res_x:
            if ("X" in i['hSeq']) or ("X" in i['lSeq']) or ("X" in i['aSeq']):
                continue

        currPdb = i['currPDB']

        print("working on classifying entry: " + currPdb)

        ChSeq = Seq(i['hLoopSeqs'])
        ClSeq = Seq(i['lLoopSeqs'])
        CaSeqs = [Seq(s) for s in i['aSeq']]
        currScore = i['minAtomSeqresFraction']
        currComplexType = i['complexType']
        currRes = i['currRes']

        if currScore < minAtomSeqresFraction:
            continue

        # check for repeats
        worseRepeat = False
        betterRepeat = False
        ii = 0

        while ii < len(acceptedAbProtStructs):

            # when considering which structures to add, compare H loop sequences, L loop seqeunces, and antigen sequences to previously accepted structures

            closestAbMatchScore = 0
            closestAgMatchScore = 0

            Ppdb = acceptedAbProtStructs[ii]['currPDB']
            PhSeq = Seq(acceptedAbProtStructs[ii]['hLoopSeqs'])
            PlSeq = Seq(acceptedAbProtStructs[ii]['lLoopSeqs'])
            PaSeqs = [Seq(s) for s in acceptedAbProtStructs[ii]['aSeq']]
            pastScore = acceptedAbProtStructs[ii]['minAtomSeqresFraction']
            pastRes = acceptedAbProtStructs[ii]['currRes']

            if len(PhSeq) > 1 and len(ChSeq) > 1:
                if use_globalSeqID:
                    hScore = seqPer(PhSeq,ChSeq,'global')
                else:
                    hScore = seqPer(PhSeq,ChSeq,'local')

                if closestAbMatchScore < hScore:
                    closestAbMatchScore = hScore

            if len(PlSeq) > 1 and len(ClSeq) > 1:
                if use_globalSeqID:
                    lScore = seqPer(PlSeq,ClSeq,'global')
                else:
                    lScore = seqPer(PlSeq,ClSeq,'local')

                if closestAbMatchScore < lScore:
                    closestAbMatchScore = lScore

            aSimBest = 0
            caLen = 0

            for a in CaSeqs:
                if len(a) == 0:
                    continue
                caLen += len(a)
                paLen = 0
                for aa in PaSeqs:
                    paLen += len(aa)

                    if use_globalSeqID:
                        aScore = seqPer(a,aa,'global')
                    else:
                        aScore = seqPer(a,aa,'local')
                    
                    if closestAgMatchScore < aScore:
                        closestAgMatchScore = aScore
            
            if (use_cutoffStrict and ((closestAgMatchScore >= percentWinDatasetCut) or (closestAbMatchScore >= percentWinDatasetCut))) or (not use_cutoffStrict and ((closestAgMatchScore >= percentWinDatasetCut) and (closestAbMatchScore >= percentWinDatasetCut))):

                knockoutWithinCount += 1

                # if something previous that it might knockout is in the whitelist, don't add it and keep the whitelisted structure. Then if the prev structure isn't whitelisted, keep whichever one has the fewest breaks in the H / L loops and antigen; if equal amounts of breaks then keep the one with the shortest antigen; if equal antigen lengths, keep the previously accepted one!
                
                if Ppdb in whiteList:
                    worseRepeat = True
                    break

                if currScore < pastScore:
                    worseRepeat = True
                    break
                elif currScore > pastScore:
                    toRemoveIfAdding.append(ii)
                    if len(toRemoveIfAdding) > 1:
                        worseRepeat = True
                        break
                    else:
                        betterRepeat = True
                        ii += 1
                        continue # if a better repeat may not be the best so keep checking
                else:
                    if paLen <= caLen:
                        worseRepeat = True
                        break
                    else:
                        toRemoveIfAdding.append(ii)
                        if len(toRemoveIfAdding) > 1:
                            worseRepeat = True
                            break
                        else:
                            betterRepeat = True
                            ii += 1
                            continue 
            else:
                ii += 1

        if worseRepeat == True:
            continue

        for ii in beforeDatePpInfoAndSeqsList:
            PhSeq = Seq(ii['hLoopSeqs'])
            PlSeq = Seq(ii['lLoopSeqs'])
            prevPDB = ii['currPDB']
            hAlighments = []
            bestHalighn = 0.0
            lAlighments = []
            bestLalighn = 0.0

            # check if the sequence for the antibody is at X% seq ID or less to a past example in the before list - so has this antibody's "binding mode" been seen before?

            pAbScore = 0.0

            if len(i['hLoopSeqs']) > 0 and len(ii['hLoopSeqs']) > 0:
                if use_globalSeqID:
                    pAbScore = seqPer(ChSeq,PhSeq,'global')
                else:
                    pAbScore = seqPer(ChSeq,PhSeq,'local')

            if len(i['lLoopSeqs']) > 0 and len(ii['lLoopSeqs']) > 0:
                if use_globalSeqID: 
                    bestLalign = seqPer(ClSeq,PlSeq,'global')
                else:
                    bestLalign = seqPer(ClSeq,PlSeq,'local')

                if bestLalign > pAbScore:
                    pAbScore = bestLalign

            if (pAbScore >= pAbSeqCut):
                prevSimilarityCount+=1
                keepI = False
                break

        # if everything's checked and there're no matches, start an accepted structure entry
        if keepI == True:

            acceptedStructInfo = i

            if len(toRemoveIfAdding) > 1: # if you'd be knocking out multiple structures to only include one more, that's NG
                continue

            acceptedAbProtStructs.append(acceptedStructInfo)

            # remove past worst duplicate, if relevant

            if len(toRemoveIfAdding) > 0:
                acceptedAbProtStructsNew = []
                iii = 0
                while iii < len(acceptedAbProtStructs):
                    if iii not in toRemoveIfAdding:
                        acceptedAbProtStructsNew.append(acceptedAbProtStructs[iii])
                    iii += 1
                acceptedAbProtStructs = acceptedAbProtStructsNew

    return acceptedAbProtStructs, prevSimilarityCount, knockoutWithinCount


