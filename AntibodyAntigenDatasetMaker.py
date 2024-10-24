# code for parsing down antibody-protein structures since a date D, removing those with matches already in SAbDab (and so the PDB) with antigens at X% seq ID or less to their antigen /and/ an antibody at Y% seq ID or less to their antibody

import argparse
import datetime
import os


from src.utils import seqPer,isScFv,tripleToSingle,seqResParser,anarciParser,abAgStructs2Seqs
from src.aadam_io import get_summary_file, create_output_files
from src.data_processing import filterbysimilarity, prefilterdataset

def unit_float(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{x} is not a valid float value")
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError(f"{x} is out of range. Must be between 0 and 1.")
    return x

def add_arguments(parser):
    parser.add_argument("-i", "--inputDb", dest="inputDb", help="path to the input database, which should be a complete download of all antibodies in complex with protein antigens, and have a summary file ending in 'summary.tsv'")
    parser.add_argument("-d", "--date", dest="splitdate", help="the date on or after which you want to make your test set; format YYYY/MM/DD")
    parser.add_argument("-cd", "--cutoffDate", dest="cutoffDate", help="a date past which you won't accept structures, useful for making training sets; format YYYY/MM/DD")
    parser.add_argument("-c1", "--abCompSeqCut", dest="abCompSeqCut", help="percent sequence ID minimum to disallow an ab-ag complex, based on past ab-ag complexes, comparing H loops to H loops and L loops to L loops")
    parser.add_argument("-c2", "--withinDatasetCut", dest="withinDatasetCut", help="percent sequence ID minimum used for knocking out complexes that are too similar w/ in the database, based on H or L or any antigen chain to any other antigen chain sequence ID")
    parser.add_argument("-m", "--methodsAllowed", dest="methodsAllowed", help="methods allowed, separated by a comma, like 'X-RAY DIFFRACTION,ELECTRON MICROSCOPY'")
    parser.add_argument("-w", "--whiteList", dest="whiteList", help="a comma separated list of pdb files that you want to definitely be included in the dataset (will still be filtered by resolution, sequence identity, method, etc. - this just ensures when structures 'knock each other out' the whitelisted ones will be given preference to remain. Useful for adding on new structures to a previously-made dataset)'")
    parser.add_argument("-o", "--outDir", dest="outDbPath", help="")
    parser.add_argument("-r", "--resCut", dest="resCut", help="defaults to 100")
    parser.add_argument("-cs", "--cutoffStrict", dest="cutoffStrict", help="make the cutoffs for sequence ID within the group strict, by applying them to the antigen and antibody loop sequences individually, instead of in combination (i.e. strict cutoffs would not allow an antibody-antigen complex with the same antigen as a previously accepted structure but different antibody in, while standard / not strict cutoffs would allow it)")
    parser.add_argument("-nx", "--nx", dest="nx", help="strip unnatural residues from the ends, and don't allow sequences with unnatural residues in the middle")
    parser.add_argument("-g", "--globalSeqID", dest="globalSeqID", help="use global sequence ID, number of matches divided by the length of the shorter sequence, to get sequence ID %. This can make the cutoff more stringent in practice, as a shorter sequence can map well to a much larger non-similar one by allowing many gaps. The default (when this flag is not given) is to do local alignment with gap penalites, then apply the sequence identity cutoffs to the best local alignment, if it's over 30 residues long.")
    parser.add_argument(
        "-mf", "--minAtomSeqresFraction", 
        dest="minAtomSeqresFraction", 
        help="Skips structures in the final dataset, where the minimum of all fractions between the respective atom and seqres sequence for the H chain, L chain and antigen chain(s), is less than the provided threshold. Values have to lie in the interval [0,1].", 
        type=unit_float, 
        default=0
    )


def main():
    parser = argparse.ArgumentParser()
    add_arguments(parser)
    arguments = parser.parse_args()


    splitdate = arguments.splitdate.split("/")
    splitdate = datetime.datetime(int(splitdate[0]), int(splitdate[1]), int(splitdate[2]))

    methodsList = arguments.methodsAllowed.split(",")

    if arguments.cutoffDate:
        cdList = arguments.cutoffDate.split("/")
        cutoffDate = datetime.datetime(int(cdList[0]), int(cdList[1]), int(cdList[2]))
    else:
        cutoffDate = None


    if not os.path.exists(arguments.outDbPath):
        os.makedirs(arguments.outDbPath)

    if not arguments.resCut:
        resCut = 100
    else:
        resCut = arguments.resCut

    if arguments.whiteList:
        whiteList = arguments.whiteList.split(",")
    else:
        whiteList = []

    print(f"Using strict cutoff: {bool(arguments.cutoffStrict)}")

    percentWinDatasetCut = float(arguments.withinDatasetCut)/100
    pAbSeqCut = float(arguments.abCompSeqCut)/100.0

    sumFileName = get_summary_file(arguments.inputDb)

    onOrAfterDateList, beforeDateProtPepComplexesList, cutoffFilterCount, methodFilterCount, qualFilterCount = prefilterdataset(sumFileName,methodsList,resCut,splitdate,cutoffDate)

    print("cutoffFilterCount:")
    print(cutoffFilterCount)
    print("qualFilterCount:")
    print(qualFilterCount)
    print("methodFilterCount:")
    print(methodFilterCount)


    # get sequences etc. from pdb files - each entry gets a dict
    onOrAfterDateInfoAndSeqsList = abAgStructs2Seqs(onOrAfterDateList,arguments.inputDb,strip_res_x=arguments.nx)

    print(f"onOrAfterDateInfoAndSeqsList len: {len(onOrAfterDateInfoAndSeqsList)}")

    if len(beforeDateProtPepComplexesList) > 0:
        beforeDatePpInfoAndSeqsList = abAgStructs2Seqs(beforeDateProtPepComplexesList,arguments.inputDb,strip_res_x=arguments.nx)
    else:
        beforeDatePpInfoAndSeqsList = []

    print(f"beforeDatePpInfoAndSeqsList len: {len(beforeDatePpInfoAndSeqsList)}")

    # filter
    acceptedAbProtStructs, prevSimilarityCount, knockoutWithinCount = filterbysimilarity(onOrAfterDateInfoAndSeqsList, percentWinDatasetCut, pAbSeqCut, beforeDatePpInfoAndSeqsList,disallow_res_x=arguments.nx, use_globalSeqID=arguments.globalSeqID, whiteList=whiteList,use_cutoffStrict=arguments.cutoffStrict, minAtomSeqresFraction=arguments.minAtomSeqresFraction)

    print("prevSimilarityCount: ")
    print(prevSimilarityCount)
    print("knockoutWithinCount: ")
    print(knockoutWithinCount)


    # save & print the accepted structures
    print("done! writing files...")

    create_output_files(acceptedAbProtStructs, outDbPath=arguments.outDbPath, inputDb=arguments.inputDb)

    # print number of structures
    print(f"base number of structure past date: {len(onOrAfterDateInfoAndSeqsList)}")
    print(f"total number of structures accepted at end: {len(acceptedAbProtStructs)}")






if __name__ == "__main__":
    main()