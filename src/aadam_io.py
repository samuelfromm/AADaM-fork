import sys
import glob
import os
import pickle
import shutil

def get_summary_file(inputDb):
    # load in summary file
    sumFilePath = inputDb + "/*.tsv"
    #print(sumFilePath)
    sumFileName = glob.glob(sumFilePath)
    if len(sumFileName) !=1:
        sys.exit("You have the wrong number of summary files, there should only be 1. Check your abag dir")
    sumFileName = sumFileName[0]

    return sumFileName


def create_output_files(acceptedAbProtStructs, outDbPath,inputDb):
    with open(outDbPath + "/fullDb.pkl", 'wb') as f:
        pickle.dump(acceptedAbProtStructs, f)

    f = open(outDbPath + "/lightDb.txt", 'w')
    f.write("pdb ID,A chain(s),H chain(s),L chain(s)\n")
    for i in acceptedAbProtStructs:
        f.write(i['currPDB'] + "," + i['currAg'] + "," + i['currAbH'] + "," + i['currAbL'] + "\n")

    structPath = outDbPath + "/structures"
    if not os.path.exists(structPath):
        os.makedirs(structPath)
    for i in acceptedAbProtStructs:
        currFilePath = inputDb + "/" + str(i['currPDB']) + ".pdb"
        shutil.copyfile(currFilePath, structPath + "/" + str(i['currPDB']) + ".pdb")

    # write complex fastas

    if not os.path.exists(outDbPath + "/complexFastas"):
        os.makedirs(outDbPath + "/complexFastas")

    for i in acceptedAbProtStructs:
        fileName = outDbPath + "/complexFastas/" + i['currPDB'] + "complex.fasta"
        f = open(fileName, 'w')
        count = 1
        for ag in i['aSeq']:
            f.write("> Antigen " + str(count) + "\n")
            f.write(ag + "\n")
            count += 1
        if i['hSeq'] == i['lSeq']:
            f.write("> Antibody H and L seq\n")
            f.write(i['hSeq'] + "\n")
        else:
            if i['hSeq'] != '':
                f.write("> Antibody H\n")
                f.write(i['hSeq'] + "\n")
            if i['lSeq'] != '':
                f.write("> Antibody L\n")
                f.write(i['lSeq'] + "\n")
        f.close()

    # write complex fastas colon seperated

    compColSepPath = outDbPath + "/complexFastasColonSep"
    if not os.path.exists(compColSepPath):
        os.makedirs(compColSepPath)

    for i in acceptedAbProtStructs:
        fileName = compColSepPath + "/" + i['currPDB'] + "complexColSep.fasta"
        f = open(fileName, 'w')
        f.write("> Complex \n")
        startBool = True
        for ag in i['aSeq']:
            if startBool:
                f.write(ag)
                startBool = False
            else:
                f.write(":" + ag)

        if i['hSeq'] == i['lSeq']:
            f.write(":" + i['hSeq'])
        else:
            if i['hSeq'] != '':
                f.write(":" + i['hSeq'])
            if i['lSeq'] != '':
                f.write(":" + i['lSeq'])
        f.close()

    # write antigen fastas

    if not os.path.exists(outDbPath + "/antigenFastas"):
        os.makedirs(outDbPath + "/antigenFastas")

    for i in acceptedAbProtStructs:
        fileName = outDbPath + "/antigenFastas/" + i['currPDB'] + "antigen.fasta"
        f = open(fileName, 'w')
        count = 1
        for ag in i['aSeq']:
            f.write("> Antigen " + str(count) + "\n")
            f.write(ag + "\n")
            count += 1
        f.close()

    # write antigen fastas colon seperated

    antiColSepPath = outDbPath + "/antigenFastasColonSep"
    if not os.path.exists(antiColSepPath):
        os.makedirs(antiColSepPath)

    for i in acceptedAbProtStructs:
        fileName = antiColSepPath + "/" + i['currPDB'] + "antigenColSep.fasta"
        f = open(fileName, 'w')
        f.write("> Antigen \n")
        startBool = True
        for ag in i['aSeq']:
            if startBool:
                f.write(ag)
                startBool = False
            else:
                f.write(":" + ag)
        f.close()

    # write antibody fastas

    if not os.path.exists(outDbPath + "/antibodyFastas"):
        os.makedirs(outDbPath + "/antibodyFastas")

    for i in acceptedAbProtStructs:
        fileName = outDbPath + "/antibodyFastas/" + i['currPDB'] + "antibody.fasta"
        f = open(fileName, 'w')
        if i['hSeq'] == i['lSeq']:
            f.write("> Antibody H and L\n")
            f.write(i['hSeq'] + "\n")
        else:
            if i['hSeq'] != '':
                f.write("> Antibody H\n")
                f.write(i['hSeq'] + "\n")
            if i['lSeq'] != '':
                f.write("> Antibody L\n")
                f.write(i['lSeq'] + "\n")
        f.close()

    # write antibody fastas colon seperated

    abColSepPath = outDbPath + "/antibodyFastasColonSep"
    if not os.path.exists(abColSepPath):
        os.makedirs(abColSepPath)

    for i in acceptedAbProtStructs:
        fileName = abColSepPath + "/" + i['currPDB'] + "antibodyColSep.fasta"
        f = open(fileName, 'w')
        f.write("> Antibody \n")
        startBool = True
        if i['hSeq'] == i['lSeq']:
            f.write(i['hSeq'])
        else:
            if i['hSeq'] != '':
                f.write(i['hSeq'])
                startBool = False
            if i['lSeq'] != '':
                if startBool:
                    f.write(i['lSeq'])
                else:
                    f.write(":" + i['lSeq'])
        f.close()
