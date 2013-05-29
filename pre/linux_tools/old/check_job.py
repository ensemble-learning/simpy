import os
import os.path

group = [
'acetylene',
'acid',
'alcohols',
'aldehydes',
'alkane',
'alkene',
'amides',
'amine',
'benzenes',
'branched_alkanes',
'ether',
'fluoro',
'glycol',
'nitriles',
'nitroalkanes',
#'noble_gas',
#'silanes',
'siloxanes',
'small_molecules',
'sulfides',
'sulfur',
'thiols',
]

def write_job(model, counter):
    if counter == 1:
        print "%30s"%model + "     complete"
    elif counter == 0:
        print "%30s"%model + "     fail for some reasons"

def check_job(folder):
    for i in group:
        currentFolder = os.path.join(folder, i)
        dirList = os.listdir(currentFolder)
        for j in dirList:
            counter = 0
            modelFolder = os.path.join(currentFolder, j)
            if os.path.isdir(modelFolder):
                mdFile = os.path.join(modelFolder, "md0.log")
                if os.path.exists(mdFile):
                    f = open(mdFile, 'r')
                    for m in f:
                        if m.startswith("Finished mdrun"):
                            counter = 1
            else:
                pass
            write_job(j, counter)

if __name__ == "__main__":
    rootPath = os.getcwd()
    check_job(rootPath)
