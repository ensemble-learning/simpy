import os

def parseFolder():
    msdfile = []
    for i in os.listdir('.'):
        if i.endswith(".msd"):
            msdfile.append(i)
    return msdfile

if __name__ == "__main__":
    o = open("checkout.dfi", 'w')
    o.write("""#DFF:DFFJOB
    $MODEL
""")
    outputmsd = parseFolder()
    for i in outputmsd:
        o.write(i + '\n')
    o.write("""	$END
	$MODEL = $SETTYPE
		RULER = STRUCTURE
	$END
	$DATABASE = $OPEN
		FILENAME = D:\DATA\DffForceField\T_amber_ver3\fit\group2\ff\2009-12-09\AMBER20.ppf
	$END
	$DATABASE = $CHECKOUT
		OUTPUT = RESULT.ppf
	$END
	$FORCEFIELD = $SETCHARGE
	$MODEL = $SAVEALL
#DFF:ENDJOB
""")
    o.close()