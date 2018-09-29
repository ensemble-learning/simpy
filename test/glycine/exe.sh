#!/bin/bash
if [ -f xmolout ] 
then
	rm xmolout
fi
if [ -f moldyn.vel ] 
then
	mv moldyn.vel vels
fi
if [ -f ffield ] 
then
	cp ffield ffieldss
fi
if [ -f fort.83 ]
then
	cp fort.83 83s
fi
if [ -f fort.4 ] 
then
	cp fort.4 4s
fi
if [ -f fort.41 ] 
then
	cp fort.41 41s
fi
if [ -f fort.91 ] 
then
	cp fort.91 91s
fi
if [ -f fort.59 ] 
then
	cp fort.59 59s
fi
if [ -f fort.58 ] 
then
	cp fort.58 58s
fi
if [ -f fort.57 ] 
then
	cp fort.57 57s
fi
if [ -f fort.71 ] 
then
	cp fort.71 71s
fi
if [ -f fort.72 ] 
then
	cp fort.72 72s
fi
if [ -f fort.73 ] 
then
	cp fort.73 73s
fi
if [ -f fort.74 ] 
then
	cp fort.74 74s
fi
if [ -f fort.76 ] 
then
	cp fort.76 76s
fi
if [ -f 13s ] 
then
	cp 13s 13s2
fi
if [ -f fort.13 ] 
then
	cp fort.13 13s
fi
if [ -f fort.99 ] 
then
	cp fort.99 99s
fi
if [ -f fort.79 ] 
then
	cp fort.79 79s
fi
if [ -f fort.7 ]  
then
	cp fort.7 7s
fi
if [ -f fort.8 ]  
then
	cp fort.8 8s
fi
if [ -f fort.81 ] 
then
	cp fort.81 81s
fi
if [ -f fort.97 ] 
then
	cp fort.97 97s
fi
if [ -f fort.90 ] 
then
	cp fort.90 90s
fi
if [ -f fort.98 ] 
then
	cp fort.98 98s
fi
if [ -f end.geo ] 
then
	rm end.geo
fi
if [ -f molfra.out ] 
then
	rm molfra.out
fi
if [ -f dipole.out ] 
then
	rm dipole.out
fi
touch fort.5a
touch moldyn.0a
touch molsav.0a
rm *.out
rm fort.*
rm moldyn.0*
rm molsav.0*
if [ ! -f control ]
then
echo 'Missing control'
exit
fi
if  [ -f geo ]
then
	cp geo fort.3
else
	echo 'Missing geo'
exit
fi
if  [ -f ffield ]  
then
	cp ffield fort.4
else
	echo 'Missing ffield'
	exit
fi
if  [ -f ranfile ]  
then
	cp ranfile fort.35
else
	echo 'Created ranfile in unit 35'
	echo '234535.1' > fort.35
fi
if [ -f iopt ] 
then
	cp iopt fort.20
else
	echo 'Created iopt in unit 20; assume normal run'
	echo '  0   0: Normal run   1: Force field optimization' > fort.20
fi
if [ -f outres ] 
then
	cp outres fort.9
else
	echo 'Touched unit 9 (outres)'
	touch fort.9
fi
if [ -f inilp ] 
then
cp inilp fort.2
fi
if [ -f params ] 
then
cp params fort.21
fi
if [ -f koppel ] 
then
cp koppel fort.22
fi
if [ -f koppel2 ] 
then
cp koppel2 fort.23
fi
if [ -f tregime ]
then
	cp tregime fort.19
fi
if [ -f restraint ] 
then
	cp restraint fort.18
fi
if [ -f restraintt ]
then
	cp restraintt fort.28
fi
if [ -f restraintv ] 
then
	cp restraintv fort.38
fi
if [ -f vels ] 
then
	cp vels moldyn.vel
fi
reac_lg > run.log
exit

