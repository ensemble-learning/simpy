#!/bin/sh

if [ $# -ne 2 ]
then
  echo "Usage : product_trace [control.file] [Mol Name]"
  exit
fi

./searchid $1 $2 > searchlist

Nline=$(wc -l searchlist)
Nline=${Nline%%searchlist}

mkdir $2
for ((i=1 ; i<Nline; i++))
do
  ID=$(head -$i searchlist | tail -1)
  echo "*Tracing $2 ($ID)"
  rxntrace $1 $ID > $2/$2_$ID.rxn
  rxn2net $2/$2_$ID.rxn $2/$2_$ID.temp 0.2
  nethighlight $1 $2/$2_$ID.temp $ID > $2/$2_$ID.net
  rm $2/$2_$ID.temp
  todos $2/$2_$ID.net
done

