#covert vasp output to raw file
echo "converting vasp output to raw file ..."
cmdall "python ~/soft/simpy/tools/deepmd/vasp2raw.py"

mkdir ../data

cat ./*/box.raw > ../data/box.raw
cat ./*/coord.raw > ../data/coord.raw
cat ./*/virial.raw > ../data/virial.raw
cat ./*/force.raw > ../data/force.raw
cat ./*/energy.raw > ../data/energy.raw

cp ./type.raw ../data/
cp ~/soft/simpy/tools/deepmd/raw_to_set.sh ../data/



