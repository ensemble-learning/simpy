"""
Csv tools

@log:
2013_10_12: initiation

"""
import os
import argparse
from utilities import parseBlock
from mytype import System, Molecule, Atom
from geo import Geo
from output_conf import toGeo, toXyz, toReaxLammps, toPdb, toGjf

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", default="csv", nargs="?", help="csv file name")
    args = parser.parse_args()

if __name__ == "__main__":
    main()
