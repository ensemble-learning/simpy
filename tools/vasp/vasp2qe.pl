#!/usr/bin/perl

#
# read POSCAR and OUTCAR from vasp run to produce force summary and/or geometry trajectory in xyz format
#

# set default values
$PrintT = 0;    # do not print trajectory
$PrintPW = 0;    # do not print trajectory
$PrintForce = 1;  # do print Force summary
$NTA=1;    # N. of replicas in A direction
$NTB=1;
$NTC=1;
# read command line
while ( $#ARGV >= 0 ) {
  $val = $ARGV[0];
  if ( $val eq "-T" ) { $PrintT = 1; }
  if ( $val eq "-TAB" ) { $PrintT = 1; $NTA=2; $NTB=2; }
  if ( $val eq "-TABC" ) { $PrintT = 1; $NTA=2; $NTB=2; $NTC=2; }
  if ( $val eq "-NoForce" ) { $PrintForce = 0; }
  if ( $val eq "-PrintPW" ) { $PrintPW = 1; }
  if ( $val eq "-PrintAXSF" ) { $PrintAXSF = 1; }
  if ( $val eq "-help" || $val eq "-h" ) {
    print "USAGE: vaspforce.pl [options]\n";
    print "       where options are:\n";
    print "              -T         Print coordinates to file trj.xyz (default=false)\n";
    print "              -TAB       Print coordinates to file trj.xyz, doubling cell size along A and B (default=false)\n";
    print "              -TABC      Print coordinates to file trj.xyz, doubling cell size along all axis (default=false)\n";
    print "              -PrintPW   Print PW style input chunk to pw.in (default=false)\n";
    print "              -PrintAXSF Print PW style AXSF to pw.axsf (default=false)\n";
    print "              -NoForce   Suppresses print of force summary (default=true)\n";
    print "              -help, -h  Print this summary and exit\n";
    exit(0);
  }
  shift @ARGV;
}

if ( $PrintAXSF == 1 ) {  # build atomic number table
  @SYMLIST = ("X ","H ","He",
            "Li","Be","B","C","N","O","F","Ne",
            "Na","Mg","Al","Si","P","S","Cl","Ar",
            "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
            "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
            "Cs","Ba","La",
                           "Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
                           "Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",);
  $NumAtomTable=87;
}

# read POSCAR to get labels and geometry constraints
$POSCARZIP=0;
if ( ! -f 'POSCAR' ) {
  if ( -f 'POSCAR.gz' ) {
    system('gunzip', 'POSCAR.gz');
    $POSCARZIP=1;
  }
}
open(IN,"<POSCAR");
$line = <IN>;   # title
$line = <IN>;   # scale factor
chomp($line); @list = split(" ",$line); 
$Scale = $list[0];
$line = <IN>;   # cell A vector
chomp($line); @list = split(" ",$line); 
$AX = $list[0]*$Scale;
$AY = $list[1]*$Scale;
$AZ = $list[2]*$Scale;
$line = <IN>;   # cell B vector
chomp($line); @list = split(" ",$line); 
$BX = $list[0]*$Scale;
$BY = $list[1]*$Scale;
$BZ = $list[2]*$Scale;
$line = <IN>;   # cell C vector
chomp($line); @list = split(" ",$line); 
$CX = $list[0]*$Scale;
$CY = $list[1]*$Scale;
$CZ = $list[2]*$Scale;
$line = <IN>;   # atomic types
chomp($line); @list = split(" ",$line); 
$NTypes = $#list + 1;
for ( $i=0; $i<$NTypes; $i++ ) { $Type[$i] = $list[$i]; }
# assign atomic numbers for AXSF
if ( $PrintAXSF ) {
  for ( $i=0; $i<$NTypes; $i++ ) {
    $t = $Type[$i];
    $val=-1;
    for ( $ia=0; $ia<$NumAtomTable; $ia++ ) {
      if ( $t eq $SYMLIST[$ia] ) { $val = $ia; }
    }
    if ( $val == -1 ) { 
      print "ERROR: could not assign atomic number to atom $t\n";
      exit(1);
    }
    $TypeNumber[$i] = $val;
  }
}
$line = <IN>;   # n. of atoms of each type
chomp($line); @list = split(" ",$line); 
$NA = 0;
for ( $it=0; $it<$NTypes; $it++ ) { 
  $NAtomType[$it] = $list[$it]; 
  $NA += $list[$it];
}
# assign atomic labels
$ia=0;
for ( $it=0; $it<$NTypes; $it++ ) {  # loop over types
  for ( $i=0; $i<$NAtomType[$it]; $i++ ) { 
    $Label[$ia] = $Type[$it]; 
    if ( $PrintAXSF ) { $AtomicNumber[$ia] = $TypeNumber[$it]; }
    $ia++;
  }
}
$line = <IN>;   # Selective dynamics?
if ( $line =~ /elective/ ) { 
  $Selective = 1; 
  $line = <IN>;   # Selective dynamics?
}
else { $Selective = 0; 
}
# build contraints array (1=free; 0=frozen)
for ( $ia=0; $ia<$NA; $ia++ ) {
  $FreezeX[$ia] = 1;
  $FreezeY[$ia] = 1;
  $FreezeZ[$ia] = 1;
} 
# will read geometry from OUTCAR
# if there are constraints (Selective Dynamics) I'll read them here, otherwise I'm done.
if ( $Selective ) {
  for ( $ia=0; $ia<$NA; $ia++ ) {
    $line = <IN>;   #  atomic coordinates
    chomp($line); @list = split(" ",$line); 
    if ( $list[3] eq 'F' ) { $FreezeX[$ia] = 0; }
    if ( $list[4] eq 'F' ) { $FreezeY[$ia] = 0; }
    if ( $list[5] eq 'F' ) { $FreezeZ[$ia] = 0; }
  }
}
close(IN);
if ( $POSCARZIP ) { system('gzip', 'POSCAR'); }
#print "Expecting $NA atoms\n";
#print "Labels will be\n";
#for ( $ia=0; $ia<$NA; $ia++ ) { print "$Label[$ia] "; }
#print "\n";
#print "Constraints\n";
#print "X: ";
#for ( $ia=0; $ia<$NA; $ia++ ) { print "$FreezeX[$ia] "; }
#print "\n";
#print "Y: ";
#for ( $ia=0; $ia<$NA; $ia++ ) { print "$FreezeY[$ia] "; }
#print "\n";
#print "Z: ";
#for ( $ia=0; $ia<$NA; $ia++ ) { print "$FreezeZ[$ia] "; }
#print "\n";

#
# read OUTCAR and extract geometries and forces
#
$OUTCARZIP=0;
if ( ! -f 'OUTCAR' ) {
  if ( -f 'OUTCAR.gz' ) {
    system('gunzip', 'OUTCAR.gz');
    $OUTCARZIP=1;
  }
}
open(IN,"<OUTCAR");
if ( $PrintT ) { open(OUT,">trj.xyz"); }
$NGEO=0; # geometries read
while ( $line = <IN> ) {
  if ( $line =~ " POSITION                                       TOTAL-FORCE" ) {
    $NGEO++;
    $line = <IN>;   # skip line of "-----"
    $Fmax = 0.0;
    for ( $ia=0; $ia<$NA; $ia++ ) {
      $line = <IN>;   # coordinates and forces
      chomp($line); @list = split(" ",$line); 
      $X[$ia] = $list[0]; 
      $Y[$ia] = $list[1]; 
      $Z[$ia] = $list[2]; 
      $FX[$ia] = $list[3]; 
      $FY[$ia] = $list[4]; 
      $FZ[$ia] = $list[5]; 
      $sum2 = $FreezeX[$ia]*$FX[$ia]*$FX[$ia]+$FreezeY[$ia]*$FY[$ia]*$FY[$ia]+$FreezeZ[$ia]*$FZ[$ia]*$FZ[$ia];
      if ( $sum2 > $Fmax ) { $Fmax = $sum2; }
    }
    # compute MAX force
    $Fmax = sqrt($Fmax);
    # read a few more lines to get the energy
    while ( $line !~ "free  energy   TOTEN" ) { $line = <IN>; }
    chomp($line); @list = split(" ",$line); 
    $Ene = $list[-2];
    if ( $PrintForce ) {
      printf("Geom. %3d: ENE = %15.8f eV  MAXF = %6.4f eV/Angst\n",$NGEO, $Ene, $Fmax);
    }
    if ( $PrintT ) {
      $NPrint = $NA*($NTA*$NTB*$NTC);
      print OUT "$NPrint\n";
      printf OUT ("Geom. %3d: ENE = %15.8f eV  MAXF = %6.4f eV/Angst\n",$NGEO, $Ene, $Fmax);
      for ( $ia=0; $ia<$NA; $ia++ ) {
        for ( $ca=0; $ca<$NTA; $ca++ ) {
          for ( $cb=0; $cb<$NTB; $cb++ ) {
            for ( $cc=0; $cc<$NTC; $cc++ ) {
              printf OUT ("%-4s  %10.5f %10.5f %10.5f\n",$Label[$ia], 
                     $X[$ia]+$ca*$AX+$cb*$BX+$cc*$CX, 
                     $Y[$ia]+$ca*$AY+$cb*$BY+$cc*$CY, 
                     $Z[$ia]+$ca*$AZ+$cb*$BZ+$cc*$CZ);
            }
          }
        }
      }
    } # end if
    if ( $PrintPW ) { #www
      open(OUTPW,">pw.in");
      $NPrint = $NA;
      print OUTPW "CELL_PARAMETERS (angstrom)\n";
      printf OUTPW (" %14.9f %14.9f %14.9f\n", $AX, $AY, $AZ);
      printf OUTPW (" %14.9f %14.9f %14.9f\n", $BX, $BY, $BZ);
      printf OUTPW (" %14.9f %14.9f %14.9f\n", $CX, $CY, $CZ);
      print OUTPW "ATOMIC_POSITIONS (angstrom)\n";
      for ( $ia=0; $ia<$NA; $ia++ ) {
        printf OUTPW ("%-4s  %14.9f %14.9f %14.9f\n",$Label[$ia], $X[$ia], $Y[$ia], $Z[$ia]);
      }
      close(OUTPW);
    } # end if
    if ( $PrintAXSF ) { #www
      open(OUTAXSF,">pw.axsf");
      $NPrint = $NA;
      print OUTAXSF "CRYSTAL\n";
      print OUTAXSF "\n";
      print OUTAXSF "PRIMVEC\n";
      printf OUTAXSF (" %14.9f %14.9f %14.9f\n", $AX, $AY, $AZ);
      printf OUTAXSF (" %14.9f %14.9f %14.9f\n", $BX, $BY, $BZ);
      printf OUTAXSF (" %14.9f %14.9f %14.9f\n", $CX, $CY, $CZ);
      print OUTAXSF "\n";
      print OUTAXSF "PRIMCOORD\n";
      print OUTAXSF "$NA 1\n";
      for ( $ia=0; $ia<$NA; $ia++ ) {
        printf OUTAXSF ("%-3d  %14.9f %14.9f %14.9f\n",$AtomicNumber[$ia], $X[$ia], $Y[$ia], $Z[$ia]);
      }
      close(OUTAXSF);
    } # end if

  }
}
close(IN);
if ( $OUTCARZIP ) { system('gzip', 'OUTCAR'); }
if ( $PrintT ) { close(OUT); }
