! Download from http://people.bath.ac.uk/td222/research/surface_area/ortho/
! see T. Duren, F. Millange, G. Ferey, K.S. Walton, R.Q. Snurr, J. Phys. Chem.
! C, 2007, 111, 15360. 

Program accessible_surface_area
!
! This program calculates the accessible surface area by using a simple Monte
! Carlo integration. To run the program the following files are needed:
! 1) an input file containing the name of the file containing the diameters of 
!    the atoms, the name of the file containing the coordinates of the 
!    framework atoms, the diameter of the probe, the number of random trials
!    around each framework atom, the lenghts of the unit cells (a, b, c) and
!    the crystal density.
! 2) A file containing the diameters of the framework atoms. The program is
!    written to read in a file that contains the coordinates in the format of 
!    a music molecule file (compare 
!    http://zeolites.cqe.northwestern.edu/Music/music.html). It can easily be
!    adapted to read in any cartesian coordinates by changing the read(10,*)
!    comments.
!
! Things that you might have to change:
! - The program uses double precision (real) variables. Depending on your 
!   compiler, you might have to change the value for realkind (E.g. 8 is the
!   value for the intel compiler, 2 is another widely used definition)
! - The maximum number of framework atoms 5000 at the moment. Change the value
!   of max_no if your unit cells contains more atoms.
!
! Changes to source code
!
! 18/05/10  Tina Duren   corrected serious bug that led to wrong results for
!                        non-cubic unit cell. When checking for overlap in the 
!                        pbc test one of the distances was divided by the wrong
!                        unit cell parameter.
! 08/08/08  Tina Duren   corrected bug that let to a wrong intialisation of 
!                        random number generator. The seed initially provided
!                        needs to be negative otherwise the same string
!                        of random numbers is produced, no matter what the
!                        seed was. Thanks to Francois-Xavier Coudert for 
!                        pointing this out.
! 08/08/08 Tina Duren    corrected slight oddity which didn't lead to wrong 
!                        results.Theta is now created between in the 
!                        interval (0,2pi) 
!
IMPLICIT NONE

Integer, parameter    :: realkind = 8
Integer, parameter    :: max_no = 1000000
Integer, parameter    :: max_cn = 48
Real, parameter       :: pi=3.14159

Integer               :: iseed
Integer               :: N, ntypes, nsample, ncount 
Integer               :: i1, i2, i3
Integer               :: i, j, k, l

Character(len = 3)    :: atom
Character(len = 10)   :: atomname(1:max_no), atomtype(1:max_no)
Character(len = 50)   :: atom_file, coord_file

Real(kind = realkind) :: sigmatype(1:max_no), atomsigma(1:max_no)
Real(kind = realkind) :: nlist_buf 
Real(kind = realkind) :: sfrac_cut_off 
Real(kind = realkind) :: nlist_cut_off
Real(kind = realkind) :: nneighbor(1:max_no)
Real(kind = realkind) :: nlist(1:max_no, 1:max_cn)
Real(kind = realkind) :: x(1:max_no), y(1:Max_no), z(1:max_no)
Real(kind = realkind) :: xl, yl, zl, rho_crys
Real(kind = realkind) :: xpoint, ypoint, zpoint
Real(kind = realkind) :: dx, dy, dz, dist2, dist
Real(kind = realkind) :: phi, theta, costheta
Real(kind = realkind) :: dprobe
Real(kind = realkind) :: stotal, sfrac, sjreal, stotalreduced
Real(kind = realkind) :: trash
Real(kind = realkind) :: ran0

Logical               :: deny, match

open(unit=99, file="frac_atom.dat")

!
! Seed for random number generator, change if you want to start from a 
! different random number.
!
iseed=-52907
!
! to properly initialise the random number generator, a negative seed 
! is needeed
!
If (iseed > 0) iseed = -1 * iseed

Read(*,*) atom_file
Read(*,*) coord_file
Open(9, file=atom_file, status='old')
Open(10, file=coord_file, status='old')


Read(*,*) dprobe ! diameter of the probe
Read(*,*) Nsample ! Number of samples per sphere
Read(*,*) XL, YL, ZL ! Cell parameters in A
Read(*,*) rho_crys   ! chrystal density in g / cm3 
Read(*,*) nlist_buf ! Neighbour list cut-off
Read(*,*) sfrac_cut_off ! sfrac cut off

ntypes = 0
DO
 Read(9,*) atom
 IF(Trim(atom) == 'EOF') EXIT
 ntypes = ntypes + 1
END DO

REWIND(9)

DO i = 1, ntypes
  read(9,*) atomtype(i), sigmatype(i)
END DO
!
! The first line of the coordinate file contains the number of molecules
!
read(10,*) N

DO i=1, N
!
! The format of the coordinate file corresponds to a music molecule file.
! Therefore there is some information that is not necessary for the surface
! area program.
! i1: continuous integer (not important for this program)
! x(i), y(i), z(i): x,y,z coordinates in Angstrom
! atomname(i): name of framework atom. Note that this is the name 
!              not the chemical symbol
! trash: a real number that doesn't play a role in this program
! i2, i3: two integer numbers that don't play a role in this program
!
! Change the following read statement if you want to use a different
! file format. But ensure that you end up with the coordinates in A and
! the name of the framework atoms!
!
   read(10,*) i1, x(i), y(i), z(i), atomname(i) , trash, i2, i3
!
! Shift unit cell so that it lies at the origin
!
   if(x(i)<0.0) x(i)=x(i)+XL
   if(x(i)>=XL) x(i)=x(i)-XL
   if(y(i)<0.0) y(i)=y(i)+YL
   if(y(i)>=YL) y(i)=y(i)-YL
   if(z(i)<0.0) z(i)=z(i)+ZL
   if(z(i)>=ZL) z(i)=z(i)-ZL
   match=.False.
   DO j = 1, ntypes
      IF(atomname(i) == atomtype(j)) Then
         atomsigma(i)=sigmatype(j)+dprobe
         match=.True.
         exit
      END IF
   END DO

   IF(.Not.match) then
      write(*,*) 'Could not find match for atom: ', i, ' ', atomname(i)
      write(*,*) 'The name is either read in incorrectly or does not exist'
      Write(*,*) 'in the list of available atoms in ', TRIM(atom_file)
      stop
   END IF
END DO
!
Do i = 1, N ! Loop over all framework atoms
   ncount=0
!  write(*,*) i, ncount
   Do j = 1, N ! Loop over all framework atoms

      IF (i == j) cycle

      dx = x(i) - x(j)
      dx = dx - XL*int(2.0*dx/XL)
      dy = y(i) - y(j)
      dy = dy - YL*int(2.0*dy/YL)
      dz = z(i) - z(j)
      dz = dz - ZL*int(2.0*dz/ZL)

      dist2 = dx*dx+dy*dy+dz*dz
      dist = sqrt(dist2)
!     write(*,*) dist, i, j
      nlist_cut_off = nlist_buf *(atomsigma(i) - dprobe)
!     write(*,*) "cutoff is ", nlist_cut_off
      IF (dist < nlist_cut_off) THEN
         ncount = ncount + 1
         nlist(i, ncount) = j
      END IF
   nneighbor(i) = ncount
!   write(*,*) i, ncount
         
  END DO
END DO

Write(*,*) 
Write(*,*) 'Calculating the accessible surface area for the following input parameters'
Write(*,*) '-------------------------------------------------------------------------'
Write(*,*) 'File with framework coordinates: ',TRIM(coord_file)
Write(*,*) 'File with atom diameters: ', atom_file
Write(*,*) 'Probe diameter in A: ',dprobe

!
! Main sampling cycle
!
stotal=0.0

DO i = 1, N ! Loop over all framework atoms

   ncount=0
   Write(*,*) "Processing", i, "out of", N

   DO j = 1, Nsample ! Number of random trials around each framework atom
!
! Generate random vector of length 1
!
! First generate phi 0:+2pi
!
      phi = ran0(iseed)*2.0*pi
!
! Generate cosTheta -1:1 to allow for even distribution of random vectors
! on the unit sphere.See http://mathworld.wolfram.com/SpherePointPicking.html
! for further explanations
!
      costheta = 1 - ran0(iseed) * 2.0
      theta = Acos(costheta)
      xpoint = sin(theta)*cos(phi)
      ypoint = sin(theta)*sin(phi)
      zpoint = costheta

! Make this vector of length (sigma+probe_diameter)/2.0 

      xpoint = xpoint*atomsigma(i)/2.0
      ypoint = ypoint*atomsigma(i)/2.0
      zpoint = zpoint*atomsigma(i)/2.0

! Translate the center of coordinate to the particle i center and apply PBC

      xpoint = xpoint + x(i)
      ypoint = ypoint + y(i)
      zpoint = zpoint + z(i)

      if(xpoint < 0.0) xpoint = xpoint + XL
      if(xpoint >= XL) xpoint = xpoint - XL
      if(ypoint < 0.0) ypoint = ypoint + YL
      if(ypoint >= YL) ypoint = ypoint - YL
      if(zpoint < 0.0) zpoint = zpoint + ZL
      if(zpoint >= ZL) zpoint = zpoint - ZL

!Test for overlap

      deny = .False.

!     write(*,*) nneighbor(i)
      DO k = 1, nneighbor(i)
         l = nlist(i, k)

         dx = xpoint - x(l)
         dx = dx - XL*int(2.0*dx/XL)
         dy = ypoint - y(l)
         dy = dy - YL*int(2.0*dy/YL)
         dz = zpoint - z(l)
         dz = dz - ZL*int(2.0*dz/ZL)

         dist2 = dx*dx+dy*dy+dz*dz

         IF(sqrt(dist2) < 0.999*atomsigma(k)/2.0) THEN
            deny=.True.
	!    if(i==3) write(*,*) k, sqrt(dist2), 0.5*atomsigma(k)
            exit
         END IF
      END DO

      IF(deny) cycle

      ncount = ncount + 1

   END DO

! Fraction of the accessible surface area for sphere i

   sfrac = real(ncount) / real(Nsample)
   Write(99, FMT="(F10.6)") sfrac 

! Surface area for sphere i in real units (A^2)

   IF (sfrac > sfrac_cut_off) THEN
      sjreal = pi*atomsigma(i)*atomsigma(i)*sfrac
      stotal = stotal + sjreal
   ENDIF
  
END DO

! Converting stotal on Surface per Volume


stotalreduced = stotal/(XL*YL*ZL)*1.E4

! Report results
     
Write(*,*)
Write(*,'(A,F12.2)') 'Accessible surface area in A^2: ', stotal
Write(*,'(A,F12.2)') 'Accessible surface area per volume in m^2/cm^3: ', stotalreduced
Write(*,'(A,F12.2)') 'Accessible surface area per mass in m^2/g: ', &
            stotalreduced / rho_crys
Write(*,*)

close(99)

END PROGRAM accessible_surface_area



!----------------FUNCTIONS-------------------------------------

FUNCTION ran0(idum)
!
! Random number generator from W.H. Press et al, Numerical Recipes in
! FORTRAN, Cambridge University Press, 1992
!
 INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
 REAL(kind = 8) ran0,AM,EPS,RNMX
 PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
   NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
 INTEGER j,k,iv(NTAB),iy
 
 SAVE iv,iy
 DATA iv /NTAB*0/, iy /0/
 if (idum.le.0.or.iy.eq.0) then
     idum=max(-idum,1)
     do 11 j=NTAB+8,1,-1
        k=idum/IQ
        idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11   continue
     iy=iv(1)
endif
k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
if (idum.lt.0) idum=idum+IM
j=1+iy/NDIV
iy=iv(j)
iv(j)=idum
ran0=min(AM*iy,RNMX)
return

END FUNCTION ran0



