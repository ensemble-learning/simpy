C
C Adapted from the program VTOTAV that was  
C distributed with VASP 4 and is of uncertain origin.
C
C V1 - Original distro version
C V2 - Various updates from JLF Da Silva and A Walsh (2007-2009)
C      http://people.bath.ac.uk/aw558 (a.walsh@bath.ac.uk)
C
C Reads in the 3D electrostatic potential (LOCPOT) from VASP
C which is generated with LVHAR= .TRUE. in the INCAR file.
C
C Outputs the 1D potential averaged along one of the lattice
C vectors, which can be used, e.g., to compute the vacuum level.
C
C Note: to obtain the Hartree + Exc contributions, the flag
C LVTOT = .TRUE. should be set instead. In VASP 4, this flag
C (confusingly)
C provided the Hatree contributions only.
C
C To do: 
C 1. Automatically align values using the calculated vacuum level.
C 2. Check if dipole corrections are required (slope of potential in the
C                                              C vacuum).
C 3. Change axis units to Angstrom (for all IDIR).
C 4. Update to read artbitrary number of ion types. 
C
      PROGRAM VTOTAV
      PARAMETER(NGXM=256,NOUTM=1024)
      CHARACTER*80 HEADER
      DIMENSION VLOCAL(NGXM*NGXM*NGXM),VAV(NOUTM)

      WRITE(*,*) '========================'
      WRITE(*,*) 'Workfunction V2 (2009)'
      WRITE(*,*) '========================'
      WRITE(*,*) ''
      WRITE(*,*) 'Which direction to average along? (1=x;2=y;3=z)'
      READ(*,*) IDIR
      IDIR=MOD(IDIR+20,3)+1

C Read in LOCPOT
      OPEN(20,FILE='LOCPOT',STATUS='OLD',ERR=1000)
      READ(20,'(A)',ERR=1000,END=1000) HEADER
      READ(20,'(A)',ERR=1000,END=1000) HEADER
      READ(20,'(A)',ERR=1000,END=1000) HEADER
      READ(20,'(A)',ERR=1000,END=1000) HEADER
      READ(20,'(A)',ERR=1000,END=1000) HEADER
      READ(20,'(A)',ERR=1000,END=1000) HEADER
      READ(20,'(A)',ERR=1000,END=1000) HEADER
      I=0; II=0; III=0; IIII=0
      READ(HEADER,*,ERR=12,END=12) I,II,III,IIII
12    NIONS=I+II+III+IIII
      READ(20,'(A)',ERR=1000,END=1000) HEADER
      WRITE(*,*) 'Number of ions:',NIONS
      DO 10 I=1,NIONS
         READ(20,*,ERR=1000,END=1000) RDUM1,RDUM2,RDUM3
   10 CONTINUE
      WRITE(*,*) '(i) Positions read...'
      READ(20,'(A)',ERR=1000,END=1000) HEADER
      READ(20,*,ERR=1000,END=1000) NGX,NGY,NGZ
      NPLWV=NGX*NGY*NGZ
      IF (IDIR.EQ.1) NOUT=NGX
      IF (IDIR.EQ.2) NOUT=NGY
      IF (IDIR.EQ.3) NOUT=NGZ
      IF (NPLWV.GT.(NGXM*NGXM*NGXM)) THEN
         WRITE(*,*) 'NPLWV .GT. NGXM**3 (',NPLWV,').'
         STOP
      ENDIF
      IF (NOUT.GT.NOUTM) THEN
         WRITE(*,*) 'NOUT .GT. NOUTM (',NOUT,').'
         STOP
      ENDIF
      READ(20,*,ERR=1000,END=1000) (VLOCAL(I),I=1,NPLWV)
      WRITE(*,*) '(ii) 3D potential read...'
      CLOSE(20)

C Average 3D Potential
      DO 20 I=1,NOUTM
   20 VAV(I)=0.
      SCALE=1./FLOAT(NPLWV/NOUT)
C      WRITE(*,*) SCALE
      IF (IDIR.EQ.1) THEN
         DO 150 IX=1,NGX
            DO 100 IZ=1,NGZ
             DO 100 IY=1,NGY
               IPL=IX+((IY-1)+(IZ-1)*NGY)*NGX
               VAV(IX)=VAV(IX)+VLOCAL(IPL)*SCALE
  100       CONTINUE
  150    CONTINUE
      ELSE IF (IDIR.EQ.2) THEN
         DO 250 IY=1,NGY
            DO 200 IZ=1,NGZ
             DO 200 IX=1,NGX
               IPL=IX+((IY-1)+(IZ-1)*NGY)*NGX
               VAV(IY)=VAV(IY)+VLOCAL(IPL)*SCALE
  200       CONTINUE
  250    CONTINUE
      ELSE IF (IDIR.EQ.3) THEN
         DO 350 IZ=1,NGZ
            DO 300 IY=1,NGY
             DO 300 IX=1,NGX
               IPL=IX+((IY-1)+(IZ-1)*NGY)*NGX
               VAV(IZ)=VAV(IZ)+VLOCAL(IPL)*SCALE
  300       CONTINUE
  350    CONTINUE
      ELSE
         WRITE(*,*) 'Wrong value of IDIR',IDIR
         STOP
      ENDIF
      Z0 = 0.d0
      OPEN(20,FILE='vplanar.txt')
      WRITE(20,*) NOUT,IDIR
C      WRITE(20,*) NOUT,IDIR, C3
      DO 500 I=1,NOUT
      WRITE(20,'(I6,2X,E18.11)') I,VAV(I)
C      WRITE(20,'(I6,2X,F10.6,E18.11)') I,(C3/(NOUT-1))*(I-1),VAV(I)
  500 CONTINUE
      WRITE(*,*) ''
      WRITE(*,*) 'Done. Always check for convergence!'
      CLOSE(20)

      STOP
 1000 WRITE(*,*) 'Error opening or reading file LOCPOT!'
      STOP
      END
