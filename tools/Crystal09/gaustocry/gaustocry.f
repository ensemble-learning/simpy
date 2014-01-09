      IMPLICIT double precision(A-H,O-Z)
      call  GAUSTOCRY
      stop
      end
      SUBROUTINE GAUSTOCRY
      IMPLICIT double precision(A-H,O-Z)
      character(len=4) mark
      character(len=2) symb,shelt
      CHARACTER*2 NDN(6)
      logical ls,lp,ld,lf,exist
      CHARACTER*2 SYMBAT(0:93)
      character*80 a
      parameter(lim015=100,lim042=100)
      dimension CHE(LIM015),EXX(LIM042,lim015),C1(LIM042,lim015),
     *C2(LIM042,lim015),LAN(LIM015),LAT(LIM015),scalef(lim015)
      dimension noccup(5,92)
      data noccup(5,1:92)/92*0.D0/
      data noccup(1:4,1:92)/
     * 1, 0, 0, 0, 2, 0, 0, 0, 3, 0, 0, 0, 4, 0, 0, 0, 4, 1, 0, 0,
     * 4, 2, 0, 0, 4, 3, 0, 0, 4, 4, 0, 0, 4, 5, 0, 0, 4, 6, 0, 0,
     * 5, 6, 0, 0, 6, 6, 0, 0, 6, 7, 0, 0, 6, 8, 0, 0, 6, 9, 0, 0,
     * 6,10, 0, 0, 6,11, 0, 0, 6,12, 0, 0, 7,12, 0, 0, 8,12, 0, 0,
     * 8,12, 1, 0, 8,12, 2, 0, 8,12, 3, 0, 8,12, 4, 0, 8,12, 5, 0,
     * 8,12, 6, 0, 8,12, 7, 0, 8,12, 8, 0, 8,12, 9, 0, 8,12,10, 0,
     * 8,13,10, 0, 8,14,10, 0, 8,15,10, 0, 8,16,10, 0, 8,17,10, 0,
     * 8,18,10, 0, 9,18,10, 0,10,18,10, 0,10,18,11, 0,10,18,12, 0,
     *10,18,13, 0,10,18,14, 0,10,18,15, 0,10,18,16, 0,10,18,17, 0,
     *10,18,18, 0,10,18,19, 0,10,18,20, 0,10,19,20, 0,10,20,20, 0,
     *10,21,20, 0,10,22,20, 0,10,23,20, 0,10,24,20, 0,11,24,20, 0,
     *12,24,20, 0,12,24,21, 0,12,24,21, 1,12,24,20, 3,12,24,20, 4,
     *12,24,20, 5,12,24,20, 6,12,24,20, 7,12,24,20, 8,12,24,20, 9,
     *12,24,20,10,12,24,20,11,12,24,20,12,12,24,20,13,12,24,20,14,
     *12,24,21,14,12,24,22,14,12,24,23,14,12,24,24,14,12,24,25,14,
     *12,24,26,14,12,24,27,14,12,24,28,14,12,24,29,14,12,24,30,14,
     *12,25,30,14,12,26,30,14,12,27,30,14,12,28,30,14,12,29,30,14,
     *12,30,30,14,13,30,30,14,14,30,30,14,14,30,31,14,14,30,32,14,
     *14,30,31,16,14,30,30,18/
      DATA SYMBAT/   'XX','H ','HE','LI','BE','B ','C ','N ',
     1'O ','F ','NE','NA','MG','AL','SI','P ','S ','CL','AR',
     2'K ','CA','SC','TI','V ','CR','MN','FE','CO','NI','CU',
     3'ZN','GA','GE','AS','SE','BR','KR','RB','SR','Y ','ZR',
     4'NB','MO','TC','RU','RH','PD','AG','CD','IN','SN','SB',
     5'TE','I ','XE','CS','BA','LA','CE','PR','ND','PM','SM',
     6'EU','GD','TB','DY','HO','ER','TM','YB','LU','HF','TA',
     7'W ','RE','OS','IR','PT','AU','HG','TL','PB','BI','PO',
     8'AT','RN','FR','RA','AC','TH','PA','U ','QC'/
      DATA NDN/'S ','SP','P ','D ','F ','G '/
      iiin=7
      idum=8
      iout=9
      INQUIRE(FILE='GAUSSIAN.BAS',EXIST=EXIST)
      IF(.NOT.EXIST)THEN
      print *,'file GAUSSIAN.BAS not found'
      stop
      ELSE
      OPEN(UNIT=iiin,FILE='GAUSSIAN.BAS',FORM='FORMATTED',STATUS='OLD')
      ENDIF
      OPEN(UNIT=iout,FILE='CRYSTAL.BAS',FORM='FORMATTED',STATUS='NEW')
c cerco la stringa BASIS per iniziare a leggere
10    read(iiin,100,end=999,err=998)mark
100   format(a4)
      if(mark.ne.'BASI')then   !if 1
      goto 10
      else                        !else if 1
11    read(iiin,1,end=21,err=998)a   
1     format(a80)
      if(a(2:3).eq.'**')then
      write(idum,*)'**'
      else
      write(idum,1)a
      endif
      goto 11
c legge il simbolo dell'atomo, fino ad un record bianco
21    rewind (idum)
      close (iiin)
20    read(idum,101,end=999)symb
101   format(1x,a2)
      if(symb.eq.'  ')then
      print *,'fine lettura basi'
      close (iout)
      return
      endif
      do iz=0,93
      if(symbat(iz).eq.symb)exit
      enddo
      print *,'leggo base gaussian per atomo ',symb,iz
c trovato il numero atomico, inizia a leggere gli shell
c per ogni atomo azzero tutto
      nat=iz
      nshell=0
30    read(idum,102)shelt,iprimg,scale
102   FORMAT(1x,A2,I3,f6.2)  
      if(shelt.eq.'**')then
      write(6,*)'****'
c scrive base per quell'atomo
      print *,'scrivo base CRYSTAL style  '
      ls=.false.
      lp=.false.
      ld=.false.
      lf=.false.
      ns=noccup(1,nat)
      np=noccup(2,nat)
      nd=noccup(3,nat)
      nf=noccup(4,nat)
      ng=noccup(5,nat)
      if(ns.ne.0)ls=.true.
      if(np.ne.0)lp=.true.
      if(nd.ne.0)ld=.true.
      if(nf.ne.0)lf=.true.
      nsc=0
      npc=0
      ndc=0
      nfc=0
      write(iout,200)nat,nshell
      write(06,202)nat,nshell,ns,np,nd,ng
200   format(2i4)
202   format(' Z = ',i2,' n. shell = ',i3,' el s ',i3,' el p ',i3,
     *' el sd',i3,' el sg',i3)
      do la=1,nshell
      ityp=lat(la)
      che(la)=0.d0
      select case(ityp)
      case(5)  !  g functions - always empty
      che(la)=0.d0
      case (0) ! orbitale s
      if(ls)then
      if(nsc+2.gt.ns)then
      che(la)=1
      ls=.false.
      else
      che(la)=2
      nsc=nsc+2
      if(nsc.eq.ns)ls=.false.
      endif
      endif
      case (1) ! sp shell - split in s AO and p AO
      nsp=0
      if(ls)then
      if((nsc+2).gt.ns)then
      nsp=1
      ls=.false.
      else
      nsp=nsp+2
      if(nsc+2.eq.ns)ls=.false.
      endif
      endif
      if(lp)then
      if(npc+6.gt.np)then
      nsp=nsp+(np-npc)
      lp=.false.
      else
      nsp=nsp+6
      npc=npc+6
      if(npc.eq.np)lp=.false.
      endif
      endif
      che(la)=nsp
      case (2)
      if(lp)then
      if(npc+6.gt.np)then
      che(la)=(np-npc)
      lp=.false.
      else
      che(la)=6.d0
      npc=npc+6
      if(npc.eq.np)lp=.false.
      endif
      endif
      case (3)
      if(ld)then
      if(ndc+10.gt.nd)then
      che(la)=(nd-ndc)
      ld=.false.
      else
      che(la)=10
      ndc=ndc+10
      if(ndc.eq.nd)ld=.false.
      endif
      endif
      case(4)
      if(lf)then
      if(nfc+14.gt.nf)then
      che(la)=(nf-nfc)
      lf=.false.
      else
      che(la)=14
      nfc=nfc+14
      if(nfc.eq.nf)lf=.false.
      endif
      endif
      end select
      write(iout,201)lat(la),lan(la),che(la),scalef(la)
      write(6,201)lat(la),lan(la),che(la),scalef(la)
      do k=1,lan(la)
      if(lat(la).eq.1)then
      write(iout,105)EXX(K,la),C1(K,la),C2(K,la)
      write(6,105)EXX(K,la),C1(K,la),C2(K,la)
      else
      write(iout,105)EXX(K,la),C1(K,la)
      write(6,105)EXX(K,la),C1(K,la)
      endif
      enddo
201   format(' 0 ',2i4,2f10.3)
      enddo
c cerca un altro atomo
      if(ls)print *,' not enough s functions - s electrons',ns
      if(lp)print *,' not enough p functions - p electrons',np
      if(ld)print *,' not enough d functions - d electrons',nd
      if(lf)print *,' not enough f functions - f electrons',nf
      goto 20
      else     
      write(6,102)shelt,iprimg,scale
c legge ancora uno shell per lo stesso atomo
      nshell=nshell+1
      lan(nshell)=iprimg 
      do iz=1,6
      if(ndn(iz).eq.shelt)exit
      enddo
      lat(nshell)=iz-1
      scalef(nshell)=scale
c  lat = iz
      lati=iz-1
      do k=1,iprimg
      if(lati.eq.1)then
      read(idum,105)EXX(K,nshell),C1(K,nshell),C2(K,nshell)
      write(6,105)EXX(K,nshell),C1(K,nshell),C2(K,nshell)
      else
      read(idum,105)EXX(K,nshell),C1(K,nshell)
      write(6,105)EXX(K,nshell),C1(K,nshell)
      endif
      enddo
      ENDIF ! fine shell dell'atomo
      goto 30
      endif ! iniziava con simbolo atomo, fine con ****
c      endif! iniziava con BASIS
      return
999   print *,'fine dati'
      return
998   print *,'errore lettura dati'
103   FORMAT(A2,' 0')
104   FORMAT(5X,A4,I4,f5.1)
105   FORMAT(1x,4G20.9)
      return
      END
