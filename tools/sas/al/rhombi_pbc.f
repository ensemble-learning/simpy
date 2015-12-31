c
c       code to produce an analysis of ONLY SURFACE atoms in terms of rhombi
c       warning : only surface atoms should be given as c(10000,3) coordinates and label(10000) elements
c
        implicit real*8 (a-h,o-z)
        dimension c(100000,3),label(100000)
        dimension neighb(100000,15),nneighb(100000),irhombi(1200000,4)
        character*2 at,label

        real*8 pbc(3)
        real*8 dx, dy, dz
        logical pbc_flag, DEBUG
        character*50 coord_file

        DEBUG = .FALSE.
        open(unit=2,file='rhombi',status='replace')

        ! read configurations
        read(*,*) coord_file
        open(9, file=coord_file, status='old')
        read(*,*) pbc_flag
        write(*,*) pbc_flag

        read (9,*) nat
c       write (6,*) ' nat=',nat
        if (nat.gt.100000) write (6,*) ' ERRORE NAT !!!!!! '

        ! read pbc if pbc_flat is true
        if (pbc_flag) then
        write(*,*) "Read pbc"
        read (9,*) pbc
        write(*,*) "a = ", pbc(1), "b= ", pbc(2), "c= ", pbc(3)
        else
        write(*,*) "no pbc"
        read (9,*)
        end if

c       write (6,*) ' '
        do i=1,nat
        read (9,*) at,x,y,z
c       write (6,4) at,x,y,z
        label(i)=at

        ! appply pbc
        if (pbc_flag) then
        x = x - floor(x/pbc(1))*pbc(1)
        y = y - floor(y/pbc(2))*pbc(1)
        z = z - floor(z/pbc(3))*pbc(1)
        end if

        c(i,1)=x
        c(i,2)=y
        c(i,3)=z
        enddo

        do i=1,nat
        jcount=0
        do j=1,nat
        if (i.eq.j) go to 2
        dx = c(i,1)-c(j,1)
        dy = c(i,2)-c(j,2)
        dz = c(i,3)-c(j,3)
        ! apply pbc
        if (pbc_flag) then
        dx = dx - pbc(1)*int(2.0*dx/pbc(1))
        dy = dy - pbc(2)*int(2.0*dy/pbc(2))
        dz = dz - pbc(3)*int(2.0*dz/pbc(3))
        end if
        dist=dx*dx + dy*dy + dz*dz
        dist=dsqrt(dist)
        if (dist.lt.3.00d0) then
        jcount=jcount+1
c        write (6,5) label(i),i,label(j),j,dist  ! c
        neighb(i,jcount)=j
        endif
 2      continue
        enddo
c        write (6,*) i,jcount  ! c
        if (jcount.gt.15) write (6,*) ' ERRORE N NEIGHBORS !!!!!! '
        nneighb(i)=jcount
        enddo

        write (6,*) ' rhombi'

        nrhombi=0
        do i=1,nat
c        write (6,*) ' i,n =',i,nneighb(i),(neighb(i,m),m=1,nneighb(i))  ! c
        do j=1,nneighb(i)
        ij=neighb(i,j)
        nj=nneighb(ij)
c        write (6,*) 'j,ij',j,ij,nj,(neighb(ij,nn),nn=1,nj)  ! c
        do k=1,nneighb(i)
        ik=neighb(i,k)
        if (k.ne.j) then
c        write (6,*) ' 1 k,ik =',k,ik  ! c
        do l=1,nneighb(ij)
        jl=neighb(ij,l)
        if (jl.ne.i.and.ik.eq.jl) then
c        write (6,*) ' 1 l,jl =',l,jl  ! c
        nk=k
        nl=l
        kij=jl
        go to 6
        endif
        enddo
        endif
        enddo
        go to 7
 6      continue
        do k=1,nneighb(i)
        ik=neighb(i,k)
        if (k.ne.j.and.k.ne.nk) then
c        write (6,*) ' 2 k,ik =',k,ik  ! c
        do l=1,nneighb(ij)
        jl=neighb(ij,l)
        if (jl.ne.i.and.ik.eq.jl) then
c        write (6,*) ' 2 l,jl =',l,jl  ! c
        nrhombi=nrhombi+1
        irhombi(nrhombi,1)=i
        irhombi(nrhombi,2)=ij
        irhombi(nrhombi,3)=kij
        irhombi(nrhombi,4)=jl
        write (2,*) i,ij,kij,jl
        endif
        enddo
        endif
        enddo
 7      continue
        enddo
        enddo
        if (nrhombi.gt.1200000) write (6,*) ' ERRORE N RHOMBI !!!!!! '

 1      format(1x,3f20.5)
 4      format(1x,a2,1x,3f20.6)
 5      format(1x,a2,1x,i5,1x,a2,1x,i5,1x,f20.8)

        stop
        end
