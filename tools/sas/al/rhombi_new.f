c
c       code to produce an analysis of ONLY SURFACE atoms in terms of rhombi
c       warning : only surface atoms should be given as c(10000,3) coordinates and label(10000) elements
c
        implicit real*8 (a-h,o-z)
        dimension c(100000,3),label(100000)
        dimension neighb(100000,15),nneighb(100000),irhombi(1200000,4)
        character*2 at,label

        read (5,*) nat
        write (6,*) ' nat=',nat
        if (nat.gt.100000) write (6,*) ' ERRORE NAT !!!!!! '
        read (5,*)
        write (6,*) ' '
        do i=1,nat
        read (5,*) at,x,y,z
        write (6,4) at,x,y,z
        label(i)=at
        c(i,1)=x
        c(i,2)=y
        c(i,3)=z
        enddo

        do i=1,nat
        jcount=0
        do j=1,nat
        if (i.eq.j) go to 2
        dist=(c(i,1)-c(j,1))**2+(c(i,2)-c(j,2))**2+(c(i,3)-c(j,3))**2
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
        write (6,*) i,ij,kij,jl
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
