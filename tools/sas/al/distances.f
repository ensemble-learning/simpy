        implicit real*8 (a-h,o-z)
        dimension c(100000,3),label(100000)
        character*2 at,label

        open(unit=2,file='surface_atoms',status='replace')

        read (5,*) nat
        write (6,*) ' nat=',nat
        dnat=dfloat(nat)
        if (nat.gt.100000) write (6,*) ' ERRORE NAT !!!!!! '
        read (5,*)
        write (6,*) ' '
        do i=1,nat
        read (5,*) at,x,y,z
        write (6,4) at,x,y,z
        c(i,1)=x
        c(i,2)=y
        c(i,3)=z
        label(i)=at
        enddo

        do i=1,nat
        jcount=0
        do j=1,nat
        if (i.eq.j) go to 2
        dist=(c(i,1)-c(j,1))**2+(c(i,2)-c(j,2))**2+(c(i,3)-c(j,3))**2
        dist=dsqrt(dist)
        if (dist.lt.3.d0) then
        jcount=jcount+1
c        write (6,5) i,j,label(i),label(j),dist
        endif
 2      continue
        enddo
c        write (6,*) i,jcount
        if (jcount.le.9) then
        write (2,1) label(i),c(i,1),c(i,2),c(i,3)
        endif
        enddo

        close(unit=2,status='keep')

 1      format(1x,a2,3f20.6)
 4      format(1x,a2,1x,3f20.5)
 5      format(1x,2i5,1x,a2,1x,a2,1x,f20.8)

        stop
        end
