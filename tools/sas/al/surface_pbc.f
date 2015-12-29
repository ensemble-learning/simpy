**********************************************************************
**********************************************************************
*                                                                    *
*     Surface analysis                                               *
*                                                                    *
*                                                                    *
**********************************************************************
        implicit real*8 (a-h,o-z)
        dimension c(100000,3),label(100000)
        real*8 pbc(3)
        real*8 d_cut ! cut_off value
        real*8 dx, dy, dz
        logical pbc_flag, DEBUG
        character*2 at,label
        character*50 coord_file

        DEBUG = .FALSE.
        ! read configurations
        read(*,*) coord_file
        open(9, file=coord_file, status='old')
        read(*,*) pbc_flag
        write(*,*) pbc_flag
        read(*,*) d_cut
        write(*,*) d_cut

        open(unit=2,file='atoms_surface',status='replace')

        read (9,*) nat
        if (DEBUG) then
        write (6,*) ' nat=',nat
        end if
        dnat=dfloat(nat)
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

        ! appply pbc
        if (pbc_flag) then
        x = x - floor(x/pbc(1))*pbc(1)
        y = y - floor(y/pbc(2))*pbc(1)
        z = z - floor(z/pbc(3))*pbc(1)
        end if
        c(i,1)=x
        c(i,2)=y
        c(i,3)=z
        label(i)=at
        enddo

        do i=1,nat
        dnormx=0.d0
        dnormy=0.d0
        dnormz=0.d0
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
        if (dist.lt.3.d0) then
        jcount=jcount+1
        dnormx=dnormx+dx
        dnormy=dnormy+dy
        dnormz=dnormz+dz
c       write (6,5) i,j,label(i),label(j),dist
c       write (6,4) label(i),c(i,1)-c(j,1),c(i,2)-c(j,2),c(i,3)-c(j,3) 
        endif
 2      continue
        enddo
        dnorm=dsqrt(dnormx**2+dnormy**2+dnormz**2)/dfloat(jcount)
c       write (6,*) ' dnorm',label(i),dnorm
        if (dnorm.ge.d_cut) then
        write (2,1) label(i),c(i,1),c(i,2),c(i,3)
        endif
        enddo

        close(unit=2,status='keep')

 1      format(1x,a2,3f20.6)
 4      format(1x,a2,1x,3f20.5)
 5      format(1x,2i5,1x,a2,1x,a2,1x,f20.8)

        stop
        end
