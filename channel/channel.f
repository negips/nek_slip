c- constants -----------------------------------------------------------

! #define tSTATSTART uparam(1) /* start time for averaging */
! #define tSTATFREQ  uparam(2) /* output frequency for statistics */

c data extraction along wall normal direction
! #define INTP_NMAX 200 /* number of sample points */
! #define XCINT 1.0     /* x coordinate of 1D line*/
! #define ZCINT 1.0     /* z coordinate of 1D line */

c mesh dimensions
! #define BETAM 2.4     /* wall normal stretching parameter */
! #define PI (4.*atan(1.))
! #define XLEN (2.*PI)
! #define ZLEN PI
! #define NUMBER_ELEMENTS_X 16
! #define NUMBER_ELEMENTS_Y 12
! #define NUMBER_ELEMENTS_Z 8

c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer e

      utrans = 1.
      udiff  = param(2)

      if (ifield .eq. 2) then
         e = gllel(ieg)
         udiff = param(8)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      ffx = 0.0 
      ffy = 0.0
      ffz = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      qvol =  0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userchk

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'SOLN'
      include 'INPUT'
      include 'TSTEP'
!      include 'TOTAL'

      real x0(3)
      data x0 /0.0, 0.0, 0.0/ 
      save x0

      integer icalld
      save    icalld
      data    icalld /0/

      real atime,timel
      save atime,timel

      integer i,j,n,m
      integer nfaces

      real ut1(lx1,ly1,lz1,lelt)
      real ut2(lx1,ly1,lz1,lelt)
      real ut3(lx1,ly1,lz1,lelt)
      integer flowtype(lelt)

      common /testvel1/ ut1,ut2,ut3,flowtype

      real ut4(lx1,ly1,lz1,lelt)
      real ut5(lx1,ly1,lz1,lelt)
      real ut6(lx1,ly1,lz1,lelt)

      common /testvel2/ ut4,ut5,ut6

      real ut7(lx1,ly1,lz1,lelt)
      real ut8(lx1,ly1,lz1,lelt)
      real ut9(lx1,ly1,lz1,lelt)

      common /testvel3/ ut7,ut8,ut9

      real ymean

      real vlsum        ! function
      real facevals(lx1,ly1)
      integer intfaces(2**ldim,lelt)
      integer intels(lelt)

      real eps

      real offst

      n     = nx1*ny1*nz1

      if (icalld.eq.0) then

        offst = 0.2

        call opzero(ut1,ut2,ut3)
     
        do i=1,nelv
          ymean = vlsum(ym1(1,1,1,i),n)/n
          if (ymean.gt.0) then
            flowtype(i)=1
          else
            flowtype(i)=2
          endif
          ymean = flowtype(i) + 0.
          call cadd(ut1(1,1,1,i),ymean,n)
        enddo  
     

        call copy(ut2,ut1,n*nelv)
        call dsavg(ut2)

!        call outpost(ut1,ut2,ut3,pr,t,'   ')


!       Mark the faces we wish to change
        nfaces = 2**ndim

        call izero(intels,nelv)
        call izero(intfaces,nelv*nfaces)

        eps = 1.0e-12

        do i=1,nelv
          if (flowtype(i).eq.1) then    
            do j=1,nfaces
              call facexs(facevals,ut2(1,1,1,i),j,0)
              ymean = vlsum(facevals,lx1*lz1)/(lx1*lz1)
              if (abs((ymean-1.5)).lt.eps) then
                intels(i) = 1
                intfaces(j,i) = 1
                write(6,*) i,j
              endif
            enddo
          endif
        enddo 


        call opzero(ut1,ut2,ut3)
        call rzero(ut3,n*nelv)
        do i=1,nelv
          if (intels(i).eq.1) then
            do j=1,nfaces
              if (intfaces(j,i).eq.1) then
                write(6,*) i,j
                call facev(ut1,i,j,offst,nx1,ny1,nz1)
              endif
            enddo  
          endif
        enddo      


!        call copy(ut3,ut1,n*nelv)
!        call dsavg(ut1)
        call outpost(ut1,ut2,ut3,pr,t,'   ')

        icalld = icalld + 1


!       Seems to be working.
!       Edge values added correctly
      endif

!      call exitt 
      
      
      if (istep.gt.0) then
!        call outpost(ut4,ut5,ut6,usrdiv,t,'ut2')
!        call outpost(vx,vy,vz,pr,t,'   ')
!        call exitt
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      if (y.lt.0) ux = -1.0
      if (y.gt.0) ux = 1.0

      uy = 0.
      uz = 0.

      return
      end
c-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer idum
      save    idum 
      data    idum / 0 /

      real C, k, kx, ky

      ux = 1.0 + (1.0e-0)*rand()
      uy = 0.0
      uz = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat   ! This routine to modify element vertices
      include 'SIZE'      ! _before_ mesh is generated, which 
      include 'TOTAL'     ! guarantees GLL mapping of mesh.

      n = nelv * 2**ldim
!      xmin = glmin(xc,n)
!      xmax = glmax(xc,n)
!      ymin = glmin(yc,n)
!      ymax = glmax(yc,n)
!      zmin = glmin(zc,n)
!      zmax = glmax(zc,n)
!
!      xscale = XLEN/(xmax-xmin)
!      yscale = 1./(ymax-ymin)
!      zscale = ZLEN/(zmax-zmin)

      pi = 4.*atan(1.0)

      write(6,*) n,nelv,lelv

      do j=1,nelv
      do i=1,2**ldim
         xc(i,j) = pi*(xc(i,j)+1.0)
!         yc(i,1) = yscale*yc(i,1)
!         yc(i,1) = tanh(BETAM*(2*yc(i,1)-1))/tanh(BETAM)
!         zc(i,1) = zscale*zc(i,1)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2   ! This routine to modify mesh coordinates


      include 'SIZE'
      include 'TOTAL'


!      call outpost(vx,vy,vz,pr,t,'   ')

!      do iel=1,nelt
!      do ifc=1,2*ndim
!         if (cbc(ifc,iel,1) .eq. 'W  ') boundaryID(ifc,iel) = 1 
!         cbc(ifc,iel,2) = cbc(ifc,iel,1) 
!         if (cbc(ifc,iel,1) .eq. 'W  ') cbc(ifc,iel,2) = 't  '
!      enddo
!      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat3
      include 'SIZE'
      include 'TOTAL'

!      param(54) = -1  ! use >0 for const flowrate or <0 bulk vel
                      ! flow direction is given by (1=x, 2=y, 3=z) 
!      param(55) = 1.0 ! flowrate/bulk-velocity 

      return
      end
c-----------------------------------------------------------------------

c automatically added by makenek
      subroutine usrdat0() 

      return
      end

c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end

c automatically added by makenek
      subroutine userqtl

      call userqtl_scig

      return
      end
