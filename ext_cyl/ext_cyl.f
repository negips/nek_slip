c-----------------------------------------------------------------------
C
C  USER SPECIFIED ROUTINES:
C
C     - boundary conditions
C     - initial conditions
C     - variable properties
C     - local acceleration for fluid (a)
C     - forcing function for passive scalar (q)
C     - general purpose routine for checking errors etc.
C
c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,eg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer e,f,eg
c     e = gllel(eg)

      udiff =0.
      utrans=0.
      return
      end
c-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,eg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer e,f,eg
c     e = gllel(eg)


c     Note: this is an acceleration term, NOT a force!
c     Thus, ffx will subsequently be multiplied by rho(x,t).


      ffx = 0.0
      ffy = 0.0
      ffz = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,eg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer e,f,eg
c     e = gllel(eg)

      qvol   = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userchk
      include 'SIZE'
      include 'TOTAL'

      real x0(3)
      save x0
      data x0 /3*0/
      
      integer bIDs(1), iobj_wall(1)

c     define objects for surface integrals
      if (istep.eq.0) then
         bIDs(1) = 1
         call create_obj(iobj_wall(1),bIDs,1)
      endif

      call slp_mark_faces() 
      
      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)
c     NOTE ::: This subroutine MAY NOT be called by every process
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      ux=1.0
      uy=0.0
      uz=0.0
      temp=0.0
      return
      end
c-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      ux=1 + 1.0*rand()
      uy=2.0*rand()
      uz=0.0
      temp=0

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat
      include 'SIZE'
      include 'TOTAL'

c     call platform_timer(0) ! not too verbose
c     call platform_timer(1) ! mxm, ping-pong, and all_reduce timer

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2
      include 'SIZE'
      include 'TOTAL'

c     param(66) = 4.   ! These give the std nek binary i/o and are 
c     param(67) = 4.   ! good default values

c     mark faces for object definition
      nface = 2*ndim
      do iel=1,nelt
         do iface = 1, nface
            if (cbc(iface,iel,1) .eq. 'W  ') then
               boundaryID(iface,iel) = 1
            endif
         enddo 
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat3
      include 'SIZE'
      include 'TOTAL'
c
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
