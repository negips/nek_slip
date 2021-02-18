!======================================================================
!     Routines for introducing third component in a 2d simulation
!     Author: Prabal S. Negi
!     Right now we assume the third component is homogeneous.
!     Later, a fourier dependence can be added for the linearized solve.
!
!====================================================================== 
!-----------------------------------------------------------------------

      subroutine init_3ds()

      implicit none

      include 'SIZE'
      include 'SOLN'

      real bfz_3ds(lx1,ly1,lz1,lelv)
      logical if3d_3ds
      common /solv_3ds/ bfz_3ds,if3d_3ds

      integer icalld
      save icalld
      data icalld /0/

      integer nxyz,ntot1
      

      nxyz  = lx1*ly1*lz1
      ntot1 = nxyz*nelv

      if3d_3ds = .true.

!     Need to initialize some variables
!     V3MASK

      call copy(v3mask,v1mask,ntot1)

!     Velocities can be initialized from useric. 


      return
      end subroutine init_3ds
!----------------------------------------------------------------------

      subroutine buildrhs_3ds

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'DEALIAS'

      integer icalld
      save icalld
      data icalld /0/

      real bfz_3ds(lx1,ly1,lz1,lelv)
      logical if3d_3ds
      common /solv_3ds/ bfz_3ds,if3d_3ds


      if (.not.if3d_3ds) then
        return
      endif


      call makef_3ds


      return
      end subroutine buildrhs_3ds

!----------------------------------------------------------------------

      subroutine makef_3ds

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'   ! ifield

      real bfz_3ds(lx1,ly1,lz1,lelv)
      logical if3d_3ds
      common /solv_3ds/ bfz_3ds,if3d_3ds


      integer nxyz,ntot1

      real ta1,ta2,ta3
      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)


      nxyz = lx1*ly1*lz1
      ntot1 = nxyz*nelv

      if (.not.if3d_3ds) then
        call rzero(bfz_3ds,ntot1)    
        return
      endif


!     Build user defined forcing for uz
      call makeuf_3ds

!     Need to see what to modify here
      ifield = 1
      if3d = .true.
      if (filterType.eq.2) call make_hpf
!     hpf field stored in ta3
      call xaddcol3(bfz_3ds,ta3,bm1,ntot1)
      if3d = .false. 

!     Put vx,vy,vz on Dealiased grid (rst form)
!     Standard nek routine. Don't need to change anything (yet) 
!      call set_convect_new(vxd,vyd,vzd,vx,vy,vz)
      call advab_3ds

!     Just leaving it here but we don't need this.
!      call admeshv      ! subroutine not defined yet

      if (iftran) call makeabf_3ds
      if ((iftran.and..not.ifchar).or.
     $    (iftran.and..not.ifnav.and.ifchar)) call makebdf_3ds


      return
      end subroutine makef_3ds

!----------------------------------------------------------------------

      subroutine makeuf_3ds

!     Compute and add: (1) user specified forcing function (FX,FY,FZ)

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      include 'PARALLEL'

      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)

      integer ntot1

      real bfz_3ds(lx1,ly1,lz1,lelv)
      logical if3d_3ds
      common /solv_3ds/ bfz_3ds,if3d_3ds

      ntot1 = lx1*ly1*lz1*nelv

      time = time-dt
      call rzero(bfz_3ds,ntot1)

      do 100 iel=1,nelv
         ielg = lglel(iel)
         do 100 k=1,lz1
         do 100 j=1,ly1
         do 100 i=1,lx1
            call nekasgn (i,j,k,iel)
            call userf   (i,j,k,ielg)
            bfz_3ds(i,j,k,iel) = ffz
 100  continue

      call col2  (bfz_3ds,bm1,ntot1)
      time = time+dt

      return
      end subroutine makeuf_3ds

!----------------------------------------------------------------------

      subroutine advab_3ds
!
!     Eulerian scheme, add convection term to forcing function 
!     at current time step.

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'

      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)

      real bfz_3ds(lx1,ly1,lz1,lelv)
      logical if3d_3ds
      common /solv_3ds/ bfz_3ds,if3d_3ds

      ntot1 = lx1*ly1*lz1*nelv

      call convop  (ta3,vz)
      call subcol3 (bfz_3ds,ta3,bm1,ntot1)

      return
      end subroutine advab_3ds
!-----------------------------------------------------------------------

      subroutine makeabf_3ds

!     Sum up contributions to kth order extrapolation scheme.

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'

      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)

      real bfz_3ds(lx1,ly1,lz1,lelv)
      logical if3d_3ds
      common /solv_3ds/ bfz_3ds,if3d_3ds

      ntot1 = lx1*ly1*lz1*nelv

      ab0 = ab(1)
      ab1 = ab(2)
      ab2 = ab(3)
      call add3s2 (ta3,abz1,abz2,ab1,ab2,ntot1)
      call copy   (abz2,abz1,ntot1)
      call copy   (abz1,bfz,ntot1)
      call add2s1 (bfz_3ds,ta3,ab0,ntot1)
      if(.not.iflomach) call col2 (bfz_3ds,vtrans,ntot1)

      return
      end subroutine makeabf_3ds

!-----------------------------------------------------------------------
      subroutine makebdf_3ds
!
!     Add contributions to F from lagged BD terms.
!
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'

      common /scrns/ ta1(lx1,ly1,lz1,lelv)
     $ ,             ta2(lx1,ly1,lz1,lelv)
     $ ,             ta3(lx1,ly1,lz1,lelv)
     $ ,             tb1(lx1,ly1,lz1,lelv)
     $ ,             tb2(lx1,ly1,lz1,lelv)
     $ ,             tb3(lx1,ly1,lz1,lelv)
     $ ,             h2 (lx1,ly1,lz1,lelv)

      real bfz_3ds(lx1,ly1,lz1,lelv)
      logical if3d_3ds
      common /solv_3ds/ bfz_3ds,if3d_3ds


      ntot1 = lx1*ly1*lz1*nelv
      const = 1./dt

      if(iflomach) then
        call cfill(h2,const,ntot1)
      else
        call cmult2(h2,vtrans(1,1,1,1,ifield),const,ntot1)
      endif

      call col3(tb3,vz,bm1,ntot1)
      call cmult(tb3,bd(2),ntot1)

      do 100 ilag=2,nbd
         if (ifgeom) then

            call col3(ta3,vzlag(1,1,1,1,ilag-1),bm1lag(1,1,1,1,ilag-1),
     $                                          ntot1)
            call cmult(ta3,bd(ilag+1),ntot1)
         else
            call col3(ta3,vzlag(1,1,1,1,ilag-1),bm1,ntot1)
            call cmult(ta3,bd(ilag+1),ntot1)

         endif
         call add2  (tb3,ta3,ntot1)
 100  continue
      call xaddcol3(bfz_3ds,tb3,h2,ntot1)     

      return
      end subroutine makebdf_3ds
!-----------------------------------------------------------------------

      subroutine lagvel_3ds

!     Keep old velocity field(s)

      implicit none 

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'

      integer ilag,ntot1

      ntot1 = lx1*ly1*lz1*nelv

      do 100 ilag=3-1,2,-1
        call copy (vzlag (1,1,1,1,ilag),vzlag (1,1,1,1,ilag-1),ntot1)
 100  continue

      call copy (vzlag,vz,ntot1)

      return
      end subroutine lagvel_3ds
!----------------------------------------------------------------------
 




!----------------------------------------------------------------------

!----------------------------------------------------------------------





