!======================================================================
!     Routines for introducing third component in a 2d simulation
!     Author: Prabal S. Negi
!     Right now we assume the third component is homogeneous.
!     Later, a fourier dependence can be added for the linearized solve.
!
!====================================================================== 
!-----------------------------------------------------------------------
      subroutine fluidp_3ds (igeom)

!     Driver for perturbation velocity

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'SOLN'

      integer igeom

      do jp=1,npert

        if (nio.eq.0.and.igeom.eq.2) write(6,1) istep,time,jp
   1    format(i9,1pe14.7,' Perturbation Solve:',i5)

        call perturbv_3ds (igeom)

      enddo

      jp=0   ! set jp to zero, for baseline flow

      return
      end subroutine fluidp_3ds
c-----------------------------------------------------------------------

      subroutine perturbv_3ds (igeom)

      implicit none

!     Solve the convection-diffusion equation for the perturbation field, 
!     with projection onto a div-free space.


      include 'SIZE'
      include 'INPUT'
      include 'EIGEN'
      include 'SOLN'
      include 'TSTEP'
      include 'MASS'

      include 'TEST'

      real resv1,resv2,resv3
      real dv1,dv2,dv3
      real h1,h2

      common /scrns/  resv1 (lx1,ly1,lz1,lelv)
     $ ,              resv2 (lx1,ly1,lz1,lelv)
     $ ,              resv3 (lx1,ly1,lz1,lelv)
     $ ,              dv1   (lx1,ly1,lz1,lelv)
     $ ,              dv2   (lx1,ly1,lz1,lelv)
     $ ,              dv3   (lx1,ly1,lz1,lelv)
      common /scrvh/  h1    (lx1,ly1,lz1,lelv)
     $ ,              h2    (lx1,ly1,lz1,lelv)

      integer intype,igeom
      integer ntot1


      ifield = 1

      ntot1 = lx1*ly1*lz1*nelv

      if (igeom.eq.1) then

!        Old geometry, old velocity

         call makefp_3ds
         call lagfieldp_3ds

!        Add third component of convective term   
!         call advab_w_3ds   

      else
c
c        New geometry, new velocity
c
         intype = -1
         call sethlm_3dsp(h1,h2,intype)
         call cresvipp_3ds(resv1,resv2,resv3,h1,h2)

         if3d = .true.   
         call ophinv   (dv1,dv2,dv3,resv1,resv2,resv3,h1,h2,tolhv,nmxv)
         if3d = .false.

         call add2(vxp(1,jp),dv1,ntot1)
         call add2(vyp(1,jp),dv2,ntot1)
         call add2(vzp(1,jp),dv3,ntot1)


!        prabal. Debugging
!--------------------------------------------------    
         call copy(tmp1,vxp(1,jp),ntot1)
         call copy(tmp2,vyp(1,jp),ntot1)
         call copy(tmp3,vzp(1,jp),ntot1)
!--------------------------------------------------    

         call incomprp_3ds (vxp(1,jp),vyp(1,jp),vzp(1,jp),prp(1,jp))

      endif

      return
      end subroutine perturbv_3ds
!-----------------------------------------------------------------------

      subroutine initp_3ds()

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'

      include '3DS'

      integer icalld
      save icalld
      data icalld /0/

      integer nxyz,ntot1
      

      if3d_3ds = .true.

      if (.not.if3d_3ds) return

      if (if3d_3ds.and.npert.ne.2) then
        write(6,'(A7,1x,I2)') 'NPERT =', npert
        write(6,*) 'Need both real and imaginary parts'
        write(6,*) 'Ensure NPERT = 2 for 3D perturbation solve'
        call exitt
      endif

      nxyz  = lx1*ly1*lz1
      ntot1 = nxyz*nelv

      k_3dsp = 1.0            ! wavenumber 
     
      call init_pertfld_3ds() 

!     Need to initialize some variables
!     V3MASK
      call copy(v3mask,v1mask,ntot1)

!     Velocities can be initialized from useric. 


      return
      end subroutine initp_3ds
!----------------------------------------------------------------------

      subroutine init_pertfld_3ds

      implicit none

      include 'SIZE'
      include 'SOLN'    ! jp
      include 'TSTEP'   ! ifield

      integer i

      do i = 1,2
        jp = i
        call nekuic

        call dsavg(vxp(1,jp))
        call dsavg(vyp(1,jp))
        call dsavg(vzp(1,jp))
      enddo
      jp = 0

      return
      end subroutine init_pertfld_3ds
!----------------------------------------------------------------------
      subroutine makefp_3ds

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'   ! ifield

      include '3DS'


      integer nxyz,ntot1

      real ta1,ta2,ta3
      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)


      nxyz = lx1*ly1*lz1
      ntot1 = nxyz*nelv

      if (.not.if3d_3ds) then
        return
      endif

!     Build user defined forcing
      call makeufp_3ds

!      if3d = .true.
!      if (filterType.eq.2) call make_hpf
!     hpf field stored in ta3
!      call xaddcol3(bfz_3ds,ta3,bm1,ntot1)
!      if3d = .false. 

      call advabp_3ds
!      if (ifnav.and.(.not.ifchar).and.(ifadj)) call advabp_adjoint_3dsp

      if (iftran) call makextp_3ds
      call makebdfp_3ds


      return
      end subroutine makefp_3ds

!----------------------------------------------------------------------

      subroutine makeufp_3ds

!     Compute and add: (1) user specified forcing function (FX,FY,FZ)

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      include 'PARALLEL'
      include 'NEKUSE'

      include '3DS'

      real ta1,ta2,ta3
      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)

      integer ntot1,iel,i,j,k,ielg
      integer ijke


      ntot1 = lx1*ly1*lz1*nelv

      time = time-dt
      call rzero(bfxp(1,jp),ntot1)
      call rzero(bfyp(1,jp),ntot1)
      call rzero(bfzp(1,jp),ntot1)

      do 100 iel=1,nelv
         ielg = lglel(iel)
         do 100 k=1,lz1
         do 100 j=1,ly1
         do 100 i=1,lx1
            call nekasgn (i,j,k,iel)
            call userf   (i,j,k,ielg)
            ijke = i+lx1*((j-1)+ly1*((k-1) + lz1*(iel-1)))
            bfxp(ijke,jp) = ffx
            bfyp(ijke,jp) = ffy
            bfzp(ijke,jp) = ffz
 100  continue

!     Not sure why we multiply by density
!      call col2(bfzp(1,jp),vtrans(1,1,1,1,ifield),nx1*ny1*nz1*nelv)


      call col2  (bfxp(1,jp),bm1,ntot1)
      call col2  (bfyp(1,jp),bm1,ntot1)
      call col2  (bfzp(1,jp),bm1,ntot1)
      time = time+dt

      return
      end subroutine makeufp_3ds

!-----------------------------------------------------------------------
      subroutine advabp_3ds

!     Eulerian scheme, add convection term to forcing function
!     at current time step.

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'

      include '3DS'

      real ta1,ta2,ta3
      real tb1,tb2,tb3
      common /scrns/ ta1 (lx1*ly1*lz1*lelv)
     $ ,             ta2 (lx1*ly1*lz1*lelv)
     $ ,             ta3 (lx1*ly1*lz1*lelv)
     $ ,             tb1 (lx1*ly1*lz1*lelv)
     $ ,             tb2 (lx1*ly1*lz1*lelv)
     $ ,             tb3 (lx1*ly1*lz1*lelv)

      integer i,ntot1
      real tmp


      ntot1 = lx1*ly1*lz1*nelv

      if (if3d.or.if3d_3ds) then
         call copy  (tb1,vx,ntot1)                   ! Save velocity
         call copy  (tb2,vy,ntot1)                   ! Save velocity
         call copy  (tb3,vz,ntot1)                   ! Save velocity

!        U <-- dU
         call copy  (vx,vxp(1,jp),ntot1)                   ! Save velocity
         call copy  (vy,vxp(1,jp),ntot1)                   ! Save velocity
         call copy  (vz,vxp(1,jp),ntot1)                   ! Save velocity

         call convop  (ta1,tb1)                                ! du.grad U
         call convop  (ta2,tb2)
         call convop  (ta3,tb3)

!        Restore velocity
         call copy  (vx,tb1,ntot1)
         call copy  (vy,tb2,ntot1)
         call copy  (vz,tb3,ntot1)

         do i=1,ntot1
            tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
            bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
            bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
            bfzp(i,jp) = bfzp(i,jp)-tmp*ta3(i)
         enddo

         call convop  (ta1,vxp(1,jp))       !  U.grad dU
         call convop  (ta2,vyp(1,jp))
         call convop  (ta3,vzp(1,jp))

         do i=1,ntot1
            tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
            bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
            bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
            bfzp(i,jp) = bfzp(i,jp)-tmp*ta3(i)
         enddo

      else  ! 2D without fourier third component

         call opcopy  (tb1,tb2,tb3,vx,vy,vz)                   ! Save velocity
         call opcopy  (vx,vy,vz,vxp(1,jp),vyp(1,jp),vzp(1,jp)) ! U <-- dU
         call convop  (ta1,tb1)                                ! du.grad U
         call convop  (ta2,tb2)
         call opcopy  (vx,vy,vz,tb1,tb2,tb3)  ! Restore velocity

         do i=1,ntot1
            tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
            bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
            bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
         enddo

         call convop  (ta1,vxp(1,jp))       !  U.grad dU
         call convop  (ta2,vyp(1,jp))

         do i=1,ntot1
            tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
            bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
            bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
         enddo

      endif

!     Add z convection for all components


      return
      end subroutine advabp_3ds
c--------------------------------------------------------------------


      subroutine makextp_3ds

!     Add extrapolation terms to perturbation source terms

!     (nek5 equivalent for velocity is "makeabf")

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'

      include '3DS'

      common /scrns/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)

      ntot1 = lx1*ly1*lz1*nelv

      ab0 = ab(1)
      ab1 = ab(2)
      ab2 = ab(3)
      call add3s2 (ta1,exx1p(1,jp),exx2p(1,jp),ab1,ab2,ntot1)
      call add3s2 (ta2,exy1p(1,jp),exy2p(1,jp),ab1,ab2,ntot1)
      call copy   (exx2p(1,jp),exx1p(1,jp),ntot1)
      call copy   (exy2p(1,jp),exy1p(1,jp),ntot1)
      call copy   (exx1p(1,jp),bfxp (1,jp),ntot1)
      call copy   (exy1p(1,jp),bfyp (1,jp),ntot1)
      call add2s1 (bfxp(1,jp),ta1,ab0,ntot1)
      call add2s1 (bfyp(1,jp),ta2,ab0,ntot1)
      if (if3d.or.if3d_3ds) then
         call add3s2 (ta3,exz1p(1,jp),exz2p(1,jp),ab1,ab2,ntot1)
         call copy   (exz2p(1,jp),exz1p(1,jp),ntot1)
         call copy   (exz1p(1,jp),bfzp (1,jp),ntot1)
         call add2s1 (bfzp(1,jp),ta3,ab0,ntot1)
      endif
c
      return
      end subroutine makextp_3ds
!-----------------------------------------------------------------------

      subroutine makebdfp_3ds

!     Add contributions to perturbation source from lagged BD terms.

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'

      include '3DS'

      real ta1,ta2,ta3,tb1,tb2,tb3,h2
      common /scrns/ ta1(lx1,ly1,lz1,lelv)
     $ ,             ta2(lx1,ly1,lz1,lelv)
     $ ,             ta3(lx1,ly1,lz1,lelv)
     $ ,             tb1(lx1,ly1,lz1,lelv)
     $ ,             tb2(lx1,ly1,lz1,lelv)
     $ ,             tb3(lx1,ly1,lz1,lelv)
     $ ,             h2 (lx1,ly1,lz1,lelv)


      integer ilag,ntot1
      real const



      ntot1 = lx1*ly1*lz1*nelv
      const = 1./dt
      call cmult2(h2,vtrans(1,1,1,1,ifield),const,ntot1)
      call opcolv3c (tb1,tb2,tb3
     $              ,vxp(1,jp),vyp(1,jp),vzp(1,jp),bm1,bd(2))

      if (if3d_3ds) then
        call col3(tb3,vzp(1,jp),bm1,ntot1)
        call cmult(tb3,bd(2),ntot1)
      endif

!     Add contribution from lag terms
      do ilag=2,nbd
         if (ifgeom) then
            call opcolv3c(ta1,ta2,ta3,vxlagp(1,ilag-1,jp),
     $                                vylagp(1,ilag-1,jp),
     $                                vzlagp(1,ilag-1,jp),
     $                                bm1lag(1,1,1,1,ilag-1),bd(ilag+1))
         else
            call opcolv3c(ta1,ta2,ta3,vxlagp(1,ilag-1,jp),
     $                                vylagp(1,ilag-1,jp),
     $                                vzlagp(1,ilag-1,jp),
     $                                bm1                   ,bd(ilag+1))
         endif
         call opadd2  (tb1,tb2,tb3,ta1,ta2,ta3)
      enddo
      call opadd2col (bfxp(1,jp),bfyp(1,jp),bfzp(1,jp),tb1,tb2,tb3,h2)

!     Add contribution of lag terms to vzp
      if (if3d_3ds) then
        do ilag=2,nbd
           if (ifgeom) then
              call col3(ta3,vzlagp(1,ilag-1,jp),bm1lag(1,1,1,1,ilag-1),
     $                                ntot1) 
              call cmult(ta3,bd(ilag+1),ntot1)           
           else
              call col3(ta3,vzlagp(1,ilag-1,jp),bm1,ntot1) 
              call cmult(ta3,bd(ilag+1),ntot1)           
           endif
           call add2(tb3,ta3,ntot1)
        enddo
        call add2col2(bfzp(1,jp),tb3,h2,ntot1)
      endif       ! if3d_3ds


      return
      end subroutine makebdfp_3ds
!-----------------------------------------------------------------------

      subroutine lagfieldp_3ds

!     Keep old velocity field(s)

      implicit none 

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'

      include '3DS'

      integer ilag,ntot1

      ntot1 = lx1*ly1*lz1*nelv

      do ilag=nbdinp-1,2,-1
         call opcopy
     $     (vxlagp(1,ilag  ,jp),vylagp(1,ilag  ,jp),vzlagp(1,ilag  ,jp)
     $     ,vxlagp(1,ilag-1,jp),vylagp(1,ilag-1,jp),vzlagp(1,ilag-1,jp))
      enddo
      call opcopy(vxlagp(1,1,jp),vylagp(1,1,jp),vzlagp(1,1,jp)
     $           ,vxp   (1,jp)  ,vyp   (1,jp)  ,vzp   (1,jp) )

!     if3d_3ds
      if (if3d_3ds) then
        do ilag=nbdinp-1,2,-1
          call copy(vzlagp(1,ilag,jp),vzlagp(1,ilag-1,jp),ntot1)
        enddo
        call copy(vzlagp(1,1,jp),vzp(1,jp),ntot1)
      endif


      return
      end subroutine lagfieldp_3ds
!----------------------------------------------------------------------
      subroutine ophx_3ds (out1,out2,out3,inp1,inp2,inp3,h1,h2)

!     OUT = (H1*A+H2*B) * INP  

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'

      include '3DS'


      real out1 (lx1,ly1,lz1,1)
      real out2 (lx1,ly1,lz1,1)
      real out3 (lx1,ly1,lz1,1)
      real inp1 (lx1,ly1,lz1,1)
      real inp2 (lx1,ly1,lz1,1)
      real inp3 (lx1,ly1,lz1,1)
      real h1   (lx1,ly1,lz1,1)
      real h2   (lx1,ly1,lz1,1)

      integer imesh,matmod
      

      imesh = 1

      if (ifstrs) then
         matmod = 0
         call axhmsf (out1,out2,out3,inp1,inp2,inp3,h1,h2,matmod)
      else

!        the numbers are only needed for axis-symmetric formulation
!        need to come back to this later. 
         call axhelm (out1,inp1,h1,h2,imesh,1)
         call axhelm (out2,inp2,h1,h2,imesh,2)
         call axhelm (out3,inp3,h1,h2,imesh,3)
      endif

      return
      end subroutine ophx_3ds
!-----------------------------------------------------------------------
      subroutine cresvipp_3ds (resv1,resv2,resv3,h1,h2)

!     Compute startresidual/right-hand-side in the velocity solver

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'    ! v?mask
      include 'MASS'    ! bm1

      include '3DS'

      real           resv1 (lx1,ly1,lz1,1)
      real           resv2 (lx1,ly1,lz1,1)
      real           resv3 (lx1,ly1,lz1,1)
      real           h1    (lx1,ly1,lz1,1)
      real           h2    (lx1,ly1,lz1,1)

      real w1,w2,w3,prextr
      common /scruz/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)

      integer igeom
      common /cgeom/ igeom

      integer ntot1,ntot2
      real const
      

      ntot1 = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*nelv

      if3d = .true.
      call bcdirvc (vxp(1,jp),vyp(1,jp),vzp(1,jp),
     $              v1mask,v2mask,v3mask)
      if3d = .false.

!     prabal. We don't care about traction conditions for now.
!     Maybe need to look at it if added stiffness terms are needed
!     Or if surface tension is needed
!      call bcneutr

      if (mod(jp,2).eq.1) then

!       Extrapolate both pressures at the same time    
        call extrapprp (prextr_3ds(1,1))
        call extrapprp (prextr_3ds(1,2))

        call opgradt (resv1,resv2,resv3,prextr_3ds(1,1))

        call mappr(resv3,prextr_3ds(1,2),w1,w2)            ! w1,w2 used as work variables
        const = k_3dsp
        call cmult(resv3,const,ntot1)
        call col2(resv3,bm1,ntot1)
      else
        call opgradt (resv1,resv2,resv3,prextr_3ds(1,2))

        call mappr(resv3,prextr_3ds(1,1),w1,w2)            ! w1,w2 used as work variables
        const = -k_3dsp
        call cmult(resv3,const,ntot1)
        call col2(resv3,bm1,ntot1)
      endif    

      call add2(resv1,bfxp(1,jp),ntot1)
      call add2(resv2,bfyp(1,jp),ntot1)
      call add2(resv3,bfzp(1,jp),ntot1)

!     prabal
      call ophx_3ds(w1,w2,w3,vxp(1,jp),vyp(1,jp),
     $              vzp(1,jp),h1,h2)
      call sub2(resv1,w1,ntot1)
      call sub2(resv2,w2,ntot1)
      call sub2(resv3,w3,ntot1)

      return
      end subroutine cresvipp_3ds
!-----------------------------------------------------------------------

      subroutine sethlm_3dsp(h1,h2,intype)

      implicit none

c     Set the variable property arrays H1 and H2
c     in the Helmholtz equation.
c     (associated with variable IFIELD)
c     INTYPE =      integration type

      include 'SIZE'
      include 'INPUT'
      include 'SOLN' 
      include 'TSTEP'   ! ifield

      include '3DS'

      real h1(1),h2(1)

      integer intype
      integer nel,ntot1

      real dtbd

      real k2

      nel   = nelfld(ifield)
      ntot1 = lx1*ly1*lz1*nel

      k2    = k_3dsp**2

      if (iftran) then
         dtbd = bd(1)/dt
         call copy  (h1,vdiff (1,1,1,1,ifield),ntot1)
         if (intype.eq.0) then
            call rzero (h2,ntot1)
         else
            call cmult2(h2,vtrans(1,1,1,1,ifield),dtbd,ntot1)

!           Add second derivative of the 3rd direction to the operator
            call add2s2(h2,vdiff(1,1,1,1,ifield),k2,ntot1)
         endif
      else
         call copy  (h1,vdiff (1,1,1,1,ifield),ntot1)
         call rzero (h2,ntot1)

!        Add second derivative of the 3rd direction to the operator
         call add2s2(h2,vdiff(1,1,1,1,ifield),k2,ntot1)
      endif


      return
      end subroutine sethlm_3dsp

!----------------------------------------------------------------------

      subroutine incomprp_3ds (ux,uy,uz,up)
c
c     Project U onto the closest incompressible field
c
c     Input:  U     := (ux,uy,uz)
c
c     Output: updated values of U, iproj, proj; and
c             up    := pressure correction req'd to impose div U = 0
c
c
c     Dependencies: ifield ==> which "density" (vtrans) is used.
c
c
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

      include '3DS'

c
      common /scrns/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)
     $ ,             dv1   (lx1,ly1,lz1,lelv)
     $ ,             dv2   (lx1,ly1,lz1,lelv)
     $ ,             dv3   (lx1,ly1,lz1,lelv)
     $ ,             dp    (lx2,ly2,lz2,lelv)
      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)
      common /scrhi/ h2inv (lx1,ly1,lz1,lelv)
      COMMON /SCRCH/ PREXTR(LX2,LY2,LZ2,LELV)
      logical ifprjp

c
      if (icalld.eq.0) tpres=0.0
      icalld=icalld+1
      npres=icalld
c
      ntot1  = lx1*ly1*lz1*nelv
      ntot2  = lx2*ly2*lz2*nelv
      intype = 1
      dtbd   = bd(1)/dt

      call rzero   (h1,ntot1)
!      call cadd    (h1,k_3dsp,ntot1)
c     call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
      call cmult2  (h2,vtrans(1,1,1,1,ifield),dtbd,ntot1)
      call invers2 (h2inv,h2,ntot1)

      call opdiv   (dp,ux,uy,uz)
      call chsign  (dp,ntot2)
      call ortho   (dp)


C******************************************************************


      ifprjp=.false.    ! project out previous pressure solutions?
      istart=param(95)  
      if (istep.ge.istart.and.istart.ne.0) ifprjp=.true.

      ! Most likely, the following can be commented out. (pff, 1/6/2010)
c     if (npert.gt.1.or.ifbase)            ifprjp=.false.
cpff  if (ifprjp)   call setrhs  (dp,h1,h2,h2inv)

                    call esolver (dp,h1,h2,h2inv,intype)

cpff  if (ifprjp)   call gensoln (dp,h1,h2,h2inv)

cNOTE:  The "cpff" comments added 11/24/17 to avoid old-style projection,
cNOTE:  which should be replaced with something more updated.

C******************************************************************

      call opgradt (w1 ,w2 ,w3 ,dp)
      call opbinv  (dv1,dv2,dv3,w1 ,w2 ,w3 ,h2inv)
      call opadd2  (ux ,uy ,uz ,dv1,dv2,dv3)

      call extrapprp(prextr)
      call lagpresp
      call add3(up,prextr,dp,ntot2)

      return
      end subroutine incomprp_3ds
c------------------------------------------------------------------------









