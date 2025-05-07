C     *** CCFULL-SC *** 
C     Program for coupled-channels equations for heavy-ion 
C     elastic/inelastic scattering
C
C     The couplings are always treated to all orders.
C     Anharmonic couplings are included in this version of CCFULL. 
C
C     The imaginary part is introduced also in the coupling Hamiltonian
C     (but in the volume nuclear part only). 
C
C     The iso-centrifugal (no-Coriolis) approximation is employed in order to 
C     reduce the dimension of the CC equations. Note that the isocentrifugal 
C     approximation works only at backward angles. 
C
C     The nuclear S-matrix for L larger than LMAX0 is computed using the 
C     interpolation method, after calculating S-matrix from LMAX0 to LMAX
C     with a step of LSTEP. 
C
C     This program computes the angular distribution of the elastic 
C     and the quasi-elastic (i.e., for inclusive process) scattering.
C     Absorption (fusion) cross section is also computed. 
c
c     In order to avoid a numerical instability, the code solves the c.c. 
c     equations from RMIN to RMAX. For the same reason, the off-diagonal 
c     components of the coupling potential are assumed to be zero for 
c     r < RCOUPCUT. These are justified under a strong absorption in 
c     the inner region. 
C
C     last modification: November 27th, 2012
C
      implicit real*8 (a-h,o-z)
      character*1 ans
      complex*16 ai,facn,facn0,fn,fc,fn0,ref0
c
C  Parameters:
C    NLEVELMAX - maximum number of intrinsic levels to be included
C
      parameter (nlevelmax=30,lmax0=300,lstep=100,nrmax=10000)
c
C Principal global variables
C
C    SIGMAF    - fusion cross section, unit mb
C
C    AP,ZP,RP  - atomic #, proton # and radius of the projectile
C    AT,ZT,RT  - those of the target
C    RMASS     - reduced mass, unit MeV/c**2
C    ANMASS    - nucleon mass, unit MeV/c**2
C    E         - bombarding energy in the center of mass frame, unit MeV
C    V0,R0,A0  - depth, range, and dissuseness parameters of uncoupled 
C                nuclear potential, which is assumed to be a Woods-Saxon form
C    IVIBROT   - option for intrinsic degree of freedom 
C                = -1; no excitation (inert)/ =0 ; vibrational coupling/ 
C                =1; rotational coupling
C    BETA      - defeormation parameter
C    OMEGA     - excitation energy of the oscillator
C    LAMBDA    - multipolarity
C    NPHONON   - the number of phonon to be included
C    BETA2     - static quadrupole deformation parameter
C    BETA4     - static hexadecapole deformation parameter
C    E2        - excitation energy of 2+ state in a rotational band
C    NROT      - the number of levels in the rotational band to be included 
C                (up to I^pi=2*NROT+ states are included)
C    L         - angular momentum of the relative motion
C    CPOT      - coupling matrix
C
C    FACN(I,L) = (2L+1)*EXP(i(sigma_L+sigma_L')-delta_{I,0})
C                         *S_L / 2i /sqrt(kk') 
C    FN0,FC = Nuclear/Coulomb scattering amplitudes for the elastic channel
C    SIGMA_L = Coulomb phase shift
C    

        dimension cpot0(nlevelmax,nlevelmax)
        dimension cpotw0(nlevelmax,nlevelmax)
        dimension facn(nlevelmax,0:3000),facn0(nlevelmax)
        dimension sigma(nlevelmax),fn(nlevelmax),ref0(nlevelmax)

	common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
	common/intr/ivibrotp,ivibrott,nproj,ntarg
        common/oscp/betap,betapn,omegap,lambdap,nphononp
	common/rotp/beta2p,beta4p,e2p,nrotp
        common/osct/betat,betatn,omegat,lambdat,nphonont
        common/osct2/betat2,betat2n,omegat2,lambdat2,nphonont2
	common/rott/beta2t,beta4t,e2t,nrott
	common/dimen/nlevel
        common/pot/v0,r0,a0,pow
        common/imagp/w0,rw,aw
        common/imagps/w0s,rws,aws
	common/dyn/e,rmass
        common/angular/l
	common/shape/rb,vb,curv
	common/grid/rmin,rmax,dr
        common/coup/cpot(nlevelmax,nlevelmax,nrmax)
        common/coupw/cpotw(nlevelmax,nlevelmax,nrmax)
        common/trans/ftr,qtrans,ntrans
	common/mut/imutual(0:nlevelmax,0:nlevelmax)
        common/eigen/ev(nlevelmax),evec(nlevelmax,nlevelmax)
        common/coulomb/rc,sigma0
        common/stability/rcoupcut
        common/legendrepol/plgndr(0:3000)
C
C Output files: 
C   OUTPUT: summary of the calculation
C   angular.dat: angular distribution for the elastic, inelastic, 
C                and quasi-elastic scattering
C   angular2.dat: contribution of each channel to angular distribution 
C                 for the inelastic scattering
C   fusion.dat: fusion and (angular averaged) total inelastic cross sections
C
	open(7,file='angular.dat',status='unknown')
	open(77,file='angular.dat2',status='unknown')
	open(8,file='fusion.dat',status='unknown')
	open(9,file='OUTPUT',status='unknown')
	open(12,file='OUTPUT2',status='unknown')
	open(13,file='OUTPUT3',status='unknown')
	open(10,file='ccfull-sc.inp',status='unknown')
	open(11,file='qel.dat',status='unknown')
C
C Define three constants used in various places
C
      hbar=197.329d0
      pi=3.141592653d0
      ai=(0.d0,1.d0)
C
C Input phase
C
      anmass=938.d0

      ap=16.
      zp=8.
      at=148.
      zt=62.
	read(10,*)ap,zp,at,zt
c  r_coup, unit fm
      rp=1.2d0*ap**(1.d0/3.d0)
      rt=1.06d0*at**(1.d0/3.d0)
	read(10,*)rp,ivibrotp,rt,ivibrott
	rp=rp*ap**(1.d0/3.d0)
	rt=rt*at**(1.d0/3.d0)

      rmass=ap*at/(ap+at)*anmass

c-----
c      ivibrott=0
c             = -1 for inert, =0 for vibration, =1 for rotation
	if(ivibrott.eq.-1) then
	      ntarg=0
              read(10,*)omega,beta,alambda,nphonon
	elseif(ivibrott.eq.0) then
c (input for target phonon)
              omegat=0.55d0
              betat=0.182d0
              lambdat=2
              nphonont=3
              read(10,*)omegat,betat,lambdat,nphonont
	      ntarg=nphonont
	else
c (input for target rotation)
	      e2t=0.d0
	      beta2t=0.d0
	      beta4t=0.d0
              nrott=0
              lambdat=2
              read(10,*)e2t,beta2t,beta4t,nrott
	      ntarg=nrott
	endif
c-----
c (input for target phonon (the second mode of excitation))
              nphonont2=3
c                     =0; no second mode in the target
              omegat2=1.16d0
              betat2=0.236d0
              lambdat2=3
              read(10,*)omegat2,betat2,lambdat2,nphonont2
	      if(ivibrott.eq.-1) nphonont2=0
c-----
c      ivibrotp=-1
c             = -1 for inert, =0 for vibration, =1 for rotation
	if(ivibrotp.eq.-1) then
	      nproj=0
              read(10,*)omegap,betap,alambdap,nphononp
	else if(ivibrotp.eq.0) then
c (input for projectile phonon)
              omegap=3.74d0
              betap=0.339d0
              lambdap=3
              nphononp=1
              read(10,*)omegap,betap,lambdap,nphononp
              nproj=nphononp
	else
c (input for projectile rotation)
              e2p=0.d0
	      beta2p=0.d0
	      beta4p=0.d0
	      nrotp=0
              read(10,*)e2p,beta2p,beta4p,nrotp
	      nproj=nrotp
	endif

      ntrans=0
c ntrans= 0 ; no transfer channel/ =1 ; with transfer channel
c ftrans(r)=ftr*dvn(r)/dr
c
      ftr=7.d0
      qtrans=5.21d0
        	read(10,*)ntrans,qtrans,ftr

c--------------- parameters for the nuclear potential
              r1=1.233d0*ap**(1.d0/3.d0)-0.98d0*ap**(-1.d0/3.d0)
              r2=1.233d0*at**(1.d0/3.d0)-0.98d0*at**(-1.d0/3.d0)
              r0=r1+r2+0.29d0
              r12=r1*r2/(r1+r2)
              a0=0.63d0
              gamma=0.95d0*(1.d0-1.8d0*(ap-2.d0*zp)/ap*(at-2.d0*zt)/at)
              v0=16.d0*pi*gamma*r12*a0
              v0=31.67*r1*r2/(r1+r2)

                v0=1551d0
                r0=0.95d0*(ap**(1.d0/3.d0)+at**(1.d0/3.d0))
                a0=1.05d0
                        read(10,*)v0,r0,a0
	                r0=r0*(ap**(1.d0/3.d0)+at**(1.d0/3.d0))
                        read(10,*)w0,rw,aw
	                rw=rw*(ap**(1.d0/3.d0)+at**(1.d0/3.d0))
                        read(10,*)w0s,rws,aws
	                rws=rws*(ap**(1.d0/3.d0)+at**(1.d0/3.d0))
                        read(10,*)rc
c	                rc=rc*(ap**(1.d0/3.d0)+at**(1.d0/3.d0))
        pow=1.d0
	write(6,*)' Modified Woods-Saxon (power of WS) (n/y)?'
	read(5,1)ans
 1	format(a1)
	     IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
	     write(6,*)'power=?'
	     read(5,*)pow
                if(pow.eq.0.d0) then
                write(6,*)'Power cannot be zero :^('
                stop
                endif
	     endif

C  the CM energy and the maximum L
      read(10,*)emin,emax,de
      read(10,*)lmax
      read(10,*)thmin,thmax,dtheta

      rmax=50.d0
      dr=0.05d0
        	read(10,*)rmax,dr
        	read(10,*)rmin,rcoupcut

C========================================== end of input phase

	if(ivibrott.ne.-1.and.ntarg.eq.0.and.nphonont2.ne.0) then
	write(6,91)
 91	format(1h ,'The target excitation should be input'
     &  /' in the first mode :^(')
	stop
	endif

        if(lmax.gt.3000) then
           write(6,*)'Too large LMAX!!! :^('
           stop
           endif

c	if(ivibrott.eq.1.and.nphonont2.ne.0) then
c 	write(6,*)'!!!!!!!!!!!!!!! Warning !!!!!!!!!!!!!!!!!!'
c	write(6,*)'    '
c	write(6,*)'The vibrational phonon in the target is not'
c        write(6,*)'a deformed phonon, but a *spherical* phonon,'
c        write(6,*)'i.e. the rotation-vibration model is not used.'
c	write(6,*)'    '
c 	write(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
c	write(6,*)'    '
c	endif

         call nucleus

C Evaluation coupling matrix elements
C------------------------------------------------------
      l=0
      call potshape
      ngrid=(rmax-rmin)/dr+1.d-6
             if(ngrid+1.gt.nrmax) then
	     write(6,*)'Too large grid !!  :^('
	     stop
	     endif
      rmax=rmin+dr*ngrid

      call cmat0
      do ir=1,ngrid+1
      r=rmin+dr*ir
      call cmat(r,cpot0,cpotw0)
      do i=1,nlevel
      do j=1,nlevel
         cpot(i,j,ir)=cpot0(i,j)
         cpotw(i,j,ir)=cpotw0(i,j)
         enddo
         enddo
         enddo

            write(7,95)
 95         format(1h '#   E (MeV)       theta_cm','     sig_el/sig_R', 
     &             '  sig_inel(mb)  sig_inel/sig_R  sig_qel/sig_R')
            write(77,96)
 96         format(1h '#     E (MeV)       theta_cm','  channel', 
     &             '  sig (mb)     sig/sig_R')

            write(12,*)'       '
            write(12,*)'Angular distribution for elastic,inel, and qel:'
            write(12,*)'       '
            write(12,80)
 80         format('    Ecm (MeV)    theta_cm (deg)  sig_el/sig_R'  
     &       ,'   sig_inel (mb)  sig_inel/sig_R  sig_qel/sig_R')
            write(12,81)
 81         format('-----------------------------------------------'
     &            ,'-----------------------------------------------')

            write(13,*)'       '
            write(13,*)'    Ecm (MeV)    Eeff (MeV)  sig_qel/sig_R'
            write(13,*)'--------------------------------------------'

c=======================================
c  start the energy loop
c-------------------------------------------------------------
	if(de.eq.0.d0) then
	iestep=0
	else
        iestep=(emax-emin)/de
	endif

        if(thmin.eq.0.) then
           thmin=max(dtheta,0.01d0)
           write(9,*)'Thmin = 0!!! ---> set Thmin=',thmin
           write(12,*)'Thmin = 0!!! ---> set Thmin=',thmin
           write(13,*)'Thmin = 0!!! ---> set Thmin=',thmin
           write(6,*)'Thmin = 0!!! ---> set Thmin=',thmin
           endif

	if(dtheta.eq.0.d0) then
	ithstep=0
	else
        ithstep=(thmax-thmin)/dtheta+1.d-6
	endif

      do 10 i=0,iestep
      e=emin+de*i

         ak=sqrt(2.d0*rmass/hbar**2*e)
         eta=(zp*zt/137.d0)*sqrt(rmass/2.d0/e)

         sigmaf=0.d0
         do ich=1,nlevelmax
         sigma(ich)=0.d0
         enddo

c  angular momentum loop
c-----------------------
      do 20 l=0,min(lmax,lmax0)

c integration of the c.c. equations
      call mnumerov(facn0,ref0,p)
      sigmaf=sigmaf+(2.d0*l+1.d0)*p*pi*hbar**2/2.d0/rmass/e*10.d0

         facn(1,l)=facn0(1)
         do ich=2,nlevel
         facn(ich,l)=facn0(ich)
         sigma(ich)=sigma(ich)+(2.d0*l+1.d0)*abs(ref0(ich))**2
     &                              *pi*hbar**2/2.d0/rmass/e*10.d0
         enddo
      if(nlevel.gt.8) write(6,*)e,l,p
c      write(20,*)e,l,p
 20   continue

        if(lmax.gt.lmax0) 
     &    call interpolation(lmax,lmax0,lstep,sigmaf,sigma,facn,ref0)

      write(9,*)'             '
      write(9,*)'********** E=',e,' (MeV) **********'
      write(9,*)'Fusion cross section =',sigmaf,' (mb)'
        do ich=2,nlevel
        write(9,33)ich,sigma(ich)
        enddo
 33     format(' Total inel sc cross section to the channel'
     &                                      ,i2,' = ',g15.5,'(mb)')
           sigmatot=0.d0
           do ich=2,nlevel
           sigmatot=sigmatot+sigma(ich)
           enddo
        if(nlevel.gt.1) write(9,*)
     &                    '  Total inclusive inel sc cross section='
     &                                  ,sigmatot,'(mb)'

      write(8,*)e,sigmaf,sigmatot
      write(6,*)e,sigmaf,sigmatot
      
C=======================================
C angular distributions

            write(9,*)'       '
            write(9,*)'Angular distribution for elastic,inel, and qel:'

      do ith=0,ithstep
         th=thmin+ith*dtheta
                  th=th/180.d0*pi

                  escale=2.d0*e*sin(th/2.d0)/(1.d0+sin(th/2.d0))
C                 Coulomb scattering amplitudes
         fc=-eta/2.d0/ak/sin(th/2.d0)**2
     &       *exp(-ai*eta*log(sin(th/2.d0)**2)+2.d0*ai*sigma0)

C       Nuclear scattering amplitudes
         call legendre(lmax,cos(th))

         do ich=1,nlevel
            fn(ich)=0.d0
               do l=0,lmax
c               fn(ich)=fn(ich)+facn(ich,l)*pn(l,cos(th))
               fn(ich)=fn(ich)+facn(ich,l)*plgndr(l)
               enddo
               enddo

            fn0=fn(1)
            sigma_el=abs(fn0+fc)**2/abs(fc)**2

                   write(77,94)e,th*180.d0/pi,'R',abs(fc)**2*10.d0,1.d0
                   write(77,93)e,th*180.d0/pi,1,
     &            abs(fn0+fc)**2*10.d0,sigma_el

                sigma_inel=0.d0
                do ich=2,nlevel
                sigma_inel=sigma_inel+abs(fn(ich))**2
                   write(77,93)e,th*180.d0/pi,ich,
     &            abs(fn(ich))**2*10.d0,abs(fn(ich))**2/abs(fc)**2

                enddo

            sigma_inel2=sigma_inel/abs(fc)**2
            sigma_inel=sigma_inel*10.d0
            sigma_qel=sigma_el+sigma_inel2

            write(7,92)e,th*180.d0/pi,sigma_el,sigma_inel,sigma_inel2,
     &                                                    sigma_qel
            write(9,92)e,th*180.d0/pi,sigma_el,sigma_inel,sigma_inel2,
     &                                                    sigma_qel
            write(12,92)e,th*180.d0/pi,sigma_el,sigma_inel,sigma_inel2,
     &                                                    sigma_qel
            write(11,92)e,escale,sigma_qel
            write(13,92)e,escale,sigma_qel
      enddo
 10   continue

 92   format(6g15.5)
 93   format(2f15.5,i5,2g15.5)
 94   format(2f15.5,a5,2g15.5)

 70   write(6,*)'The program terminated normally :^)'	

      stop
      end

c*************************************************************************
      subroutine interpolation(lmax,lmax0,lstep,sigmaf,sigma,facn,ref0)
C
C Interpolation of nuclear S matrix. 
C
C************************************************************************
      implicit real*8 (a-h,o-z)
      complex*16 ai,facn,facn0,ref0,ref

      parameter (nlevelmax=30)

      dimension facn(nlevelmax,0:3000),facn0(nlevelmax)
      dimension sigma(nlevelmax),ref0(nlevelmax)
      dimension ref(nlevelmax,0:3000)
      dimension yr(0:10),yi(0:10),dyr(0:10,0:10),dyi(0:10,0:10)
      dimension cr(nlevelmax,0:10),ci(nlevelmax,0:10)
      dimension fcw(0:3000),gcw(0:3000),fpcw(0:3000),gpcw(0:3000)
      dimension sigmad(0:3000),iexp(0:3000)

	common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
	common/dimen/nlevel
	common/dyn/e,rmass
        common/angular/l
        common/coup1/cpot1(nlevelmax)

      ai=(0.d0,1.d0)

      lstep0=min((lmax-lmax0)/5.d0+1.d-6,lstep*1.d0)
      lstep0=max(1,lstep0)

      n=(lmax-lmax0)/lstep0+1

      do 20 l0=0,n
      l=lmax0+l0*lstep0
      call mnumerov(facn0,ref0,p)
      do ich=1,nlevel
      ref(ich,l)=ref0(ich)
      enddo
 20   continue

c----------------------------------------- interpolation
      i0=0
      do 10 l=lmax0+1,lmax
         i00=i0
           i0=max(5,(l-lmax0)/lstep0)
                   if(l.gt.lmax0+lstep0*i0) i0=i0+1
           if(lmax0+i0*lstep0.gt.3000) i0=i0-1

c----
        if(i00.ne.i0) then
        lmin=lmax0+(i0-5)*lstep0
        do ich=1,nlevel
           yr(0)=ref(ich,lmin)
           yr(1)=ref(ich,lmin+lstep0)
           yr(2)=ref(ich,lmin+2*lstep0)
           yr(3)=ref(ich,lmin+3*lstep0)
           yr(4)=ref(ich,lmin+4*lstep0)
           yr(5)=ref(ich,lmin+5*lstep0)

           yi(0)=(ref(ich,lmin)-yr(0))/ai
           yi(1)=(ref(ich,lmin+lstep0)-yr(1))/ai
           yi(2)=(ref(ich,lmin+2*lstep0)-yr(2))/ai
           yi(3)=(ref(ich,lmin+3*lstep0)-yr(3))/ai
           yi(4)=(ref(ich,lmin+4*lstep0)-yr(4))/ai
           yi(5)=(ref(ich,lmin+5*lstep0)-yr(5))/ai

        nstep=5
        do j=0,nstep-1
           dyr(1,j)=yr(j+1)-yr(j)
           dyi(1,j)=yi(j+1)-yi(j)
           enddo

        do k=2,nstep
        do j=0,nstep-k
           dyr(k,j)=dyr(k-1,j+1)-dyr(k-1,j)
           dyi(k,j)=dyi(k-1,j+1)-dyi(k-1,j)
           enddo
           enddo

        i=1
        cr(ich,0)=yr(0)
        ci(ich,0)=yi(0)
        do k=1,nstep
           i=i*k
           cr(ich,k)=dyr(k,0)/(i*1.d0)/(lstep0*1.d0)**k
           ci(ich,k)=dyi(k,0)/(i*1.d0)/(lstep0*1.d0)**k
           enddo
        enddo
        endif

c----
        p=0.d0
        do ich=1,nlevel
           xr=1.d0
           xi=1.d0
           ref0(ich)=cr(ich,0)+ai*ci(ich,0)
           do k=1,nstep
              xr=xr*(l-(lmin+(k-1.d0)*lstep0))
              xi=xi*(l-(lmin+(k-1.d0)*lstep0))
              ref0(ich)=ref0(ich)+cr(ich,k)*xr+ai*ci(ich,k)*xi
              enddo

        sigma(ich)=sigma(ich)+(2.d0*l+1.d0)*abs(ref0(ich))**2
     &                              *pi*hbar**2/2.d0/rmass/e*10.d0
        p=p+abs(ref0(ich))**2
        enddo
        p=1.d0-p
            if(p.lt.0.d0) p=0.d0

        sigmaf=sigmaf+(2.d0*l+1.d0)*p*pi*hbar**2/2.d0/rmass/e*10.d0

        do ich=1,nlevel
         ech=e-cpot1(ich)
         ak=sqrt((2.d0*rmass/hbar**2*ech))
         rho=ak*1.d0
         eta=(zp*zt/137.d0)*sqrt(rmass/2.d0/ech)
         call dfcoul(eta,rho,fcw,fpcw,gcw,gpcw,sigmad,l,iexp)
           if(ich.eq.1) then
              sigma00=sigmad(l)
              ak0=ak
              facn(1,l)=(2.d0*l+1.d0)*exp(2.d0*ai*sigma00)
     &                                   *(ref0(1)-1.d0)/2.d0/ai/ak0
           else
              facn(ich,l)=(2.d0*l+1.d0)*exp(ai*(sigma00+sigmad(l)))
     &                       *ref0(ich)/2.d0/ai/ak0
           endif
        enddo

 10   continue

      return
      end

C*************************************************************
      subroutine nucleus
C
C Subroutine to record several input parameters in the 
C 'OUTPUT' file
C
C*************************************************************
      implicit real*8(a-h,o-z)
      parameter (nlevelmax=30)
      character*1 ans
      character*2 text,pro,targ
      dimension text(109)
	common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
	common/intr/ivibrotp,ivibrott,nproj,ntarg
        common/oscp/betap,betapn,omegap,lambdap,nphononp
	common/rotp/beta2p,beta4p,e2p,nrotp
        common/osct/betat,betatn,omegat,lambdat,nphonont
        common/osct2/betat2,betat2n,omegat2,lambdat2,nphonont2
	common/rott/beta2t,beta4t,e2t,nrott
	common/rotpn/beta2pn,beta4pn
	common/rottn/beta2tn,beta4tn
	common/phonon/nphonon
	common/dimen/nlevel
        common/pot/v0,r0,a0,pow
        common/imagp/w0,rw,aw
        common/imagps/w0s,rws,aws
	common/dyn/e,rmass
	common/shape/rb,vb,curv
	common/grid/rmin,rmax,dr
        common/trans/ftr,qtrans,ntrans
	common/mut/imutual(0:nlevelmax,0:nlevelmax)
        common/ahv/betnahv(0:nlevelmax,0:nlevelmax),
     &             betcahv(0:nlevelmax,0:nlevelmax),
     &             omeahv(0:nlevelmax)
        common/ahv2/betnahv2(0:nlevelmax,0:nlevelmax),
     &             betcahv2(0:nlevelmax,0:nlevelmax),
     &             omeahv2(0:nlevelmax)
        common/ahvp/betnahvp(0:nlevelmax,0:nlevelmax),
     &             betcahvp(0:nlevelmax,0:nlevelmax),
     &             omeahvp(0:nlevelmax)
        common/coup1/cpot1(nlevelmax)

      data text/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     &'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ',
     &'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
     &'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',
     &'Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm',
     &'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re',
     &'Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra',
     &'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md',
     &'No','Lr','XX','X1','X2','X3','X4','04'/
 
        nlevel=(nproj+1)*(ntarg+1)
	np=zp
	npp=ap
	pro=text(np)
	nt=zt
	ntt=at
	targ=text(nt)
	write(9,*)npp,pro,'  +  ',ntt,targ, '     Reaction'
	write(12,*)npp,pro,'  +  ',ntt,targ, '     Reaction'
	write(13,*)npp,pro,'  +  ',ntt,targ, '     Reaction'
	write(6,*)npp,pro,'  +  ',ntt,targ, '     Reaction'

	write(9,30)
	write(12,30)
	write(13,30)
	write(6,30)

	r0p=rp/ap**(1.d0/3.d0)
	r0t=rt/at**(1.d0/3.d0)

	if(ntarg.ne.0) then
	if(ivibrott.eq.0) then
	write(6,21)betat,omegat,lambdat,nphonont 
 21	format('Phonon Excitation in the targ.: beta=',F6.3,
     &    ', omega=',F5.2, '(MeV), Lambda=',I2,', Nph=',I2)

	betatn=betat
	write(6,*)'    '
	write(6,*)' Different beta_N from beta_C for this mode(n/y)?'
	read(5,1)ans
	     if(ans.eq.'Y' .or. ans.eq.'y') then
	     write(6,*)'beta_N=?'
	     read(5,*)betatn
	     endif
	write(9,211)betatn,betat,r0t,omegat,lambdat,nphonont 
	write(12,211)betatn,betat,r0t,omegat,lambdat,nphonont 
	write(13,211)betatn,betat,r0t,omegat,lambdat,nphonont 
 211	format('Phonon Excitation in the targ.: beta_N=',F6.3,
     &    ', beta_C=',F6.3,', r0=',F5.2,'(fm),',
     &     /32X,'omega=',F5.2, '(MeV), Lambda=',I2,', Nph=',I2)
	endif

	if(ivibrott.eq.1) then

	write(9,22)beta2t,beta4t,r0t,e2t,nrott
	write(12,22)beta2t,beta4t,r0t,e2t,nrott
	write(13,22)beta2t,beta4t,r0t,e2t,nrott
	write(6,22)beta2t,beta4t,r0t,e2t,nrott
 22     format('Rotational Excitation in the targ.: beta2=',F6.3,
     &    ', beta4=',F6.3,', r0=',F5.2,'(fm),',
     &      /36X,'E2=',F5.2,'(MeV), Nrot=',I2)

	beta2tn=beta2t
	beta4tn=beta4t
	write(6,*)'    '
	write(6,*)' Different beta_N from beta_C (n/y)?'
	read(5,1)ans
	     if(ans.eq.'Y' .or. ans.eq.'y') then
	     write(6,*)'beta2_N,beta4_N=?'
	     read(5,*)beta2tn,beta4tn
	     endif

	write(9,222)beta2tn,beta4tn
	write(12,222)beta2tn,beta4tn
	write(13,222)beta2tn,beta4tn
	write(6,222)beta2tn,beta4tn
 222     format('beta2_N=',F6.3,', beta4_N=',F6.3)

	endif
	endif

	if(nphonont2.ne.0) then
	write(6,23)betat2,omegat2,lambdat2,nphonont2 
 23	format('Phonon Excitation in the targ.: beta=',F6.3,
     &    ', omega=',F5.2, '(MeV), Lambda=',I2,', Nph=',I2)

	betat2n=betat2
	write(6,*)'    '
	write(6,*)' Different beta_N from beta_C for this mode(n/y)?'
	read(5,1)ans
	     if(ans.eq.'Y' .or. ans.eq.'y') then
	     write(6,*)'beta_N=?'
	     read(5,*)betat2n
	     endif

	write(9,233)betat2n,betat2,r0t,omegat2,lambdat2,nphonont2 
	write(12,233)betat2n,betat2,r0t,omegat2,lambdat2,nphonont2 
	write(13,233)betat2n,betat2,r0t,omegat2,lambdat2,nphonont2 
 233	format('Phonon Excitation in the targ.: beta_N=',F6.3,
     &    ', beta_C=',F6.3,', r0=',F5.2,'(fm),',
     &    /32X, 'omega=',F5.2, '(MeV), Lambda=',I2,', Nph=',I2)

	call mutual
	endif

	if(nproj.ne.0) then
	if(ivibrotp.eq.0) then
	write(6,24)betap,omegap,lambdap,nphononp 
 24	format('Phonon Excitation in the proj.: beta=',F6.3,
     &    ', omega=',F5.2, '(MeV), Lambda=',I2,', Nph=',I2)

	betapn=betap
	write(6,*)'    '
	write(6,*)' Different beta_N from beta_C for this mode(n/y)?'
	read(5,1)ans
	     if(ans.eq.'Y' .or. ans.eq.'y') then
	     write(6,*)'beta_N=?'
	     read(5,*)betapn
	     endif
	write(9,244)betapn,betap,r0p,omegap,lambdap,nphononp 
	write(12,244)betapn,betap,r0p,omegap,lambdap,nphononp 
	write(13,244)betapn,betap,r0p,omegap,lambdap,nphononp 
 244	format('Phonon Excitation in the proj.: beta_N=',F6.3,
     &    ', beta_C=',F6.3,', r0=',F5.2,'(fm),',
     &    /32X, 'omega=',F5.2, '(MeV), Lambda=',I2,', Nph=',I2)
	endif

	if(ivibrotp.eq.1) then

	write(9,25)beta2p,beta4p,r0p,e2p,nrotp
	write(12,25)beta2p,beta4p,r0p,e2p,nrotp
	write(13,25)beta2p,beta4p,r0p,e2p,nrotp
	write(6,25)beta2p,beta4p,r0p,e2p,nrotp
 25     format('Rotational Excitation in the proj.: beta2=',F6.3,
     &    ', beta4=',F6.3,', r0=',F5.2,'(fm),',
     &      /36X,'E2=',F5.2,'(MeV), Nrot=',I2)

	beta2pn=beta2p
	beta4pn=beta4p
	write(6,*)'    '
	write(6,*)' Different beta_N from beta_C (n/y)?'
	read(5,1)ans
	     if(ans.eq.'Y' .or. ans.eq.'y') then
	     write(6,*)'beta2_N,beta4_N=?'
	     read(5,*)beta2pn,beta4pn
	     endif

	write(9,222)beta2pn,beta4pn
	write(12,222)beta2pn,beta4pn
	write(13,222)beta2pn,beta4pn
	write(6,222)beta2pn,beta4pn

	endif
	endif

        call anharmonicity

	if(ntrans.ne.0) then
	write(9,30)
	write(9,26)ftr, qtrans
	write(12,30)
	write(12,26)ftr, qtrans
	write(13,30)
	write(13,26)ftr, qtrans
	write(6,26)ftr, qtrans
 26	format('Transfer channel: Strength=',F5.2,
     &    ', Q=',F5.2, '(MeV)')
	endif

	write(9,30)
	write(12,30)
	write(13,30)

	r00=r0/(ap**(1.d0/3.d0)+at**(1.d0/3.d0))
	write(9,27)v0,r00,a0,pow
	write(12,27)v0,r00,a0,pow
	write(13,27)v0,r00,a0,pow
 27	format('Potential parameters: V0=', F8.2, '(MeV), r0=', 
     &        F5.2, '(fm), a=', F5.2, '(fm), power=',F5.2)
		call potshape
		write(9,28)rb,vb,curv
		write(12,28)rb,vb,curv
		write(13,28)rb,vb,curv
 28 	format('   Uncoupled barrier: Rb=', F5.2, '(fm), Vb=',F8.2, 
     &        '(MeV), Curv=',F5.2, '(MeV)')
	write(9,30)
	write(12,30)
	write(13,30)

	r00=rw/(ap**(1.d0/3.d0)+at**(1.d0/3.d0))
	write(9,29)w0,r00,aw
	write(12,29)w0,r00,aw
	write(13,29)w0,r00,aw
 29	format('Potential parameters (volume imaginary): W0=', F8.2
     &        , '(MeV), rw=', 
     &        F5.2, '(fm), aw=', F5.2, '(fm)')
	r00=rws/(ap**(1.d0/3.d0)+at**(1.d0/3.d0))
	write(9,291)w0s,r00,aws
	write(12,291)w0s,r00,aws
	write(13,291)w0s,r00,aws
 291	format('Potential parameters (surface imaginary): W0s=', F8.2
     &        , '(MeV), rws=', 
     &        F5.2, '(fm), aws=', F5.2, '(fm)')
	write(9,30)
	write(12,30)
	write(13,30)

        write(9,*)'Channels:'
        write(12,*)'Channels:'
        write(13,*)'Channels:'
        nlevel=nlevel+ntrans
	 if(nlevel.gt.nlevelmax) then
	 write(6,*)'Too many channels :^('
	 stop
	 endif
         write(6,*)'         '
         write(6,*) 'nlevel=',nlevel

        call cmat1

	   ich=0
           do 31 ip=0,nproj
           do 32 it=0,ntarg
           do 33 it2=0,nphonont2
			if(nphonont2.ne.0) then
			if(imutual(it,it2).eq.0) goto 33
			endif
	   ich=ich+1

           write(9,34)ich,it,it2,ip,cpot1(ich)
           write(12,34)ich,it,it2,ip,cpot1(ich)
           write(13,34)ich,it,it2,ip,cpot1(ich)
           write(6,34)ich,it,it2,ip,cpot1(ich)
 34        format(i4,' Ch: IT=',I3,', IT2=',I3,', IP=',I3,',    E*='
     &                                                ,f8.3,' MeV')

 33        continue
 32        continue
 31        continue
           if(ntrans.ne.0)write(9,36)nlevel,cpot1(nlevel)
           if(ntrans.ne.0)write(12,36)nlevel,cpot1(nlevel)
           if(ntrans.ne.0)write(13,36)nlevel,cpot1(nlevel)
           if(ntrans.ne.0)write(6,36)nlevel,cpot1(nlevel)
 36        format(i4,' Ch: Transfer    -Q=',f12.3,'  MeV')
         write(6,*)'         '

	write(9,30)
	write(12,30)
	write(13,30)
 30	format('-------------------------------------------------')

 1	format(a1)
      return
      end

c*************************************************************
      subroutine anharmonicity
C
C Anharmonic couplings for vibrations
C
C*************************************************************
      implicit real*8(a-h,o-z)
      parameter (nlevelmax=30)
      character*1 ans,sign
	common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
	common/intr/ivibrotp,ivibrott,nproj,ntarg
        common/oscp/betap,betapn,omegap,lambdap,nphononp
	common/rotp/beta2p,beta4p,e2p,nrotp
        common/osct/betat,betatn,omegat,lambdat,nphonont
        common/osct2/betat2,betat2n,omegat2,lambdat2,nphonont2
	common/rott/beta2t,beta4t,e2t,nrott
	common/phonon/nphonon
	common/dimen/nlevel
        common/pot/v0,r0,a0,pow
	common/dyn/e,rmass
	common/shape/rb,vb,curv
	common/grid/rmin,rmax,dr
        common/trans/ftr,qtrans,ntrans
	common/mut/imutual(0:nlevelmax,0:nlevelmax)
        common/ahv/betnahv(0:nlevelmax,0:nlevelmax),
     &             betcahv(0:nlevelmax,0:nlevelmax),
     &             omeahv(0:nlevelmax)
        common/ahv2/betnahv2(0:nlevelmax,0:nlevelmax),
     &             betcahv2(0:nlevelmax,0:nlevelmax),
     &             omeahv2(0:nlevelmax)
        common/ahvp/betnahvp(0:nlevelmax,0:nlevelmax),
     &             betcahvp(0:nlevelmax,0:nlevelmax),
     &             omeahvp(0:nlevelmax)

C-----------------------------
c AHV coupling for the first mode in the target
           do it=0,ntarg
           do jt=0,ntarg

           betnahv(it,jt)=0.d0
           betcahv(it,jt)=0.d0

           if(it.eq.jt-1) then
              betnahv(it,jt)=betatn*sqrt(jt*1.d0)
              betcahv(it,jt)=betat*sqrt(jt*1.d0)
           else if(jt.eq.it-1) then
              betnahv(it,jt)=betatn*sqrt(it*1.d0)
              betcahv(it,jt)=betat*sqrt(it*1.d0)
           endif

           enddo
              omeahv(it)=it*omegat
           enddo

        if(ivibrott.eq.0.and.ntarg.gt.1) then

	write(6,90)'AHV couplings for the first mode ',
     &             'in the target phonon (n/y)?'
 90	format(2a)
	read(5,1)ans
        if(ans.eq.'Y' .or. ans.eq.'y') then
           write(9,*)'            '
           write(9,90)'**** AHV couplings in the target',
     &                ' (the 1st mode) ****'
           write(9,*)'             '
           write(12,*)'            '
           write(12,90)'**** AHV couplings in the target',
     &                ' (the 1st mode) ****'
           write(12,*)'             '
           write(13,*)'            '
           write(13,90)'**** AHV couplings in the target',
     &                ' (the 1st mode) ****'
           write(13,*)'             '
	if((-1.d0)**lambdat.eq.1.d0) sign='+'
	if((-1.d0)**lambdat.eq.-1.d0) sign='-'

           do it=0,ntarg
           do 44 jt=0,ntarg
           if(it.gt.jt) goto 44
           
           if(it.eq.0.and.jt.eq.0) goto 44
           if(it.eq.0.and.jt.eq.1) goto 44

        write(6,*)'             '
	write(6,94)'Transition from the',lambdat,sign,'^',it,
     &                  ' to the',lambdat,sign,'^',jt,'state:'
        write(6,95)'   beta_N and beta_C in the HO limit=',
     &             betnahv(it,jt),betcahv(it,jt)
        write(6,90)'    Would you like to change these ',
     &                 'beta_N a/o beta_C (n/y)?'
	read(5,1)ans

        if(ans.eq.'Y' .or. ans.eq.'y') then
	write(9,94)'Transition from the',lambdat,sign,'^',it,
     &                  ' to the',lambdat,sign,'^',jt,'state:'
        write(9,95)'   beta_N and beta_C in the HO limit=',
     &             betnahv(it,jt),betcahv(it,jt)
	write(12,94)'Transition from the',lambdat,sign,'^',it,
     &                  ' to the',lambdat,sign,'^',jt,'state:'
        write(12,95)'   beta_N and beta_C in the HO limit=',
     &             betnahv(it,jt),betcahv(it,jt)
	write(13,94)'Transition from the',lambdat,sign,'^',it,
     &                  ' to the',lambdat,sign,'^',jt,'state:'
        write(13,95)'   beta_N and beta_C in the HO limit=',
     &             betnahv(it,jt),betcahv(it,jt)
        write(6,*)'   beta_N and beta_C =?'
        read(5,*)betnahv(it,jt),betcahv(it,jt)
        write(9,97)betnahv(it,jt),betcahv(it,jt)
        write(12,97)betnahv(it,jt),betcahv(it,jt)
        write(13,97)betnahv(it,jt),betcahv(it,jt)
        endif

 94     format(1h ,a,i2,2a,i2,a,i2,2a,i2,a)
 95     format(1h ,a,2g15.5)
 97     format('    beta_N=',G15.5,',   beta_C=',G15.5)

 44        enddo
           enddo

              write(6,*)'           '
              write(6,*)'Excitation energy for the first mode:'
              do i=2,ntarg
              write(6,98)'    Energy of the',i,
     &              '-phonon state in the HO=',i*omegat
              write(6,*)'   Would you like to change this energy(n/y)?'
              read(5,1)ans
              if(ans.eq.'Y' .or. ans.eq.'y') then
              write(6,*)'   Energy=?'
              read(5,*)omeahv(i)
              write(6,*)'       '
              write(9,61)'Energy of the',i,'-phonon state=',omeahv(i),
     &                         '(HO: ',i*omegat,')'
              write(12,61)'Energy of the',i,'-phonon state=',omeahv(i),
     &                         '(HO: ',i*omegat,')'
              write(13,61)'Energy of the',i,'-phonon state=',omeahv(i),
     &                         '(HO: ',i*omegat,')'
              endif
              enddo
              write(9,*)'          '              
              write(12,*)'          '              
              write(13,*)'          '              
              write(6,*)'          '              
 98           format(a,i2,a,g15.5)
 61           format(1h ,a,i2,a,g15.5,a,g12.5,a)
            endif
            endif
c-----------------------------
c AHV coupling for the second mode in the target
           do it=0,nphonont2
           do jt=0,nphonont2

           betnahv2(it,jt)=0.d0
           betcahv2(it,jt)=0.d0

           if(it.eq.jt-1) then
              betnahv2(it,jt)=betat2n*sqrt(jt*1.d0)
              betcahv2(it,jt)=betat2*sqrt(jt*1.d0)
           else if(jt.eq.it-1) then
              betnahv2(it,jt)=betat2n*sqrt(it*1.d0)
              betcahv2(it,jt)=betat2*sqrt(it*1.d0)
           endif

           enddo
              omeahv2(it)=it*omegat2
           enddo

        if(nphonont2.gt.1) then

	write(6,90)'AHV couplings for the second mode ',
     &             'in the target phonon (n/y)?'
	read(5,1)ans
        if(ans.eq.'Y' .or. ans.eq.'y') then
           write(9,*)'            '
           write(9,90)'**** AHV Couplings in the target',
     &                ' (the 2nd mode) ****'
           write(9,*)'             '
           write(12,*)'            '
           write(12,90)'**** AHV Couplings in the target',
     &                ' (the 2nd mode) ****'
           write(12,*)'             '
           write(13,*)'            '
           write(13,90)'**** AHV Couplings in the target',
     &                ' (the 2nd mode) ****'
           write(13,*)'             '

	if((-1.d0)**lambdat2.eq.1.d0) sign='+'
	if((-1.d0)**lambdat2.eq.-1.d0) sign='-'

           do it=0,nphonont2
           do 54 jt=0,nphonont2
           if(it.gt.jt) goto 54
           
           if(it.eq.0.and.jt.eq.0) goto 54
           if(it.eq.0.and.jt.eq.1) goto 54

        write(6,*)'             '
	write(6,94)'Transition from the',lambdat2,sign,'^',it,
     &                  ' to the',lambdat2,sign,'^',jt,'state:'
        write(6,95)'    beta_N and beta_C in the HO limit=',
     &             betnahv2(it,jt),betcahv2(it,jt)
        write(6,90)'     Would you like to change these ',
     &                 'beta_N a/o beta_C (n/y)?'
	read(5,1)ans

        if(ans.eq.'Y' .or. ans.eq.'y') then
	write(9,94)'Transition from the',lambdat2,sign,'^',it,
     &                  ' to the',lambdat2,sign,'^',jt,'state:'
        write(9,95)'   beta_N and beta_C in the HO limit=',
     &             betnahv(it,jt),betcahv(it,jt)
	write(12,94)'Transition from the',lambdat2,sign,'^',it,
     &                  ' to the',lambdat2,sign,'^',jt,'state:'
        write(12,95)'   beta_N and beta_C in the HO limit=',
     &             betnahv(it,jt),betcahv(it,jt)
	write(13,94)'Transition from the',lambdat2,sign,'^',it,
     &                  ' to the',lambdat2,sign,'^',jt,'state:'
        write(13,95)'   beta_N and beta_C in the HO limit=',
     &             betnahv(it,jt),betcahv(it,jt)
        write(6,*)'    beta_N and beta_C =?'
        read(5,*)betnahv2(it,jt),betcahv2(it,jt)
        write(9,97)betnahv2(it,jt),betcahv2(it,jt)
        write(12,97)betnahv2(it,jt),betcahv2(it,jt)
        write(13,97)betnahv2(it,jt),betcahv2(it,jt)
        endif

 54        enddo
           enddo

              write(6,*)'           '
              write(6,*)'Excitation energy for the second mode:'
              do i=2,nphonont2
              write(6,98)'    Energy of the',i,
     &              '-phonon state in the HO=',i*omegat2
            write(6,*)'   Would you like to change this energy(n/y)?'
              read(5,1)ans
              if(ans.eq.'Y' .or. ans.eq.'y') then
              write(6,*)'   Energy=?'
              read(5,*)omeahv2(i)
              write(9,61)'Energy of the',i,'-phonon state=',omeahv2(i),
     &                         '(HO: ',i*omegat2,')'
              write(12,61)'Energy of the',i,'-phonon state=',omeahv2(i),
     &                         '(HO: ',i*omegat2,')'
              write(13,61)'Energy of the',i,'-phonon state=',omeahv2(i),
     &                         '(HO: ',i*omegat2,')'
              endif
              enddo
              write(9,*)'          '              
              write(12,*)'          '              
              write(13,*)'          '              
              write(6,*)'          '              
            endif
            endif

C-----------------------------
c AHV coupling for the projectile phonon 
           do it=0,nproj
           do jt=0,nproj

           betnahvp(it,jt)=0.d0
           betcahvp(it,jt)=0.d0

           if(it.eq.jt-1) then
              betnahvp(it,jt)=betapn*sqrt(jt*1.d0)
              betcahvp(it,jt)=betap*sqrt(jt*1.d0)
           else if(jt.eq.it-1) then
              betnahvp(it,jt)=betapn*sqrt(it*1.d0)
              betcahvp(it,jt)=betap*sqrt(it*1.d0)
           endif

           enddo
              omeahvp(it)=it*omegap
           enddo

        if(ivibrotp.eq.0.and.nproj.gt.1) then

	write(6,72)
 72     format('AHV couplings for the projectile phonon (n/y)?')
	read(5,1)ans
        if(ans.eq.'Y' .or. ans.eq.'y') then
           write(9,*)'            '
           write(9,73)
           write(12,*)'            '
           write(12,73)
           write(13,*)'            '
           write(13,73)
 73        format(1h ,'**** AHV couplings in the projectile ****') 
           write(9,*)'            '
           write(12,*)'            '
           write(13,*)'            '
	if((-1.d0)**lambdap.eq.1.d0) sign='+'
	if((-1.d0)**lambdap.eq.-1.d0) sign='-'

           do it=0,nproj
           do 64 jt=0,nproj
           if(it.gt.jt) goto 64
           
           if(it.eq.0.and.jt.eq.0) goto 64
           if(it.eq.0.and.jt.eq.1) goto 64

        write(6,*)'             '
	write(6,94)'Transition from the',lambdap,sign,'^',it,
     &                  ' to the',lambdap,sign,'^',jt,'state:'
        write(6,95)'  beta_N and beta_C in the HO limit=',
     &             betnahvp(it,jt),betcahvp(it,jt)
        write(6,90)'   Would you like to change these ',
     &                 'beta_N a/o beta_C (n/y)?'
	read(5,1)ans

        if(ans.eq.'Y' .or. ans.eq.'y') then
	write(9,94)'Transition from the',lambdap,sign,'^',it,
     &                  ' to the',lambdap,sign,'^',jt,'state:'
        write(9,95)'   beta_N and beta_C in the HO limit=',
     &             betnahvp(it,jt),betcahvp(it,jt)
	write(12,94)'Transition from the',lambdap,sign,'^',it,
     &                  ' to the',lambdap,sign,'^',jt,'state:'
        write(12,95)'   beta_N and beta_C in the HO limit=',
     &             betnahvp(it,jt),betcahvp(it,jt)
	write(13,94)'Transition from the',lambdap,sign,'^',it,
     &                  ' to the',lambdap,sign,'^',jt,'state:'
        write(13,95)'   beta_N and beta_C in the HO limit=',
     &             betnahvp(it,jt),betcahvp(it,jt)
        write(6,*)'  beta_N and beta_C =?'
        read(5,*)betnahvp(it,jt),betcahvp(it,jt)
        write(9,97)betnahvp(it,jt),betcahvp(it,jt)
        write(12,97)betnahvp(it,jt),betcahvp(it,jt)
        write(13,97)betnahvp(it,jt),betcahvp(it,jt)
        endif

 64        enddo
           enddo

              write(6,*)'           '
           write(6,*)'Excitation energy for the projectile phonon:'
           do i=2,nproj
              write(6,98)'    Energy of the',i,
     &              '-phonon state in the HO=',i*omegap
            write(6,*)'   Would you like to change this energy(n/y)?'
              read(5,1)ans
              if(ans.eq.'Y' .or. ans.eq.'y') then
              write(6,*)'   Energy=?'
              read(5,*)omeahvp(i)
              write(9,61)'Energy of the',i,'-phonon state=',omeahvp(i),
     &                         '(HO: ',i*omegap,')'
              write(12,61)'Energy of the',i,'-phonon state=',omeahvp(i),
     &                         '(HO: ',i*omegap,')'
              write(13,61)'Energy of the',i,'-phonon state=',omeahvp(i),
     &                         '(HO: ',i*omegap,')'
              endif
              enddo
              write(9,*)'          '              
              write(12,*)'          '              
              write(13,*)'          '              
              write(6,*)'          '              
            endif
            endif

 1	format(a1)

      return
      end

C*****************************************************************************
      subroutine mutual
C
C Subroutine to select levels to be included in the c.c. calculations
C
C*****************************************************************************
        implicit real*8(a-h,m,o-z)
        parameter (nlevelmax=30)
	character*1 ans,sign1,sign2
	common/intr/ivibrotp,ivibrott,nproj,ntarg
        common/osct/betat,betatn,omegat,lambdat,nphonont
        common/osct2/betat2,betat2n,omegat2,lambdat2,nphonont2
	common/rott/beta2t,beta4t,e2t,nrott
	common/dimen/nlevel
        common/trans/ftr,qtrans,ntrans
	common/mut/imutual(0:nlevelmax,0:nlevelmax)

	do 10 i=0,nlevelmax
	do 20 j=0,nlevelmax
	imutual(i,j)=0
 20	continue
 10	continue

	write(6,*)'    '
	write(6,*)'Mutual excitations in the *target* nucleus'

	write(6,*)'  Include the mutual excitations (y/n)?'
	read(5,1)ans
 1	format(a1)
	     if(ans.eq.'N' .or. ans.eq.'n') then
	     imut=0
	     nlevel=(ntarg+nphonont2+1)*(nproj+1)
c	     write(9,99)
c 99	     FORMAT(1H ,'No mutual excitations in the target are included.')
	     do 80 i=0,ntarg
	     do 81 j=0,nphonont2
	     if(i.eq.0.or.j.eq.0) imutual(i,j)=1
 81          continue
 80	     continue
	     goto 70
	     endif

	write(6,*)'  All the possible mutual excitation channels (n/y)?'
	read(5,1)ans
	     if(ans.eq.'Y' .or. ans.eq.'y') then
	     imut=1
	     nlevel=(ntarg+1)*(nphonont2+1)*(nproj+1)
c	     write(9,90)
c 90	     FORMAT(1H ,'All the possible mutual excitation channels'
c     &                  1H ,'in the target are included.')
	     do 82 i=0,ntarg
	     do 83 j=0,nphonont2
	     imutual(i,j)=1
 83          continue
 82	     continue
	     goto 70
	     endif

	imut=2
	if((-1.d0)**lambdat.eq.1.d0) sign1='+'
	if((-1.d0)**lambdat.eq.-1.d0) sign1='-'
	if((-1.d0)**lambdat2.eq.1.d0) sign2='+'
	if((-1.d0)**lambdat2.eq.-1.d0) sign2='-'
	do 30 i=0,ntarg
	do 40 j=0,nphonont2
		if(i.eq.0.or.j.eq.0) then 
		imutual(i,j)=1
		goto 40
		endif
	write(6,94)'  Include (',lambdat,sign1,'^',i,',',lambdat2,sign2,
     &            '^',j,') state ? (y/n)'
	read(5,1)ans
	     if(ans.eq.'N' .or. ans.eq.'n') then
	     imutual(i,j)=0
	     goto 40
	     endif
	imutual(i,j)=1
 40	continue
 30	continue

 70	write(9,*)'    '
	write(9,93)
 	write(12,*)'    '
	write(12,93)
 	write(13,*)'    '
	write(13,93)
	write(6,93)
 93	format(1h ,'  Excited states in the target to be included: ')

	nlevel=0
	if((-1.d0)**lambdat.eq.1.d0) sign1='+'
	if((-1.d0)**lambdat.eq.-1.d0) sign1='-'
	if((-1.d0)**lambdat2.eq.1.d0) sign2='+'
	if((-1.d0)**lambdat2.eq.-1.d0) sign2='-'
	do 50 i=0,ntarg
	do 60 j=0,nphonont2
		if(imutual(i,j).eq.0) goto 60
	write(9,94)'    (',lambdat,sign1,'^',i,',',lambdat2,sign2,
     &                   '^',j,') state'
	write(12,94)'    (',lambdat,sign1,'^',i,',',lambdat2,sign2,
     &                   '^',j,') state'
	write(13,94)'    (',lambdat,sign1,'^',i,',',lambdat2,sign2,
     &                   '^',j,') state'
	write(6,94)'    (',lambdat,sign1,'^',i,',',lambdat2,sign2,
     &                   '^',j,') state'
 94	format(1h ,a,i2,2a,i2,a,i2,2a,i2,a)
	nlevel=nlevel+1
 60	continue
 50	continue
        nlevel=nlevel*(nproj+1)

	return
	end
C*****************************************************************************
      subroutine potshape
C
C shape of the bare Coulomb barrier, i.e. the barrier position, the curvature,
C and the position of the Coulomb pocket
C
C*****************************************************************************
      implicit real*8(a-h,m,o-z)
      common/shape/rb,vb,curv
      common/grid/rmin,rmax,dr
      common/dyn/e,rmass
      common/const/hbar,pi
      external dv,v

        r=50.5d0
        u0=dv(r)
 10	r=r-1.d0

           if(r.lt.0.d0) then
              rb=-5.d0
	      rmin0=-5.d0
              return
              endif

	u1=dv(r)
	if(u0*u1.gt.0.d0) goto 10
		ra=r+1.d0
		rb=r
	        tolk=1.d-6
                n=log10(abs(rb-ra)/tolk)/log10(2.d0)+0.5d0
                do 20 i=1,n
                r=(ra+rb)/2.d0
                u=dv(r)
                if(u0*u.lt.0.d0) then
                rb=r
                else
                ra=r
                endif
 20             continue

	rb=r

        ddv00=(dv(rb+1.d-5)-dv(rb-1.d-5))/2.d-5

	if(ddv00.gt.0.d0) then
	write(6,*)'something is strange :^('
	write(9,*)'something is strange :^('
	write(12,*)'something is strange :^('
	write(13,*)'something is strange :^('
	stop
	endif

	curv=hbar*sqrt(abs(ddv00)/rmass)
	vb=v(rb)

		ra=rb-0.5d0
		u0=dv(ra)
		rbb=0.1d0
		u1=dv(rbb)
	        tolk=1.d-6
                n=log10(abs(rb-ra)/tolk)/log10(2.d0)+0.5d0
                do 21 i=1,n
                r=(ra+rbb)/2.d0
                u=dv(r)
                if(u0*u.lt.0.d0) then
                rbb=r
                else
                ra=r
                endif
 21             continue

	rmin0=r

c	write(6,*)'rmin=',rmin,v(rmin)
c       write(6,*)'rb=',rb,v(rb)
c       write(6,*)'curv=',curv00

        rmin0=0.d0

	return
	end

c**************************************************************
      subroutine mnumerov(facn,ref,p)
C
C Subroutine for integration of the c.c. eqs. by modified 
C Numerov method
C
C**************************************************************
      implicit real*8 (a-h,o-z)
      external v,w
      parameter (nlevelmax=30,nrmax=10000)

      dimension psi(nlevelmax),psi0(nlevelmax),psi1(nlevelmax)
      dimension xi(nlevelmax,nlevelmax),xi0(nlevelmax,nlevelmax),
     &          xi1(nlevelmax,nlevelmax)
      dimension phi0(nlevelmax)
      dimension bb(nlevelmax,nlevelmax),
     &          bin(nlevelmax,nlevelmax)
      dimension cc(nlevelmax,nlevelmax),
     &          cin(nlevelmax,nlevelmax)
      dimension bb0(nlevelmax,nlevelmax),bb20(nlevelmax,nlevelmax)
      dimension bb2(nlevelmax,nlevelmax),ref(nlevelmax)
      dimension dd0(nlevelmax,nlevelmax),dd1(nlevelmax,nlevelmax)
      dimension dd(nlevelmax)
      dimension fcw(0:3000),gcw(0:3000),fpcw(0:3000),gpcw(0:3000)
      dimension sigmad(0:3000),iexp(0:3000)
      dimension ech(nlevelmax)
      dimension facn(nlevelmax),sigmas(nlevelmax)
 
	common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
	common/dimen/nlevel
	common/dyn/e,rmass
        common/angular/l
	common/grid/rmin,rmax,dr
        common/coup/cpot(nlevelmax,nlevelmax,nrmax)
        common/coupw/cpotw(nlevelmax,nlevelmax,nrmax)
        common/coup1/cpot1(nlevelmax)
        common/coulomb/rc,sigma0
	common/shape/rb,vb,curv

      complex*16 k,kk
      complex*16 psi,psi0,psi1,phi0
      complex*16 xi,xi0,xi1
      complex*16 ai
      complex*16 cwup0,cwdown0,cwup1,cwdown1
      complex*16 bb,bin,cc,cin,dd,dd0,dd1
      complex*16 bb2,ref,bb20,bb0
      complex*16 facn
      complex*16 aa(nlevelmax,nlevelmax)

      call cmat1
      ai=(0.d0,1.d0)
c define two constants used in this subroutine
      fac=dr**2*(2.d0*rmass/hbar**2) 
      ngrid=(rmax-rmin)/dr+1.d-6
      rmax=rmin+ngrid*dr
      ibarrier=(rb-rmin)/dr+1.d-6
      istab=1.d0/dr+1.d-6

	do 22 lc=0,3000
	fcw(lc)=0.d0
	gcw(lc)=0.d0
	fpcw(lc)=0.d0
	gpcw(lc)=0.d0
	sigmad(lc)=0.d0
	iexp(lc)=0d0
 22	continue

      do 51 i=1,nlevelmax
      do 52 j=1,nlevelmax
      bb(i,j)=0.d0
      bin(i,j)=0.d0
      cc(i,j)=0.d0
      cin(i,j)=0.d0
      aa(i,j)=0.d0
 52   continue
      dd(i)=0.d0
 51   continue

      do i=1,nlevel
         aa(i,i)=1.d0
         enddo

c initial conditions
c==========================================================
      do 15 io=1,nlevel

      do 200 j1=1,nlevel
      psi(j1)=0.d0
      psi0(j1)=0.d0
      psi1(j1)=0.d0
      phi0(j1)=0.d0
 200  continue

      psi1(io)=0.01d0

      do 91 i0=1,nlevel
      xi0(i0,io)=0.d0
      xi1(i0,io)=(1.d0-fac/12.d0*(v(rmin+dr)-ai*w(rmin+dr)-e))*psi1(i0)
      do 92 ic=1,nlevel
      xi1(i0,io)=xi1(i0,io)-fac/12.d0*(cpot(i0,ic,1)-ai*cpotw(i0,ic,1))
     &                                                       *psi1(ic)
 92   continue
 91   continue
 15   continue

c===============
      do 29 ir=2,ngrid+1
c----------------------------------------------  ngridions start
      r=rmin+dr*ir
      r0=rmin+dr*(ir-2)
      r1=rmin+dr*(ir-1)

      do 71 i0=1,nlevel
      do 72 ic=1,nlevel
       dd0(i0,ic)=fac/sqrt(12.d0)*(cpot(i0,ic,ir-1)
     &                                    -ai*cpotw(i0,ic,ir-1))
      if(i0.eq.ic) dd0(i0,ic)=dd0(i0,ic)+fac/sqrt(12.d0)
     &                        *(v(r1)-ai*w(r1)-e)+sqrt(3.d0)
 72   continue
 71   continue

      do 73 i0=1,nlevel
      do 74 ic=1,nlevel
      dd1(i0,ic)=0.d0
      if(i0.eq.ic) dd1(i0,ic)=dd1(i0,ic)-1.d0
      do 75 ik=1,nlevel
         dd1(i0,ic)=dd1(i0,ic)+dd0(i0,ik)*dd0(ik,ic)
 75   continue
 74   continue
 73   continue
	
      do 16 ich=1,nlevel
      do 54 i0=1,nlevel
         xi(i0,ich)=-xi0(i0,ich)
      do 53 ic=1,nlevel
         xi(i0,ich)=xi(i0,ich)+dd1(i0,ic)*xi1(ic,ich)
 53   continue
 54   continue
 16   continue

      if(ir.eq.ngrid+1) goto 66
c      if(ir.eq.ibarrier) call stabilize(xi1,xi,aa,ir)
      if(mod(ir,istab).eq.0.and.r.lt.15.d0) 
     &                                 call stabilize(xi1,xi,aa,ir)

      do ich=1,nlevel

            dummy=0.d0
            do i0=1,nlevel
               if(abs(xi(i0,ich)).gt.dummy) dummy=abs(xi(i0,ich))
            enddo

      do i0=1,nlevel
         xi0(i0,ich)=xi1(i0,ich)/dummy
         xi1(i0,ich)=xi(i0,ich)/dummy
      enddo
      enddo

 29   continue
c--------------------------------------------------------------
c  matching to the coulomb wave function at rmax
 66   do 17 io=1,nlevel

      do 56 i0=1,nlevel
      do 55 ic=1,nlevel
      cc(i0,ic)=-fac/12.d0*(cpot(i0,ic,ngrid-1)
     &                               -ai*cpotw(i0,ic,ngrid-1))
      if(i0.eq.ic) cc(i0,ic)=cc(i0,ic)
     &          -fac/12.d0*(v(rmax-dr)-ai*w(rmax-dr)-e)+1.d0
 55   continue
 56   continue
      call matinv(nlevel,cc,cin)
      do 94 i0=1,nlevel
         psi0(i0)=0.d0
      do 93 ic=1,nlevel
         psi0(i0)=psi0(i0)+cin(i0,ic)*xi0(ic,io)
 93   continue
 94   continue

      do 556 i0=1,nlevel
      do 555 ic=1,nlevel
      cc(i0,ic)=-fac/12.d0*(cpot(i0,ic,ngrid+1)
     &                                      -ai*cpotw(i0,ic,ngrid+1))
      if(i0.eq.ic) cc(i0,ic)=cc(i0,ic)
     &            -fac/12.d0*(v(rmax+dr)-ai*w(rmax+dr)-e)+1.d0
 555   continue
 556   continue
      call matinv(nlevel,cc,cin)
      do 944 i0=1,nlevel
         psi(i0)=0.d0
      do 933 ic=1,nlevel
         psi(i0)=psi(i0)+cin(i0,ic)*xi(ic,io)
 933   continue
 944   continue

      do 60 ii=1,nlevel
c   coulomb wave function
      ec=e-cpot1(ii)
      rho=sqrt(2.d0*rmass*ec)/hbar*(rmax-dr)
      eta=(zp*zt/137.d0)*sqrt(rmass/2.d0/ec)
      call dfcoul(eta,rho,fcw,fpcw,gcw,gpcw,sigmad,l,iexp)
      cwup0=(gcw(l)+ai*fcw(l))
      cwdown0=gcw(l)-ai*fcw(l)
            if(l.eq.0.and.ii.eq.1) sigma0=sigmad(l)
            
               sigmas(ii)=sigmad(l)

      ec=e-cpot1(ii)
      rho=sqrt(2.d0*rmass*ec)/hbar*(rmax+dr)
      eta=(zp*zt/137.d0)*sqrt(rmass/2.d0/ec)
      call dfcoul(eta,rho,fcw,fpcw,gcw,gpcw,sigmad,l,iexp)
      cwup1=(gcw(l)+ai*fcw(l))
      cwdown1=gcw(l)-ai*fcw(l)

      bb0(ii,io)=(cwup0*psi(ii)-cwup1*psi0(ii))
     &                         /(cwup0*cwdown1-cwup1*cwdown0)
      bb20(ii,io)=-(cwdown0*psi(ii)-cwdown1*psi0(ii))
     &                         /(cwup0*cwdown1-cwup1*cwdown0)
 60   continue
c---------------------------------------------------------------
 17   continue

c===============================================================
c                                                       s-matrix
      do i=1,nlevel
      do j=1,nlevel
         bb(i,j)=0.d0
         bb2(i,j)=0.d0
         do j2=1,nlevel
            bb(i,j)=bb(i,j)+bb0(i,j2)*aa(j2,j)
            bb2(i,j)=bb2(i,j)+bb20(i,j2)*aa(j2,j)
            enddo
            enddo
            enddo

c  calculate the inversion of the matrix bb
      call matinv(nlevel,bb,bin)
c
      do i=1,nlevel
         ref(i)=0.d0
         do io=1,nlevel
         ref(i)=ref(i)+bb2(i,io)*bin(io,1)
         enddo
         enddo

         kk=sqrt(2.d0*rmass/hbar**2*e)
         do 85 io=1,nlevel
         ech(io)=e-cpot1(io)
         if(ech(io).lt.0.d0) ech(io)=0.d0
         k=sqrt((2.d0*rmass/hbar**2*ech(io)))
         ref(io)=-ref(io)*sqrt(k/kk)
 85   continue

c                   (partial fusion cross section)
      p=0.d0
      do io=1,nlevel
      p=p+(abs(ref(io)))**2
      enddo
      p=1.d0-p

c             if(p.lt.0d0) p=0.d0

      sigma00=sigmas(1)
      do ich=2,nlevel
         ak=sqrt((2.d0*rmass/hbar**2*ech(ich)))
         facn(ich)=(2.d0*l+1.d0)*exp(ai*(sigma00+sigmas(ich)))
     &                       *ref(ich)/2.d0/ai/sqrt(kk*ak)
         facn(ich)=facn(ich)*sqrt(ak/kk)
      enddo

      facn(1)=(2.d0*l+1.d0)*exp(2.d0*ai*sigma00)
     &                                   *(ref(1)-1.d0)/2.d0/ai/kk

      return
      end

c**************************************************************
      subroutine stabilize(xi1,xi,aa,ir1)
C
C A subroutine to stabilize the solutions by imposing the linear 
C independence. 
C
C Refs.
C   Z.H. Levine, PRA30('84) 1120
C   T.N. Rescigno and A.E. Orel, PRA25('82)2402
C
C**************************************************************
      implicit real*8 (a-h,o-z)
      parameter (nlevelmax=30,nrmax=10000)

      dimension psi(nlevelmax,nlevelmax)
      dimension xi(nlevelmax,nlevelmax),xi1(nlevelmax,nlevelmax)
      dimension cc(nlevelmax,nlevelmax),
     &          cin(nlevelmax,nlevelmax)

	common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
	common/dimen/nlevel
	common/dyn/e,rmass
        common/angular/l
	common/grid/rmin,rmax,dr
        common/coup/cpot(nlevelmax,nlevelmax,nrmax)
        common/coupw/cpotw(nlevelmax,nlevelmax,nrmax)

      complex*16 psi,xi,xi1
      complex*16 ai
      complex*16 cc,cin
      complex*16 aa(nlevelmax,nlevelmax)
      complex*16 aa0(nlevelmax,nlevelmax)
      complex*16 xid(nlevelmax,nlevelmax),xid1(nlevelmax,nlevelmax)

      do i=1,nlevel
      do j=1,nlevel
         aa0(i,j)=aa(i,j)
         xid(i,j)=xi(i,j)
         xid1(i,j)=xi1(i,j)
         enddo
         enddo

      ai=(0.d0,1.d0)
      fac=dr**2*(2.d0*rmass/hbar**2) 

c     transforme to the physical wave functions psi from xi 
      r=rmin+ir1*dr
         do i0=1,nlevel
         do ic=1,nlevel
          cc(i0,ic)=-fac/12.d0*(cpot(i0,ic,ir1)-ai*cpotw(i0,ic,ir1))
          if(i0.eq.ic) cc(i0,ic)=cc(i0,ic)
     &                        -fac/12.d0*(v(r)-ai*w(r)-e)+1.d0
          enddo
          enddo
        call matinv(nlevel,cc,cin)
      do ich=1,nlevel
        do i0=1,nlevel
          psi(i0,ich)=0.d0
          do ic=1,nlevel
           psi(i0,ich)=psi(i0,ich)+cin(i0,ic)*xi(ic,ich)
           enddo
           enddo
           enddo

         do i=1,nlevel
         do j=1,nlevel
              aa(i,j)=0.d0
              do k=1,nlevel
                 aa(i,j)=aa(i,j)+psi(i,k)*aa0(k,j)
                 enddo
                 enddo
                 enddo

        call matinv(nlevel,psi,cin)
        
        do i=1,nlevel
        do j=1,nlevel
           xi(i,j)=0.d0
           xi1(i,j)=0.d0
           do k=1,nlevel
              xi(i,j)=xi(i,j)+xid(i,k)*cin(k,j)
              xi1(i,j)=xi1(i,j)+xid1(i,k)*cin(k,j)
              enddo
              enddo
              enddo

      return
      end

c******************************************************************
	subroutine cmat0
C
C Preparation for the nuclear coupling matrix
C
C*****************************************************************
      implicit real*8 (a-h,o-z)
      parameter (nlevelmax=30)
      dimension a(nlevelmax,nlevelmax)
      common/hion/ap,zp,rp,at,zt,rt
      common/const/hbar,pi
      common/intr/ivibrotp,ivibrott,nproj,ntarg
      common/oscp/betap,betapn,omegap,lambdap,nphononp
      common/osct/betat,betatn,omegat,lambdat,nphonont
      common/osct2/betat2,betat2n,omegat2,lambdat2,nphonont2
      common/rotp/beta2p,beta4p,e2p,nrotp
      common/rott/beta2t,beta4t,e2t,nrott
      common/rotpn/beta2pn,beta4pn
      common/rottn/beta2tn,beta4tn
      common/dimen/nlevel
      common/mut/imutual(0:nlevelmax,0:nlevelmax)
      common/eigen/ev(nlevelmax),evec(nlevelmax,nlevelmax)
      common/trans/ftr,qtrans,ntrans
        common/ahv/betnahv(0:nlevelmax,0:nlevelmax),
     &             betcahv(0:nlevelmax,0:nlevelmax),
     &             omeahv(0:nlevelmax)
        common/ahv2/betnahv2(0:nlevelmax,0:nlevelmax),
     &             betcahv2(0:nlevelmax,0:nlevelmax),
     &             omeahv2(0:nlevelmax)
        common/ahvp/betnahvp(0:nlevelmax,0:nlevelmax),
     &             betcahvp(0:nlevelmax,0:nlevelmax),
     &             omeahvp(0:nlevelmax)
      external cg

        fvibt=rt*betatn/sqrt(4.d0*pi)
	frot2t=rt*beta2tn
	frot4t=rt*beta4tn
        fvibt2=rt*betat2n/sqrt(4.d0*pi)
        fvibp=rp*betapn/sqrt(4.d0*pi)
	frot2p=rp*beta2pn
	frot4p=rp*beta4pn

	do 10 i=1,nlevelmax
	do 20 j=1,nlevelmax
	a(i,j)=0.d0
 20	continue
 10	continue

	ndim=nlevel-ntrans
	if(ndim.eq.1) return
	
	   i=0
           do 31 ip=0,nproj
           do 32 it=0,ntarg
           do 33 it2=0,nphonont2
			if(nphonont2.ne.0) then
			if(imutual(it,it2).eq.0) goto 33
			endif
	   i=i+1
           j=0
           do 41 jp=0,nproj
           do 42 jt=0,ntarg
           do 43 jt2=0,nphonont2
			if(nphonont2.ne.0) then
			if(imutual(jt,jt2).eq.0) goto 43
			endif
	   j=j+1

              if(i.gt.j) then
                 a(i,j)=a(j,i)
                 goto 43
                 endif

	   c=0.d0
	   if(ip.eq.jp.and.it2.eq.jt2) then
             if(ivibrott.eq.0) then
c               if(it.eq.jt-1) c=sqrt(jt*1.d0)*fvibt
c               if(jt.eq.it-1) c=sqrt(it*1.d0)*fvibt

               c=rt*betnahv(it,jt)/sqrt(4.d0*pi)

	     else
		c=frot2t*sqrt((2.d0*2*it+1.d0)*5.d0*(2.d0*2*jt+1.d0)
     &    /4.d0/pi)*cg(2*it,0,2,0,2*jt,0)**2/(2.d0*2*jt+1.d0)
		c=c+frot4t*sqrt((2.d0*2*it+1.d0)*9.d0*(2.d0*2*jt+1.d0)
     &     /4.d0/pi)*cg(2*it,0,4,0,2*jt,0)**2/(2.d0*2*jt+1.d0)
            endif
	   endif
	   if(ip.eq.jp.and.it.eq.jt) then
c               if(it2.eq.jt2-1) c=c+sqrt(jt2*1.d0)*fvibt2
c               if(jt2.eq.it2-1) c=c+sqrt(it2*1.d0)*fvibt2

               c=c+rt*betnahv2(it2,jt2)/sqrt(4.d0*pi)

	   endif

c-------
	   if(it.eq.jt.and.it2.eq.jt2) then
             if(ivibrotp.eq.0) then
c               if(ip.eq.jp-1) c=c+sqrt(jp*1.d0)*fvibp
c               if(jp.eq.ip-1) c=c+sqrt(ip*1.d0)*fvibp

               c=c+rp*betnahvp(ip,jp)/sqrt(4.d0*pi)

	     else
		c=c+frot2p*sqrt((2.d0*2*ip+1.d0)*5.d0*(2.d0*2*jp+1.d0)
     &     /4.d0/pi)*cg(2*ip,0,2,0,2*jp,0)**2/(2.d0*2*jp+1.d0)
		c=c+frot4p*sqrt((2.d0*2*ip+1.d0)*9.d0*(2.d0*2*jp+1.d0)
     &      /4.d0/pi)*cg(2*ip,0,4,0,2*jp,0)**2/(2.d0*2*jp+1.d0)
             endif
            endif


	a(i,j)=c

 43              continue
 42              continue
 41              continue
 33              continue
 32              continue
 31              continue

	call mdiag(a,ndim)

      return
      end

c******************************************************************
	subroutine cmat(r,cpot,cpotw)
C
C Coupling matrix
C
C*****************************************************************
      implicit real*8 (a-h,o-z)
      parameter (nlevelmax=30)
      dimension cpot(nlevelmax,nlevelmax)
      dimension cpotw(nlevelmax,nlevelmax)
      common/hion/ap,zp,rp,at,zt,rt
      common/const/hbar,pi
      common/intr/ivibrotp,ivibrott,nproj,ntarg
      common/oscp/betap,betapn,omegap,lambdap,nphononp
      common/osct/betat,betatn,omegat,lambdat,nphonont
      common/osct2/betat2,betat2n,omegat2,lambdat2,nphonont2
      common/rotp/beta2p,beta4p,e2p,nrotp
      common/rott/beta2t,beta4t,e2t,nrott
      common/dimen/nlevel
      common/mut/imutual(0:nlevelmax,0:nlevelmax)
      common/eigen/ev(nlevelmax),evec(nlevelmax,nlevelmax)
      common/trans/ftr,qtrans,ntrans
        common/ahv/betnahv(0:nlevelmax,0:nlevelmax),
     &             betcahv(0:nlevelmax,0:nlevelmax),
     &             omeahv(0:nlevelmax)
        common/ahv2/betnahv2(0:nlevelmax,0:nlevelmax),
     &             betcahv2(0:nlevelmax,0:nlevelmax),
     &             omeahv2(0:nlevelmax)
        common/ahvp/betnahvp(0:nlevelmax,0:nlevelmax),
     &             betcahvp(0:nlevelmax,0:nlevelmax),
     &             omeahvp(0:nlevelmax)
        common/stability/rcoupcut
      external fct,fcp,fact,vncc,cg,fct2,fct4,fcp2,fcp4,fctt
      external fct2v,fcp2v
      external ftrans,wcc,w,ws

	do 10 i=1,nlevelmax
	do 20 j=1,nlevelmax
	cpot(i,j)=0.d0
	cpotw(i,j)=0.d0
 20	continue
 10	continue

	ndim=nlevel-ntrans
	if(nlevel.eq.1) return
	if(ndim.eq.1) goto 55
	
	   i=0
           do 31 ip=0,nproj
           do 32 it=0,ntarg
           do 33 it2=0,nphonont2
			if(nphonont2.ne.0) then
			if(imutual(it,it2).eq.0) goto 33
			endif
	   i=i+1

           j=0
           do 41 jp=0,nproj
           do 42 jt=0,ntarg
           do 43 jt2=0,nphonont2
			if(nphonont2.ne.0) then
			if(imutual(jt,jt2).eq.0) goto 43
			endif
	   j=j+1

              if(i.gt.j) then
                 cpot(i,j)=cpot(j,i)
                 cpotw(i,j)=cpotw(j,i)
                 goto 43
                 endif
	c=0.d0
	cw=0.d0
c nuclear coupling

	do 50 k=1,ndim
	c=c+vncc(r,ev(k))*evec(i,k)*evec(j,k)
	cw=cw+wcc(r,ev(k))*evec(i,k)*evec(j,k)
 50	continue
        cpotw(i,j)=cw

c coulomb coupling 

	   if(ip.eq.jp.and.it2.eq.jt2) then
             if(ivibrott.eq.0) then
c               if(it.eq.jt-1) c=c+sqrt(jt*1.d0)*fct(r)
c               if(jt.eq.it-1) c=c+sqrt(it*1.d0)*fct(r)

               if(it.ne.jt) then
               cc=betcahv(it,jt)*fct(r)
               if(betat.ne.0.d0) cc=cc/betat
               else
               cc=betcahv(it,jt)*fct2v(r)/sqrt(4.d0*pi)
               endif

               c=c+cc

	     else
		c=c+sqrt((2.d0*2*it+1.d0)*5.d0*(2.d0*2*jt+1.d0)
     &     /4.d0/pi)*cg(2*it,0,2,0,2*jt,0)**2/(2.d0*2*jt+1.d0)*fct2(r)
     &		   +sqrt((2.d0*2*it+1.d0)*9.d0*(2.d0*2*jt+1.d0)
     &     /4.d0/pi)*cg(2*it,0,4,0,2*jt,0)**2/(2.d0*2*jt+1.d0)*fct4(r)
            endif
	   endif

	   if(it.eq.jt.and.it2.eq.jt2) then
             if(ivibrotp.eq.0) then
c               if(ip.eq.jp-1) c=c+sqrt(jp*1.d0)*fcp(r)
c               if(jp.eq.ip-1) c=c+sqrt(ip*1.d0)*fcp(r)

               if(ip.ne.jp) then
               cc=betcahvp(ip,jp)*fcp(r)
               if(betap.ne.0.d0) cc=cc/betap
               else
               cc=betcahvp(ip,jp)*fcp2v(r)/sqrt(4.d0*pi)
               endif
               c=c+cc

	     else
		c=c+sqrt((2.d0*2*ip+1.d0)*5.d0*(2.d0*2*jp+1.d0)
     &     /4.d0/pi)*cg(2*ip,0,2,0,2*jp,0)**2/(2.d0*2*jp+1.d0)*fcp2(r)
     &            +sqrt((2.d0*2*ip+1.d0)*9.d0*(2.d0*2*jp+1.d0)
     &     /4.d0/pi)*cg(2*ip,0,4,0,2*jp,0)**2/(2.d0*2*jp+1.d0)*fcp4(r)
             endif
            endif

	   if(ip.eq.jp.and.it.eq.jt) then
c               if(it2.eq.jt2-1) c=c+sqrt(jt2*1.d0)*fctt(r)
c               if(jt2.eq.it2-1) c=c+sqrt(it2*1.d0)*fctt(r)

               if(it2.ne.jt2) then
               cc=betcahv2(it2,jt2)*fctt(r)
               if(betat2.ne.0.d0) cc=cc/betat2
               else
               cc=betcahv2(it2,jt2)*fct2v(r)/sqrt(4.d0*pi)
               endif
               c=c+cc

	   endif

c excitation energy
           if(it.eq.jt.and.ip.eq.jp.and.it2.eq.jt2) then
		if(ivibrott.eq.0) then
c		c=c+it*omegat
		c=c+omeahv(it)
		else
		c=c+2.d0*it*(2.d0*it+1.d0)/6.d0*e2t
		endif
		if(ivibrotp.eq.0) then
c		c=c+ip*omegap
		c=c+omeahvp(ip)
		else
		c=c+2.d0*ip*(2.d0*ip+1.d0)/6.d0*e2p
		endif
c	        c=c+it2*omegat2
	        c=c+omeahv2(it2)
 	   endif

	cpot(i,j)=c

 43              continue
 42              continue
 41              continue
 33              continue
 32              continue
 31              continue

      c0=cpot(1,1)
      do 51 i=1,ndim
          cpot(i,i)=cpot(i,i)-vn(r)
          cpotw(i,i)=cpotw(i,i)-(w(r)-ws(r))
c         cpot(i,i)=cpot(i,i)-c0
 51      continue

c    transfer channel
 55      if(ntrans.eq.1) then
            cpot(1,nlevel)=ftrans(r)
            cpot(nlevel,1)=ftrans(r)
	    cpot(nlevel,nlevel)=cpot(1,1)-qtrans
            endif

      do i=1,nlevel
      do j=1,nlevel
        if(i.ne.j.and.r.lt.rcoupcut) then 
         cpot(i,j)=0.d0
         cpotw(i,j)=0.d0
        endif
        enddo
        enddo

      return
      end

c******************************************************************
	subroutine cmat1
C
C Coupling matrix
C
C*****************************************************************
      implicit real*8 (a-h,o-z)
      parameter (nlevelmax=30)

      common/hion/ap,zp,rp,at,zt,rt
      common/const/hbar,pi
      common/intr/ivibrotp,ivibrott,nproj,ntarg
      common/oscp/betap,betapn,omegap,lambdap,nphononp
      common/osct/betat,betatn,omegat,lambdat,nphonont
      common/osct2/betat2,betat2n,omegat2,lambdat2,nphonont2
      common/rotp/beta2p,beta4p,e2p,nrotp
      common/rott/beta2t,beta4t,e2t,nrott
      common/dimen/nlevel
      common/mut/imutual(0:nlevelmax,0:nlevelmax)
      common/eigen/ev(nlevelmax),evec(nlevelmax,nlevelmax)
      common/trans/ftr,qtrans,ntrans
        common/ahv/betnahv(0:nlevelmax,0:nlevelmax),
     &             betcahv(0:nlevelmax,0:nlevelmax),
     &             omeahv(0:nlevelmax)
        common/ahv2/betnahv2(0:nlevelmax,0:nlevelmax),
     &             betcahv2(0:nlevelmax,0:nlevelmax),
     &             omeahv2(0:nlevelmax)
        common/ahvp/betnahvp(0:nlevelmax,0:nlevelmax),
     &             betcahvp(0:nlevelmax,0:nlevelmax),
     &             omeahvp(0:nlevelmax)
        common/coup1/cpot1(nlevelmax)

	do 10 i=1,nlevelmax
	cpot1(i)=0.d0
 10	continue

	ndim=nlevel-ntrans
	if(nlevel.eq.1) return
	if(ndim.eq.1) goto 55
	
	   i=0
           do 31 ip=0,nproj
           do 32 it=0,ntarg
           do 33 it2=0,nphonont2
			if(nphonont2.ne.0) then
			if(imutual(it,it2).eq.0) goto 33
			endif
	   i=i+1
             
           c=0.d0
c excitation energy
		if(ivibrott.eq.0) then
c		c=c+it*omegat
		c=c+omeahv(it)
		else
		c=c+2.d0*it*(2.d0*it+1.d0)/6.d0*e2t
		endif
		if(ivibrotp.eq.0) then
c		c=c+ip*omegap
		c=c+omeahvp(ip)
		else
		c=c+2.d0*ip*(2.d0*ip+1.d0)/6.d0*e2p
		endif
c	        c=c+it2*omegat2
	        c=c+omeahv2(it2)

	cpot1(i)=c

 33              continue
 32              continue
 31              continue

c    transfer channel
 55      if(ntrans.eq.1) then
	    cpot1(nlevel)=-qtrans
            endif

      return
      end

c*****************************************************************
	function v(r)
C
C  potential for the relative motion
C
C*****************************************************************
	implicit real*8(a-h,o-z)
	external vc,vn,vcent
	v=vc(r)+vn(r)+vcent(r)
	return
	end
c*****************************************************************
	function vc(r)
c
c coulomb potential
c
c*****************************************************************
	implicit real*8(a-h,o-z)
	common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
        common/coulomb/rc,sigma0
c  coulomb radius
	if(r.gt.rc) then
	vc=zp*zt/r*hbar/137.d0
	else
	vc=zp*zt*hbar/137.d0*(3.d0*rc**2-r**2)/2.d0/rc**3
	endif
	return
	end
c*****************************************************************
	function dvc(r)
c
c first derivative of the coulomb potential
c
c*****************************************************************
	implicit real*8(a-h,o-z)
	common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
c  coulomb radius
c	rc=1.2d0*(ap**(1.d0/3.d0)+at**(1.d0/3.d0))
c	if(r.gt.rc) then
	dvc=-zp*zt/r**2*hbar/137.d0
c	else
c	dvc=zp*zt*hbar/137.d0*(3.d0*rc**2-r**2)/2.d0/rc**3
c	endif
	return
	end
c****************************************************************
	function vn(r)
c
c (power of) woods-saxon potential for nuclear potential
c
c****************************************************************
	implicit real*8 (a-h,o-z)
	common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
        common/pot/v0,r0,a0,pow

	if(r.lt.50.d0) then
        vn=-v0/(1.d0+exp((r-r0)/a0/pow))**pow
	else
	vn=0.d0
        endif

	return
	end
c****************************************************************
	function w(r)
c
c imaginary part of the optical potential
c
c****************************************************************
	implicit real*8 (a-h,o-z)
	common/hion/ap,zp,rp,at,zt,rt
        common/imagp/w0,rw,aw
        external ws

	if(r.lt.50.d0) then
        w=w0/(1.d0+exp((r-rw)/aw))
	else
	w=0.d0
        endif

        w=w+ws(r)

	return
	end

c****************************************************************
	function wcc(r,xt)
c
c imaginary part of the optical potential
c
c****************************************************************
	implicit real*8 (a-h,o-z)
	common/hion/ap,zp,rp,at,zt,rt
        common/imagp/w0,rw,aw

	if(r.lt.50.d0) then
        wcc=w0/(1.d0+exp((r-rw-xt)/aw))
	else
	wcc=0.d0
        endif

	return
	end

c****************************************************************
	function ws(r)
c
c imaginary part of the optical potential (surface absorption)
c
c****************************************************************
	implicit real*8 (a-h,o-z)
	common/hion/ap,zp,rp,at,zt,rt
        common/imagps/w0s,rws,aws

	if(r.lt.50.d0) then
        ws=w0s/aws*exp((r-rws)/aws)/(1.d0+exp((r-rws)/aws))**2
	else
	ws=0.d0
        endif

	return
	end

c****************************************************************
	function vncc(r,xt)
c
c (power of) woods-saxon potential for nuclear potential 
c
c****************************************************************
	implicit real*8 (a-h,o-z)

	common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
        common/pot/v0,r0,a0,pow

	if(r.lt.50.d0) then
           rr=r-r0-xt
        vncc=-v0/(1.d0+exp(rr/a0/pow))**pow
	else
	vncc=0.d0
        endif

	return
	end

c***************************************************************
      function dvn(r)
c     
c the first derivative of vn(r)
c
c**************************************************************
      implicit real*8(a-h,o-z)
      common/pot/v0,r0,a,pow
      external vn

      if(r.lt.50.d0) then
c      dvn=v0/a*exp((r-r0)/a)/(1.d0+exp((r-r0)/a))**2
      dvn=(vn(r+1.d-5)-vn(r-1.d-5))/2.d-5
      else
      dvn=0.d0
      endif

      return
      end
c************************************************************
	function vcent(r)
c
c centrifugal potential
c
c************************************************************
	implicit real*8(a-h,o-z)
	common/const/hbar,pi
	common/angular/l
	common/dyn/e,rmass
	vcent=l*(l+1.d0)*hbar**2/2.d0/rmass/r**2
	return
	end
c************************************************************
	function dvcent(r)
c
c the first derivative of the centrifugal potential
c
c************************************************************
	implicit real*8(a-h,o-z)
	common/const/hbar,pi
	common/angular/l
	common/dyn/e,rmass
	dvcent=-2.d0*l*(l+1.d0)*hbar**2/2.d0/rmass/r**3
	return
	end

c***************************************************************
	function dv(r)
c
c the first derivative of v(r)
c
c***************************************************************
	implicit real*8(a-h,o-z)
	external dvn,dvc,dvcent

	dv=dvn(r)+dvc(r)+dvcent(r)

	return
	end
c*****************************************************************
	function fcp(r)
c
c coulomb coupling form factor for projectile phonon excitation
c
c*****************************************************************
	implicit real*8(a-h,o-z)
        common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
        common/oscp/betap,betapn,omega,lambda,nphonon

	rc=rp
	if(r.gt.rc) then
	fcp=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/r
     &                                         *(rp/r)**lambda
	else
	fcp=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/rc
     &                                         *(r/rp)**lambda
	endif

        fcp=fcp*betap/sqrt(4.d0*pi)

	return
	end

c*****************************************************************
	function fcp2v(r)
c
c coulomb coupling form factor for projectile excitation
c (for anharmonic vibration)
c
c*****************************************************************
	implicit real*8(a-h,o-z)
        common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi

	lambda=2
	rc=rp
	if(r.gt.rc) then
	fcp2v=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/r
     &                                         *(rp/r)**lambda
	else
	fcp2v=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/rc
     &                                         *(r/rp)**lambda
	endif

	return
	end

c*****************************************************************
	function fcp2(r)
c
c coulomb coupling form factor for projectile excitation
c (rotational e2 coupling)
c
c*****************************************************************
	implicit real*8(a-h,o-z)
        common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
	common/rotp/beta2p,beta4p,e2p,nrotp

	lambda=2
	rc=rp
	if(r.gt.rc) then
	fcp2=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/r
     &                                         *(rp/r)**lambda
	else
	fcp2=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/rc
     &                                         *(r/rp)**lambda
	endif

	fcp2=fcp2*(beta2p+2.d0*sqrt(5.d0/pi)*beta2p**2/7.d0)

	return
	end
c*****************************************************************
	function fcp4(r)
c
c coulomb coupling form factor for projectile excitation
c (rotational e4 coupling)
c
c*****************************************************************
	implicit real*8(a-h,o-z)
        common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
	common/rotp/beta2p,beta4p,e2p,nrotp

	lambda=4
	rc=rp
	if(r.gt.rc) then
	fcp4=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/r
     &                                         *(rp/r)**lambda
	else
	fcp4=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/rc
     &                                         *(r/rp)**lambda
	endif

	fcp4=fcp4*(beta4p+9.d0*beta2p**2/7.d0/sqrt(pi))

	return
	end
c*****************************************************************
	function fct(r)
c
c coulomb coupling form factor for the target phonon excitation
c
c*****************************************************************
	implicit real*8(a-h,o-z)
        common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
        common/osct/betat,betatn,omega,lambda,nphonon

	rc=rt
	if(r.gt.rc) then
	fct=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/r
     &                                         *(rt/r)**lambda
	else
	fct=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/rc
     &                                         *(r/rt)**lambda
	endif

        fct=fct*betat/sqrt(4.d0*pi)

	return
	end

c*****************************************************************
	function fctt(r)
c
c coulomb coupling form factor for the second target phonon excitation
c
c*****************************************************************
	implicit real*8(a-h,o-z)
        common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
        common/osct2/betat2,betat2n,omega2,lambda2,nphonon2

	rc=rt
	if(r.gt.rc) then
	fctt=3.d0/(2.d0*lambda2+1.d0)*zp*zt/137.d0*hbar/r
     &                                         *(rt/r)**lambda2
	else
	fctt=3.d0/(2.d0*lambda2+1.d0)*zp*zt/137.d0*hbar/rc
     &                                         *(r/rt)**lambda2
	endif

        fctt=fctt*betat2/sqrt(4.d0*pi)

	return
	end

c*****************************************************************
	function fct2v(r)
c
c coulomb coupling form factor for target excitation
c (for anharmonic vibration)
c
c*****************************************************************
	implicit real*8(a-h,o-z)
        common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi

	lambda=2
	rc=rt
	if(r.gt.rc) then
	fct2v=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/r
     &                                         *(rt/r)**lambda
	else
	fct2v=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/rc
     &                                         *(r/rt)**lambda
	endif

	return
	end


c*****************************************************************
	function fct2(r)
c
c coulomb coupling form factor for target excitation
c (rotational e2 coupling)
c
c*****************************************************************
	implicit real*8(a-h,o-z)
        common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
	common/rott/beta2t,beta4t,e2t,nrott

	lambda=2
	rc=rt
	if(r.gt.rc) then
	fct2=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/r
     &                                         *(rt/r)**lambda
	else
	fct2=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/rc
     &                                         *(r/rt)**lambda
	endif

	fct2=fct2*(beta2t+2.d0*sqrt(5.d0/pi)*beta2t**2/7.d0)

	return
	end

c*****************************************************************
	function fct4(r)
c
c coulomb coupling form factor for target excitation
c (rotational e4 coupling)
c
c*****************************************************************
	implicit real*8(a-h,o-z)
        common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
	common/rott/beta2t,beta4t,e2t,nrott

	lambda=4
	rc=rt
	if(r.gt.rc) then
	fct4=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/r
     &                                         *(rt/r)**lambda
	else
	fct4=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/rc
     &                                         *(r/rt)**lambda
	endif

	fct4=fct4*(beta4t+9.d0*beta2t**2/7.d0/sqrt(pi))

	return
	end

c***************************************************************
      function ddvn(r)
c     
c the second derivative of vn(r)
c
c**************************************************************
      implicit real*8(a-h,o-z)
      common/pot/v0,r0,a,pow
      external dvn

      if(r.lt.50.d0) then
c      ddvn=v0/a**2*exp((r-r0)/a)/(1.d0+exp((r-r0)/a))**2
c     &    -2.d0*v0/a**2*exp(2.d0*(r-r0)/a)
c     &                           /(1.d0+exp((r-r0)/a))**3
      ddvn=(dvn(r+1.d-5)-dvn(r-1.d-5))/2.d-5
      else
      ddvn=0.d0
      endif

      return
      end
c****************************************************************
	function ftrans(r)
c
c coupling form factor for transfer reactions
c
c***************************************************************
	implicit real*8(a-h,o-z)
	common/const/hbar,pi
        common/trans/ftr,qtrans,ntrans
	external dvn

	ftrans=ftr*dvn(r)
	return
	end

c**************************************************************
                function pn(n,x)
c
c legendre polynomials
c
c**************************************************************
                implicit real*8(a-h,o-z)
                if (n.eq.0) then
                pn=1.d0
                go to 10
                elseif (n.eq.1) then
                   pn=x
                   go to 10
                   else
                      pm=1.d0
                      pz=x
                      do 20 i=1,n-1
                         pp=((2.d0*i+1.d0)*x*pz-i*pm)/(i+1.d0)
                         pm=pz
                         pz=pp
 20                      continue
                         pn=pz
                         endif
 10                      return
                         end


c**************************************************************
                subroutine legendre(lmax,x)
c
c legendre polynomials
c
c**************************************************************
                implicit real*8(a-h,o-z)
                common/legendrepol/plgndr(0:3000)

                plgndr(0)=1.d0
                plgndr(1)=x

                      pm=1.d0
                      pz=x
                      do 20 i=1,lmax-1
                         pp=((2.d0*i+1.d0)*x*pz-i*pm)/(i+1.d0)
                         plgndr(i+1)=pp

                         pm=pz
                         pz=pp
 20                      continue

                         return
                         end


c*************************************************************
	subroutine mdiag(a,n)
C
C The Jacobi method to compute all eigenvalues and eigenvectors 
C of a real symmetric matrix A, which is of size N by N, stored 
C in a physical NP by NP array. On output, the elements of A 
C above the diagonal are destroyed. D returns the eigenvalues of 
C A in its first N elements. V is a matrix with the same logical 
C and physical dimensions as A whose columns contain, on output, 
C the normalized eigenvectors of A. NROT returns the number of 
C Jacobi rotations which were required. 
C
C N.B.; The I-th component of the eigenvector for the K-th 
C       eigenvalue is given by V(I,K). 
C
C************************************************************
	implicit real*8(a-h,o-z)
	parameter(nmax=100)
        parameter (nlevelmax=30)
        dimension a(nlevelmax,nlevelmax)
	dimension b(nmax),z(nmax)
        common/eigen/d(nlevelmax),v(nlevelmax,nlevelmax)

	do 12 ip=1,n
	do 11 iq=1,n
	v(ip,iq)=0.d0
 11	continue
	v(ip,ip)=1.d0
 12	continue

	do 13 ip=1,n
	b(ip)=a(ip,ip)
	d(ip)=b(ip)
	z(ip)=0.d0
 13	continue

	nrot=0

	do 24 i=1,50
	sm=0.d0

	do 15 ip=1,n-1
	do 14 iq=ip+1,n
	sm=sm+abs(a(ip,iq))
 14	continue
 15     continue

	if(sm.eq.0.d0) return

	if(i.lt.4) then
	tresh=0.2d0*sm/n**2
	else
	tresh=0.d0
	endif

	do 22 ip=1,n-1
	do 21 iq=ip+1,n
	g=100.d0*abs(a(ip,iq))
	
	if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip)))
     &     .and.(abs(d(iq))+g.eq.abs(d(iq)))) then
	   a(ip,iq)=0.d0
	elseif(abs(a(ip,iq)).gt.tresh) then
	   h=d(iq)-d(ip)

	   if(abs(h)+g.eq.abs(h)) then
	     t=a(ip,iq)/h
	   else
	     theta=0.5d0*h/a(ip,iq)
	     t=1.d0/(abs(theta)+sqrt(1.d0+theta**2))
	     if(theta.lt.0.d0) t=-t
	   endif

	c=1.d0/sqrt(1.d0+t**2)
	s=t*c
	tau=s/(1.d0+c)
	h=t*a(ip,iq)
	z(ip)=z(ip)-h
	z(iq)=z(iq)+h
	d(ip)=d(ip)-h
	d(iq)=d(iq)+h
	a(ip,iq)=0.d0

	do 16 j=1,ip-1
	g=a(j,ip)
	h=a(j,iq)
	a(j,ip)=g-s*(h+g*tau)
	a(j,iq)=h+s*(g-h*tau)
 16	continue

	do 17 j=ip+1,iq-1
	g=a(ip,j)
	h=a(j,iq)
	a(ip,j)=g-s*(h+g*tau)
	a(j,iq)=h+s*(g-h*tau)
 17	continue

	do 18 j=iq+1,n
	g=a(ip,j)
	h=a(iq,j)
	a(ip,j)=g-s*(h+g*tau)
	a(iq,j)=h+s*(g-h*tau)
 18	continue

	do 19 j=1,n
	g=v(j,ip)
	h=v(j,iq)
	v(j,ip)=g-s*(h+g*tau)
	v(j,iq)=h+s*(g-h*tau)
 19	continue
	
	nrot=nrot+1
	endif

 21	continue
 22	continue

	do 23 ip=1,n
	b(ip)=b(ip)+z(ip)
	d(ip)=b(ip)
	z(ip)=0.d0

 23	continue
 24	continue

	pause '50 iterations should never happen'

	return
	end

C****************************************************************
      function cg(j1,m1,j2,m2,j3,m3)
C
C  Clebsh-Gordan Coefficient <j1 m1 j2 m2 | j3 m3>
C****************************************************************
      implicit real*8(a-h,o-z)
	external fact

	if(m1+m2.ne.m3) then
	cg=0.d0
	return
	endif

	if(j3.lt.abs(j1-j2)) then
	cg=0.d0
	return
	endif

	if(j3.gt.j1+j2) then
	cg=0.d0
	return
	endif

      ka=j1+j2-j3
      kb=j3+j1-j2
      kc=j2+j3-j1
      kd=j1+j2+j3+1
      del=sqrt(fact(ka)*fact(kb)*fact(kc)/fact(kd))

      cg=0.d0

      do 10 n=0,max(j1+j2-j3,j1-m1,j2+m2)

      ka1=j1+j2-j3-n
      if(ka1.lt.0.d0) go to 10
      ka2=j3-j2+m1+n
      if(ka2.lt.0.d0) go to 10
      ka3=j3-j1-m2+n
      if(ka3.lt.0.d0) go to 10
      ka4=j1-m1-n
      if(ka4.lt.0.d0) go to 10
      ka5=j2+m2-n
      if(ka5.lt.0.d0) go to 10
      an=n*1.d0

      cg=cg+(-1.d0)**n
     &  /(fact(n)*fact(ka1)*fact(ka2)*fact(ka3)*fact(ka4)*fact(ka5))

 10   continue

 20   cg=cg*sqrt(fact(j1+m1)*fact(j1-m1))
      cg=cg*sqrt(fact(j2+m2)*fact(j2-m2))
      cg=cg*sqrt(fact(j3+m3)*fact(j3-m3))

      cg=cg*sqrt(2.d0*j3+1.d0)*del

      return
      end
C***************************************************************
      subroutine matinv(nmax,c,d)
C     
C subroutine for calculating the inversion of a matrix C, 
C whose dimension is NMAX
C
C***************************************************************
      implicit complex*16 (a-h,o-z)
      parameter (nlevelmax=30)
      dimension u(nlevelmax),v(nlevelmax),
     &          c(nlevelmax,nlevelmax),d(nlevelmax,nlevelmax)
      deter=1.0
      do 1 m=1,nmax
      do 1 n=1,nmax
      d(n,m)=0.0
      if(n.eq.m) d(n,m)=1.0
 1    continue
      do 10 n=1,nmax
      t=c(n,n)
      if (abs(t).lt.1.e-20) go to3
      go to 6
 3    j=n
 4    j=j+1
      if (j.gt.nmax) go to 11
      t=c(n,j)
      deter=-deter
      if (abs(t).le.1.e-30) go to 4
      do 5 k=1,nmax
      u(k)=c(n,k)
      v(k)=d(n,k)
      c(n,k)=c(j,k)
      d(n,k)=d(j,k)
      c(j,k)=u(k)
 5    d(j,k)=v(k)
 6    do 7 k=1,nmax
      if (k.eq.n) go to 7
      a=c(k,n)/c(n,n)
      do 8 l=1,nmax
      c(k,l)=c(k,l)-a*c(n,l)                              
      d(k,l)=d(k,l)-a*d(n,l)
 8    continue
 7    continue
      b=c(n,n)
      deter=b*deter
      do 10 m=1,nmax
      c(n,m)=c(n,m)/b
      d(n,m)=d(n,m)/b
 10   continue
      return
 11   print 12
 12   format('matrix not invertible')
      return
      end

                         function fact(n)
                         implicit real*8(a-h,o-z)
			 if(n.lt.0) then
			 fact=0.d0
			 return
			 endif
                         if (n.eq.0) then
                            fact=1.d0
                            go to 10
                            endif
                            fact=1.d0
                            do 20 i=1,n
                               fact=fact*i*1.d0
 20                            continue
 10                            return
                               end

C*****************************************************************
C Subroutine for the Coulomb wave function
C*****************************************************************
      SUBROUTINE JFLGAM(XD,YD,PAR,PAI,NBCHIF)
      DOUBLE PRECISION XD,YD,PAR,PAI,TEST,C,HLO2PI,PI,PIS2,PIS4
      DOUBLE PRECISION X,Y,U,V,TRA,TRA1,ALO2PI,RAC2,COSI,SINI
      DOUBLE PRECISION COS2I,SIN2I,ZMOD,DEPI
      DOUBLE PRECISION XX
      DOUBLE PRECISION ALOPI
      DOUBLE PRECISION SUPINT
      DIMENSION TEST(7),C(6)
      DATA TEST/2.9152D+7,2.2958D+3,1.4124D+2,3.9522D+1,19.6611D0,
     &12.791D0,-10.D0/
      DATA C/8.333333333333333D-2,-2.777777777777777D-3,
     &7.936507936507937D-4,-5.952380952380952D-4,8.417508417508418D-4,
     &-1.917526917526918D-3/
      DATA HLO2PI/0.918938533204672/,PI/3.141592653589793/
      DATA PIS2/1.570796326794897/,PIS4/0.785398163397448/
      DATA ALO2PI/1.837877066409345/,RAC2/0.3465735902799726/
      DATA DEPI/6.283185307179586/,ALOPI/1.1447298858494002/
      DATA SUPINT/2147483647.D0/
      NBCHIF=15
      X=DABS(XD)
      XX=X
      IF(YD)1,2,1
    1 Y=DABS(YD)
      KR=1
      I=DMOD(10.99D0-X,SUPINT)
C     TRANSLATION
      IF(I)3,3,4
    4 TRA=I
      X=X+TRA
C     LOGARITHME(X+IY) (X,Y POSITIFS)
    3 IF(X-Y)5,6,7
    5 TRA1=X/Y
      IF(TRA1)8,8,9
    8 U=DLOG(Y)
      V=PIS2
      GO TO 10
    6 U=RAC2+DLOG(X)
      V=PIS4
      GO TO 10
    9 TRA=Y*DSQRT(1.D0+TRA1*TRA1)
      TRA1=Y/X
   11 U=DLOG(TRA)
      V=DATAN(TRA1)
   10 GO TO (12,19,23),KR
    7 TRA1=Y/X
      TRA=X*DSQRT(1.D0+TRA1*TRA1)
      GO TO 11
   12 KR=2
C     DEVELOPPEMENT ASYMPTOTIQUE ( X SUPERIEUR A 10 )
      TRA=X-0.5D0
      PAI=V*TRA+Y*(U-1.D0)
      PAR=-X+HLO2PI+U*TRA-V*Y
      ZMOD=X+Y
      IF(ZMOD-TEST(1))13,13,14
   13 TRA=X*X+Y*Y
      COSI=X/TRA
      SINI=Y/TRA
      SIN2I=(SINI*COSI)+(SINI*COSI)
      COS2I=(COSI+SINI)*(COSI-SINI)
      K=1
      GO TO 15
   16 TRA=COSI*COS2I-SINI*SIN2I
      SINI=SINI*COS2I+COSI*SIN2I
      COSI=TRA
   15 PAR=PAR+C(K)*COSI
      PAI=PAI-C(K)*SINI
      K=K+1
      IF(ZMOD-TEST(K))16,16,14
C     TRANSLATION INVERSE
   17 I=I-1
      X=I
      X=XX+X
      GO TO 3
   19 PAR=PAR-U
      PAI=PAI-V
   14 IF(I-1)18,60,17
   60 IF(XD)17,61,17
C     CONTROLE DU QUADRANT
   18 IF(XD)20,61,21
   61 TRA=PI*Y
      IF(TRA-1.D-2)300,300,301
  300 PAR= TRA*(2.D0+TRA*(-2.D0+TRA*(1.333333333333333D0+TRA*(
     &-0.6666666666666666D0+TRA*(0.2666666666666666D0+TRA*(
     &-0.08888888888888888D0+TRA*0.02539682539682540D0))))))
      TRA1=-DLOG(Y)-DLOG(PAR)
      GO TO 302
  301 PAR=1.D0-DEXP(-TRA-TRA)
      TRA1=-DLOG(Y*PAR)
  302 PAR=0.5D0*(ALO2PI-TRA+TRA1)
      PAI=PAI-PIS2
   21 IF(YD)28,100,100
C     X+IY CHANGE EN -X-IY
   20 TRA=PI*Y
      PAR=ALO2PI-U-PAR-TRA
      PAI=PI-V-PAI
      TRA=DEXP(-TRA-TRA)
      X=PI*DMOD(X,2.D0)
      SINI=(1.D0-TRA)*DCOS(X)
      COSI=(1.D0+TRA)*DSIN(X)
      KR=3
      X=DABS(COSI)
      Y=DABS(SINI)
      GO TO 3
   23 IF(COSI)24,25,25
   24 V=PI-DSIGN(V,SINI)
      GO TO 26
   25 IF(SINI)27,26,26
   27 V=-V
   26 PAR=PAR-U
      PAI=PAI-V
      IF(YD)100,100,28
   28 PAI=-PAI
C     ARGUMENT DANS -PI,PI
  100 TRA=DABS(PAI/DEPI)
      IF(TRA-1.D+15)203,204,204
  204 NBCHIF=0
      PAI=0.D0
      GO TO 29
  203 IF(TRA-1.D0)205,205,206
  206 NBCHIF=DLOG10(TRA)
      NBCHIF=14-NBCHIF
      TRA=DMOD(TRA,SUPINT)
      PAI=DMOD(TRA,1.D0)*DSIGN(DEPI,PAI)
  205 IF(DABS(PAI)-PI)29,29,207
  207 PAI=PAI-DSIGN(DEPI,PAI)
      GO TO 29
C     JFLGAM REEL
    2 PAI=0.D0
      IF(XD)31,32,33
C     CONDITIONS D EXISTENCE
   32 WRITE (6,1000)
 1000 FORMAT (21H JFLGAM(0) EST INFINI)
      GO TO 50
   31 IF(X-4503599627370496.D0)34,35,35
   35 WRITE (6,1001)
 1001 FORMAT (30H ARGUMENT DE JFLGAM TROP GRAND)
      GO TO 50
   34 Y=DMOD(X,SUPINT)
      IF(Y)400,36,400
  400 IF(Y-0.99D0)33,33,405
  405 TRA=IDINT(Y+0.1D0)
      IF(DABS(Y-TRA)-5.D-15)36,36,33
   36 WRITE (6,1002)
 1002 FORMAT (28H JFLGAM (-ENTIER) EST INFINI)
   50 PAR=1.D+74
      NBCHIF=0
      GO TO 29
C     TRANSLATION
   33 I=DMOD(10.99D0-X,SUPINT)
      IF(I)37,37,38
   38 TRA=I
      X=X+TRA
C     DEVELOPPEMENT ASYMPTOTIQUE
   37 Y=DLOG(X)
      PAR=-X+HLO2PI+Y*(X-0.5D0)
      IF(X-TEST(1))39,39,43
   39 COSI=1.D0/X
      COS2I=COSI*COSI
      K=1
      GO TO 41
   42 COSI=COSI*COS2I
   41 PAR=PAR+C(K)*COSI
      K=K+1
      IF(X-TEST(K))42,42,40
C     TRANSLATION INVERSE
   44 X=X-1.D0
   48 Y=DLOG(X)
      PAR=PAR-Y
      I=I-1
   40 IF(I-1)43,49,44
   49 X=DABS(XD)
      GO TO 48
C     X NEGATIF
   43 IF(XD)45,29,29
   45 PAR=ALOPI-PAR-Y
      Y=PI*DMOD(X,2.D0)
      Y=-DSIN(Y)
      IF(Y)46,36,47
   46 Y=-Y
      PAI=PI
   47 PAR=PAR-DLOG(Y)
      ENTRY JFLGV1
   29 RETURN
      END

      SUBROUTINE YFCLEN(ETA,RO,U,UP,V,VP,SIGMA0,IDIV,NN)
      IMPLICIT COMPLEX*16(A-D,F-H),REAL*8(E,P-Z)
C
      IF(NN.EQ.1)GO TO 20
C
      ETA2=ETA*ETA
      FA=DCMPLX(1.D0,ETA)
      M=0.25*ETA+4.D1
C
C          POLYNOMES DE TCHEBICHEV JUSQU'AU RANG M
C
      K=M+2
      X=(ETA+ETA)/RO
      XX=X+X-1.D0
      T0=1.D0
      T1=XX
      XX=XX+XX
      DO 6 J=2,K
      TJ=XX*T1-T0
      T0=T1
      T1=TJ
    6 CONTINUE
      TM=T1
      TL=T0
C
C          INITIALISATION
C
      AM=(0.D0,0.D0)
      AL=(1.D-40,1.D-40)
      BN=(0.D0,0.D0)
      BM=(1.D-40,1.D-40)
      BL=(0.D0,0.D0)
      BK=DCMPLX(4.D0*DFLOAT(M+1),0.D0)*AL+BM
      F=(0.D0,0.D0)
      G =(0.D0,0.D0)
      GP=(0.D0,0.D0)
C
C          RECURRENCE DESCENDANTE
C
      K=M
   10 R=K
      TK=XX*TL-TM
      TM=TL
      TL=TK
      HK=DCMPLX(TK,0.D0)
      C1=DCMPLX(R*(R+1.D0)-ETA2,ETA*(R+R+1.D0))
      C2=(4.D0,0.D0)*DCMPLX(R+1.D0,0.D0)
      C2=C2*DCMPLX(-R-1.D0,ETA*3.D0)
      C3=FA*DCMPLX(-R-R-4.D0,ETA)
      C4=DCMPLX((7.D0*R+5.D0)/4.D0,0.D0)
      C5=DCMPLX(R+R+2.D0,0.D0)
      C6=DCMPLX((R+3.D0)/4.D0,0.D0)
      AK=(C2*AL+C3*AM-C4*BL-C5*BM-C6*BN)/C1
      J=K/2
      J=K-J-J
      IF(J)1,2,1
    1 F=F-AK
      GO TO 3
    2 F=F+AK
    3 Z=ABS(AK)
      G=G+HK*AK
      GP=GP+HK*BK
C
C          F=A0/2-A1+A2-A3+A4-A5+...
C
C          CONGRUENCE MODULO 10**60
C
      IF(Z-1.D60)4,5,5
    5 D=(1.D60,0.D0)
      AK=AK/D
      AL=AL/D
      AM=AM/D
      BK=BK/D
      BL=BL/D
      BM=BM/D
      BN=BN/D
      F=F/D
      G=G/D
      GP=GP/D
    4 IF(K)8,8,9
    9 D=DCMPLX(4.D0*R,0.D0)
      BJ=D*AK+BL
      AM=AL
      AL=AK
      BN=BM
      BM=BL
      BL=BK
      BK=BJ
      K=K-1
      GO TO 10
C
C          NORMALISATION ET CALCUL DE Z(RO)
C
    8 D=(0.5D0,0.D0)*AK
      F=F-D
      G=G-D
      GP=GP-(0.5D0,0.D0)*BK
      D=DCMPLX(0.D0,-ETA*DLOG(2.D0)+SIGMA0)
      AXPO=EXP(D)
      F=F/AXPO
      G=G/F
      GP=GP/F
C
C          CALCUL DE F ET G
C
      D=DCMPLX(0.D0,RO-ETA*DLOG(RO))
      AXPO=EXP(D)
      D=G*AXPO
      GP=AXPO*(DCMPLX(0.D0,1.D0-ETA/RO)*G-DCMPLX(X/RO,0.D0)*GP)
      V=D
      D=(0.D0,-1.D0)*D
      U=D
      VP=GP
      GP=(0.D0,-1.D0)*GP
      UP=GP
      IDIV=0
      RETURN
C
C          SERIE ORIGINE
C
   20 PI=3.141592653589793D0
      XA=0.577215664901533D0
      RO2=RO*RO
      ETAP=ETA+ETA
      PIETA=PI*ETA
      Z=138.15510557964276D0
      IDIV=0
      IF(DABS(PIETA)-Z)21,21,22
   22 INDG=IDINT(PIETA/Z)
      IDIV=60*INDG
      IF(ETA.LT.0) IDIV=0
      PIETA=PIETA-Z*DFLOAT(INDG)
   21 PIETA2=0.5D0*PIETA
      P=DEXP(PIETA2)*DSQRT(DSINH(PIETA)/PIETA)
      CALL JFDELG(1.D0,ETA,PAR,PAI,NB)
      Z1=ETAP*(XA+XA+DLOG(2.D0)-1.D0+PAR)
      U0=0.D0
      U1=RO
      V0=1.D0
      V1=Z1*RO
      U=U0+U1
      V=V0+V1
      UP=1.D0
      VP=Z1
      XN=2.D0
      DO 104N=2,10000
      XN1=XN*(XN-1.D0)
      U2=(ETAP*RO*U1-RO2*U0)/XN1
      U=U+U2
      V2=(ETAP*RO*V1-RO2*V0-ETAP*(XN+XN-1.D0)*U2)/XN1
      V=V+V2
      UP=UP+XN*U2/RO
      VP=VP+XN*V2/RO
      IF(DABS(U2/U).GT.1.D-14)GOTO18
      IF(DABS(V2/V).LE.1.D-14)GOTO19
   18 U0=U1
      U1=U2
      V0=V1
      V1=V2
      XN=XN+1.D0
  104 CONTINUE
   19 PP=V+ETAP*U*DLOG(RO)
      W=U/P
      WP=UP/P
      V=P*PP
      VP=P*(VP+ETAP*(UP*DLOG(RO)+U/RO))
      U=W
      UP=WP
      RETURN
      END

      SUBROUTINE YFASYM(ETA,RAU,FO,FPO,GO,GPO,SIGO,IEXP)
      IMPLICIT REAL*8 (A-H,O-Z)
      IEXP=0
      TRB=0.D0
      RAU2=RAU+RAU
      ETAC=ETA*ETA
      CALL JFLGAM(1.D0,ETA,TRA,SIGO,NTRUC)
   40 N=0
      PS=1.D0
      GS=0.D0
      PT=0.D0
      GT=1.D0-ETA/RAU
      SF=PS
      SG=GS
      SPF=PT
      SPG=GT
   45 DENOM= DFLOAT(N+1)*RAU2
      AN= DFLOAT(N+N+1)*ETA/DENOM
      BN= (ETAC-DFLOAT(N*(N+1)))/DENOM
      PS1=AN*PS-BN*PT
      GS1=AN*GS-BN*GT-PS1/RAU
      PT1=AN*PT+BN*PS
      GT1=AN*GT+BN*GS-PT1/RAU
   42 SF=SF+PS1
      SG=SG+GS1
      SPF=SPF+PT1
      SPG=SPG+GT1
      N=N+1
      IF(N-17)46,48,44
   48 TRA=PS*PS+PT*PT
   44 TRB=PS1*PS1+PT1*PT1
      TEST=TRA-TRB
      IF(TEST)47,47,46
   46 PS=PS1
      GS=GS1
      PT=PT1
      GT=GT1
      TRA=TRB
      GOTO 45
   47 TETAO= RAU-ETA*DLOG (RAU2)+SIGO
      TRA= DSIN(TETAO)
      TRB=DCOS(TETAO)
      GO=SF*TRB-SPF*TRA
      GPO=SG*TRB-SPG*TRA
      FO=SPF*TRB+SF*TRA
      FPO=SPG*TRB+SG*TRA
      RETURN
      END

      SUBROUTINE DFCOUL(ETA,RO,F,FP,G,GP,SIGMA,L,IEXP)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(*),FP(*),G(*),GP(*),IEXP(*),SIGMA(*)
      DATA DEPI/6.283185307179586D0/
      ETAC=ETA*ETA
      L1=L+1
      CALL DFCZ0(ETA,RO,F0,FP0,G0,GP0,S,I)
      F(1)=F0
      FP(1)=FP0
      G(1)=G0
      GP(1)=GP0
      IEXP(1)=I
      SIGMA(1)=S
      IF(L)1,1,2
    1 RETURN
    2 LINF=0
      IND=0
      IF((ETA.GT.0).AND.(RO.LT.(ETA+ETA)))GO TO 21
      Z=ETA+DSQRT(ETAC+DFLOAT(L*(L+1)))
      IF(RO.GE.Z)GO TO 20
    7 ROINF=ETA+DSQRT(ETAC+DFLOAT(LINF*(LINF+1)))
      IF(RO-ROINF)3,4,4
    4 IF(LINF-L)5,6,6
    5 LINF=LINF+1
      GO TO 7
    3 IND=1
    6 LIN=LINF+1
   20 XM=1.D0
      IF(IND.EQ.0)LIN=L1
      DO 8 J=2,LIN
      ZIG=(DSQRT(ETAC+XM*XM))/XM
      ZAG=ETA/XM+XM/RO
      F(J)=(ZAG*F(J-1)-FP(J-1))/ZIG
      FP(J)=ZIG*F(J-1)-ZAG*F(J)
      G(J)=(ZAG*G(J-1)-GP(J-1))/ZIG
      GP(J)=ZIG*G(J-1)-ZAG*G(J)
      IEXP(J)=I
      SIG=SIGMA(J-1)+DATAN(ETA/(J-1))
      IPI=SIG/DEPI
      SIG=SIG-IPI*DEPI
      IF(SIG)60,50,70
   60 IF(SIG.LT.(-DEPI/2.D0))SIG=SIG+DEPI
      GO TO 50
   70 IF(SIG.GT.(DEPI/2.D0))SIG=SIG-DEPI
   50 SIGMA(J)=SIG
    8 XM=XM+1.D0
      IF(IND.EQ.0)RETURN
      GO TO 22
   21 LIN=1
   22 FTEST=F(LIN)
      FPTEST=FP(LIN)
      LMAX=LINF+25+IDINT(5.D0*DABS(ETA))
      IF(LMAX-L)9,10,10
    9 LMAX=L
   10 FI=1.D0
      FPI=1.D0
   13 XM=LMAX+1
      ZIG=(DSQRT(ETAC+XM*XM))/XM
      ZAG=ETA/XM+XM/RO
      FL=(ZAG*FI+FPI)/ZIG
      FPL=ZAG*FL-ZIG*FI
      IF(DABS(FL)-1.D15)26,27,27
   26 IF(DABS(FPL)-1.D15)28,27,27
   27 FL=FL*1.D-15
      FPL=FPL*1.D-15
   28 FI=FL
      FPI=FPL
      IF(LMAX-L)11,11,12
   12 LMAX=LMAX-1
      GO TO 13
   11 F(LMAX+1)=FL
      FP(LMAX+1)=FPL
      IF(LMAX-LINF)15,15,14
   14 GO TO 12
   15 FACT=FTEST/F(LIN)
      FACTP=FPTEST/FP(LIN)
      INDICE=I/60
      XM=LINF
      DO 16 J=LIN,L1
      F(J)=F(J)*FACT
      FP(J)=FP(J)*FACTP
   25 IF(J.EQ.1)GO TO 16
      ZIG=(DSQRT(ETAC+XM*XM))/XM
      ZAG=ETA/XM+XM/RO
      G(J)=(ZAG*G(J-1)-GP(J-1))/ZIG
      GP(J)=ZIG*G(J-1)-ZAG*G(J)
      IF(DABS(G(J))-1.D60)17,18,18
   17 IF(DABS(GP(J))-1.D60)19,18,18
   18 G(J)=G(J)/1.D60
      GP(J)=GP(J)/1.D60
      INDICE=INDICE+1
   19 IEXP(J)=INDICE*60
      A=FP(J)*G(J)
      B=-F(J)*GP(J)
      IF(A-1.D0)29,30,30
   29 I1=IDINT(DLOG10(A))
      I2=IDINT(DLOG10(B))
      GO TO 31
   30 I1=IDINT(DLOG10(A))+1
      I2=IDINT(DLOG10(B))+1
   31 F(J)=F(J)*1.D1**(-I2)
      FP(J)=FP(J)*1.D1**(-I1)
      SIG=SIGMA(J-1)+DATAN(ETA/(J-1))
      IPI=SIG/DEPI
      SIG=SIG-IPI*DEPI
      IF(SIG)61,51,71
   61 IF(SIG.LT.(-DEPI/2.D0))SIG=SIG+DEPI
      GO TO 51
   71 IF(SIG.GT.(DEPI/2.D0))SIG=SIG-DEPI
   51 SIGMA(J)=SIG
   16 XM=XM+1.D0
      RETURN
      END

      SUBROUTINE YFIREG(ETA,RO,G0,GP0)
      IMPLICIT REAL*8(A-H,O-Z)
      IF(ETA.LE.0.D0)GOTO250
      IF(ETA.LE.3.D0)GOTO251
      IF(ETA.LE.1.D1)GOTO252
      IF(ETA.LE.18.D0)GOTO253
      IF(ETA.LE.22.D0)GOTO254
      IF(RO.LE.0.3D0+(3.D1-ETA)/8.D1)GOTO200
C   SERIE DE TAYLOR DEPART RAU0
  300 CONTINUE
      RAU0=1.666666666666667D0*DABS(ETA)+7.5D0
      CALL YFASYM(ETA,RAU0,F0,FP0,G0,GP0,SIGMA0,IEXP)
      X=RAU0-RO
      X2=X*X
      X3=X*X2
      UNR=1.D0/RAU0
      ETR0=1.D0-2.D0*ETA*UNR
      U0=G0
      U1=-X*GP0
      U2=-0.5D0*ETR0*X2*U0
      S=U0+U1+U2
      V1=U1/X
      V2=2.D0*U2/X
      T=V1+V2
      XN=3.D0
      DO10N=3,10000
C     N=N
      XN1=XN-1.D0
      XN1=XN*XN1
      U3=X*U2*UNR*(1.D0-2.D0/XN)-ETR0*U1*X2/XN1+X3*U0*UNR/XN1
      S=S+U3
      V3=XN*U3/X
      T=T+V3
   16 IF(DABS(U3/S).GT.1.D-11)GO TO 11
      IF(DABS(V3/T).LE.1.D-11)GO TO 20
   11 U0=U1
      U1=U2
      U2=U3
      XN=XN+1.D0
   10 CONTINUE
   20 G0=S
      GP0=-T
      RETURN
C   SERIE  ORIGINE
  200 CONTINUE
      PI=3.141592653589793D0
      GA=0.577215664901533D0
      ETA2=ETA*ETA
      RO2=RO*RO
      ETAP=ETA+ETA
      PIETA=PI*ETA
      PIETA2=0.5D0*PIETA
      B=DEXP(PIETA2)*DSQRT(DSINH(PIETA)/PIETA)
      CALL JFDELG(1.D0,ETA,PAR,PAI,NB)
      C1=ETAP*(GA+GA+DLOG(2.D0)-1.D0+PAR)
      U0=0.D0
      U1=RO
      V0=1.D0
      V1=C1*RO
      U=U0+U1
      V=V0+V1
      UP=1.D0
      VP=C1
      XN=2.D0
      DO 104N=2,10000
      XN1=XN*(XN-1.D0)
      U2=(ETAP*RO*U1-RO2*U0)/XN1
      U=U+U2
      V2=(ETAP*RO*V1-RO2*V0-ETAP*(XN+XN-1.D0)*U2)/XN1
      V=V+V2
      UP=UP+XN*U2/RO
      VP=VP+XN*V2/RO
   17 IF(DABS(U2/U).GT.1.D-14)GOTO18
      IF(DABS(V2/V).LE.1.D-14)GOTO19
   18 U0=U1
      U1=U2
      V0=V1
      V1=V2
      XN=XN+1.D0
  104 CONTINUE
   19 GP=V+ETAP*U*DLOG(RO)
      G0=B*GP
      GP0=B*(VP+ETAP*(UP*DLOG(RO)+U/RO))
      RETURN
  250 IF(RO.LE.0.5D0*ETA+9.D0)GOTO200
      GOTO 300
  251 IF(RO.LE.2.25D0+7.35D0*(3.D0-ETA))GOTO200
      GOTO 300
  252 IF(RO.LE.1.2D0+1.5D-1*(1.D1-ETA))GOTO200
      GOTO 300
  253 IF(RO.LE.0.6D0+0.75D-1*(18.D0-ETA))GOTO200
      GOTO 300
  254 IF(RO.LE.0.4D0+0.5D-1*(22.D0-ETA))GOTO200
      GOTO 300
      END

      SUBROUTINE YFRICA(ETA,RAU,FO,FPO,GO,GPO,SIGMA0,IDIV)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Q(5),QP(5)
C
C        COEFFICIENTS RICCATI
C
      DATA G61,G62,G63,G64,G65,G66,G67,G68,G69,G610,
     &G611/ 0.1159057617187498D-01,0.3863525390624998D-01,
     &0.4660034179687498D-01,0.4858398437499998D-01,
     &0.1156514485677080D 01,0.5687475585937496D 01,
     &0.1323888288225445D 02,0.1713083224826384D 02,
     &0.1269003295898436D 02,0.5055236816406248D 01,
     &0.8425394694010415D 00/
      DATA G81,G82,G83,G84,G85,G86,G87,G88,G89,G810,G811,G812,G813,G814,
     &G815/ 0.1851092066083633D-01,0.8638429641723630D-01,
     &0.1564757823944092D 00,0.1430139541625977D 00,
     &0.1924622058868408D 00,0.8500803152720129D 01,
     &0.7265429720878595D 02,0.3057942376817972D 03,
     &0.7699689544836672D 03,0.1254157054424285D 04,
     &0.1361719536066055D 04,0.9831831171035763D 03,
     &0.4547869927883148D 03,0.1222640538215636D 03,
     &0.1455524450256709D 02/
      DATA GP61,GP62,GP63,GP64,GP65,GP66/ 0.2897644042968748D-01,
     &0.2318115234375000D 00,0.8056640625000000D 00,
     &0.1601562499999998D 01,0.3046875000000000D 00,
     &0.5624999999999998D 01/
      DATA GP81,GP82,GP83,GP84,GP85,GP86,GP87,
     &GP88/ 0.6478822231292720D-01,0.6910743713378906D 00,
     &0.3322952270507811D 01,0.9483032226562498D 01,
     &0.1769653320312499D 02,0.3478710937499998D 02,
     &0.5020312499999999D 02,0.7874999999999999D 02/
      DATA Q /0.4959570165D-1,0.8888888889D-2,0.2455199181D-2,
     &0.9108958061D-3,0.2534684115D-3/
      DATA QP /0.1728260369D0,0.3174603174D-3,0.3581214850D-2,
     &0.3117824680D-3,0.9073966427D-3/
      CALL JFLGAM(1.D0,ETA,TRA,SIGMA0,IND)
      RAU2=RAU+RAU
      RAUC=RAU*RAU
      ETAC=ETA*ETA
      ETA2=ETA+ETA
      ETARO=ETA*RAU
      ETARO2=ETARO+ETARO
      PIETA=3.141592653589793*ETA
      IND=0
      JND=0
      IG=0
      IF(ETA)20,20,21
   20 IF(-ETARO-14.0625D0)32,22,22
   22 INDICE=1
 
C             RICCATI 3
      IDIV=0
      GO TO 2
   21 IF(DABS(RAU-ETA2).LE.1.D-9)GO TO 18
      IF(RAU-ETA2)30,18,31
   31 IF(RAU-ETA2-2.D1*(ETA**0.25D0))34,33,33
   33 INDICE=0
C             RICCATI  2
      IDIV=0
      GO TO 2
   32 NN=1
      GO TO 35
   34 NN=0
   35 CALL YFCLEN(ETA,RAU,FO,FPO,GO,GPO,SIGMA0,IDIV,NN)
      RETURN
   30 IF(ETARO-12.D0)32,32,23
   23 TRA=ETA2-6.75D0*(ETA**0.4D0)
      IF(RAU-TRA)6,6,24
   24 IND=1
      JND=1
      RO=RAU
      RAU=TRA
      RAU0=TRA
C             RICCATI  1
 
    6 X=RAU/ETA2
      U= (1.D0-X)/X
      X2=X*X
      RU= DSQRT(U)
      RX= DSQRT(X)
      TRE= 1.D0/(U*RU*ETA2)
      TRB=TRE*TRE
      FI= (DSQRT((1.D0-X)*X)+DASIN(RX)-1.570796326794897D0)*ETA2
      TR1= -0.25D0*DLOG(U)
  602 TR2= -((9.D0*U+6.D0)*U+5.D0)/48.D0
  603 TR3= ((((-3.D0*U-4.D0)*U+6.D0)*U+12.D0)*U+5.D0)/64.D0
  604 TR4=- ((((((U+2.D0)*945.D0*U+1395.D0)*U+12300.D0)*U+25191.D0)
     &*U+19890.D0)*U+5525.D0)/46080.D0
  605 TR5= ((((((((-27.D0*U-72.D0)*U-68.D0)*U+360.D0)*U+2190.D0)
     &*U+4808.D0)*U+5148.D0)*U+2712.D0)*U+565.D0)/2048.D0
  606 TR6=- (((((((((G61*U+G62)*U+G63)*U+G64)*U+G65)*U+G66)*U+G67)
     &*U+G68)*U+G69)*U+G610)*U+G611
  607 TR7= ((((((((((((-81.*U-324.)*U-486.)*U-404.)*U+4509.)*U+52344.)
     &*U+233436.)*U+567864.)*U+838521.)*U+775884.)*U+441450.)
     &*U+141660.)*U+19675.) /6144.
  608 TR8= (((((((((((((G81*U+G82)*U+G83)*U+G84)*U+G85)*U+G86)*U+G87)
     &*U+G88)*U+G89)*U+G810)*U+G811)*U+G812)*U+G813)*U+G814)*U+G815
      PSIP=PSIP+TRA
      XXX=138.1551055796428D0
      FI= FI+TRE*(TR2+TRB*(TR4+TRB*(TR6+TRB*TR8)))
      PSI=-FI
      INDG=IDINT(PSI/XXX)
      IDIV=60*INDG
      TRA= TR1+TRB*(TR3+TRB*(TR5+TRB*TR7))
      FI=FI+TRA
      PSI=PSI+TRA
 
      FIP=RU*ETA2
      TRA=1.D0/(X2*U)
      TR1=0.25D0
      TRE= TRE/(X2*X2*U)
      TRB=TRB/(X2*X2)
  702 TR2= -(8.D0*X-3.D0)/32.D0
  703 TR3= ((24.D0*X-12.D0)*X+3.D0)/64.D0
  704 TR4= (((-1536.D0*X+704.D0)*X-336.D0)*X+63.D0)/2048.D0
  705 TR5= ((((1920.D0*X-576.D0)*X+504.D0)*X-180.D0)*X+27.D0)/1024.D0
  706 TR6 = ((((-GP66*X+GP65)*X-GP64)*X+GP63)*X-GP62)*X+GP61
  707 TR7= - ((((((-40320.D0*X-10560.D0)*X-13248.D0)*X+7560.D0)
     &*X-3132.D0)*X+756.D0)*X-81.D0) / 2048.D0
  708 TR8 =- ((((((GP88*X+GP87)*X+GP86)*X-GP85)*X+GP84)*X-GP83)*X+GP82)
     &*X-GP81
      FIP=FIP +TRE*(TR2+TRB*(TR4+TRB*(TR6+TRB*TR8)))
      TRA= TRA*(TR1+TRB*(TR3+TRB*(TR5+TRB*TR7)))
      FIP=FIP+TRA
      PSIP=-FIP
      IF(INDG.EQ.0)GO TO 8
      PSI=PSI-XXX*DFLOAT(INDG)
      FI =FI +XXX*DFLOAT(INDG)
    8 FO=0.5D0*DEXP(FI)
      GO= DEXP(PSI)
      FPO= FO*FIP/ETA2
      GPO=GO*PSIP/ETA2
      IF(JND.EQ.0)RETURN
      RAU=RO
      GO=FO
      GPO=FPO
   27 X=RAU0-RO
      X2=X*X
      X3=X*X2
      UNR=1.D0/RAU0
      ETR0=1.D0-2.D0*ETA*UNR
      U0=GO
      U1=-X*GPO
      U2=-0.5D0*ETR0*X2*U0
      S=U0+U1+U2
      V1=U1/X
      V2=2.D0*U2/X
      T=V1+V2
      XN=3.D0
 
      DO10N=3,10000
C     N=N
      XN1=XN-1.D0
      XN1=XN*XN1
      U3=X*U2*UNR*(1.D0-2.D0/XN)-ETR0*U1*X2/XN1+X3*U0*UNR/XN1
      S=S+U3
      V3=XN*U3/X
      T=T+V3
   16 IF(DABS(U3/S).GT.1.D-10)GO TO 11
      IF(DABS(V3/T).LE.1.D-10)GO TO 17
   11 U0=U1
      U1=U2
      U2=U3
      XN=XN+1.D0
   10 CONTINUE
   17 IF(IG)25,26,25
   25 GO=S
      GPO=-T
      FO=HO
      FPO=HPO
      RETURN
   26 HO=S
      HPO=-T
   18 ET0=ETA**(0.166666666666667)
      ETAD=ETAC*ETAC
      ET=ETA**(0.6666666666666667)
      ET1=ET*ET
      ET2=ET1*ET1
      ET3=ET2*ET
      ET4=ETAD*ET
      ET5=ET4*ET
      FO=1.D0-Q(1)/ET1-Q(2)/ETAC-Q(3)/ET3-Q(4)/ETAD-Q(5)/ET5
      GO=1.D0+Q(1)/ET1-Q(2)/ETAC+Q(3)/ET3-Q(4)/ETAD+Q(5)/ET5
      FPO=1.D0+QP(1)/ET+QP(2)/ETAC+QP(3)/ET2+QP(4)/ETAD+QP(5)/ET4
      GPO=1.D0-QP(1)/ET+QP(2)/ETAC-QP(3)/ET2+QP(4)/ETAD-QP(5)/ET4
      FO=0.7063326373D0*ET0*FO
      GO=1.223404016D0*ET0*GO
      FPO=0.4086957323D0*FPO/ET0
      GPO=-0.7078817734D0*GPO/ET0
      IDIV=0
      IF(IND.EQ.0)RETURN
      IG=1
      RAU0=ETA2
      GO TO 27
    2 X=ETA2/RAU
      X2=X*X
      U=1.D0-X
      RU= DSQRT(U)
      U3=U*U*U
      TRD= 1.D0/(U3*ETA2*ETA2)
      TRC=X2*TRD
      TRE=1.D0/(U*RU*ETA2)
      FI= -0.25D0*DLOG(U)
      TRB=TRD/64.D0
      TR3= (((3.D0*U-4.D0)*U-6.D0)*U+12.D0)*U-5.D0
  501 TR5= ((((((((-27.D0*U+72.D0)*U-68.D0)*U-360.D0)*U+2190.D0)
     &*U-4808.D0)*U+5148.D0)*U-2712.D0)*U+565.D0)/32.D0
  502 TR7= ((((((((((((81.D0*U-324.D0)*U+486.D0)*U-404.D0)*U-4509.D0)
     &*U+52344.D0)*U-233436.D0)*U+567864.D0)*U-838521.D0)*U+775884.D0)
     &*U-441450.D0)*U+141660.D0)*U-19675.D0)/96.D0
      FI= FI+TRB*(TR3+TRD*(TR5+TRD*TR7))
 
      FIP=0.25D0/U
      TRB=3.D0*TRC/(64.D0*U)
      TR3= (X-4.D0)*X+8.D0
  511 TR5= ((((9.D0*X-60.D0)*X+168.D0)*X-192.D0)*X+640.D0)/16.D0
  512 TR7= ((((((-27.D0*X+252.D0)*X-1044.D0)*X+2520.D0)*X-4416.D0)
     &*X-3520.D0)*X-13440.D0)/32.D0
      FIP =FIP+TRB*(TR3+TRC*(TR5+TRC*TR7))
      TRA= DABS((RU-1.D0)/(RU+1.D0))
      PSI= (0.5D0*DLOG(TRA)+RU/X)*ETA2+0.785398163397448D0
      TR2= -((9.D0*U-6.D0)*U+5.D0)/48.D0
  521 TR4= ((((((U-2.D0)*945.D0*U+1395.D0)*U-12300.D0)*U+25191.D0)
     &*U-19890.D0)*U+5525.D0)/46080.D0
  522 TR6 = (((((((((-G61*U+G62)*U-G63)*U+G64)*U-G65)*U+G66)*U-G67)
     &*U+G68)*U-G69)*U+G610)*U-G611
  523 TR8= (((((((((((((G81*U-G82)*U+G83)*U-G84)*U+G85)*U-G86)*U+G87)
     &*U-G88)*U+G89)*U-G810)*U+G811)*U-G812)*U+G813)*U-G814)*U+G815
      PSI= PSI+TRE*(TR2+TRD*(TR4+TRD*(TR6+TRD*TR8)))
      PSIP = -RU*ETA2/X2
      TRB=TRE*X/U
      TR2= (3.D0*X-8.D0)/32.D0
  531 TR4= - (((63.D0*X-336.D0)*X+704.D0)*X-1536.D0)/2048.D0
  532 TR6 = ((((GP61*X-GP62)*X+GP63)*X-GP64)*X+GP65)*X-GP66
  533 TR8 = ((((((-GP81*X+GP82)*X-GP83)*X+GP84)*X-GP85)*X+GP86)*X+GP87)
     &*X+GP88
      PSIP =PSIP+ TRB*(TR2+TRC*(TR4+TRC*(TR6+TRC*TR8)))
      TRA= DEXP(FI)
      FO= TRA*DSIN(PSI)
      GO= TRA*DCOS(PSI)
      IF(INDICE)535,536,535
  535 TRA=FO
      FO=GO
      GO=-TRA
  536 TRA=-ETA2/RAUC
      FPO=(FIP*FO+PSIP*GO)*TRA
      GPO=(FIP*GO-PSIP*FO)*TRA
      RETURN
      END

      SUBROUTINE DFCZ0(ETA,RO,F0,FP0,G0,GP0,SIGMA0,IEXP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A1(110),A2(110),B1(110),B2(110)
      IF(RO)2,2,1
    2 WRITE (6,1000)
 1000 FORMAT (21H RO NEGATIF OU NUL **)
      RETURN
    1 IF(ETA-30.D0)3,3,4
    4 IF(DABS(ETA)-5.D2)28,28,29
   28 CALL YFRICA(ETA,RO,F0,FP0,G0,GP0,SIGMA0,IEXP)
      RETURN
   29 WRITE (6,1001)
 1001 FORMAT (42H VALEUR ABSOLUE DE ETA SUPE-&EU-E A 500 **)
      RETURN
    3 IF(ETA+8.D0)4,5,5
    5 IF(ETA)6,7,6
    7 F0=DSIN(RO)
      G0=DCOS(RO)
      FP0=G0
      GP0=-F0
      IEXP=0
      SIGMA0=0.D0
      RETURN
    6 BORNE=1.666666666666667D0*DABS(ETA)+7.5D0
      IF(RO-BORNE)8,9,9
    9 CALL YFASYM(ETA,RO,F0,FP0,G0,GP0,SIGMA0,IEXP)
      RETURN
    8 IF(ETA-10.D0)10,11,11
   10 IF(ETA)12,12,13
   13 IF(RO-2.D0)14,12,12
   11 IF(ETA-(5.D0*RO+6.D1)/7.D0)12,12,14
   12 CALL YFASYM(ETA,BORNE,F0,FP0,G0,GP0,SIGMA0,IEXP)
      H=BORNE
      DH=F0/H
      IF(ETA)20,21,21
   20 N=-0.5D0*ETA+5.D0
      GO TO 22
   21 N=ETA/5.D0+5.D0
   22 N=5*(N+1)
      Z=4.D0/H
      Y=1.D0-(ETA+ETA)*Z
      A1(N+2)=1.D-55
      A1(N+3)=0.D0
      A1(N+4)=1.D-64
      B1(N+3)=1.D-50
      B1(N+4)=1.D-68
      A2(N+2)=0.D0
      A2(N+3)=1.D-74
      A2(N+4)=1.D-53
      B2(N+3)=0.D0
      B2(N+4)=1.D-66
      M=N+2
      DI=N
      DO 23 II=2,M
      I=M-II+2
      B1(I)=B1(I+2)+Z*(DI+1.D0)*A1(I+1)
      S=A1(I+2)+Y*(A1(I+1)-A1(I))
      Q=(DI+2.D0)*B1(I)+(DI-1.D0)*B1(I+1)
      A1(I-1)=S-Z*Q
      B2(I)=B2(I+2)+Z*(DI+1.D0)*A2(I+1)
      S=A2(I+2)+Y*(A2(I+1)-A2(I))
      Q=(DI+2.D0)*B2(I)+(DI-1.D0)*B2(I+1)
      A2(I-1)=S-Z*Q
      IF(I.GE.N)GO TO 23
      D=-(B2(I+2)+B2(I))/(B1(I+2)+B1(I))
      DO 24 J=I,M
      A2(J)=A2(J)+D*A1(J)
      B2(J)=B2(J)+D*B1(J)
   24 CONTINUE
      A2(I-1)=A2(I-1)+D*A1(I-1)
   23 DI=DI-1.D0
      Q=A1(3)-A1(1)
      C=A2(3)-A2(1)
      C=Q/C
      X1=0.5D0*(A1(2)-C*A2(2))
      DO 25 I=3,M
      X1=X1+A1(I)-C*A2(I)
   25 CONTINUE
      X1=DH/X1
      X2=-C*X1
      DO 26 I=2,M
      B1(I)=X1*B1(I)+X2*B2(I)
      A1(I)=X1*A1(I)+X2*A2(I)
   26 CONTINUE
      A1(1)=X1*A1(1)+X2*A2(1)
      B1(1)=0.D0
      X=RO/H
      Y=2.D0*X-1.D0
      T1=1.D0
      T2=Y
      RESULT=0.5D0*A1(2)+Y*A1(3)
      DERIVE=0.5D0*B1(2)+Y*B1(3)
      DO 27 I=2,N
      TI=2.D0*Y*T2-T1
      T1=T2
      T2=TI
      RESULT=RESULT+TI*A1(I+2)
      DERIVE=DERIVE+TI*B1(I+2)
   27 CONTINUE
      F0=RESULT*RO
      FP0=DERIVE*RO+RESULT
      GO TO 30
C   SERIE ORIGINE REGULIERE
   14 PI=3.141592653589793D0
      CALL JFLGAM(1.D0,ETA,TRA,SIGMA0,NTRUC)
      IEXP=0
      RO2=RO*RO
      ETAP=ETA+ETA
      PIETA=PI*ETA
      PIETA2=0.5D0*PIETA
      B=DEXP(PIETA2)*DSQRT(DSINH(PIETA)/PIETA)
      U0=0.D0
      U1=RO
      U=U0+U1
      UP=1.D0
      XN=2.D0
      DO 15 N=2,10000
      XN1=XN*(XN-1.D0)
      U2=(ETAP*RO*U1-RO2*U0)/XN1
      U=U+U2
      UP=UP+XN*U2/RO
   17 IF(DABS(U2/U).LT.1.D-10)GO TO 19
   18 U0=U1
      U1=U2
      XN=XN+1.D0
   15 CONTINUE
   19 F0=U/B
      FP0=UP/B
   30 CALL YFIREG(ETA,RO,G0,GP0)
      RETURN
      END

      SUBROUTINE JFDELG (XD,YD,PAR,PAI,NBCHIF)
      DOUBLE PRECISION XD,YD,PAR,PAI,TEST,C,PI
      DOUBLE PRECISION X,Y,U,V,TRA,TRA1,COSI,SINI
      DOUBLE PRECISION COS2I,SIN2I,ZMOD,DEPI
      DOUBLE PRECISION TRB,XX
      DOUBLE PRECISION RAC2,PIS4
      DOUBLE PRECISION SUPINT
      DIMENSION TEST(7),C(6)
      DATA TEST/2.9152D+7,2.2958D+3,1.4124D+2,3.9522D+1,19.6611D0,
     &12.791D0,-10.D0/
      DATA RAC2/0.3465735902799726/,PIS4/0.785398163397448/
      DATA C/8.333333333333333D-2,-8.33333333333333D-3,
     &3.968253968253968D-3,-4.166666666666667D-3,7.575757575757576D-3,
     &-2.109279609279609D-2/
      DATA PI/3.141592653589793/
      DATA DEPI/6.283185307179586/
      DATA SUPINT/2147483647.D0/
      X=DABS(XD)
      XX=X
      NBCHIF=15
      IF(YD)1,2,1
    1 Y=DABS(YD)
      KR=1
      I=DMOD(10.99D0-X,SUPINT)
C     TRANSLATION
      IF(I)3,3,4
    4 TRA=I
      X=X+TRA
C     LOGARITHME(X+IY) (X,Y POSITIFS)
    3 IF(X-Y)5,6,7
    5 TRA1=X/Y
      TRB=1.D0+TRA1*TRA1
      TRA=Y*DSQRT(TRB)
      SINI=1./(TRB*Y)
      COSI=SINI*TRA1
      TRA1=Y/X
      GO TO 11
    6 U=RAC2+DLOG(X)
      V=PIS4
      SINI=0.5D0/X
      COSI=SINI
      GO TO 10
    7 TRA1=Y/X
      TRB=1.D0+TRA1*TRA1
      TRA=X*DSQRT(TRB)
      COSI=1./(TRB*X)
      SINI=COSI*TRA1
   11 U=DLOG(TRA)
      V=DATAN(TRA1)
C     DEVELOPPEMENT ASYMPTOTIQUE ( X SUPERIEUR A 10 )
   10 PAR=U-0.5*COSI
      PAI=V+0.5*SINI
      ZMOD=X+Y
      IF(ZMOD-TEST(1))13,13,14
   13 SIN2I=(SINI*COSI)+(SINI*COSI)
      COS2I=(COSI+SINI)*(COSI-SINI)
      SINI=SIN2I
      COSI=COS2I
      K=1
      GO TO 15
   16 TRA=COSI*COS2I-SINI*SIN2I
      SINI=SINI*COS2I+COSI*SIN2I
      COSI=TRA
   15 PAR=PAR-C(K)*COSI
      PAI=PAI+C(K)*SINI
      K=K+1
      IF(ZMOD-TEST(K))16,16,14
C     TRANSLATION INVERSE
   17 I=I-1
      X=I
      X=XX+X
   56 IF(X-Y)55,55,57
   55 TRA1=X/Y
      TRB=X*TRA1+Y
      SINI=1.D0/TRB
      COSI=TRA1/TRB
      GO TO 19
   57 TRA1=Y/X
      TRB=X+Y*TRA1
      COSI=1.D0/TRB
      SINI=TRA1/TRB
   19 PAR=PAR-COSI
      PAI=PAI+SINI
   14 IF(I)18,18,17
 
C     CONTROLE DU QUADRANT
   18 IF(XD)20,61,21
   61 TRA=PI*Y
      IF(TRA-1.D-2)300,300,301
  300 TRB= TRA*(2.D0+TRA*(-2.D0+TRA*(1.333333333333333D0+TRA*(
     &-0.6666666666666666D0+TRA*(0.2666666666666666D0+TRA*(
     &-0.08888888888888888D0+TRA*0.02539682539682540D0))))))
      TRB=(2.D0-TRB)/TRB
      GO TO 302
  301 TRB= DEXP(-TRA-TRA)
      TRB=(1.D0+TRB)/(1.D0-TRB)
  302 PAI=0.5D0*(1.D0/Y+PI*TRB)
   21 IF(YD)28,100,100
C     X+IY CHANGE EN -X-IY
   20 TRA=DEXP(-DEPI*Y)
      TRB=TRA*TRA
      COS2I=DEPI*DMOD(X,1.D0)
      SIN2I=-2.D0*TRA*DCOS(COS2I)+1.D0+TRB
      PAR=PAR+COSI+DEPI*TRA*DSIN(COS2I)/SIN2I
      PAI=PAI-SINI+PI*(TRB-1.D0)/SIN2I
      IF(YD)100,100,28
   28 PAI=-PAI
   35 WRITE (6,1001)
 1001 FORMAT (30H ARGUMENT DE JFDELG TROP GRAND)
      GO TO 50
   34 Y=DMOD(X,SUPINT)
      IF(Y) 400,36,400
  400 IF(Y-0.99D0) 33,33,405
  405 TRA= IDINT(Y+0.1D0)
      IF(DABS(Y-TRA)-5.D-15)36,36,33
   31 IF(X-4503599627370496.D0)34,35,35
 
C     ARGUMENT DANS -PI,PI
  100 TRA=DABS(PAI/DEPI)
      IF(TRA-1.D+15)203,204,204
  204 NBCHIF=0
      PAI=0.D0
      GO TO 29
  203 IF(TRA-1.D0)205,205,206
  206 NBCHIF=DLOG10(TRA)
      NBCHIF=14-NBCHIF
      TRA=DMOD(TRA,SUPINT)
      PAI=DMOD(TRA,1.D0)*DSIGN(DEPI,PAI)
  205 IF(DABS(PAI)-PI)29,29,207
  207 PAI=PAI-DSIGN(DEPI,PAI)
      GO TO 29
C        DELGAMMA REEL
    2 PAI=0.D0
      IF(XD)31,32,33
C     CONDITIONS D EXISTENCE
   32 WRITE (6,1000)
 1000 FORMAT (21H JFDELG(0) EST INFINI)
      GO TO 50
   36 WRITE (6,1002)
 1002 FORMAT (28H JFDELG (-ENTIER) EST INFINI)
   50 PAR=1.D+74
      NBCHIF=0
      GO TO 29
C     TRANSLATION
   33 I=DMOD(10.99D0-X,SUPINT)
      IF(I)37,37,38
   38 TRA=I
      X=X+TRA
C     DEVELOPPEMENT ASYMPTOTIQUE
   37 Y=DLOG(X)
      PAR=Y-0.5D0/X
      IF(X-TEST(1))39,39,43
   39 COS2I=1.D0/(X*X)
      COSI=COS2I
      K=1
      GO TO 41
   42 COSI=COSI*COS2I
   41 PAR=PAR-C(K)*COSI
      K=K+1
      IF(X-TEST(K))42,42,40
C     TRANSLATION INVERSE
   44 I=I-1
      X=I
      X=XX+X
      PAR=PAR-1.D0/X
   40 IF(I)43,43,44
C     X NEGATIF
   43 IF(XD)45,29,29
   45 PAR=PAR+1.D0/X
      Y=PI*DMOD(X,2.D0)
      PAR=PAR+PI*DCOS(Y)/DSIN(Y)
      ENTRY JFDEV1
   29 RETURN
      END



