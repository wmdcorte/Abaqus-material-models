c CDP_DCM : Concrete Damaged Plasticity (CDP) model De Corte/Miarka 
c -- UMAT for ABAQUS/Standard
c Developed by Wouter De Corte (Ghent University) and Petr Miarka 
c (Institute of Physics, CAS - Brno)
c 
c The code may be freely used in academic courses and educational 
c settings.
c  
c If you use this subroutine, or code derived from it, in work that 
c you publish, please cite the paper in which it is described and
c verified:
c
c W. De Corte, P. Miarka, "A user-material subroutine for the Abaqus 
c concrete damaged plasticity model: implementation and systematic 
c verifications", Appl. Sci. 2026, 16(21), 10798;
c https://doi.org/10.3390/app162110798 
c
c DISCLAIMER
c
c This is research software, provided "as is" without warranty of any
c kind. It was verified against the native Abaqus CDP model only on
c the cases reported in the paper cited above. The authors accept no
c responsibility for the correctness of results, including numerical
c errors, convergence failures or use outside the verified range.
c Verifying the subroutine for the user's own parameters, elements
c and loading, and any use of the results, is the user's
c responsibility. 
c
c This is a from-scratch UMAT reimplementation of the yield surface,
c flow rule and damage/hardening laws that ABAQUS documents for its
c BUILT-IN "*CONCRETE DAMAGED PLASTICITY" material (Abaqus Analysis
c User's Guide, "Concrete damaged plasticity"; Abaqus Theory Guide).
c That native model is not exposed to users as callable/portable
c source code, so this file re-derives the algorithm from the
c published equations. It is NOT a copy of SIMULIA's internal code.
c
c The model is the Lubliner/Lee-Fenves plastic-damage formulation:
c   - Yield surface: Lubliner et al. (1989), modified by Lee & Fenves
c     (1998), function of the effective stress invariants p (mean),
c     q (Mises) and the algebraically largest principal effective
c     stress, with two independent hardening/softening variables
c     (tension, compression).
c   - Flow potential: non-associated Drucker-Prager hyperbolic
c     potential in the (p,q) plane, dilation angle psi, eccentricity
c     ecc.
c   - Damage: two independent scalar damage variables dt (tension),
c     dc (compression), applied through ABAQUS's multiplicative
c     stiffness-recovery combination with tension/compression
c     recovery factors wt, wc.
c   - The return mapping algorithm is inspired by the CDPM2 model 
c     implemented by Seungwook Seok available as VUMAT on Github. 
c     The repository contains the Abaqus/Explicit user-material 
c     (VUMAT) of the concrete damage-plasticity model 2 (CDPM2).
c
c References:
c 1) J. Lubliner, J. Oliver, S. Oller, E. Onate: "A plastic-damage
c    model for concrete", Int. J. Solids Structures 25(3), 1989.
c 2) J. Lee, G.L. Fenves: "Plastic-Damage Model for Cyclic Loading of
c    Concrete Structures", J. Eng. Mech. 124(8), 1998.
c 3) J. Lee, G.L. Fenves, "A return-mapping algorithm for plastic-damage 
c    models: 3-D and plane stress formulation". International Journal for
c    Numerical Methods in Engineering. Volume 50, pp. 487-506, 2001.
c 4) P. Grassl, D. Xenos, U. Nyström, R. Rempling, K. Gylltoft.: 
c    "CDPM2: A damage-plasticity approach to modelling the failure of 
c    concrete". International Journal of Solids and Structures. 
c    Volume 50, Issue 24, pp. 3805-3816, 2013.
c 5) seungwookseok/ABAQUS-version-CDPM2 [Computer Software]. (2019)
c    Retrieved from https://github.com/seungwookseok/ABAQUS-version-CDPM2
c 6) Abaqus Analysis User's Guide, section "Concrete damaged
c    plasticity", and Abaqus Theory Guide (SIMULIA/Dassault Systemes).
c ----------------------------------------------------------------------
c User-defined material properties (props)
c
c props(1)  E        props(2)  nu       props(3)  psi [deg]
c props(4)  ecc (default 0.1)           props(5)  fb0/fc0 (default 1.16)
c props(6)  Kc (default 2/3)            props(7)  mu (viscosity, >=0)
c props(8)  wt (default 0)              props(9)  wc (default 1)
c props(10) itype : 0 = tension tables in cracking STRAIN,
c                   1 = in crack-opening DISPLACEMENT (divided by CELENT)
c props(11) ntt   rows tension stiffening   (>=1)
c props(12) ntc   rows compression hardening (>=1)
c props(13) ntdt  rows tension damage        (0 = no tension damage)
c props(14) ntdc  rows compression damage    (0 = no compr. damage)
c then the four 2-column tables in this order (sigma, eck) /
c (sigma, ein) / (dt, eck) / (dc, ein).
c
c OPTIONAL extra constant after the tables (only read if present):
c props(15+2*(ntt+ntc+ntdt+ntdc)) itang : 0 = analytic tangent with
c           finite-difference fallback (default), 1 = always finite
c           differences.
c
c Declare  *USER MATERIAL, CONSTANTS=<NPROPS>, UNSYMM   and Depvar=27.
c
c State variables:
c  1-2   kappa_t, kappa_c          3-8   plastic strain
c  9-10  d_t, d_c                  11-16 total strain
c  17-18 viscous d_t, d_c          19    combined degradation d
c  20    stress-weight factor r    21-26 effective stress
c  27    number of sub-steps used in the last increment (1 = none)
c For plane stress statev(13) holds the total out-of-plane strain
c and statev(5) the out-of-plane plastic strain (internal ordering).
c
c Internal ordering (11,22,33,12,23,13), engineering shear strains;
c Abaqus 3-D order (11,22,33,12,13,23) -> components 5,6 swapped;
c plane strain/axisymmetric (11,22,33,12); plane stress (11,22,12).
c ======================================================================
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     & rpl,ddsddt,drplde,drpldt,
     & stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     & ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     & celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)

      include 'aba_param.inc'

      character*80 cmname
      dimension stress(ntens),statev(nstatv),
     & ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     & stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     & props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      integer MAXTAB
      parameter (MAXTAB=400)
      real*8 cst(20),tab(MAXTAB,10)
      integer ntab(4)

      real*8 deps(6),depsP(6),stateO(27),stateN(27),stateP(27)
      real*8 sigeff(6),sigeffP(6),sigeffM(6)
      real*8 sig(6),sigP(6),sigM(6),sigUM(6),sigUMP(6),sigUMM(6)
      real*8 evtr(3),evec(3,3),xsol(2),evtrP(3),evecP(3,3),xsolP(2)
      real*8 dtv,dcv,dtvP,dcvP,dtvM,dcvM,dtvold,dcvold,vfac,vfacP
      real*8 elmat(6,6),dd6(6,6),dsmax,pert,dcomb,epsel(6),depl(6)
      real*8 e33,s33,d33,tol33,e33old,s33old,scal
      integer imap(6),ipFail,ipFailP,ipFailM,k1,k2,jj,nsub,nsubP
      integer iplas,iplasP,itang,ierr,iana,ips,it33,ii,ij
      integer ifirst
      save ifirst
      data ifirst /0/

c ---------------------------------------------------------------------
c Material constants and tables (local, per call -> thread safe)
c ---------------------------------------------------------------------
      call cdp_build(props,nprops,celent,cst,tab,ntab,itang,ierr)
      if (ierr .ne. 0) call xit

      if (ifirst .eq. 0) then
         ifirst = 1
         call cdp_print(cst,ntab,itang)
      end if

      if (ntens .lt. 3) then
         write(*,*) 'CDP UMAT: NTENS < 3 not supported.'
         call xit
      end if
      ips = 0
      if (ntens .eq. 3) ips = 1

c ---------------------------------------------------------------------
c Map strain increment to internal ordering
c ---------------------------------------------------------------------
      do k1 = 1,6
         deps(k1) = 0.d0
         imap(k1) = 0
      end do
      if (ips .eq. 1) then
         deps(1) = dstran(1)
         deps(2) = dstran(2)
         deps(4) = dstran(3)
         imap(1) = 1
         imap(2) = 2
         imap(4) = 3
      else
         deps(1) = dstran(1)
         deps(2) = dstran(2)
         deps(3) = dstran(3)
         deps(4) = dstran(4)
         imap(1) = 1
         imap(2) = 2
         imap(3) = 3
         imap(4) = 4
         if (ntens .eq. 6) then
            deps(5) = dstran(6)
            deps(6) = dstran(5)
            imap(5) = 6
            imap(6) = 5
         end if
      end if

      do k1 = 1,27
         stateO(k1) = 0.d0
      end do
      do k1 = 1,min(nstatv,27)
         stateO(k1) = statev(k1)
      end do
      dtvold = stateO(17)
      dcvold = stateO(18)

c ---------------------------------------------------------------------
c Stress update. Plane stress: Newton on deps(3) so that sigma_33 = 0.
c ---------------------------------------------------------------------
      if (ips .eq. 1) then
c        elastic estimate of the out-of-plane strain increment
         deps(3) = -cst(2)/(1.d0-cst(2))*(deps(1)+deps(2))
         tol33 = 1.d-9*max(cst(11),1.d-30)
         e33old = 0.d0
         s33old = 0.d0
         do it33 = 1,40
            call cdp_update(cst,tab,ntab,deps,stateO,stateN,sigeff,
     $           evtr,evec,xsol,iplas,nsub,ipFail)
            if (ipFail .eq. 1) goto 700
            call cdp_combine(cst,sigeff,stateN(9),stateN(10),
     $           stateN(20),dtvold,dcvold,dtime,sig,dtv,dcv,vfac)
            s33 = sig(3)
            scal = max(abs(sig(1)),abs(sig(2)),abs(sig(4)),cst(11))
            if (abs(s33) .le. 1.d-11*scal) goto 350
c           derivative d sigma_33 / d eps_33 from the analytic tangent
            iana = 0
            if (nsub .eq. 1 .and. itang .eq. 0) then
               call cdp_tangent(cst,tab,ntab,evtr,evec,xsol,iplas,
     $              stateO,stateN,sigeff,dtvold,dcvold,dtime,vfac,
     $              dd6,iana)
            end if
            if (iana .eq. 1) then
               d33 = dd6(3,3)
            else if (it33 .gt. 1 .and. abs(deps(3)-e33old)
     $              .gt. 1.d-300) then
               d33 = (s33-s33old)/(deps(3)-e33old)
            else
               d33 = (1.d0-stateN(19))*cst(1)*(1.d0-cst(2))/
     $              ((1.d0+cst(2))*(1.d0-2.d0*cst(2)))
            end if
            if (abs(d33) .lt. 1.d-12*cst(1)) d33 = 1.d-12*cst(1)
            e33old = deps(3)
            s33old = s33
            deps(3) = deps(3) - s33/d33
         end do
c        no convergence of the plane-stress iteration
         goto 700
 350     continue
      else
         call cdp_update(cst,tab,ntab,deps,stateO,stateN,sigeff,
     $        evtr,evec,xsol,iplas,nsub,ipFail)
         if (ipFail .eq. 1) goto 700
         call cdp_combine(cst,sigeff,stateN(9),stateN(10),
     $        stateN(20),dtvold,dcvold,dtime,sig,dtv,dcv,vfac)
      end if

      stateN(17) = dtv
      stateN(18) = dcv
      dcomb = 1.d0 - (1.d0-(1.d0-cst(9)*stateN(20))*dcv)*
     $     (1.d0-(1.d0-cst(10)*(1.d0-stateN(20)))*dtv)
      stateN(19) = dcomb
      stateN(27) = dble(nsub)

c ---------------------------------------------------------------------
c Energies (informational): SSE elastic, SPD effective plastic work
c ---------------------------------------------------------------------
      do k1 = 1,6
         epsel(k1) = stateN(10+k1) - stateN(2+k1)
         depl(k1) = stateN(2+k1) - stateO(2+k1)
      end do
      sse = 0.d0
      do k1 = 1,6
         sse = sse + 0.5d0*(1.d0-dcomb)*sigeff(k1)*epsel(k1)
         spd = spd + sigeff(k1)*depl(k1)
      end do

      do k1 = 1,min(nstatv,27)
         statev(k1) = stateN(k1)
      end do
      call cdp_mapsig(sig,sigUM,ips)

c ---------------------------------------------------------------------
c Consistent tangent: analytic (6x6 internal) unless sub-incremented
c or disabled; otherwise central finite differences.
c ---------------------------------------------------------------------
      iana = 0
      if (nsub .eq. 1 .and. itang .eq. 0) then
         call cdp_tangent(cst,tab,ntab,evtr,evec,xsol,iplas,
     $        stateO,stateN,sigeff,dtvold,dcvold,dtime,vfac,dd6,iana)
      end if
      if (iana .eq. 1) then
         if (ips .eq. 1) then
c           static condensation of eps_33 (internal index 3)
            do k1 = 1,6
               do k2 = 1,6
                  if (k1.ne.3 .and. k2.ne.3) then
                     dd6(k1,k2) = dd6(k1,k2) -
     $                    dd6(k1,3)*dd6(3,k2)/dd6(3,3)
                  end if
               end do
            end do
         end if
         do k1 = 1,6
            do k2 = 1,6
               if (imap(k1).ne.0 .and. imap(k2).ne.0) then
                  ddsdde(imap(k1),imap(k2)) = dd6(k1,k2)
               end if
            end do
         end do
         goto 900
      end if

c ---- finite-difference tangent (fallback) ----
      dsmax = 0.d0
      do k1 = 1,ntens
         if (abs(dstran(k1)) .gt. dsmax) dsmax = abs(dstran(k1))
      end do
      pert = 1.0d-4*cst(11)/cst(1)
      if (pert .lt. 1.d-9) pert = 1.d-9
      if (pert .gt. 1.d-2*dsmax .and. dsmax .gt. 0.d0) then
         if (1.d-2*dsmax .gt. 1.d-9) pert = 1.d-2*dsmax
      end if

      do jj = 1,ntens
         do k1 = 1,6
            depsP(k1) = deps(k1)
            if (imap(k1) .eq. jj) depsP(k1) = depsP(k1) + pert
         end do
         call cdp_fdupdate(cst,tab,ntab,depsP,stateO,ips,dtvold,
     $        dcvold,dtime,itang,sigP,ipFailP)
         if (ipFailP .eq. 1) goto 800
         do k1 = 1,6
            depsP(k1) = deps(k1)
            if (imap(k1) .eq. jj) depsP(k1) = depsP(k1) - pert
         end do
         call cdp_fdupdate(cst,tab,ntab,depsP,stateO,ips,dtvold,
     $        dcvold,dtime,itang,sigM,ipFailM)
         if (ipFailM .eq. 1) goto 800
         call cdp_mapsig(sigP,sigUMP,ips)
         call cdp_mapsig(sigM,sigUMM,ips)
         do k1 = 1,ntens
            ddsdde(k1,jj) = (sigUMP(k1)-sigUMM(k1))/(2.d0*pert)
         end do
      end do
      goto 900

c ---- update failed: ask for a cutback, return elastic stiffness ----
 700  continue
      pnewdt = 0.25d0
      call cdp_elasticmatrix(cst,elmat)
      call cdp_fillel(elmat,imap,ntens,ips,1.d0,ddsdde)
      return

c ---- secant fallback (perturbed update failed) ----
 800  continue
      call cdp_elasticmatrix(cst,elmat)
      dsmax = 1.d0 - dcomb
      if (dsmax .lt. 1.d-4) dsmax = 1.d-4
      call cdp_fillel(elmat,imap,ntens,ips,dsmax,ddsdde)

 900  continue
      do k1 = 1,ntens
         stress(k1) = sigUM(k1)
      end do
      return
      end

c ======================================================================
c Fill DDSDDE with fac * elastic matrix (condensed for plane stress)
c ======================================================================
      subroutine cdp_fillel(elmat,imap,ntens,ips,fac,ddsdde)
      integer imap(6),ntens,ips,k1,k2
      real*8 elmat(6,6),fac,ddsdde(ntens,ntens),e2(6,6)
      do k1 = 1,6
         do k2 = 1,6
            e2(k1,k2) = elmat(k1,k2)
         end do
      end do
      if (ips .eq. 1) then
         do k1 = 1,6
            do k2 = 1,6
               if (k1.ne.3 .and. k2.ne.3) then
                  e2(k1,k2) = e2(k1,k2) - e2(k1,3)*e2(3,k2)/e2(3,3)
               end if
            end do
         end do
      end if
      do k1 = 1,ntens
         do k2 = 1,ntens
            ddsdde(k1,k2) = 0.d0
         end do
      end do
      do k1 = 1,6
         do k2 = 1,6
            if (imap(k1).ne.0 .and. imap(k2).ne.0) then
               ddsdde(imap(k1),imap(k2)) = fac*e2(k1,k2)
            end if
         end do
      end do
      return
      end

c ======================================================================
c Full update (plasticity + damage + viscous combination), with the
c plane-stress iteration when ips=1. Used by the finite-difference
c tangent so that each perturbed column sees the same algorithm.
c ======================================================================
      subroutine cdp_fdupdate(cst,tab,ntab,deps,stateO,ips,dtvold,
     $     dcvold,dtime,itang,sig,ipFail)
      integer MAXTAB
      parameter (MAXTAB=400)
      real*8 cst(20),tab(MAXTAB,10)
      integer ntab(4),ips,itang,ipFail
      real*8 deps(6),stateO(27),dtvold,dcvold,dtime,sig(6)

      real*8 dp(6),stateN(27),sigeff(6),evtr(3),evec(3,3),xsol(2)
      real*8 dtv,dcv,vfac,dd6(6,6),s33,d33,e33old,s33old,scal
      integer k,iplas,nsub,iana,it33

      do k = 1,6
         dp(k) = deps(k)
      end do
      ipFail = 0
      if (ips .eq. 0) then
         call cdp_update(cst,tab,ntab,dp,stateO,stateN,sigeff,
     $        evtr,evec,xsol,iplas,nsub,ipFail)
         if (ipFail .eq. 1) return
         call cdp_combine(cst,sigeff,stateN(9),stateN(10),
     $        stateN(20),dtvold,dcvold,dtime,sig,dtv,dcv,vfac)
         return
      end if
      dp(3) = -cst(2)/(1.d0-cst(2))*(dp(1)+dp(2))
      e33old = 0.d0
      s33old = 0.d0
      do it33 = 1,40
         call cdp_update(cst,tab,ntab,dp,stateO,stateN,sigeff,
     $        evtr,evec,xsol,iplas,nsub,ipFail)
         if (ipFail .eq. 1) return
         call cdp_combine(cst,sigeff,stateN(9),stateN(10),
     $        stateN(20),dtvold,dcvold,dtime,sig,dtv,dcv,vfac)
         s33 = sig(3)
         scal = max(abs(sig(1)),abs(sig(2)),abs(sig(4)),cst(11))
         if (abs(s33) .le. 1.d-11*scal) return
         iana = 0
         if (nsub .eq. 1 .and. itang .eq. 0) then
            call cdp_tangent(cst,tab,ntab,evtr,evec,xsol,iplas,
     $           stateO,stateN,sigeff,dtvold,dcvold,dtime,vfac,
     $           dd6,iana)
         end if
         if (iana .eq. 1) then
            d33 = dd6(3,3)
         else if (it33 .gt. 1 .and. abs(dp(3)-e33old)
     $           .gt. 1.d-300) then
            d33 = (s33-s33old)/(dp(3)-e33old)
         else
            d33 = (1.d0-stateN(19))*cst(1)*(1.d0-cst(2))/
     $           ((1.d0+cst(2))*(1.d0-2.d0*cst(2)))
         end if
         if (abs(d33) .lt. 1.d-12*cst(1)) d33 = 1.d-12*cst(1)
         e33old = dp(3)
         s33old = s33
         dp(3) = dp(3) - s33/d33
      end do
      ipFail = 1
      return
      end

c ======================================================================
c Parse props() -> constants cst(), preprocessed tables tab(:,1..10),
c row counts ntab(4), tangent flag itang.
c tab columns: 1 e  2 sigma  3 d  4 kappa  at the merged nodes,
c   tension (ntab(1) rows); 5-8 the same for compression (ntab(2)).
c cst: 1 E 2 nu 3 psi 4 tanpsi 5 ecc 6 alpha 7 gamma 8 mu 9 wt 10 wc
c      11 sigt0 12 K 13 G 14 ndens 15 itype
c ======================================================================
      subroutine cdp_build(props,nprops,celent,cst,tab,ntab,itang,
     $     ierr)
      integer MAXTAB
      parameter (MAXTAB=400)
      integer nprops,ntab(4),itang,ierr
      real*8 props(nprops),celent,cst(20),tab(MAXTAB,10)

      real*8 psideg,fb0fc0,rkc,pi,dval,ym,x0,x1,y0,y1,xx,fr
      real*8 sraw(MAXTAB),eraw(MAXTAB),craw(MAXTAB),ecraw(MAXTAB)
      real*8 dtx0(MAXTAB),dty0(MAXTAB),dcx0(MAXTAB),dcy0(MAXTAB)
      real*8 xs(MAXTAB),ys(MAXTAB)
      integer i,j,k,ioff,ntt,ntc,ntdt,ntdc,itype,ndens,n,ns,nxt

      pi = 3.14159265358979d0
      ierr = 0
      do i = 1,20
         cst(i) = 0.d0
      end do

      ym = props(1)
      cst(1) = ym
      cst(2) = props(2)
      psideg = props(3)
      if (psideg .lt. 0.d0) psideg = 0.d0
      if (psideg .gt. 89.9d0) psideg = 89.9d0
      cst(3) = psideg*pi/180.d0
      cst(4) = tan(cst(3))
      cst(5) = props(4)
      if (cst(5) .le. 0.d0 .or. cst(5) .gt. 1.d0) cst(5) = 0.1d0
      fb0fc0 = props(5)
      if (fb0fc0 .le. 1.d0) fb0fc0 = 1.16d0
      rkc = props(6)
      if (rkc .le. 0.5d0 .or. rkc .gt. 1.d0) rkc = 2.d0/3.d0
      cst(6) = (fb0fc0-1.d0)/(2.d0*fb0fc0-1.d0)
      cst(7) = 3.d0*(1.d0-rkc)/(2.d0*rkc-1.d0)
      cst(8) = props(7)
      if (cst(8) .lt. 0.d0) cst(8) = 0.d0
      cst(9) = props(8)
      if (cst(9) .lt. 0.d0 .or. cst(9) .gt. 1.d0) cst(9) = 0.d0
      cst(10) = props(9)
      if (cst(10) .lt. 0.d0 .or. cst(10) .gt. 1.d0) cst(10) = 1.d0
      i = int(props(10) + 0.5d0)
      if (i .lt. 0) i = 0
      itype = mod(i,10)
      ndens = 1
c     (props(10) >= 10 is accepted; only mod(props(10),10) is used)
      cst(14) = dble(ndens)
      cst(15) = dble(itype)
      cst(12) = ym/(3.d0*(1.d0-2.d0*cst(2)))
      cst(13) = ym/(2.d0*(1.d0+cst(2)))

      ntt = int(props(11) + 0.5d0)
      ntc = int(props(12) + 0.5d0)
      ntdt = int(props(13) + 0.5d0)
      ntdc = int(props(14) + 0.5d0)
      if (ntt .lt. 1) ntt = 1
      if (ntc .lt. 1) ntc = 1
      if (ntdt .lt. 0) ntdt = 0
      if (ntdc .lt. 0) ntdc = 0
      if (ntt.gt.MAXTAB .or. ntc.gt.MAXTAB .or. ntdt.gt.MAXTAB
     $     .or. ntdc.gt.MAXTAB) then
         write(*,*) 'CDP UMAT: a table exceeds MAXTAB =',MAXTAB
         ierr = 1
         return
      end if
      if (nprops .lt. 14+2*(ntt+ntc+ntdt+ntdc)) then
         write(*,*) 'CDP UMAT: NPROPS too small for the tables.'
         ierr = 1
         return
      end if

      ioff = 14
      do i = 1,ntt
         sraw(i) = props(ioff+2*i-1)
         eraw(i) = props(ioff+2*i)
      end do
      ioff = ioff + 2*ntt
      do i = 1,ntc
         craw(i) = props(ioff+2*i-1)
         ecraw(i) = props(ioff+2*i)
      end do
      ioff = ioff + 2*ntc
      if (ntdt .gt. 0) then
         do i = 1,ntdt
            dty0(i) = props(ioff+2*i-1)
            dtx0(i) = props(ioff+2*i)
         end do
      else
         ntdt = 2
         dtx0(1) = 0.d0
         dty0(1) = 0.d0
         dtx0(2) = 1.d10
         dty0(2) = 0.d0
      end if
      ioff = ioff + 2*int(props(13)+0.5d0)
      if (ntdc .gt. 0) then
         do i = 1,ntdc
            dcy0(i) = props(ioff+2*i-1)
            dcx0(i) = props(ioff+2*i)
         end do
      else
         ntdc = 2
         dcx0(1) = 0.d0
         dcy0(1) = 0.d0
         dcx0(2) = 1.d10
         dcy0(2) = 0.d0
      end if
      ioff = ioff + 2*int(props(14)+0.5d0)
c     optional tangent flag
      itang = 0
      if (nprops .ge. ioff+1) then
         if (props(ioff+1) .gt. 0.5d0) itang = 1
      end if

c ---- TYPE=DISPLACEMENT: divide the tension abscissae by CELENT ----
      if (itype .eq. 1 .and. celent .gt. 1.d-12) then
         do i = 1,ntt
            eraw(i) = eraw(i)/celent
         end do
         if (int(props(13)+0.5d0) .gt. 0) then
            do i = 1,ntdt
               dtx0(i) = dtx0(i)/celent
            end do
         end if
      end if

c ---- merged node sets: sigma nodes + damage nodes inside the range,
c      so that sigma and d are both linear on every segment ----
      call cdp_mergenodes(sraw,eraw,ntt,dtx0,dty0,ntdt,ym,MAXTAB,
     $     tab(1,1),tab(1,2),tab(1,3),tab(1,4),ntab(1))
      call cdp_mergenodes(craw,ecraw,ntc,dcx0,dcy0,ntdc,ym,MAXTAB,
     $     tab(1,5),tab(1,6),tab(1,7),tab(1,8),ntab(2))
      ntab(3) = ntdt
      ntab(4) = ntdc

      cst(11) = sraw(1)
      if (cst(11) .le. 0.d0) cst(11) = 1.d-3*ym
      return
      end

c ======================================================================
c Build the merged node table of one sense: nodes e(k) = union of the
c sigma-table abscissae and the damage-table abscissae inside the
c range (sorted, unique); sg(k), dd(k) = sigma and d at the nodes
c (linear in the input tables); kp(k) = kappa at the nodes,
c kappa = e - d/(1-d) sigma/E, forced non-decreasing.
c ======================================================================
      subroutine cdp_mergenodes(sig,x,n,dx,dy,nd,ym,maxn,e,sg,dd,kp,m)
      integer n,nd,maxn,m
      real*8 sig(*),x(*),dx(*),dy(*),ym,e(*),sg(*),dd(*),kp(*)
      real*8 xn(2000),tmp,dval
      integer i,j,k
      m = 0
      do i = 1,n
         m = m+1
         xn(m) = x(i)
      end do
      do i = 1,nd
         if (dx(i).gt.x(1) .and. dx(i).lt.x(n)) then
            m = m+1
            xn(m) = dx(i)
         end if
      end do
      do i = 2,m
         tmp = xn(i)
         j = i-1
 10      if (j.ge.1) then
            if (xn(j) .gt. tmp) then
               xn(j+1) = xn(j)
               j = j-1
               goto 10
            end if
         end if
         xn(j+1) = tmp
      end do
      k = 1
      do i = 2,m
         if (xn(i) .gt. xn(k)*(1.d0+1.d-12)+1.d-300) then
            k = k+1
            xn(k) = xn(i)
         end if
      end do
      m = k
      if (m .gt. maxn) then
         write(*,*) 'CDP UMAT: merged table exceeds MAXTAB'
         call xit
      end if
      do i = 1,m
         e(i) = xn(i)
         call cdp_interp(e(i),x,sig,n,sg(i))
         call cdp_interp(e(i),dx,dy,nd,dval)
         if (dval .lt. 0.d0) dval = 0.d0
         if (dval .gt. 0.99d0) dval = 0.99d0
         dd(i) = dval
         kp(i) = e(i) - dval/(1.d0-dval)*sg(i)/ym
      end do
      do i = 2,m
         if (kp(i) .lt. kp(i-1)) kp(i) = kp(i-1)
      end do
      return
      end

c ======================================================================
c Exact table lookup for one sense (tension: cols 1-4, compression:
c cols 5-8 of tab). Given kappa, find the segment k with
c kp(k) <= kappa <= kp(k+1), solve kappa(e) = kappa for e in
c [e(k), e(k+1)] (sigma and d linear there; safeguarded Newton), and
c return sigma_bar = sigma/(1-d), d(sigma_bar)/d(kappa), d and
c d(d)/d(kappa). Flat extrapolation outside the table.
c ======================================================================
      subroutine cdp_lookup(cst,tab,ntab,isense,kappa,sbar,dsbar,
     $     dval,ddval)
      integer MAXTAB
      parameter (MAXTAB=400)
      real*8 cst(20),tab(MAXTAB,10)
      integer ntab(4),isense
      real*8 kappa,sbar,dsbar,dval,ddval

      real*8 ym,e0,e1,s0,s1,d0,d1,k0,k1,e,sg,dd,kp,dkde,ds,ddd,g,omd
      real*8 lo,hi,glo,sgp,ddp
      integer m,k,ic,it

      ym = cst(1)
      ic = 1
      if (isense .eq. 2) ic = 5
      m = ntab(isense)
      dsbar = 0.d0
      ddval = 0.d0
c     below the first node / above the last: flat
      if (kappa .le. tab(1,ic+3)) then
         dval = tab(1,ic+2)
         sbar = tab(1,ic+1)/(1.d0-dval)
         return
      end if
      if (kappa .ge. tab(m,ic+3)) then
         dval = tab(m,ic+2)
         sbar = tab(m,ic+1)/(1.d0-dval)
         return
      end if
c     segment search (kp non-decreasing)
      k = 1
      do it = 1,m-1
         if (kappa .ge. tab(it,ic+3) .and. kappa .le. tab(it+1,ic+3))
     $        then
            k = it
            goto 20
         end if
      end do
 20   continue
      e0 = tab(k,ic)
      e1 = tab(k+1,ic)
      s0 = tab(k,ic+1)
      s1 = tab(k+1,ic+1)
      d0 = tab(k,ic+2)
      d1 = tab(k+1,ic+2)
      k0 = tab(k,ic+3)
      k1 = tab(k+1,ic+3)
      if (k1 .le. k0 .or. e1 .le. e0) then
         dval = d0
         sbar = s0/(1.d0-d0)
         return
      end if
      ds = (s1-s0)/(e1-e0)
      ddd = (d1-d0)/(e1-e0)
c     Newton on g(e) = e - d(e) sg(e)/((1-d(e)) E) - kappa, bracketed
      lo = e0
      hi = e1
      e = e0 + (e1-e0)*(kappa-k0)/(k1-k0)
      glo = -1.d0
      do it = 1,60
         sg = s0 + ds*(e-e0)
         dd = d0 + ddd*(e-e0)
         omd = 1.d0-dd
         g = e - dd*sg/(omd*ym) - kappa
         if (abs(g) .le. 1.d-14*max(abs(kappa),e1)) goto 30
         dkde = 1.d0 - (ddd*sg*omd + dd*ds*omd + dd*sg*ddd)/
     $        (omd*omd*ym)
         if (g .gt. 0.d0) then
            hi = e
         else
            lo = e
         end if
         if (dkde .gt. 1.d-12) then
            e = e - g/dkde
         else
            e = 0.5d0*(lo+hi)
         end if
         if (e .le. lo .or. e .ge. hi) e = 0.5d0*(lo+hi)
      end do
 30   continue
      sg = s0 + ds*(e-e0)
      dd = d0 + ddd*(e-e0)
      omd = 1.d0-dd
      dkde = 1.d0 - (ddd*sg*omd + dd*ds*omd + dd*sg*ddd)/
     $     (omd*omd*ym)
      dval = dd
      sbar = sg/omd
      if (dkde .gt. 1.d-12) then
         dsbar = (ds*omd + sg*ddd)/(omd*omd)/dkde
         ddval = ddd/dkde
      end if
      return
      end

c ======================================================================
c Piecewise-linear interpolation, flat extrapolation
c ======================================================================
      subroutine cdp_interp(x,xtab,ytab,n,y)
      real*8 x,xtab(*),ytab(*),y
      integer n,i
      if (n .le. 1) then
         y = ytab(1)
         return
      end if
      if (x .le. xtab(1)) then
         y = ytab(1)
         return
      end if
      if (x .ge. xtab(n)) then
         y = ytab(n)
         return
      end if
      do i = 1,n-1
         if (x .ge. xtab(i) .and. x .le. xtab(i+1)) then
            if (xtab(i+1) .eq. xtab(i)) then
               y = ytab(i)
            else
               y = ytab(i) + (ytab(i+1)-ytab(i))*
     $              (x-xtab(i))/(xtab(i+1)-xtab(i))
            end if
            return
         end if
      end do
      y = ytab(n)
      return
      end

c ======================================================================
c Same, also returning the local slope (0 outside the table)
c ======================================================================
      subroutine cdp_interp2(x,xtab,ytab,n,y,dydx)
      real*8 x,xtab(*),ytab(*),y,dydx
      integer n,i
      dydx = 0.d0
      if (n .le. 1) then
         y = ytab(1)
         return
      end if
      if (x .le. xtab(1)) then
         y = ytab(1)
         return
      end if
      if (x .ge. xtab(n)) then
         y = ytab(n)
         return
      end if
      do i = 1,n-1
         if (x .ge. xtab(i) .and. x .le. xtab(i+1)) then
            if (xtab(i+1) .eq. xtab(i)) then
               y = ytab(i)
            else
               dydx = (ytab(i+1)-ytab(i))/(xtab(i+1)-xtab(i))
               y = ytab(i) + dydx*(x-xtab(i))
            end if
            return
         end if
      end do
      y = ytab(n)
      return
      end

c ======================================================================
c Lee-Fenves / Lubliner yield function
c   F = [ q - 3 alpha p + beta <smax> - gamma <-smax> ]/(1-alpha)
c       - sigma_bar_c(kappa_c),    p = -tr(sigma)/3
c ======================================================================
      subroutine cdp_yieldf(cst,tab,ntab,p,q,smax,kt,kc,fval)
      integer MAXTAB
      parameter (MAXTAB=400)
      real*8 cst(20),tab(MAXTAB,10)
      integer ntab(4)
      real*8 p,q,smax,kt,kc,fval
      real*8 sbt,sbc,beta,alpha,dum1,dum2,dum3
      alpha = cst(6)
      call cdp_lookup(cst,tab,ntab,1,kt,sbt,dum1,dum2,dum3)
      call cdp_lookup(cst,tab,ntab,2,kc,sbc,dum1,dum2,dum3)
      if (sbt .lt. 1.d-8) sbt = 1.d-8
      beta = sbc/sbt*(1.d0-alpha) - (1.d0+alpha)
      fval = (q - 3.d0*alpha*p + beta*max(smax,0.d0)
     $     - cst(7)*max(-smax,0.d0))/(1.d0-alpha) - sbc
      return
      end

c ======================================================================
c Stress-weight factor r = sum<si>/sum|si|
c ======================================================================
      subroutine cdp_rfactor(s1,s2,s3,rval)
      real*8 s1,s2,s3,rval,num,den
      num = max(s1,0.d0)+max(s2,0.d0)+max(s3,0.d0)
      den = abs(s1)+abs(s2)+abs(s3)
      if (den .lt. 1.d-10) then
         rval = 0.d0
      else
         rval = num/den
      end if
      return
      end

c ======================================================================
c Cyclic Jacobi eigen-decomposition of a symmetric 3x3 matrix
c ======================================================================
      subroutine cdp_jacobi3(ain,evals,evecs)
      real*8 ain(3,3),evals(3),evecs(3,3)
      real*8 a(3,3),v(3,3)
      real*8 off,theta,t,c,s,tau,apq,app,aqq,akp,akq,vkp,vkq
      integer i,j,k,p,q,sweep
      integer pidx(3),qidx(3)
      data pidx /1,1,2/
      data qidx /2,3,3/
      do i = 1,3
         do j = 1,3
            a(i,j) = ain(i,j)
            v(i,j) = 0.d0
         end do
         v(i,i) = 1.d0
      end do
      do sweep = 1,60
         off = abs(a(1,2))+abs(a(1,3))+abs(a(2,3))
         if (off .lt. 1.d-13*max(abs(a(1,1)),abs(a(2,2)),
     $        abs(a(3,3)),1.d-300)) goto 200
         do k = 1,3
            p = pidx(k)
            q = qidx(k)
            if (abs(a(p,q)) .gt. 1.d-300) then
               theta = (a(q,q)-a(p,p))/(2.d0*a(p,q))
               if (theta .ge. 0.d0) then
                  t = 1.d0/(theta+sqrt(theta*theta+1.d0))
               else
                  t = -1.d0/(-theta+sqrt(theta*theta+1.d0))
               end if
               c = 1.d0/sqrt(t*t+1.d0)
               s = t*c
               tau = s/(1.d0+c)
               apq = a(p,q)
               app = a(p,p)
               aqq = a(q,q)
               a(p,p) = app - t*apq
               a(q,q) = aqq + t*apq
               a(p,q) = 0.d0
               a(q,p) = 0.d0
               do i = 1,3
                  if (i .ne. p .and. i .ne. q) then
                     akp = a(i,p)
                     akq = a(i,q)
                     a(i,p) = akp - s*(akq+tau*akp)
                     a(p,i) = a(i,p)
                     a(i,q) = akq + s*(akp-tau*akq)
                     a(q,i) = a(i,q)
                  end if
               end do
               do i = 1,3
                  vkp = v(i,p)
                  vkq = v(i,q)
                  v(i,p) = vkp - s*(vkq+tau*vkp)
                  v(i,q) = vkq + s*(vkp-tau*vkq)
               end do
            end if
         end do
      end do
 200  continue
      evals(1) = a(1,1)
      evals(2) = a(2,2)
      evals(3) = a(3,3)
      do i = 1,3
         do j = 1,3
            evecs(i,j) = v(i,j)
         end do
      end do
      return
      end

      subroutine cdp_sort3desc(evals,evecs)
      real*8 evals(3),evecs(3,3),tmp,tv(3)
      integer i,j,imax
      do i = 1,2
         imax = i
         do j = i+1,3
            if (evals(j) .gt. evals(imax)) imax = j
         end do
         if (imax .ne. i) then
            tmp = evals(i)
            evals(i) = evals(imax)
            evals(imax) = tmp
            do j = 1,3
               tv(j) = evecs(j,i)
               evecs(j,i) = evecs(j,imax)
               evecs(j,imax) = tv(j)
            end do
         end if
      end do
      return
      end

c ======================================================================
c Material-point update: elastic predictor, eigen-decomposition, yield
c check, implicit return mapping, adaptive sub-incrementation.
c Extra outputs for the analytic tangent: trial principal stresses
c evtr(3) and eigenvectors evec(3,3) of the LAST sub-step, the solved
c xsol = (dlam, q) of that sub-step, iplas (1 = plastic), and the
c number of accepted sub-steps nsub.
c ======================================================================
      subroutine cdp_update(cst,tab,ntab,deps,stateOld,stateNew,
     $     sigeff,evtr,evec,xsol,iplas,nsub,ipFail)
      integer MAXTAB
      parameter (MAXTAB=400)
      real*8 cst(20),tab(MAXTAB,10)
      integer ntab(4)
      real*8 deps(6),stateOld(27),stateNew(27),sigeff(6)
      real*8 evtr(3),evec(3,3),xsol(2)
      integer iplas,nsub,ipFail

      real*8 epsplold(6),epstotold(6),epstotcur(6),epsel(6)
      real*8 sigtrialv(6),Amat(3,3),evals(3),evecs(3,3)
      real*8 ptrial,qtrial,mean,d1,d2,d3
      real*8 kappatcur,kappaccur,kappatn,kappacn
      real*8 ftrial,smaxtrial,ftol
      real*8 dlam,q,p
      real*8 epsplcur(6)
      real*8 sigeffprinc(3),nmat(3,3),dtval,dcval,rstar,dum
      real*8 fracdone,fracstep,fracmin
      integer i,j,converged,attempt,maxattempt

      do i = 1,6
         epsplold(i) = stateOld(i+2)
         epstotold(i) = stateOld(i+10)
         epsplcur(i) = epsplold(i)
      end do
      kappatcur = stateOld(1)
      kappaccur = stateOld(2)

      ipFail = 0
      nsub = 0
      iplas = 0
      xsol(1) = 0.d0
      xsol(2) = 0.d0
      fracdone = 0.d0
      fracstep = 1.d0
      fracmin = 1.d-4
      attempt = 0
      maxattempt = 60
      ftol = 1.d-9

 100  continue
      attempt = attempt+1
      if (attempt .gt. maxattempt) then
         ipFail = 1
         return
      end if
      if (fracstep .gt. 1.d0-fracdone) fracstep = 1.d0-fracdone

      do i = 1,6
         epstotcur(i) = epstotold(i) + (fracdone+fracstep)*deps(i)
         epsel(i) = epstotcur(i) - epsplcur(i)
      end do
      call cdp_stressfromstrain(cst,epsel,sigtrialv)

      Amat(1,1) = sigtrialv(1)
      Amat(2,2) = sigtrialv(2)
      Amat(3,3) = sigtrialv(3)
      Amat(1,2) = sigtrialv(4)
      Amat(2,1) = sigtrialv(4)
      Amat(2,3) = sigtrialv(5)
      Amat(3,2) = sigtrialv(5)
      Amat(1,3) = sigtrialv(6)
      Amat(3,1) = sigtrialv(6)
      call cdp_jacobi3(Amat,evals,evecs)
      call cdp_sort3desc(evals,evecs)

      mean = (evals(1)+evals(2)+evals(3))/3.d0
      ptrial = -mean
      d1 = evals(1)-mean
      d2 = evals(2)-mean
      d3 = evals(3)-mean
      qtrial = sqrt(1.5d0*(d1*d1+d2*d2+d3*d3))
      smaxtrial = evals(1)

      call cdp_yieldf(cst,tab,ntab,ptrial,qtrial,smaxtrial,
     $     kappatcur,kappaccur,ftrial)

      if (ftrial .le. ftol*max(abs(evals(1)),abs(evals(3)),1.d0))
     $     then
         sigeffprinc(1) = evals(1)
         sigeffprinc(2) = evals(2)
         sigeffprinc(3) = evals(3)
         kappatn = kappatcur
         kappacn = kappaccur
         converged = 1
         iplas = 0
         xsol(1) = 0.d0
         xsol(2) = qtrial
      else
         call cdp_returnmap(cst,tab,ntab,evals,kappatcur,kappaccur,
     $        dlam,q,p,kappatn,kappacn,converged)
         if (converged .eq. 1) then
            if (qtrial .gt. 1.d-10) then
               sigeffprinc(1) = -p + (q/qtrial)*d1
               sigeffprinc(2) = -p + (q/qtrial)*d2
               sigeffprinc(3) = -p + (q/qtrial)*d3
            else
               sigeffprinc(1) = -p
               sigeffprinc(2) = -p
               sigeffprinc(3) = -p
            end if
            iplas = 1
            xsol(1) = dlam
            xsol(2) = q
         end if
      end if

      if (converged .ne. 1) then
         fracstep = fracstep*0.5d0
         if (fracstep .lt. fracmin) then
            ipFail = 1
            return
         end if
         goto 100
      end if

c     accept the sub-step
      do i = 1,3
         evtr(i) = evals(i)
         do j = 1,3
            evec(i,j) = evecs(i,j)
         end do
      end do
      do i = 1,3
         do j = 1,3
            nmat(i,j) = evecs(i,1)*evecs(j,1)*sigeffprinc(1)
     $                + evecs(i,2)*evecs(j,2)*sigeffprinc(2)
     $                + evecs(i,3)*evecs(j,3)*sigeffprinc(3)
         end do
      end do
      sigeff(1) = nmat(1,1)
      sigeff(2) = nmat(2,2)
      sigeff(3) = nmat(3,3)
      sigeff(4) = nmat(1,2)
      sigeff(5) = nmat(2,3)
      sigeff(6) = nmat(1,3)

      call cdp_strainfromstress(cst,sigeff,epsel)
      do i = 1,6
         epsplcur(i) = epstotcur(i) - epsel(i)
      end do
      kappatcur = kappatn
      kappaccur = kappacn
      nsub = nsub+1

      fracdone = fracdone + fracstep
      if (fracdone .lt. 0.999999999d0) goto 100

c     damage and r* from the converged end state
      call cdp_dmgt(cst,tab,ntab,kappatcur,dtval,dum)
      call cdp_dmgc(cst,tab,ntab,kappaccur,dcval,dum)
      call cdp_rfactor(sigeffprinc(1),sigeffprinc(2),sigeffprinc(3),
     $     rstar)

      do i = 1,27
         stateNew(i) = stateOld(i)
      end do
      stateNew(1) = kappatcur
      stateNew(2) = kappaccur
      do i = 1,6
         stateNew(i+2) = epsplcur(i)
         stateNew(i+10) = epstotcur(i)
         stateNew(i+20) = sigeff(i)
      end do
      stateNew(9) = dtval
      stateNew(10) = dcval
      stateNew(20) = rstar
      return
      end

c ======================================================================
c Return mapping: 2x2 system in x = (dlam, q)
c   R1 = q (1 + 3 G dlam / Gt(q)) - qtr = 0
c   R2 = F(p(dlam), q, smax(dlam,q), kt(dlam,q), kc(dlam,q)) = 0
c Phase 1: Newton with analytical Jacobian + Armijo backtracking.
c Phase 2: bisection on dlam with the unique q(dlam) from R1 = 0.
c ev(3) = sorted trial principal effective stresses.
c ======================================================================
      subroutine cdp_returnmap(cst,tab,ntab,ev,kto,kco,dlam,q,p,
     $     ktn,kcn,converged)
      integer MAXTAB
      parameter (MAXTAB=400)
      real*8 cst(20),tab(MAXTAB,10)
      integer ntab(4)
      real*8 ev(3),kto,kco,dlam,q,p,ktn,kcn
      integer converged

      real*8 x(2),xn(2),R(2),Rn(2),dR(2,5),dx(2),dkt(5),dkc(5)
      real*8 s(3),ds(3,5),kt,kc,ptr,qtr,mean,d1,d2,d3
      real*8 detJ,lam,normR,normRn,tol
      real*8 glo,ghi,gm,lo,hi,mid,qlo,qhi,qm
      integer it,maxit,ls,tries

      mean = (ev(1)+ev(2)+ev(3))/3.d0
      ptr = -mean
      d1 = ev(1)-mean
      d2 = ev(2)-mean
      d3 = ev(3)-mean
      qtr = sqrt(1.5d0*(d1*d1+d2*d2+d3*d3))

      x(1) = 1.d-10
      x(2) = 0.5d0*qtr
      maxit = 60
      tol = max(abs(ptr),abs(qtr),1.d0)*1.d-11
      converged = 0

c ---- Phase 1 ----
      do it = 1,maxit
         call cdp_resid(cst,tab,ntab,x,ev,kto,kco,R,dR,kt,kc,
     $        dkt,dkc,s,ds,p)
         normR = sqrt(R(1)*R(1)+R(2)*R(2))
         if (normR .lt. tol) then
            converged = 1
            dlam = x(1)
            q = x(2)
            ktn = kt
            kcn = kc
            return
         end if
         detJ = dR(1,1)*dR(2,2)-dR(1,2)*dR(2,1)
         if (abs(detJ) .lt. 1.d-30) goto 200
         dx(1) = -( dR(2,2)*R(1)-dR(1,2)*R(2))/detJ
         dx(2) = -(-dR(2,1)*R(1)+dR(1,1)*R(2))/detJ
         lam = 1.d0
         do ls = 1,30
            xn(1) = max(x(1)+lam*dx(1),0.d0)
            xn(2) = max(x(2)+lam*dx(2),0.d0)
            call cdp_resid(cst,tab,ntab,xn,ev,kto,kco,Rn,dR,kt,kc,
     $           dkt,dkc,s,ds,p)
            normRn = sqrt(Rn(1)*Rn(1)+Rn(2)*Rn(2))
            if (normRn .lt. (1.d0-1.d-4*lam)*normR) goto 100
            lam = lam*0.5d0
         end do
         goto 200
 100     continue
         x(1) = xn(1)
         x(2) = xn(2)
      end do

c ---- Phase 2 ----
 200  continue
      lo = 0.d0
      call cdp_solve_q(cst,tab,ntab,lo,ev,kto,kco,qlo)
      x(1) = lo
      x(2) = qlo
      call cdp_resid(cst,tab,ntab,x,ev,kto,kco,R,dR,kt,kc,
     $     dkt,dkc,s,ds,p)
      glo = R(2)
      hi = 1.d-6
      tries = 0
 210  continue
      call cdp_solve_q(cst,tab,ntab,hi,ev,kto,kco,qhi)
      x(1) = hi
      x(2) = qhi
      call cdp_resid(cst,tab,ntab,x,ev,kto,kco,R,dR,kt,kc,
     $     dkt,dkc,s,ds,p)
      ghi = R(2)
      if (glo*ghi .gt. 0.d0 .and. tries .lt. 80 .and.
     $     hi .le. 10.d0) then
         hi = hi*2.d0
         tries = tries+1
         goto 210
      end if
      if (glo*ghi .gt. 0.d0) then
         converged = 0
         return
      end if
      do it = 1,200
         mid = 0.5d0*(lo+hi)
         call cdp_solve_q(cst,tab,ntab,mid,ev,kto,kco,qm)
         x(1) = mid
         x(2) = qm
         call cdp_resid(cst,tab,ntab,x,ev,kto,kco,R,dR,kt,kc,
     $        dkt,dkc,s,ds,p)
         gm = R(2)
         if (abs(gm) .lt. tol .or.
     $        (hi-lo) .lt. 1.d-14*max(hi,1.d0)) goto 300
         if (glo*gm .le. 0.d0) then
            hi = mid
            ghi = gm
         else
            lo = mid
            glo = gm
         end if
      end do
 300  continue
      dlam = mid
      q = qm
      ktn = kt
      kcn = kc
      converged = 1
      return
      end

c ======================================================================
c Residual of the return mapping and ALL its first derivatives.
c Variables v = 1: dlam, 2: q, 3..5: trial principal stresses ev(j).
c Returns R(2), dR(2,5), kt, kc, dkt(5), dkc(5), the end-of-step
c principal effective stresses s(3) with ds(3,5), and p.
c Columns 1-2 of dR are the Newton Jacobian; columns 3-5 feed the
c implicit-function derivatives of the consistent tangent.
c Every piecewise term (Macaulay brackets, r, table slopes, kappa
c clamps) is differentiated on the branch of the current iterate.
c ======================================================================
      subroutine cdp_resid(cst,tab,ntab,x,ev,kto,kco,R,dR,kt,kc,
     $     dkt,dkc,s,ds,p)
      integer MAXTAB
      parameter (MAXTAB=400)
      real*8 cst(20),tab(MAXTAB,10)
      integer ntab(4)
      real*8 x(2),ev(3),kto,kco,R(2),dR(2,5),kt,kc,dkt(5),dkc(5)
      real*8 s(3),ds(3,5),p

      real*8 dlam,q,mean,ptr,qtr,d(3),dd(3,5),dp(5),dqtr(5)
      real*8 asq,Gt,dGt(5),ratio,dratio(5)
      real*8 Dmax,Dmin,dDmax(5),dDmin(5),demax,demin,ddemax(5)
      real*8 ddemin(5),num,den,dnum(5),dden(5),rw,drw(5)
      real*8 hdkt,hdkc,sbt,sbc,slt,slc,dsbt(5),dsbc(5)
      real*8 beta,dbeta(5),macp,macm,dmacp(5),dmacm(5)
      real*8 alpha,gam,bulkK,shearG,tanpsi,tiny,dum1,dum2
      integer i,j,v

      alpha = cst(6)
      gam = cst(7)
      bulkK = cst(12)
      shearG = cst(13)
      tanpsi = cst(4)
      tiny = 1.d-10

      dlam = x(1)
      q = x(2)
      mean = (ev(1)+ev(2)+ev(3))/3.d0
      ptr = -mean
      do i = 1,3
         d(i) = ev(i)-mean
      end do
      qtr = sqrt(1.5d0*(d(1)**2+d(2)**2+d(3)**2))

      do v = 1,5
         dp(v) = 0.d0
         dqtr(v) = 0.d0
         dGt(v) = 0.d0
         dratio(v) = 0.d0
         dDmax(v) = 0.d0
         dDmin(v) = 0.d0
         drw(v) = 0.d0
         dkt(v) = 0.d0
         dkc(v) = 0.d0
         dmacp(v) = 0.d0
         dmacm(v) = 0.d0
         do i = 1,3
            dd(i,v) = 0.d0
         end do
      end do
      do j = 1,3
         dp(2+j) = -1.d0/3.d0
         do i = 1,3
            dd(i,2+j) = -1.d0/3.d0
         end do
         dd(j,2+j) = dd(j,2+j) + 1.d0
      end do

      p = ptr + bulkK*dlam*tanpsi
      dp(1) = bulkK*tanpsi

      asq = (cst(5)*cst(11)*tanpsi)**2
      Gt = sqrt(asq + q*q)
      if (Gt .lt. 1.d-12) Gt = 1.d-12
      dGt(2) = q/Gt

      R(1) = q*(1.d0+3.d0*shearG*dlam/Gt) - qtr
      do v = 1,5
         dR(1,v) = 0.d0
      end do
      dR(1,1) = 3.d0*shearG*q/Gt
      dR(1,2) = 1.d0 + 3.d0*shearG*dlam*asq/Gt**3

      if (qtr .gt. tiny) then
         do j = 1,3
            dqtr(2+j) = 1.5d0*d(j)/qtr
            dR(1,2+j) = -dqtr(2+j)
         end do
         ratio = q/qtr
         dratio(2) = 1.d0/qtr
         do j = 1,3
            dratio(2+j) = -ratio*dqtr(2+j)/qtr
         end do
      else
         ratio = 0.d0
      end if

      do i = 1,3
         s(i) = -p + ratio*d(i)
         do v = 1,5
            ds(i,v) = -dp(v) + dratio(v)*d(i) + ratio*dd(i,v)
         end do
      end do

      if (qtr .gt. tiny) then
         Dmax = 1.5d0*ratio*d(1)/Gt + tanpsi/3.d0
         Dmin = 1.5d0*ratio*d(3)/Gt + tanpsi/3.d0
         do v = 1,5
            dDmax(v) = 1.5d0*(dratio(v)*d(1)/Gt + ratio*dd(1,v)/Gt
     $           - ratio*d(1)*dGt(v)/Gt**2)
            dDmin(v) = 1.5d0*(dratio(v)*d(3)/Gt + ratio*dd(3,v)/Gt
     $           - ratio*d(3)*dGt(v)/Gt**2)
         end do
      else
         Dmax = tanpsi/3.d0
         Dmin = tanpsi/3.d0
      end if
      demax = dlam*Dmax
      demin = dlam*Dmin
      do v = 1,5
         ddemax(v) = dlam*dDmax(v)
         ddemin(v) = dlam*dDmin(v)
      end do
      ddemax(1) = ddemax(1) + Dmax
      ddemin(1) = ddemin(1) + Dmin

c     stress-weight factor rw and its derivatives
      den = abs(s(1))+abs(s(2))+abs(s(3))
      if (den .lt. tiny) then
         rw = 0.d0
      else
         num = max(s(1),0.d0)+max(s(2),0.d0)+max(s(3),0.d0)
         rw = num/den
         do v = 1,5
            dnum(v) = 0.d0
            dden(v) = 0.d0
            do i = 1,3
               if (s(i) .gt. 0.d0) then
                  dnum(v) = dnum(v) + ds(i,v)
                  dden(v) = dden(v) + ds(i,v)
               else if (s(i) .lt. 0.d0) then
                  dden(v) = dden(v) - ds(i,v)
               end if
            end do
            drw(v) = (dnum(v)*den - num*dden(v))/den**2
         end do
      end if

c     hardening increments (clamped non-negative)
      hdkt = rw*demax
      if (hdkt .gt. 0.d0) then
         do v = 1,5
            dkt(v) = drw(v)*demax + rw*ddemax(v)
         end do
      else
         hdkt = 0.d0
      end if
      hdkc = (rw-1.d0)*demin
      if (hdkc .gt. 0.d0) then
         do v = 1,5
            dkc(v) = drw(v)*demin + (rw-1.d0)*ddemin(v)
         end do
      else
         hdkc = 0.d0
      end if
      kt = kto + hdkt
      kc = kco + hdkc

c     cohesions and beta
      call cdp_lookup(cst,tab,ntab,1,kt,sbt,slt,dum1,dum2)
      call cdp_lookup(cst,tab,ntab,2,kc,sbc,slc,dum1,dum2)
      if (sbt .lt. 1.d-8) then
         sbt = 1.d-8
         slt = 0.d0
      end if
      do v = 1,5
         dsbt(v) = slt*dkt(v)
         dsbc(v) = slc*dkc(v)
      end do
      beta = sbc/sbt*(1.d0-alpha) - (1.d0+alpha)
      do v = 1,5
         dbeta(v) = (1.d0-alpha)*(dsbc(v)*sbt-sbc*dsbt(v))/sbt**2
      end do

c     Macaulay brackets of smax = s(1)
      macp = max(s(1),0.d0)
      macm = max(-s(1),0.d0)
      if (s(1) .gt. 0.d0) then
         do v = 1,5
            dmacp(v) = ds(1,v)
         end do
      else if (s(1) .lt. 0.d0) then
         do v = 1,5
            dmacm(v) = -ds(1,v)
         end do
      end if

      R(2) = (q - 3.d0*alpha*p + beta*macp - gam*macm)/(1.d0-alpha)
     $     - sbc
      do v = 1,5
         dR(2,v) = (-3.d0*alpha*dp(v) + dbeta(v)*macp
     $        + beta*dmacp(v) - gam*dmacm(v))/(1.d0-alpha) - dsbc(v)
      end do
      dR(2,2) = dR(2,2) + 1.d0/(1.d0-alpha)
      return
      end

c ======================================================================
c Unique q >= 0 root of R1(dlam, q) = 0 for fixed dlam (dR1/dq >= 1)
c ======================================================================
      subroutine cdp_solve_q(cst,tab,ntab,dlam,ev,kto,kco,qsol)
      integer MAXTAB
      parameter (MAXTAB=400)
      real*8 cst(20),tab(MAXTAB,10)
      integer ntab(4)
      real*8 dlam,ev(3),kto,kco,qsol
      real*8 xx(2),R(2),dR(2,5),kt,kc,dkt(5),dkc(5),s(3),ds(3,5),p
      real*8 lo,hi,flo,fhi,mid,fm,qtr,mean
      integer tries,it

      mean = (ev(1)+ev(2)+ev(3))/3.d0
      qtr = sqrt(1.5d0*((ev(1)-mean)**2+(ev(2)-mean)**2
     $     +(ev(3)-mean)**2))
      lo = 0.d0
      xx(1) = dlam
      xx(2) = lo
      call cdp_resid(cst,tab,ntab,xx,ev,kto,kco,R,dR,kt,kc,
     $     dkt,dkc,s,ds,p)
      flo = R(1)
      if (abs(flo) .lt. 1.d-13) then
         qsol = lo
         return
      end if
      hi = max(qtr,1.d0)
      xx(2) = hi
      call cdp_resid(cst,tab,ntab,xx,ev,kto,kco,R,dR,kt,kc,
     $     dkt,dkc,s,ds,p)
      fhi = R(1)
      tries = 0
 10   continue
      if (flo*fhi .gt. 0.d0 .and. tries .lt. 60) then
         hi = hi*2.d0
         xx(2) = hi
         call cdp_resid(cst,tab,ntab,xx,ev,kto,kco,R,dR,kt,kc,
     $        dkt,dkc,s,ds,p)
         fhi = R(1)
         tries = tries+1
         goto 10
      end if
      mid = 0.5d0*(lo+hi)
      do it = 1,100
         mid = 0.5d0*(lo+hi)
         xx(2) = mid
         call cdp_resid(cst,tab,ntab,xx,ev,kto,kco,R,dR,kt,kc,
     $        dkt,dkc,s,ds,p)
         fm = R(1)
         if (abs(fm) .lt. 1.d-13 .or.
     $        (hi-lo) .lt. 1.d-13*max(hi,1.d0)) goto 20
         if (flo*fm .le. 0.d0) then
            hi = mid
            fhi = fm
         else
            lo = mid
            flo = fm
         end if
      end do
 20   continue
      qsol = mid
      return
      end

c ======================================================================
c Damage lookups kappa -> eck -> d, with the chained slope dd/dkappa
c ======================================================================
      subroutine cdp_dmgt(cst,tab,ntab,kt,dtval,ddt)
      integer MAXTAB
      parameter (MAXTAB=400)
      real*8 cst(20),tab(MAXTAB,10)
      integer ntab(4)
      real*8 kt,dtval,ddt,sb,dsb
      call cdp_lookup(cst,tab,ntab,1,kt,sb,dsb,dtval,ddt)
      return
      end

      subroutine cdp_dmgc(cst,tab,ntab,kc,dcval,ddc)
      integer MAXTAB
      parameter (MAXTAB=400)
      real*8 cst(20),tab(MAXTAB,10)
      integer ntab(4)
      real*8 kc,dcval,ddc,sb,dsb
      call cdp_lookup(cst,tab,ntab,2,kc,sb,dsb,dcval,ddc)
      return
      end

c ======================================================================
c Isotropic elasticity, internal ordering, engineering shear
c ======================================================================
      subroutine cdp_stressfromstrain(cst,strain,stress)
      real*8 cst(20),strain(6),stress(6),factor,pr
      pr = cst(2)
      factor = cst(1)/((1.d0+pr)*(1.d0-2.d0*pr))
      stress(1) = factor*((1.d0-pr)*strain(1)+pr*strain(2)+
     $     pr*strain(3))
      stress(2) = factor*(pr*strain(1)+(1.d0-pr)*strain(2)+
     $     pr*strain(3))
      stress(3) = factor*(pr*strain(1)+pr*strain(2)+
     $     (1.d0-pr)*strain(3))
      stress(4) = factor*((1.d0-2.d0*pr)/2.d0)*strain(4)
      stress(5) = factor*((1.d0-2.d0*pr)/2.d0)*strain(5)
      stress(6) = factor*((1.d0-2.d0*pr)/2.d0)*strain(6)
      return
      end

      subroutine cdp_strainfromstress(cst,stress,strain)
      real*8 cst(20),strain(6),stress(6),ym,pr
      ym = cst(1)
      pr = cst(2)
      strain(1) = (stress(1)-pr*stress(2)-pr*stress(3))/ym
      strain(2) = (-pr*stress(1)+stress(2)-pr*stress(3))/ym
      strain(3) = (-pr*stress(1)-pr*stress(2)+stress(3))/ym
      strain(4) = 2.d0*(1.d0+pr)*stress(4)/ym
      strain(5) = 2.d0*(1.d0+pr)*stress(5)/ym
      strain(6) = 2.d0*(1.d0+pr)*stress(6)/ym
      return
      end

      subroutine cdp_elasticmatrix(cst,elmat)
      real*8 cst(20),elmat(6,6),factor,pr
      integer i,j
      pr = cst(2)
      do i = 1,6
         do j = 1,6
            elmat(i,j) = 0.d0
         end do
      end do
      factor = cst(1)/((1.d0+pr)*(1.d0-2.d0*pr))
      do i = 1,3
         do j = 1,3
            elmat(i,j) = factor*pr
         end do
         elmat(i,i) = factor*(1.d0-pr)
         elmat(i+3,i+3) = factor*(1.d0-2.d0*pr)/2.d0
      end do
      return
      end

c ======================================================================
c Internal (11,22,33,12,23,13) -> Abaqus ordering
c ======================================================================
      subroutine cdp_mapsig(sig,sigum,ips)
      real*8 sig(6),sigum(6)
      integer ips,k
      do k = 1,6
         sigum(k) = 0.d0
      end do
      if (ips .eq. 1) then
         sigum(1) = sig(1)
         sigum(2) = sig(2)
         sigum(3) = sig(4)
      else
         sigum(1) = sig(1)
         sigum(2) = sig(2)
         sigum(3) = sig(3)
         sigum(4) = sig(4)
         sigum(5) = sig(6)
         sigum(6) = sig(5)
      end if
      return
      end

c ======================================================================
c Viscous (Duvaut-Lions) regularisation of dt, dc and combination
c   (1-d) = (1 - st dc_v)(1 - sc dt_v), st = 1-wt r, sc = 1-wc(1-r)
c vfac = dt/(mu+dt) (1 if mu = 0) is returned for the tangent.
c ======================================================================
      subroutine cdp_combine(cst,sigeff,dtinv,dcinv,rstar,dtvold,
     $     dcvold,dtime,signom,dtv,dcv,vfac)
      real*8 cst(20),sigeff(6),dtinv,dcinv,rstar,dtvold,dcvold,dtime
      real*8 signom(6),dtv,dcv,vfac
      real*8 st,sc,dcomb,visc
      integer k
      visc = cst(8)
      if (visc .gt. 0.d0) then
         vfac = dtime/(visc+dtime)
         dtv = (1.d0-vfac)*dtvold + vfac*dtinv
         dcv = (1.d0-vfac)*dcvold + vfac*dcinv
         if (dtv .lt. dtvold) dtv = dtvold
         if (dcv .lt. dcvold) dcv = dcvold
      else
         vfac = 1.d0
         dtv = dtinv
         dcv = dcinv
      end if
      st = 1.d0 - cst(9)*rstar
      sc = 1.d0 - cst(10)*(1.d0-rstar)
      dcomb = 1.d0 - (1.d0-st*dcv)*(1.d0-sc*dtv)
      do k = 1,6
         signom(k) = (1.d0-dcomb)*sigeff(k)
      end do
      return
      end

c ======================================================================
c Analytic consistent tangent dd6(6,6) = d(sigma)/d(eps), internal
c Voigt ordering with engineering shear strains.
c
c   sigma = (1-d) sigma_bar,   sigma_bar = sum_i s_i(lambda) n_i x n_i
c with lambda = trial principal stresses (functions of eps through
c the elastic predictor) and s_i from the return mapping. For an
c isotropic tensor function (Mandel notation):
c   d sigma_bar / d sigma_tr = sum_ij S_ij M_i M_j^T
c                             + sum_{i<j} 2 theta_ij E_ij E_ij^T
c   S_ij = D s_i / D lambda_j  (implicit function theorem on R = 0)
c   theta_ij = (s_i - s_j)/(lambda_i - lambda_j), or its limit
c              (S_ii - S_ij) when lambda_i = lambda_j.
c The damage d depends on kappa_t, kappa_c (through dt, dc and the
c viscous factor) and on r(s). iana = 0 if the tangent could not be
c formed (singular Jacobian) -> caller falls back to finite differences.
c ======================================================================
      subroutine cdp_tangent(cst,tab,ntab,evtr,evec,xsol,iplas,
     $     stateO,stateN,sigeff,dtvold,dcvold,dtime,vfac,dd6,iana)
      integer MAXTAB
      parameter (MAXTAB=400)
      real*8 cst(20),tab(MAXTAB,10)
      integer ntab(4),iplas,iana
      real*8 evtr(3),evec(3,3),xsol(2),stateO(27),stateN(27)
      real*8 sigeff(6),dtvold,dcvold,dtime,vfac,dd6(6,6)

      real*8 R(2),dR(2,5),kt,kc,dkt(5),dkc(5),s(3),ds(3,5),p
      real*8 SS(3,3),tkt(3),tkc(3),dx(2,3),detJ,rhs1,rhs2
      real*8 dt,dc,ddt,ddc,dtv,dcv,fdt,fdc,rstar,drds(3),Drw(3)
      real*8 num,den,st,sc,A,B,dcomb,pdt,pdc,pr,Dd(3)
      real*8 PM(6,3),PE(6,3),th(3),L(6,6),elmat(6,6),C(6,6),g(6)
      real*8 ddeps(6),sq2,scal,ni(3),nj(3)
      integer i,j,k,m,ip,jp,pairs(2,3)
      data pairs /1,2, 1,3, 2,3/

      iana = 0
      sq2 = sqrt(2.d0)

c ---- principal-space derivatives S_ij, Dkt, Dkc ----
      if (iplas .eq. 1) then
         call cdp_resid(cst,tab,ntab,xsol,evtr,stateO(1),stateO(2),
     $        R,dR,kt,kc,dkt,dkc,s,ds,p)
         detJ = dR(1,1)*dR(2,2)-dR(1,2)*dR(2,1)
         if (abs(detJ) .lt. 1.d-30) return
         do j = 1,3
            rhs1 = -dR(1,2+j)
            rhs2 = -dR(2,2+j)
            dx(1,j) = ( dR(2,2)*rhs1-dR(1,2)*rhs2)/detJ
            dx(2,j) = (-dR(2,1)*rhs1+dR(1,1)*rhs2)/detJ
         end do
         do i = 1,3
            do j = 1,3
               SS(i,j) = ds(i,2+j) + ds(i,1)*dx(1,j) + ds(i,2)*dx(2,j)
            end do
         end do
         do j = 1,3
            tkt(j) = dkt(2+j) + dkt(1)*dx(1,j) + dkt(2)*dx(2,j)
            tkc(j) = dkc(2+j) + dkc(1)*dx(1,j) + dkc(2)*dx(2,j)
         end do
      else
         do i = 1,3
            s(i) = evtr(i)
            tkt(i) = 0.d0
            tkc(i) = 0.d0
            do j = 1,3
               SS(i,j) = 0.d0
            end do
            SS(i,i) = 1.d0
         end do
      end if

c ---- degenerate trial eigenvalues: the sorted-principal functions
c      are not smooth under a swap of equal eigenvalues (max/min
c      assignment). Use the symmetric generalized derivative: average
c      S, tkt, tkc over the permutations that leave lambda unchanged.
c      This is what a central-difference tangent converges to and what
c      keeps the tangent symmetric under the degenerate directions.
      scal = max(abs(evtr(1)),abs(evtr(3)),1.d-300)
      call cdp_symdeg(evtr,scal,SS,tkt,tkc)

c ---- damage chain ----
      call cdp_dmgt(cst,tab,ntab,stateN(1),dt,ddt)
      call cdp_dmgc(cst,tab,ntab,stateN(2),dc,ddc)
      if (cst(8) .gt. 0.d0) then
         dtv = (1.d0-vfac)*dtvold + vfac*dt
         dcv = (1.d0-vfac)*dcvold + vfac*dc
         fdt = vfac
         fdc = vfac
         if (dtv .lt. dtvold) then
            dtv = dtvold
            fdt = 0.d0
         end if
         if (dcv .lt. dcvold) then
            dcv = dcvold
            fdc = 0.d0
         end if
      else
         dtv = dt
         dcv = dc
         fdt = 1.d0
         fdc = 1.d0
      end if
c     r and dr/ds
      den = abs(s(1))+abs(s(2))+abs(s(3))
      do i = 1,3
         drds(i) = 0.d0
      end do
      if (den .lt. 1.d-10) then
         rstar = 0.d0
      else
         num = max(s(1),0.d0)+max(s(2),0.d0)+max(s(3),0.d0)
         rstar = num/den
         do i = 1,3
            if (s(i) .gt. 0.d0) then
               drds(i) = (den - num)/den**2
            else if (s(i) .lt. 0.d0) then
               drds(i) = num/den**2
            end if
         end do
      end if
      do j = 1,3
         Drw(j) = 0.d0
         do i = 1,3
            Drw(j) = Drw(j) + drds(i)*SS(i,j)
         end do
      end do
      st = 1.d0 - cst(9)*rstar
      sc = 1.d0 - cst(10)*(1.d0-rstar)
      A = 1.d0 - st*dcv
      B = 1.d0 - sc*dtv
      dcomb = 1.d0 - A*B
      pdt = A*sc
      pdc = st*B
      pr = -(cst(9)*dcv*B - cst(10)*dtv*A)
      do j = 1,3
         Dd(j) = pdt*fdt*ddt*tkt(j) + pdc*fdc*ddc*tkc(j) + pr*Drw(j)
      end do

c ---- Mandel eigenprojections ----
      do i = 1,3
         do k = 1,3
            ni(k) = evec(k,i)
         end do
         PM(1,i) = ni(1)*ni(1)
         PM(2,i) = ni(2)*ni(2)
         PM(3,i) = ni(3)*ni(3)
         PM(4,i) = sq2*ni(1)*ni(2)
         PM(5,i) = sq2*ni(2)*ni(3)
         PM(6,i) = sq2*ni(1)*ni(3)
      end do
      scal = max(abs(evtr(1)),abs(evtr(3)),1.d-300)
      do m = 1,3
         ip = pairs(1,m)
         jp = pairs(2,m)
         do k = 1,3
            ni(k) = evec(k,ip)
            nj(k) = evec(k,jp)
         end do
         PE(1,m) = ni(1)*nj(1)
         PE(2,m) = ni(2)*nj(2)
         PE(3,m) = ni(3)*nj(3)
         PE(4,m) = sq2*0.5d0*(ni(1)*nj(2)+ni(2)*nj(1))
         PE(5,m) = sq2*0.5d0*(ni(2)*nj(3)+ni(3)*nj(2))
         PE(6,m) = sq2*0.5d0*(ni(1)*nj(3)+ni(3)*nj(1))
         if (abs(evtr(ip)-evtr(jp)) .gt. 1.d-7*scal) then
            th(m) = (s(ip)-s(jp))/(evtr(ip)-evtr(jp))
         else
            th(m) = 0.5d0*(SS(ip,ip)-SS(ip,jp)+SS(jp,jp)-SS(jp,ip))
         end if
      end do

c ---- L (Mandel 6x6) ----
      do i = 1,6
         do j = 1,6
            L(i,j) = 0.d0
            do ip = 1,3
               do jp = 1,3
                  L(i,j) = L(i,j) + SS(ip,jp)*PM(i,ip)*PM(j,jp)
               end do
            end do
            do m = 1,3
               L(i,j) = L(i,j) + 2.d0*th(m)*PE(i,m)*PE(j,m)
            end do
         end do
      end do

c ---- d sigma_bar_V / d eps_V = T^-1 L T elmat ----
      call cdp_elasticmatrix(cst,elmat)
      do i = 1,6
         do j = 1,6
            C(i,j) = L(i,j)
            if (j .ge. 4) C(i,j) = C(i,j)*sq2
            if (i .ge. 4) C(i,j) = C(i,j)/sq2
         end do
      end do
      do i = 1,6
         do j = 1,6
            L(i,j) = 0.d0
            do k = 1,6
               L(i,j) = L(i,j) + C(i,k)*elmat(k,j)
            end do
         end do
      end do
c ---- d d / d eps_V ----
      do m = 1,6
         g(m) = 0.d0
         do j = 1,3
            g(m) = g(m) + Dd(j)*PM(m,j)
         end do
         if (m .ge. 4) g(m) = g(m)*sq2
      end do
      do k = 1,6
         ddeps(k) = 0.d0
         do m = 1,6
            ddeps(k) = ddeps(k) + g(m)*elmat(m,k)
         end do
      end do
c ---- assemble ----
      do i = 1,6
         do j = 1,6
            dd6(i,j) = (1.d0-dcomb)*L(i,j) - sigeff(i)*ddeps(j)
         end do
      end do
      iana = 1
      return
      end

c ======================================================================
c Average S(3,3), v1(3), v2(3) over the permutations of the indices
c of (numerically) equal eigenvalues ev(3) (sorted descending).
c ======================================================================
      subroutine cdp_symdeg(ev,scal,S,v1,v2)
      real*8 ev(3),scal,S(3,3),v1(3),v2(3)
      real*8 Sa(3,3),w1(3),w2(3),tol
      integer perm(3,6),np,k,a,b,i,j
      logical d12,d23
      data perm /1,2,3, 2,1,3, 1,3,2, 3,2,1, 2,3,1, 3,1,2/
      tol = 1.d-7*scal
      d12 = abs(ev(1)-ev(2)) .le. tol
      d23 = abs(ev(2)-ev(3)) .le. tol
      if (.not.d12 .and. .not.d23) return
      if (d12 .and. d23) then
         np = 6
      else if (d12) then
         np = 2
      else
         np = 3
      end if
      do a = 1,3
         w1(a) = 0.d0
         w2(a) = 0.d0
         do b = 1,3
            Sa(a,b) = 0.d0
         end do
      end do
      do k = 1,6
         if (np.eq.6 .or. (np.eq.2 .and. (k.eq.1.or.k.eq.2)) .or.
     $        (np.eq.3 .and. (k.eq.1.or.k.eq.3))) then
            do a = 1,3
               i = perm(a,k)
               w1(a) = w1(a) + v1(i)
               w2(a) = w2(a) + v2(i)
               do b = 1,3
                  j = perm(b,k)
                  Sa(a,b) = Sa(a,b) + S(i,j)
               end do
            end do
         end if
      end do
      if (np .eq. 3) np = 2
      do a = 1,3
         v1(a) = w1(a)/dble(np)
         v2(a) = w2(a)/dble(np)
         do b = 1,3
            S(a,b) = Sa(a,b)/dble(np)
         end do
      end do
      return
      end

c ======================================================================
c One-time diagnostic printout
c ======================================================================
      subroutine cdp_print(cst,ntab,itang)
      real*8 cst(20)
      integer ntab(4),itang
      write(*,*) '****************************************'
      write(*,*) ' ABAQUS is using the user material CDP_DCM'
      write(*,*) ' (Concrete Damaged Plasticity) UMAT'
      write(*,*) '----------------------------------------'
      write(*,2) '  E ..................... = ',cst(1)
      write(*,2) '  nu .................... = ',cst(2)
      write(*,2) '  dilation angle (rad) .. = ',cst(3)
      write(*,2) '  eccentricity .......... = ',cst(5)
      write(*,2) '  alpha ................. = ',cst(6)
      write(*,2) '  gamma ................. = ',cst(7)
      write(*,2) '  viscosity mu .......... = ',cst(8)
      write(*,2) '  wt .................... = ',cst(9)
      write(*,2) '  wc .................... = ',cst(10)
      write(*,2) '  sigma_t0 .............. = ',cst(11)
      write(*,*) '  tension type (0 strain/1 disp) = ',int(cst(15))
      write(*,*) '  table rows t/c/dt/dc .. = ',ntab(1),ntab(2),
     $     ntab(3),ntab(4)
      if (itang .eq. 0) then
         write(*,*) '  tangent: analytic (FD fallback)'
      else
         write(*,*) '  tangent: finite differences'
      end if
      write(*,*) '****************************************'
 2    format(1x,A,1pE12.5)
      return
      end
