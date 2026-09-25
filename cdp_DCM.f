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
c verifications", Applied Sciences (under review)
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
c 3) Abaqus Analysis User's Guide, section "Concrete damaged
c    plasticity", and Abaqus Theory Guide (SIMULIA/Dassault Systemes).
c 4) P. Grassl, D. Xenos, U. Nyström, R. Rempling, K. Gylltoft.: 
c    "CDPM2: A damage-plasticity approach to modelling the failure of 
c    concrete". International Journal of Solids and Structures. 
c    Volume 50, Issue 24, pp. 3805-3816, 2013.
c 5) seungwookseok/ABAQUS-version-CDPM2 [Computer Software]. (2019)
c    Retrieved from https://github.com/seungwookseok/ABAQUS-version-CDPM2
c
c
c
c
c User-defined material properties (props):
c
c props(1)  E    : Young's modulus
c props(2)  nu   : Poisson's ratio
c props(3)  psi  : Dilation angle, degrees (measured in the p-q plane
c                  at high confining pressure; typical 25-42 deg)
c props(4)  ecc  : Flow potential eccentricity (default 0.1 if <=0)
c props(5)  fb0/fc0 : ratio of biaxial to uniaxial compressive stress
c                  at first yield (default 1.16 if <=1 given)
c props(6)  Kc   : ratio of the 2nd stress invariant on the tensile
c                  meridian to that on the compressive meridian, at
c                  any p with max principal stress negative
c                  (default 2/3 if given outside (0.5,1])
c props(7)  mu   : viscosity parameter (>=0). Practical Duvaut-Lions-
c                  type regularization of the damage variables dt,dc
c                  (helps convergence once softening/localization
c                  starts). 0 = off = rate-independent model. This is
c                  an engineering approximation of Abaqus's internal
c                  viscoplastic regularization, not a bit-identical
c                  reproduction of it.
c props(8)  wt   : tension stiffness recovery factor, 0<=wt<=1
c                  (default 0 = no recovery of tensile stiffness when
c                  the stress state turns compressive)
c props(9)  wc   : compression stiffness recovery factor, 0<=wc<=1
c                  (default 1 = full recovery of compressive stiffness
c                  on crack closure, i.e. load reversal into
c                  compression)
c props(10) itype: tension-table independent-variable flag.
c                  0 = values are cracking STRAIN directly (matches
c                      *CONCRETE TENSION STIFFENING, TYPE=STRAIN)
c                  1 = values are crack-opening DISPLACEMENT; this
c                      UMAT converts them to strain by dividing by the
c                      element characteristic length CELENT (crack-
c                      band regularization), matching TYPE=DISPLACEMENT
c                  Also governs the independent-variable convention of
c                  the tension damage table (props block 4 below).
c                  Compression tables are always strain-based.
c props(11) ntt  : number of rows in the tension stiffening table (>=1)
c props(12) ntc  : number of rows in the compression hardening table
c                  (>=1)
c props(13) ntdt : number of rows in the tension damage table
c                  (0 = no tension damage: dt = 0 always)
c props(14) ntdc : number of rows in the compression damage table
c                  (0 = no compression damage: dc = 0 always)
c
c Followed by four variable-length tables, each row = 2 properties,
c in this fixed order (only as many rows as ntt/ntc/ntdt/ntdc say):
c
c   Block 1 (ntt rows): tension stiffening, (sigma_t, cracking strain
c            or displacement per props(10)). Row 1 must be at cracking
c            strain/displacement = 0; rows in increasing order.
c   Block 2 (ntc rows): compression hardening, (sigma_c, inelastic
c            strain). Row 1 at inelastic strain = 0; increasing order.
c   Block 3 (ntdt rows, only if ntdt>0): tension damage, (dt, cracking
c            strain or displacement per props(10)).
c   Block 4 (ntdc rows, only if ntdc>0): compression damage,
c            (dc, inelastic strain).
c
c Total NPROPS = 14 + 2*ntt + 2*ntc + 2*ntdt + 2*ntdc.
c Declare  *USER MATERIAL, CONSTANTS=<NPROPS>, UNSYMM
c (UNSYMM because the flow rule is non-associated, so DDSDDE is
c genuinely unsymmetric -- see notes near the tangent computation).
c
c ----------------------------------------------------------------------
c State variables (statev), Depvar = 27:
c
c statev(1)      kappa_t  : tension hardening/cracking variable
c                           (equivalent tensile plastic strain measure)
c statev(2)      kappa_c  : compression hardening/crushing variable
c statev(3:8)    plastic strain (11,22,33,12,23,13; eng. shear;
c                internal component ordering, see below)
c statev(9)      dt       : tension damage (inviscid/rate-independent)
c statev(10)     dc       : compression damage (inviscid)
c statev(11:16)  total strain (11,22,33,12,23,13; internal ordering)
c statev(17)     dt_v     : viscous (regularized) tension damage
c                           actually applied to the reported stress
c                           (= statev(9) when props(7)=mu=0)
c statev(18)     dc_v     : viscous (regularized) compression damage
c                           actually applied
c statev(19)     d_comb   : combined scalar stiffness degradation
c                           actually applied this increment
c                           (1-d_comb) = (1-st*dc_v)(1-sc*dt_v)
c statev(20)     r*       : multiaxial stress-weight factor at the
c                           converged end-of-increment effective
c                           stress (1 = all tension, 0 = all
c                           compression)
c statev(21:26)  effective (undamaged) stress (11,22,33,12,23,13;
c                internal ordering) -- diagnostic
c statev(27)     status flag, reserved (always 1.0)
c
c Set Depvar = 27 in the material definition.
c
c ----------------------------------------------------------------------
c Conventions:
c Internally stresses/strains are stored as 6-vectors in the order
c (11,22,33,12,23,13) with ENGINEERING shear strains. ABAQUS hands
c UMAT strains in the order (11,22,33,12,13,23); components 5 and 6
c are swapped at the interface (same convention as many published
c ABAQUS UMATs).
c
c Notes/limitations:
c   - Small-strain: state tensors are not rotated with DROT; with
c     NLGEOM=YES use only where rotations remain small over an
c     increment (standard UMAT caveat).
c   - NTENS=6 (3D) and NTENS=4 (plane strain, axisymmetric) are
c     supported. Plane STRESS (NTENS=3) is NOT supported and the
c     routine stops with an error message if it is detected, because
c     the strain/stress component ordering for plane stress
c     (11,22,12; no independent component 3) is fundamentally
c     different and is not handled by the mapping used here.
c   - As documented for the native ABAQUS CDP model itself, this
c     formulation has no compressive "cap": very high hydrostatic
c     compression does not saturate. This is a property of the base
c     Lubliner/Lee-Fenves surface, not an implementation shortcut.
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

      real*8 ym,pr,psi,tanpsi,ecc,alpha,gammap,fb0fc0,rkc,visc,wt,wc,
     1     sigt0,bulkK,shearG
      common/cdpc/ym,pr,psi,tanpsi,ecc,alpha,gammap,fb0fc0,rkc,visc,
     1     wt,wc,sigt0,bulkK,shearG

      real*8 deps(6),depsP(6),stateO(27),stateN(27),stateP(27)
      real*8 sigeff(6),sigeffP(6),sigeffM(6)
      real*8 sig(6),sigP(6),sigM(6),sigUM(6),sigUMP(6),sigUMM(6)
      real*8 dtinv,dcinv,dcinvP,dtinvP,dcinvM,dtinvM
      real*8 rstar,rstarP,rstarM,dtvold,dcvold,dtv,dcv,dtvP,dcvP,
     1     dtvM,dcvM
      real*8 elmat(6,6),dsmax,pert
      integer imap(6),ipFail,ipFailP,ipFailM,k1,k2,jj
      integer ifirst
      save ifirst
      data ifirst /0/

      if (ntens .lt. 4) then
         write(*,*) 'CDP UMAT: NTENS < 4 (plane stress) is not',
     $        ' supported by this UMAT. Use 3D solid or plane',
     $        ' strain/axisymmetric elements.'
         call xit
      end if

c ---------------------------------------------------------------------
c Build material constants + preprocessed hardening/damage tables into
c common blocks 
c ---------------------------------------------------------------------
      call cdp_buildtables(props,nprops,celent)

      if (ifirst .eq. 0) then
         call cdp_printinput()
         ifirst = 1
      end if

c ---------------------------------------------------------------------
c Map UMAT strain increment (order 11,22,33,12,13,23) to the internal
c 6-vector (order 11,22,33,12,23,13). imap(k) = UMAT index feeding
c internal component k (0 = none, cannot occur since ntens>=4 here).
c ---------------------------------------------------------------------
      do k1 = 1,6
         deps(k1) = 0.d0
         imap(k1) = 0
      end do
      deps(1) = dstran(1)
      deps(2) = dstran(2)
      deps(3) = dstran(3)
      imap(1) = 1
      imap(2) = 2
      imap(3) = 3
      deps(4) = dstran(4)
      imap(4) = 4
      if (ntens .eq. 6) then
         deps(5) = dstran(6)
         deps(6) = dstran(5)
         imap(5) = 6
         imap(6) = 5
      end if

c Old state (zero-pad if fewer than 27 state variables are declared)
      do k1 = 1,27
         stateO(k1) = 0.d0
      end do
      do k1 = 1,min(nstatv,27)
         stateO(k1) = statev(k1)
      end do
      dtvold = stateO(17)
      dcvold = stateO(18)

c ---------------------------------------------------------------------
c Stress update at this material point (plasticity + inviscid damage)
c ---------------------------------------------------------------------
      call cdp_umatupdate(deps,stateO,stateN,sigeff,ipFail)

      if (ipFail .eq. 1) then
c        Return mapping did not converge even after sub-incrementing:
c        ask ABAQUS for a smaller time increment and hand back the
c        undamaged elastic stiffness.
         pnewdt = 0.25d0
         call cdp_elasticmatrix(elmat)
         do k1 = 1,ntens
            do k2 = 1,ntens
               ddsdde(k1,k2) = 0.d0
            end do
         end do
         do k1 = 1,6
            do k2 = 1,6
               if (imap(k1).ne.0 .and. imap(k2).ne.0) then
                  ddsdde(imap(k1),imap(k2)) = elmat(k1,k2)
               end if
            end do
         end do
         return
      end if

      dtinv = stateN(9)
      dcinv = stateN(10)
      rstar = stateN(20)

c ---------------------------------------------------------------------
c Viscous regularization of the damage variables + combination into
c the nominal (damaged) stress.
c ---------------------------------------------------------------------
      call cdp_combine(sigeff,dtinv,dcinv,rstar,dtvold,dcvold,dtime,
     $     sig,dtv,dcv)

      stateN(17) = dtv
      stateN(18) = dcv
      stateN(19) = 1.d0 - (1.d0-(1.d0-wt*rstar)*dcv)*
     $     (1.d0-(1.d0-wc*(1.d0-rstar))*dtv)
      stateN(27) = 1.d0

c Store new state
      do k1 = 1,min(nstatv,27)
         statev(k1) = stateN(k1)
      end do

c Map internal stress back to UMAT ordering
      call cdp_mapsig(sig,sigUM)

c Specific work (approximate, informational only -- not split into
c elastic/plastic/damage parts)
      do k1 = 1,ntens
         sse = sse + 0.5d0*(stress(k1)+sigUM(k1))*dstran(k1)
      end do

c ---------------------------------------------------------------------
c Numerical consistent tangent: central-difference perturbation of
c each strain-increment component, re-running the full update
c (plasticity + damage + viscous combination) from the same old state.
c ---------------------------------------------------------------------
      dsmax = 0.d0
      do k1 = 1,ntens
         if (abs(dstran(k1)) .gt. dsmax) dsmax = abs(dstran(k1))
      end do
      pert = 1.0d-4*sigt0/ym
      if (pert .lt. 1.d-9) pert = 1.d-9
      if (pert .gt. 1.d-2*dsmax .and. dsmax .gt. 0.d0) then
         if (1.d-2*dsmax .gt. 1.d-9) pert = 1.d-2*dsmax
      end if

      do jj = 1,ntens
c        --- forward perturbation ---
         do k1 = 1,6
            depsP(k1) = deps(k1)
            if (imap(k1) .eq. jj) depsP(k1) = depsP(k1) + pert
         end do
         call cdp_umatupdate(depsP,stateO,stateP,sigeffP,ipFailP)
         if (ipFailP .eq. 1) goto 800
         dtinvP = stateP(9)
         dcinvP = stateP(10)
         rstarP = stateP(20)
         call cdp_combine(sigeffP,dtinvP,dcinvP,rstarP,dtvold,dcvold,
     $        dtime,sigP,dtvP,dcvP)
         call cdp_mapsig(sigP,sigUMP)
c        --- backward perturbation ---
         do k1 = 1,6
            depsP(k1) = deps(k1)
            if (imap(k1) .eq. jj) depsP(k1) = depsP(k1) - pert
         end do
         call cdp_umatupdate(depsP,stateO,stateP,sigeffM,ipFailM)
         if (ipFailM .eq. 1) goto 800
         dtinvM = stateP(9)
         dcinvM = stateP(10)
         rstarM = stateP(20)
         call cdp_combine(sigeffM,dtinvM,dcinvM,rstarM,dtvold,dcvold,
     $        dtime,sigM,dtvM,dcvM)
         call cdp_mapsig(sigM,sigUMM)
c        --- central difference ---
         do k1 = 1,ntens
            ddsdde(k1,jj) = (sigUMP(k1)-sigUMM(k1))/(2.d0*pert)
         end do
      end do
      goto 900

c ---------------------------------------------------------------------
c Fallback tangent: secant (damaged) elastic stiffness, used only if a
c perturbed state itself failed to converge (rare; the base state
c already succeeded above).
c ---------------------------------------------------------------------
 800  continue
      call cdp_elasticmatrix(elmat)
      dsmax = 1.d0 - max(dtv,dcv)
      if (dsmax .lt. 1.d-4) dsmax = 1.d-4
      do k1 = 1,6
         do k2 = 1,6
            if (imap(k1).ne.0 .and. imap(k2).ne.0) then
               ddsdde(imap(k1),imap(k2)) = dsmax*elmat(k1,k2)
            end if
         end do
      end do

 900  continue

c Update stress
      do k1 = 1,ntens
         stress(k1) = sigUM(k1)
      end do

      return
      end
c ======================================================================
c Parse props() into material constants + preprocessed (kappa,sigma_bar)
c hardening tables, stored in common blocks cdpc/cdptab/cdpn.
c ======================================================================
      subroutine cdp_buildtables(props,nprops,celent)

      integer nprops
      real*8 props(nprops),celent

      real*8 ym,pr,psi,tanpsi,ecc,alpha,gammap,fb0fc0,rkc,visc,wt,wc,
     1     sigt0,bulkK,shearG
      common/cdpc/ym,pr,psi,tanpsi,ecc,alpha,gammap,fb0fc0,rkc,visc,
     1     wt,wc,sigt0,bulkK,shearG

      real*8 kttab(60),sbttab(60),ecktab(60),
     1     kctab(60),sbctab(60),ecktabc(60),
     2     dtx(60),dty(60),dcx(60),dcy(60)
      common/cdptab/kttab,sbttab,ecktab,kctab,sbctab,ecktabc,
     1     dtx,dty,dcx,dcy
      integer ntt,ntc,ntdt,ntdc
      common/cdpn/ntt,ntc,ntdt,ntdc

      real*8 psideg,fb0fc0in,rkcin,dval,pi
      real*8 sigtraw(60),ecktraw(60),sigcraw(60),ecktrawc(60)
      integer i,ioff,ntdt0,ntdc0,itype

      pi = 3.14159265358979d0

      ym = props(1)
      pr = props(2)
      psideg = props(3)
      if (psideg .lt. 0.d0) psideg = 0.d0
      if (psideg .gt. 89.9d0) psideg = 89.9d0
      ecc = props(4)
      if (ecc .le. 0.d0 .or. ecc .gt. 1.d0) ecc = 0.1d0
      fb0fc0in = props(5)
      if (fb0fc0in .le. 1.d0) fb0fc0in = 1.16d0
      rkcin = props(6)
      if (rkcin .le. 0.5d0 .or. rkcin .gt. 1.d0) rkcin = 2.d0/3.d0
      visc = props(7)
      if (visc .lt. 0.d0) visc = 0.d0
      wt = props(8)
      if (wt .lt. 0.d0 .or. wt .gt. 1.d0) wt = 0.d0
      wc = props(9)
      if (wc .lt. 0.d0 .or. wc .gt. 1.d0) wc = 1.d0
      itype = int(props(10) + 0.5d0)
      ntt = int(props(11) + 0.5d0)
      ntc = int(props(12) + 0.5d0)
      ntdt0 = int(props(13) + 0.5d0)
      ntdc0 = int(props(14) + 0.5d0)

c ---- MAXTAB = 60: every kttab/sbttab/.../dcy array above is fixed at
c      this size (Fortran 77 fixed-size arrays, shared via common
c      /cdptab/). A table row count above 60 must ABORT here rather
c      than be silently clamped: clamping ntt/ntc/etc without also
c      re-deriving the ioff step below would misalign every table read
c      after the truncated one against the real props() layout Abaqus
c      built from constants=. If you need more than 60 rows per table,
c      raise every "(60)" in this file's real*8 xxxtab/xxxraw
c      declarations (5 places) and the 4 limits just below, together,
c      and recompile.
      if (ntt .lt. 1) ntt = 1
      if (ntt .gt. 60) then
         write(*,*) 'CDP UMAT: tension stiffening table has',ntt,
     $        ' rows, exceeds the compiled-in limit of 60. Increase',
     $        ' MAXTAB (every real*8 xxxtab(60) declaration in this',
     $        ' file) and recompile, or use a coarser table.'
         call xit
      end if
      if (ntc .lt. 1) ntc = 1
      if (ntc .gt. 60) then
         write(*,*) 'CDP UMAT: compression hardening table has',ntc,
     $        ' rows, exceeds the compiled-in limit of 60. Increase',
     $        ' MAXTAB (every real*8 xxxtab(60) declaration in this',
     $        ' file) and recompile, or use a coarser table.'
         call xit
      end if
      if (ntdt0 .lt. 0) ntdt0 = 0
      if (ntdt0 .gt. 60) then
         write(*,*) 'CDP UMAT: tension damage table has',ntdt0,
     $        ' rows, exceeds the compiled-in limit of 60. Increase',
     $        ' MAXTAB (every real*8 xxxtab(60) declaration in this',
     $        ' file) and recompile, or use a coarser table.'
         call xit
      end if
      if (ntdc0 .lt. 0) ntdc0 = 0
      if (ntdc0 .gt. 60) then
         write(*,*) 'CDP UMAT: compression damage table has',ntdc0,
     $        ' rows, exceeds the compiled-in limit of 60. Increase',
     $        ' MAXTAB (every real*8 xxxtab(60) declaration in this',
     $        ' file) and recompile, or use a coarser table.'
         call xit
      end if

      bulkK = ym/(3.d0*(1.d0-2.d0*pr))
      shearG = ym/(2.d0*(1.d0+pr))
      psi = psideg*pi/180.d0
      tanpsi = tan(psi)
      alpha = (fb0fc0in-1.d0)/(2.d0*fb0fc0in-1.d0)
      gammap = 3.d0*(1.d0-rkcin)/(2.d0*rkcin-1.d0)

c ---- read the four tables from props ----
      ioff = 14
      do i = 1,ntt
         sigtraw(i) = props(ioff+2*i-1)
         ecktraw(i) = props(ioff+2*i)
      end do
      ioff = ioff + 2*ntt
      do i = 1,ntc
         sigcraw(i) = props(ioff+2*i-1)
         ecktrawc(i) = props(ioff+2*i)
      end do
      ioff = ioff + 2*ntc
      if (ntdt0 .gt. 0) then
         do i = 1,ntdt0
            dty(i) = props(ioff+2*i-1)
            dtx(i) = props(ioff+2*i)
         end do
         ntdt = ntdt0
      else
         ntdt = 2
         dtx(1) = 0.d0
         dty(1) = 0.d0
         dtx(2) = 1.d10
         dty(2) = 0.d0
      end if
      ioff = ioff + 2*ntdt0
      if (ntdc0 .gt. 0) then
         do i = 1,ntdc0
            dcy(i) = props(ioff+2*i-1)
            dcx(i) = props(ioff+2*i)
         end do
         ntdc = ntdc0
      else
         ntdc = 2
         dcx(1) = 0.d0
         dcy(1) = 0.d0
         dcx(2) = 1.d10
         dcy(2) = 0.d0
      end if

c ---- TYPE=DISPLACEMENT conversion of the tension side, using the
c      element characteristic length (crack-band regularization) ----
      if (itype .eq. 1 .and. celent .gt. 1.d-12) then
         do i = 1,ntt
            ecktraw(i) = ecktraw(i)/celent
         end do
         if (ntdt0 .gt. 0) then
            do i = 1,ntdt0
               dtx(i) = dtx(i)/celent
            end do
         end if
      end if

c ---- build (kappa, sigma_bar) tables, converting cracking/inelastic
c      strain to plastic strain via the damage at that same strain
c      point: kappa = eck - (d/(1-d))*(sigma/E)  ----
      do i = 1,ntt
         call cdp_interp(ecktraw(i),dtx,dty,ntdt,dval)
         if (dval .lt. 0.d0) dval = 0.d0
         if (dval .gt. 0.99d0) dval = 0.99d0
         sbttab(i) = sigtraw(i)/(1.d0-dval)
         kttab(i) = ecktraw(i) - (dval/(1.d0-dval))*(sigtraw(i)/ym)
         ecktab(i) = ecktraw(i)
      end do
      do i = 2,ntt
         if (kttab(i) .lt. kttab(i-1)) kttab(i) = kttab(i-1)
      end do

      do i = 1,ntc
         call cdp_interp(ecktrawc(i),dcx,dcy,ntdc,dval)
         if (dval .lt. 0.d0) dval = 0.d0
         if (dval .gt. 0.99d0) dval = 0.99d0
         sbctab(i) = sigcraw(i)/(1.d0-dval)
         kctab(i) = ecktrawc(i) - (dval/(1.d0-dval))*(sigcraw(i)/ym)
         ecktabc(i) = ecktrawc(i)
      end do
      do i = 2,ntc
         if (kctab(i) .lt. kctab(i-1)) kctab(i) = kctab(i-1)
      end do

      sigt0 = sigtraw(1)
      if (sigt0 .le. 0.d0) sigt0 = 1.d-3*ym

      return
      end
c ======================================================================
c Piecewise-linear interpolation with flat extrapolation beyond the
c table's range
c xtab/ytab must be non-decreasing in x for i=1..n.
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
c Piecewise-linear interpolation, ALSO returning the local segment
c slope dydx (flat extrapolation outside the table -> slope 0 there).
c Same value convention as cdp_interp so the two never disagree; 
c used by cdp_rmresid2 for the analytical Jacobian.
c ======================================================================
      subroutine cdp_interp2(x,xtab,ytab,n,y,dydx)
      real*8 x,xtab(*),ytab(*),y,dydx
      integer n,i

      if (n .le. 1) then
         y = ytab(1)
         dydx = 0.d0
         return
      end if
      if (x .le. xtab(1)) then
         y = ytab(1)
         dydx = 0.d0
         return
      end if
      if (x .ge. xtab(n)) then
         y = ytab(n)
         dydx = 0.d0
         return
      end if
      do i = 1,n-1
         if (x .ge. xtab(i) .and. x .le. xtab(i+1)) then
            if (xtab(i+1) .eq. xtab(i)) then
               y = ytab(i)
               dydx = 0.d0
            else
               dydx = (ytab(i+1)-ytab(i))/(xtab(i+1)-xtab(i))
               y = ytab(i) + dydx*(x-xtab(i))
            end if
            return
         end if
      end do
      y = ytab(n)
      dydx = 0.d0
      return
      end


c ======================================================================
c Lee-Fenves yield function:
c   F = [ q - 3*alpha*p + beta*<smax> - gamma*<-smax> ]/(1-alpha)
c       - sigma_bar_c(kappac)
c p = -1/3 trace(sigma_eff)  (positive in compression)
c q = Mises equivalent of sigma_eff (>=0)
c smax = algebraically largest principal effective stress
c ======================================================================
      subroutine cdp_yieldf(p,q,smax,kappat,kappac,fval)
      real*8 p,q,smax,kappat,kappac,fval

      real*8 ym,pr,psi,tanpsi,ecc,alpha,gammap,fb0fc0,rkc,visc,wt,wc,
     1     sigt0,bulkK,shearG
      common/cdpc/ym,pr,psi,tanpsi,ecc,alpha,gammap,fb0fc0,rkc,visc,
     1     wt,wc,sigt0,bulkK,shearG

      real*8 kttab(60),sbttab(60),ecktab(60),
     1     kctab(60),sbctab(60),ecktabc(60),
     2     dtx(60),dty(60),dcx(60),dcy(60)
      common/cdptab/kttab,sbttab,ecktab,kctab,sbctab,ecktabc,
     1     dtx,dty,dcx,dcy
      integer ntt,ntc,ntdt,ntdc
      common/cdpn/ntt,ntc,ntdt,ntdc

      real*8 sbt,sbc,beta,macp,macm

      call cdp_interp(kappat,kttab,sbttab,ntt,sbt)
      call cdp_interp(kappac,kctab,sbctab,ntc,sbc)
      if (sbt .lt. 1.d-8) sbt = 1.d-8

      beta = sbc/sbt*(1.d0-alpha) - (1.d0+alpha)
      macp = max(smax,0.d0)
      macm = max(-smax,0.d0)

      fval = (q - 3.d0*alpha*p + beta*macp - gammap*macm)/
     $     (1.d0-alpha) - sbc

      return
      end


c ======================================================================
c Multiaxial stress-weight factor r(sigma_hat) = sum<si>/sum|si|,
c from three principal (effective) stresses, any order.
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
c Cyclic Jacobi eigenvalue algorithm for a symmetric 3x3 matrix.
c Eigenvalues returned UNSORTED in evals(3); evecs columns are the
c corresponding (orthonormal) eigenvectors. Call cdp_sort3desc after
c this to get descending order.
c ======================================================================
      subroutine cdp_jacobi3(ain,evals,evecs)
      real*8 ain(3,3),evals(3),evecs(3,3)
      real*8 a(3,3),v(3,3)
      real*8 off,theta,t,c,s,tau,apq,app,aqq,akp,akq,vkp,vkq
      integer i,j,k,p,q,sweep,maxsweep
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

      maxsweep = 60
      do sweep = 1,maxsweep
         off = abs(a(1,2))+abs(a(1,3))+abs(a(2,3))
         if (off .lt. 1.d-13) goto 200
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


c ======================================================================
c Sort the 3 eigenpairs from cdp_jacobi3 into descending eigenvalue
c order (evals(1) >= evals(2) >= evals(3)), permuting evecs columns
c to match.
c ======================================================================
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
c Material-point update: elastic predictor, eigen-decomposition,
c yield check, implicit (backward-Euler) return mapping with adaptive
c strain sub-incrementation on Newton non-convergence.
c
c Input : deps(6)      strain increment, internal ordering
c         stateOld(27) state at start of increment
c Output: stateNew(1),(2)     kappat, kappac
c         stateNew(3:8)       plastic strain
c         stateNew(9),(10)    dt, dc (inviscid damage)
c         stateNew(11:16)     total strain
c         stateNew(20)        r* (multiaxial stress-weight factor)
c         stateNew(21:26)     effective (undamaged) stress
c         sigeff(6)           effective (undamaged) stress (same as
c                              stateNew(21:26), returned separately
c                              for convenience)
c         ipFail              =1 if the return mapping failed to
c                              converge even after sub-incrementing
c                              (caller should request a smaller time
c                              increment via PNEWDT)
c
c Common blocks cdpc/cdptab/cdpn (see cdp_buildtables) must already be
c populated for the CURRENT material before this is called.
c ======================================================================
      subroutine cdp_umatupdate(deps,stateOld,stateNew,sigeff,ipFail)

      real*8 deps(6),stateOld(27),stateNew(27),sigeff(6)
      integer ipFail

      real*8 epsplold(6),epstotold(6),epstotcur(6),epsel(6)
      real*8 sigtrialv(6),Amat(3,3),evals(3),evecs(3,3)
      real*8 d1,d2,d3,ptrial,qtrial,mean
      real*8 kappatcur,kappaccur,kappatold,kappacold,kappatn,kappacn
      real*8 ftrial,smaxtrial,ftol
      real*8 dlam,q,p
      real*8 epsplcur(6),epsplnext(6)
      real*8 sigeffprinc(3),nmat(3,3),dtval,dcval,rstar
      real*8 fracdone,fracstep,fracmin
      integer i,j,converged,attempt,maxattempt

      do i = 1,6
         epsplold(i) = stateOld(i+2)
         epstotold(i) = stateOld(i+10)
      end do
      kappatold = stateOld(1)
      kappacold = stateOld(2)

      ipFail = 0

      do i = 1,6
         epsplcur(i) = epsplold(i)
      end do
      kappatcur = kappatold
      kappaccur = kappacold

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

      call cdp_stressfromstrain(epsel,sigtrialv)

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

      call cdp_yieldf(ptrial,qtrial,smaxtrial,kappatcur,kappaccur,
     $     ftrial)

      if (ftrial .le. ftol*max(abs(evals(1)),abs(evals(3)),1.d0)) then
         sigeffprinc(1) = evals(1)
         sigeffprinc(2) = evals(2)
         sigeffprinc(3) = evals(3)
         kappatn = kappatcur
         kappacn = kappaccur
         converged = 1
      else
         call cdp_returnmap(ptrial,qtrial,d1,d2,d3,kappatcur,
     $        kappaccur,dlam,q,p,kappatn,kappacn,converged)
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

c     --- accept this sub-step: reconstruct the stress tensor from
c         principal values + trial eigenvectors, update state ---
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

      call cdp_strainfromstress(sigeff,epsel)
      do i = 1,6
         epsplnext(i) = epstotcur(i) - epsel(i)
      end do
      do i = 1,6
         epsplcur(i) = epsplnext(i)
      end do
      kappatcur = kappatn
      kappaccur = kappacn

      fracdone = fracdone + fracstep
      if (fracdone .lt. 0.999999999d0) then
         goto 100
      end if

c ---- final damage evaluation + r* from the converged end state ----
      call cdp_dmgt(kappatcur,dtval)
      call cdp_dmgc(kappaccur,dcval)
      call cdp_rfactor(sigeffprinc(1),sigeffprinc(2),sigeffprinc(3),
     $     rstar)

      stateNew(1) = kappatcur
      stateNew(2) = kappaccur
      do i = 1,6
         stateNew(i+2) = epsplcur(i)
      end do
      stateNew(9) = dtval
      stateNew(10) = dcval
      do i = 1,6
         stateNew(i+10) = epstotcur(i)
      end do
      stateNew(20) = rstar
      do i = 1,6
         stateNew(i+20) = sigeff(i)
      end do

      return
      end
c ======================================================================
c Implicit (backward-Euler) return mapping: solves the 2x2 system
c   R1(dlam,q) = q*(1 + 3*G*dlam/Gt(q)) - qtrial      = 0
c   R2(dlam,q) = F(p(dlam), q, smax(dlam,q), kt, kc)  = 0
c with p(dlam) = ptrial + K*dlam*tan(psi),
c      Gt(q)   = sqrt((ecc*sigt0*tan(psi))^2 + q^2),
c      smax, kt, kc as in cdp_rmresid. (Derivation notes: because the
c flow potential G depends only on p and q, not the third invariant,
c the plastic flow direction is coaxial with the trial effective
c stress -- the deviatoric principal stresses simply scale by q/qtrial
c and p shifts by K*dlam*tan(psi), reducing the return map to this
c 2-unknown system.)
c
c Phase 1: Newton with the ANALYTICAL Jacobian from cdp_rmresid2
c (closed-form, branch-consistent -- see that routine's header) instead
c of a finite-difference one, damped by a real Armijo backtracking line
c search that requires ||R|| to actually decrease (not merely
c dlam,q>=0 as before).
c Phase 2 (only if phase 1 fails to converge or hits a singular
c Jacobian): R1(dlam,q)=0 has a UNIQUE q>=0 root for any fixed dlam>=0
c (cdp_solve_q's header proves dR1/dq>=1 everywhere), so the 2D problem
c reduces to a 1D bisection on dlam of g(dlam)=R2(dlam,q_of(dlam)).
c Bisection needs only a sign-changing bracket (found by doubling) and
c then cannot fail to converge, unlike a 2D Newton step -- the
c equivalent, for this yield surface, of cdpm2umat.f's dedicated
c bisection-based vertex/degenerate-case fallback.
c
c ======================================================================
      subroutine cdp_returnmap(ptrial,qtrial,d1,d2,d3,kappatold,
     $     kappacold,dlam,q,p,kappatnew,kappacnew,converged)

      real*8 ptrial,qtrial,d1,d2,d3,kappatold,kappacold
      real*8 dlam,q,p,kappatnew,kappacnew
      integer converged

      real*8 x(2),xn(2),R(2),Rn(2),Jac(2,2),dx(2)
      real*8 detJ,lam,normR,normRn,tol,kt,kc,smax
      real*8 glo,ghi,gm,lo,hi,mid,qlo,qhi,qm
      integer it,maxit,ls,tries

      x(1) = 1.d-10
      if (qtrial .gt. 0.d0) then
         x(2) = 0.5d0*qtrial
      else
         x(2) = 0.d0
      end if

      maxit = 60
      tol = max(abs(ptrial),abs(qtrial),1.d0)*1.d-11

c ---- Phase 1: analytical-Jacobian Newton, Armijo backtracking ----
      do it = 1,maxit
         call cdp_rmresid2(x,ptrial,qtrial,d1,d2,d3,kappatold,
     $        kappacold,R,Jac,kt,kc,smax,p)
         normR = sqrt(R(1)*R(1)+R(2)*R(2))
         if (normR .lt. tol) then
            converged = 1
            dlam = x(1)
            q = x(2)
            kappatnew = kt
            kappacnew = kc
            return
         end if

         detJ = Jac(1,1)*Jac(2,2)-Jac(1,2)*Jac(2,1)
         if (abs(detJ) .lt. 1.d-30) goto 200
         dx(1) = -( Jac(2,2)*R(1)-Jac(1,2)*R(2))/detJ
         dx(2) = -(-Jac(2,1)*R(1)+Jac(1,1)*R(2))/detJ

         lam = 1.d0
         do ls = 1,30
            xn(1) = max(x(1)+lam*dx(1),0.d0)
            xn(2) = max(x(2)+lam*dx(2),0.d0)
            call cdp_rmresid(xn,ptrial,qtrial,d1,d2,d3,kappatold,
     $           kappacold,Rn,kt,kc,smax,p)
            normRn = sqrt(Rn(1)*Rn(1)+Rn(2)*Rn(2))
            if (normRn .lt. (1.d0-1.d-4*lam)*normR) goto 100
            lam = lam*0.5d0
         end do
         goto 200
 100     continue
         x(1) = xn(1)
         x(2) = xn(2)
      end do

c ---- Phase 2: dimension-reduced bisection fallback (see header) ----
 200  continue
      lo = 0.d0
      call cdp_solve_q(lo,ptrial,qtrial,d1,d2,d3,kappatold,
     $     kappacold,qlo)
      x(1) = lo
      x(2) = qlo
      call cdp_rmresid(x,ptrial,qtrial,d1,d2,d3,kappatold,
     $     kappacold,R,kt,kc,smax,p)
      glo = R(2)

      hi = 1.d-6
      tries = 0
 210  continue
      call cdp_solve_q(hi,ptrial,qtrial,d1,d2,d3,kappatold,
     $     kappacold,qhi)
      x(1) = hi
      x(2) = qhi
      call cdp_rmresid(x,ptrial,qtrial,d1,d2,d3,kappatold,
     $     kappacold,R,kt,kc,smax,p)
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
         call cdp_solve_q(mid,ptrial,qtrial,d1,d2,d3,kappatold,
     $        kappacold,qm)
         x(1) = mid
         x(2) = qm
         call cdp_rmresid(x,ptrial,qtrial,d1,d2,d3,kappatold,
     $        kappacold,R,kt,kc,smax,p)
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
      kappatnew = kt
      kappacnew = kc
      converged = 1

      return
      end


c ======================================================================
c Residual evaluator for cdp_returnmap. Also returns the implied
c hardening variables kt,kc, the largest principal effective stress
c smax (needed by the yield function), and the volumetric effective
c stress p at x=(dlam,q).
c ======================================================================
      subroutine cdp_rmresid(x,ptrial,qtrial,d1,d2,d3,kappatold,
     $     kappacold,R,kt,kc,smax,p)

      real*8 x(2),ptrial,qtrial,d1,d2,d3,kappatold,kappacold
      real*8 R(2),kt,kc,smax,p

      real*8 ym,pr,psi,tanpsi,ecc,alpha,gammap,fb0fc0,rkc,visc,wt,wc,
     1     sigt0,bulkK,shearG
      common/cdpc/ym,pr,psi,tanpsi,ecc,alpha,gammap,fb0fc0,rkc,visc,
     1     wt,wc,sigt0,bulkK,shearG

      real*8 dlam,q,Gt,ratio,smid,smin
      real*8 depsmax,depsmin,rfac,dkt,dkc,fval

      dlam = x(1)
      q = x(2)

      p = ptrial + bulkK*dlam*tanpsi
      Gt = sqrt((ecc*sigt0*tanpsi)**2 + q*q)
      if (Gt .lt. 1.d-12) Gt = 1.d-12

      R(1) = q*(1.d0+3.d0*shearG*dlam/Gt) - qtrial

      if (qtrial .gt. 1.d-10) then
         ratio = q/qtrial
      else
         ratio = 0.d0
      end if
      smax = -p + ratio*d1
      smid = -p + ratio*d2
      smin = -p + ratio*d3

      if (qtrial .gt. 1.d-10) then
         depsmax = dlam*(3.d0*ratio*d1/(2.d0*Gt) + tanpsi/3.d0)
         depsmin = dlam*(3.d0*ratio*d3/(2.d0*Gt) + tanpsi/3.d0)
      else
         depsmax = dlam*tanpsi/3.d0
         depsmin = dlam*tanpsi/3.d0
      end if

      call cdp_rfactor(smax,smid,smin,rfac)

      dkt = rfac*depsmax
      if (dkt .lt. 0.d0) dkt = 0.d0
      dkc = -(1.d0-rfac)*depsmin
      if (dkc .lt. 0.d0) dkc = 0.d0

      kt = kappatold+dkt
      kc = kappacold+dkc

      call cdp_yieldf(p,q,smax,kt,kc,fval)
      R(2) = fval

      return
      end
c ======================================================================
c Analytical (closed-form) residual + Jacobian for the SAME 2x2 system
c as cdp_rmresid (R1: consistency of q with the trial deviatoric
c stress and the flow rule; R2: the Lee-Fenves yield condition F=0).
c ======================================================================
      subroutine cdp_rmresid2(x,ptrial,qtrial,d1,d2,d3,kappatold,
     $     kappacold,R,Jac,kt,kc,smax,p)

      real*8 x(2),ptrial,qtrial,d1,d2,d3,kappatold,kappacold
      real*8 R(2),Jac(2,2),kt,kc,smax,p

      real*8 ym,pr,psi,tanpsi,ecc,alpha,gammap,fb0fc0,rkc,visc,wt,wc,
     1     sigt0,bulkK,shearG
      common/cdpc/ym,pr,psi,tanpsi,ecc,alpha,gammap,fb0fc0,rkc,visc,
     1     wt,wc,sigt0,bulkK,shearG

      real*8 kttab(60),sbttab(60),ecktab(60),
     1     kctab(60),sbctab(60),ecktabc(60),
     2     dtx(60),dty(60),dcx(60),dcy(60)
      common/cdptab/kttab,sbttab,ecktab,kctab,sbctab,ecktabc,
     1     dtx,dty,dcx,dcy
      integer ntt,ntc,ntdt,ntdc
      common/cdpn/ntt,ntc,ntdt,ntdc

      real*8 dlam,q,Gt,asq,ratio,smid,smin
      real*8 dp_ddlam,dR1_ddlam,dR1_dq,dratio_dq
      real*8 dsmax_ddlam,dsmax_dq,dsmid_ddlam,dsmid_dq
      real*8 dsmin_ddlam,dsmin_dq
      real*8 depsmax,depsmin,Dmaxv,Dminv
      real*8 ddepsmax_ddlam,ddepsmax_dq,ddepsmin_ddlam,ddepsmin_dq
      real*8 dDmax_dq,dDmin_dq
      real*8 rfac,num,den,dnum_ddlam,dnum_dq,dden_ddlam,dden_dq
      real*8 drfac_ddlam,drfac_dq
      real*8 hdkt,hdkc,ddkt_ddlam,ddkt_dq,ddkc_ddlam,ddkc_dq
      real*8 sbt,sbc,slopet,slopec
      real*8 dsbt_ddlam,dsbt_dq,dsbc_ddlam,dsbc_dq
      real*8 beta,dbeta_ddlam,dbeta_dq
      real*8 macp,macm,dmacp_ddlam,dmacp_dq,dmacm_ddlam,dmacm_dq

      dlam = x(1)
      q = x(2)

      p = ptrial + bulkK*dlam*tanpsi
      dp_ddlam = bulkK*tanpsi

      asq = (ecc*sigt0*tanpsi)**2
      Gt = sqrt(asq + q*q)
      if (Gt .lt. 1.d-12) Gt = 1.d-12

      R(1) = q*(1.d0+3.d0*shearG*dlam/Gt) - qtrial
      dR1_ddlam = 3.d0*shearG*q/Gt
      dR1_dq = 1.d0 + 3.d0*shearG*dlam*asq/Gt**3

      if (qtrial .gt. 1.d-10) then
         ratio = q/qtrial
         dratio_dq = 1.d0/qtrial
      else
         ratio = 0.d0
         dratio_dq = 0.d0
      end if

      smax = -p + ratio*d1
      smid = -p + ratio*d2
      smin = -p + ratio*d3
      dsmax_ddlam = -dp_ddlam
      dsmid_ddlam = -dp_ddlam
      dsmin_ddlam = -dp_ddlam
      dsmax_dq = d1*dratio_dq
      dsmid_dq = d2*dratio_dq
      dsmin_dq = d3*dratio_dq

      if (qtrial .gt. 1.d-10) then
         Dmaxv = 3.d0*ratio*d1/(2.d0*Gt) + tanpsi/3.d0
         Dminv = 3.d0*ratio*d3/(2.d0*Gt) + tanpsi/3.d0
         depsmax = dlam*Dmaxv
         depsmin = dlam*Dminv
         ddepsmax_ddlam = Dmaxv
         ddepsmin_ddlam = Dminv
         dDmax_dq = 3.d0*d1*asq/(2.d0*qtrial*Gt**3)
         dDmin_dq = 3.d0*d3*asq/(2.d0*qtrial*Gt**3)
         ddepsmax_dq = dlam*dDmax_dq
         ddepsmin_dq = dlam*dDmin_dq
      else
         depsmax = dlam*tanpsi/3.d0
         depsmin = dlam*tanpsi/3.d0
         ddepsmax_ddlam = tanpsi/3.d0
         ddepsmin_ddlam = tanpsi/3.d0
         ddepsmax_dq = 0.d0
         ddepsmin_dq = 0.d0
      end if

      den = abs(smax)+abs(smid)+abs(smin)
      if (den .lt. 1.d-10) then
         rfac = 0.d0
         drfac_ddlam = 0.d0
         drfac_dq = 0.d0
      else
         num = max(smax,0.d0)+max(smid,0.d0)+max(smin,0.d0)
         rfac = num/den
         dnum_ddlam = 0.d0
         dnum_dq = 0.d0
         dden_ddlam = 0.d0
         dden_dq = 0.d0
         if (smax .gt. 0.d0) then
            dnum_ddlam = dnum_ddlam + dsmax_ddlam
            dnum_dq = dnum_dq + dsmax_dq
            dden_ddlam = dden_ddlam + dsmax_ddlam
            dden_dq = dden_dq + dsmax_dq
         else if (smax .lt. 0.d0) then
            dden_ddlam = dden_ddlam - dsmax_ddlam
            dden_dq = dden_dq - dsmax_dq
         end if
         if (smid .gt. 0.d0) then
            dnum_ddlam = dnum_ddlam + dsmid_ddlam
            dnum_dq = dnum_dq + dsmid_dq
            dden_ddlam = dden_ddlam + dsmid_ddlam
            dden_dq = dden_dq + dsmid_dq
         else if (smid .lt. 0.d0) then
            dden_ddlam = dden_ddlam - dsmid_ddlam
            dden_dq = dden_dq - dsmid_dq
         end if
         if (smin .gt. 0.d0) then
            dnum_ddlam = dnum_ddlam + dsmin_ddlam
            dnum_dq = dnum_dq + dsmin_dq
            dden_ddlam = dden_ddlam + dsmin_ddlam
            dden_dq = dden_dq + dsmin_dq
         else if (smin .lt. 0.d0) then
            dden_ddlam = dden_ddlam - dsmin_ddlam
            dden_dq = dden_dq - dsmin_dq
         end if
         drfac_ddlam = (dnum_ddlam*den-num*dden_ddlam)/den**2
         drfac_dq = (dnum_dq*den-num*dden_dq)/den**2
      end if

      hdkt = rfac*depsmax
      if (hdkt .gt. 0.d0) then
         ddkt_ddlam = drfac_ddlam*depsmax + rfac*ddepsmax_ddlam
         ddkt_dq = drfac_dq*depsmax + rfac*ddepsmax_dq
      else
         hdkt = 0.d0
         ddkt_ddlam = 0.d0
         ddkt_dq = 0.d0
      end if

      hdkc = (rfac-1.d0)*depsmin
      if (hdkc .gt. 0.d0) then
         ddkc_ddlam = drfac_ddlam*depsmin + (rfac-1.d0)*
     $        ddepsmin_ddlam
         ddkc_dq = drfac_dq*depsmin + (rfac-1.d0)*ddepsmin_dq
      else
         hdkc = 0.d0
         ddkc_ddlam = 0.d0
         ddkc_dq = 0.d0
      end if

      kt = kappatold + hdkt
      kc = kappacold + hdkc

      call cdp_interp2(kt,kttab,sbttab,ntt,sbt,slopet)
      call cdp_interp2(kc,kctab,sbctab,ntc,sbc,slopec)
      if (sbt .lt. 1.d-8) sbt = 1.d-8
      dsbt_ddlam = slopet*ddkt_ddlam
      dsbt_dq = slopet*ddkt_dq
      dsbc_ddlam = slopec*ddkc_ddlam
      dsbc_dq = slopec*ddkc_dq

      beta = sbc/sbt*(1.d0-alpha) - (1.d0+alpha)
      dbeta_ddlam = (1.d0-alpha)*(dsbc_ddlam*sbt-sbc*dsbt_ddlam)/
     $     sbt**2
      dbeta_dq = (1.d0-alpha)*(dsbc_dq*sbt-sbc*dsbt_dq)/sbt**2

      macp = max(smax,0.d0)
      macm = max(-smax,0.d0)
      if (smax .gt. 0.d0) then
         dmacp_ddlam = dsmax_ddlam
         dmacp_dq = dsmax_dq
         dmacm_ddlam = 0.d0
         dmacm_dq = 0.d0
      else if (smax .lt. 0.d0) then
         dmacp_ddlam = 0.d0
         dmacp_dq = 0.d0
         dmacm_ddlam = -dsmax_ddlam
         dmacm_dq = -dsmax_dq
      else
         dmacp_ddlam = 0.d0
         dmacp_dq = 0.d0
         dmacm_ddlam = 0.d0
         dmacm_dq = 0.d0
      end if

      R(2) = (q - 3.d0*alpha*p + beta*macp - gammap*macm)/
     $     (1.d0-alpha) - sbc

      Jac(1,1) = dR1_ddlam
      Jac(1,2) = dR1_dq
      Jac(2,1) = (-3.d0*alpha*dp_ddlam + dbeta_ddlam*macp +
     $     beta*dmacp_ddlam - gammap*dmacm_ddlam)/(1.d0-alpha) -
     $     dsbc_ddlam
      Jac(2,2) = (1.d0 + dbeta_dq*macp + beta*dmacp_dq -
     $     gammap*dmacm_dq)/(1.d0-alpha) - dsbc_dq

      return
      end
c ======================================================================
c Finds the unique q>=0 root of R1(dlam,q)=0 for a FIXED dlam>=0.
c dR1/dq = 1 + 3*shearG*dlam*(ecc*sigt0*tanpsi)^2/Gt(q)^3 >= 1 always
c (every term is a product/ratio of non-negative quantities), so R1 is
c strictly increasing in q on q>=0 -- bisection on q therefore cannot
c fail once a bracket is found (unconditionally safe, unlike a 2D
c Newton step). Used only by cdp_returnmap's fallback phase.
c ======================================================================
      subroutine cdp_solve_q(dlam,ptrial,qtrial,d1,d2,d3,kappatold,
     $     kappacold,qsol)

      real*8 dlam,ptrial,qtrial,d1,d2,d3,kappatold,kappacold,qsol

      real*8 xx(2),R(2),kt,kc,smax,p
      real*8 lo,hi,flo,fhi,mid,fm
      integer tries,it

      lo = 0.d0
      xx(1) = dlam
      xx(2) = lo
      call cdp_rmresid(xx,ptrial,qtrial,d1,d2,d3,kappatold,
     $     kappacold,R,kt,kc,smax,p)
      flo = R(1)
      if (abs(flo) .lt. 1.d-13) then
         qsol = lo
         return
      end if

      hi = max(qtrial,1.d0)
      xx(2) = hi
      call cdp_rmresid(xx,ptrial,qtrial,d1,d2,d3,kappatold,
     $     kappacold,R,kt,kc,smax,p)
      fhi = R(1)
      tries = 0
 10   continue
      if (flo*fhi .gt. 0.d0 .and. tries .lt. 60) then
         hi = hi*2.d0
         xx(2) = hi
         call cdp_rmresid(xx,ptrial,qtrial,d1,d2,d3,kappatold,
     $        kappacold,R,kt,kc,smax,p)
         fhi = R(1)
         tries = tries+1
         goto 10
      end if

      do it = 1,100
         mid = 0.5d0*(lo+hi)
         xx(2) = mid
         call cdp_rmresid(xx,ptrial,qtrial,d1,d2,d3,kappatold,
     $        kappacold,R,kt,kc,smax,p)
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
c Damage lookups: kappa -> cracking/inelastic strain (via the
c preprocessed hardening table) -> damage (via the user damage table).
c ======================================================================
      subroutine cdp_dmgt(kappat,dtval)
      real*8 kappat,dtval

      real*8 kttab(60),sbttab(60),ecktab(60),
     1     kctab(60),sbctab(60),ecktabc(60),
     2     dtx(60),dty(60),dcx(60),dcy(60)
      common/cdptab/kttab,sbttab,ecktab,kctab,sbctab,ecktabc,
     1     dtx,dty,dcx,dcy
      integer ntt,ntc,ntdt,ntdc
      common/cdpn/ntt,ntc,ntdt,ntdc

      real*8 eck
      call cdp_interp(kappat,kttab,ecktab,ntt,eck)
      call cdp_interp(eck,dtx,dty,ntdt,dtval)
      if (dtval .lt. 0.d0) dtval = 0.d0
      if (dtval .gt. 0.99d0) dtval = 0.99d0
      return
      end

      subroutine cdp_dmgc(kappac,dcval)
      real*8 kappac,dcval

      real*8 kttab(60),sbttab(60),ecktab(60),
     1     kctab(60),sbctab(60),ecktabc(60),
     2     dtx(60),dty(60),dcx(60),dcy(60)
      common/cdptab/kttab,sbttab,ecktab,kctab,sbctab,ecktabc,
     1     dtx,dty,dcx,dcy
      integer ntt,ntc,ntdt,ntdc
      common/cdpn/ntt,ntc,ntdt,ntdc

      real*8 eck
      call cdp_interp(kappac,kctab,ecktabc,ntc,eck)
      call cdp_interp(eck,dcx,dcy,ntdc,dcval)
      if (dcval .lt. 0.d0) dcval = 0.d0
      if (dcval .gt. 0.99d0) dcval = 0.99d0
      return
      end


c ======================================================================
c Isotropic elastic stress<->strain relations and stiffness, internal
c ordering (11,22,33,12,23,13), engineering shear.
c ======================================================================
      subroutine cdp_stressfromstrain(strain,stress)
      real*8 strain(6),stress(6),factor

      real*8 ym,pr,psi,tanpsi,ecc,alpha,gammap,fb0fc0,rkc,visc,wt,wc,
     1     sigt0,bulkK,shearG
      common/cdpc/ym,pr,psi,tanpsi,ecc,alpha,gammap,fb0fc0,rkc,visc,
     1     wt,wc,sigt0,bulkK,shearG

      factor = ym/((1.d0+pr)*(1.d0-2.d0*pr))
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

      subroutine cdp_strainfromstress(stress,strain)
      real*8 strain(6),stress(6)

      real*8 ym,pr,psi,tanpsi,ecc,alpha,gammap,fb0fc0,rkc,visc,wt,wc,
     1     sigt0,bulkK,shearG
      common/cdpc/ym,pr,psi,tanpsi,ecc,alpha,gammap,fb0fc0,rkc,visc,
     1     wt,wc,sigt0,bulkK,shearG

      strain(1) = (stress(1)-pr*stress(2)-pr*stress(3))/ym
      strain(2) = (-pr*stress(1)+stress(2)-pr*stress(3))/ym
      strain(3) = (-pr*stress(1)-pr*stress(2)+stress(3))/ym
      strain(4) = 2.d0*(1.d0+pr)*stress(4)/ym
      strain(5) = 2.d0*(1.d0+pr)*stress(5)/ym
      strain(6) = 2.d0*(1.d0+pr)*stress(6)/ym
      return
      end

      subroutine cdp_elasticmatrix(elmat)
      real*8 elmat(6,6),factor

      real*8 ym,pr,psi,tanpsi,ecc,alpha,gammap,fb0fc0,rkc,visc,wt,wc,
     1     sigt0,bulkK,shearG
      common/cdpc/ym,pr,psi,tanpsi,ecc,alpha,gammap,fb0fc0,rkc,visc,
     1     wt,wc,sigt0,bulkK,shearG

      integer i,j
      do i = 1,6
         do j = 1,6
            elmat(i,j) = 0.d0
         end do
      end do
      factor = ym/((1.d0+pr)*(1.d0-2.d0*pr))
      elmat(1,1) = factor*(1.d0-pr)
      elmat(2,2) = factor*(1.d0-pr)
      elmat(3,3) = factor*(1.d0-pr)
      elmat(1,2) = factor*pr
      elmat(1,3) = factor*pr
      elmat(2,3) = factor*pr
      elmat(2,1) = factor*pr
      elmat(3,1) = factor*pr
      elmat(3,2) = factor*pr
      elmat(4,4) = factor*(1.d0-2.d0*pr)/2.d0
      elmat(5,5) = factor*(1.d0-2.d0*pr)/2.d0
      elmat(6,6) = factor*(1.d0-2.d0*pr)/2.d0
      return
      end


c ======================================================================
c Map internal 6-vector (11,22,33,12,23,13) to UMAT ordering
c (11,22,33,12,13,23).
c ======================================================================
      subroutine cdp_mapsig(sig,sigum)
      real*8 sig(6),sigum(6)
      sigum(1) = sig(1)
      sigum(2) = sig(2)
      sigum(3) = sig(3)
      sigum(4) = sig(4)
      sigum(5) = sig(6)
      sigum(6) = sig(5)
      return
      end


c ======================================================================
c Viscous (Duvaut-Lions type) regularization of the damage variables
c and combination into the nominal (damaged) stress:
c   (1-d) = (1 - st*dc_v)(1 - sc*dt_v),  st = 1-wt*r*, sc = 1-wc*(1-r*)
c With visc (=mu, from common cdpc) = 0 this reduces exactly to the
c rate-independent combination using the inviscid dt,dc.
c ======================================================================
      subroutine cdp_combine(sigeff,dtinv,dcinv,rstar,dtvold,dcvold,
     $     dtime,signom,dtv,dcv)

      real*8 sigeff(6),dtinv,dcinv,rstar,dtvold,dcvold,dtime
      real*8 signom(6),dtv,dcv

      real*8 ym,pr,psi,tanpsi,ecc,alpha,gammap,fb0fc0,rkc,visc,wt,wc,
     1     sigt0,bulkK,shearG
      common/cdpc/ym,pr,psi,tanpsi,ecc,alpha,gammap,fb0fc0,rkc,visc,
     1     wt,wc,sigt0,bulkK,shearG

      real*8 vfac,st,sc,dcomb
      integer k

      if (visc .gt. 0.d0) then
         vfac = dtime/(visc+dtime)
         dtv = (1.d0-vfac)*dtvold + vfac*dtinv
         dcv = (1.d0-vfac)*dcvold + vfac*dcinv
         if (dtv .lt. dtvold) dtv = dtvold
         if (dcv .lt. dcvold) dcv = dcvold
      else
         dtv = dtinv
         dcv = dcinv
      end if

      st = 1.d0 - wt*rstar
      sc = 1.d0 - wc*(1.d0-rstar)
      dcomb = 1.d0 - (1.d0-st*dcv)*(1.d0-sc*dtv)

      do k = 1,6
         signom(k) = (1.d0-dcomb)*sigeff(k)
      end do

      return
      end


c ======================================================================
c One-time diagnostic printout of the parsed material constants
c ======================================================================
      subroutine cdp_printinput()

      real*8 ym,pr,psi,tanpsi,ecc,alpha,gammap,fb0fc0,rkc,visc,wt,wc,
     1     sigt0,bulkK,shearG
      common/cdpc/ym,pr,psi,tanpsi,ecc,alpha,gammap,fb0fc0,rkc,visc,
     1     wt,wc,sigt0,bulkK,shearG
      integer ntt,ntc,ntdt,ntdc
      common/cdpn/ntt,ntc,ntdt,ntdc

      write(*,*) '****************************************'
      write(*,*) ' ABAQUS is using the user material CDP   '
      write(*,*) ' (Concrete Damaged Plasticity) UMAT       '
      write(*,*) '----------------------------------------  '
      write(*,2) '  E ..................... = ',ym
      write(*,2) '  nu .................... = ',pr
      write(*,2) '  dilation angle (rad) ... = ',psi
      write(*,2) '  eccentricity ecc ....... = ',ecc
      write(*,2) '  alpha ................ = ',alpha
      write(*,2) '  gamma ................. = ',gammap
      write(*,2) '  viscosity mu .......... = ',visc
      write(*,2) '  wt (tension recovery) . = ',wt
      write(*,2) '  wc (compr. recovery) .. = ',wc
      write(*,2) '  sigma_t0 (from table) . = ',sigt0
      write(*,2) '  bulk modulus K ........ = ',bulkK
      write(*,2) '  shear modulus G ....... = ',shearG
      write(*,*) '  tension table rows .... = ',ntt
      write(*,*) '  compression table rows  = ',ntc
      write(*,*) '  tension damage rows .... = ',ntdt
      write(*,*) '  compression damage rows  = ',ntdc
      write(*,*) '****************************************'
 2    format(1x,A,1pE12.5)

      return
      end
