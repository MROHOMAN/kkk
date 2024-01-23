c
c=====================================================================72
c
c  UMAT subroutine for crystal plasticity - SNL & MSU
c
c=====================================================================72
c
      subroutine umat ( stress,  statev,  ddsdde,  sse,     spd,
     &                  scd,     rpl,     ddsddt,  drplde,  drpldt,
     &                  strain,  dstrain, time,    dtime,   temp,
     &                  dtemp,   predef,  dpred,   cmname,  ndi,
     &                  nshr,    ntens,   nstatv,  props,   nprops,
     &                  coords,  drot,    pnewdt,  celent,  dfgrd0,
     &                  dfgrd1,  noel,    npt,     layer,   kspt,
     &                  kstep,   kinc )

      include 'ABA_PARAM.INC'

      character*8 cmname
c
c---------------------------------------------------------------------72
c---- Dimension arrays passed into the UMAT sub
c---------------------------------------------------------------------72
c
      dimension
     &  coords(3),      ! Coordinates of Gauss pt. being evaluated
     &  ddsdde(ntens,ntens), ! Tangent Stiffness Matrix
     &  ddsddt(ntens),  ! Change in stress per change in temperature
     &  dfgrd0(3,3),    ! Deformation gradient at t_n
     &  dfgrd1(3,3),    ! Deformation gradient at t_(n+1)
     &  dpred(1),       ! Change in predefined state variables
     &  drplde(ntens),  ! Change in heat generation per change in strain
     &  drot(3,3),      ! Incremental rotation matrix
     &  dstrain(ntens), ! Strain increment tensor (vector form)
     &  predef(1),      ! Predefined state vars dependent on field variables
     &  props(nprops),  ! Material properties
     &  statev(nstatv), ! State variables
     &  strain(ntens),  ! Strain tensor (vector form)
     &  stress(ntens),  ! Cauchy stress (vector form)
     &  time(2)         ! Time Step and Total Time
c
c---------------------------------------------------------------------72
c---- Contents of the props array (nprops = 21)
c---------------------------------------------------------------------72
c
c------- crystal type and # of grain per ip
c          props(1)  crystalID : 1-FCC, 2-BCC, 3HCP
c          props(2)  numgrn
c
c------- crystal elasticity
c          props(3)  elastID : 1-ISO, 2-ANISO
c                    ISO      ANISO
c                          CUBIC  HCP
c          props(4)  eMod   C11   C11
c          props(5)  eMu    C12   C12
c          props(6)         C44   C13
c          props(7)               C33
c          props(8)               C44
c
c------- slip system kinetic equation: power law
c          props( 9)  xm     : strain-rate sensitivity (m)
c          props(10)  gam0   : reference shear rate (d_0)
c
c------- slip system hardening law : Voce's type
c          props(11)  h0     : initial work hardening rate (Theta_0)
c          props(12)  tausi  : slip system strength (tau_0)
c          props(13)  taus0  : saturation threshold strength (tau_s0)
c          props(14)  xms    : strain-rate sensitivity - (m')
c          props(15)  gamss0 : ref deformation rate - g_s^star
c          props(16)  kappa0 : initial slip system strength / slip
c                                  taus0 .ge. kappa0
c------- Lattice orientation codes
c          props(17)  kODF
c          props(18)  kODFOut
c
c------- iterations/tolerance data for state and Newton's method
c          props(19)  maxIterState
c          props(20)  tolerState
c          props(21)  maxIterNewt
c          props(22)  tolerNewt
c
c---------------------------------------------------------------------72
c---- Some local variables
c---------------------------------------------------------------------72
c
      real*8 InnerProductVec

      character*160 filePath
      common /XtalWorkDir/ filePath

      character*160 rootName
      common /XtalRootFile/ rootName

      data  filePath
     & /'C:/PBC_6x6x6/'/
      data  rootName / 'test' /

      data  numel  /1/, 
     &      numqpt /1/, 
     &      iprint /1/            ! dummy variables

      integer numel_aba, numqpt_aba, kelem
      common /data_aba/ numel_aba, numqpt_aba

      data kelem /0/
      save kelem
	  
	  integer   GIndex 
      common  /newAddGIndex/ GIndex
      save /newAddGIndex/
	    
      real*8 time_end
      common /timeEnd/ time_end
      data   time_end  /1.5d0/
      
      character*160 filename
      integer length1, I, flag
	  
	  integer NincrF, NincrG
      save    NincrF, NincrG
	  
      
      numel_aba  = nint (props(1))
      numqpt_aba = nint (props(2))
c---------------------------------------------------------------------72
      call UpdateCoordinates(coords, npt, noel) 
c---------------------------------------------------------------------72
c---------------------------------------------------------------------72  
      GIndex = noel
      if (npt .eq. 1) kelem = kelem + 1

      if ( (time(2) .eq. 0.0) ) then
c      if ( (time(2) .eq. 0.0) .and. (kinc .eq. 0) ) then
c         if ( (noel .eq. 1) .and. (npt .eq. 1) ) then
c         if ( (kelem .eq. 1) .and. (npt .eq. 1) ) then
            NincrF = 0
            NincrG = 0
            call Initializegrain(noel)
            call SetUpCrystalProps2(props, nprops, numqpt, numel)
c         endif
      endif

      if ( (time(2) .eq. 0.0) .and. (kinc .gt. 0) ) then
         call SetUpInitialState2(numqpt, numel, kelem, npt,
     &                          nstatv, statev)
      endif
c
c---------------------------------------------------------------------72
c---- Evolve material state - compute stresses and algorithmic moduli
c---------------------------------------------------------------------72
c
      if ( (time(2) .ge. 0.0) .and. (kinc .ge. 1) ) then
c         if (time(2) .eq. 0.0 .and.
c     &      InnerProductVec(dstrain,dstrain,ntens) .eq. 0.d0 ) then
c            call EffectiveElastStiffness(ddsdde, ntens, ndi)
c         else
            call EvolveMatlState2(stress, statev, ddsdde, strain,
     &         dstrain, time, dtime, temp, dtemp, ! ndi, nshr, ntens,
     &         nstatv, drot, dfgrd0, dfgrd1, kelem, npt,
     &         kinc, numel, numqpt, iprint, pnewdt, coords)
c         endif
      endif

      if ( (kelem.eq.numel_aba) .and. (npt.eq.numqpt_aba) ) kelem = 0
      
      
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetUpCrystalProps(
     &   props, mprops, numqpt, numel
     &   )
     
      implicit none
      include 'params_xtal.inc'

      integer numel, numqpt, mprops
      real*8  props(mprops)

      character*160 filePath, fileRoot
      common /XtalWorkDir/ filePath  /XtalRootFile/ fileRoot
c
      integer numgrn, numslip, numvtx, kODFout
      real*8  kappa0(NKAPP), fCeDevVol(NVECS), fCeVol
      real*8  matProp(NPROPS, MAX_SLIP), tauSlip(MAX_SLIP)
      real*8  fCeDev(NVECS, NVECS), fCeiDev(NVECS, NVECS)
      real*8  zBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0(DIMS, DIMS, MAX_SLIP), qBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0Vec(NVECS, MAX_SLIP), qBar0Vec(DIMS, MAX_SLIP)
      real*8  ppTBar0(NVECS, NVECS, MAX_SLIP)
      real*8  gcrot0(DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  sigfs(5, MAX_VTX)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)

      common /XtalNumGrn/ numgrn, numslip
      common /XtalNumVtx/ numvtx
      common /XtalODFOut/ kODFout
      common /XtalState0/ kappa0
      common /XtalProps / fCeDev, fCeiDev, matProp,
     &                    fCeDevVol, fCeVol, hardmtx
      common /XtalSlipG / tauSlip, zBar0, pBar0, qBar0,
     &                    pBar0Vec, qBar0Vec, ppTBar0
      common /XtalCrot0 / gcrot0
      common /XtalStrVtx/ sigfs

      real*8  overstress(MAX_SLIP)
      common /XtalOverS/ overstress

      integer maxIterState, maxIterNewt
      real*8  tolerState, tolerNewt

      common /nlcsolve1/ maxIterState, maxIterNewt
      common /nlcsolve2/ tolerState, tolerNewt
c
      real*8  gstress    (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstress_n  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran    (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran_n  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa     (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa_n   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev    (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev_n  (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geqvalues  (NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot    (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot      (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot_n    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot      (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot_n    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps       (9, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps_n     (9, MAX_GRN, NUMQPT_T, NUMEL_T)

      common  /XtalVars  / gstress, gestran, gkappa, gstatev, geqvalues,
     &                     gcrot, grrot, ggamdot, geps
      common  /XtalVars_n/ gstress_n, gestran_n, gkappa_n, gstatev_n,
     &                     gcrot_n, grrot_n, geps_n
c
      integer kODF
      common /XtalODFCode/ kODF

      real*8  angles(DIMS, MAX_ORIEN)
      real*8  euler(DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      common /EulerAngles/ angles, euler

      integer numor, seed
      common  /XtalOrients/ numor, seed
      
      real*8  gshear(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamtar(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gswitch(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      common  /XtalStep/ gshear, ggamtar, gswitch
      
      real*8  tauc(NKAPP, MAX_GRN), gammat(MAX_SLIP, MAX_GRN)
      common  /TauCrit/ tauc, gammat
      
      real*8  grainsize
      common  /XtalGrainSize/ grainsize  
      
      real*8  ggndden    (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden_n  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      common  /Dislocation/  ggndden, ggndden_n
c
c---------------------------------------------------------------------72
c
      call CrystalInitialize(props, mprops, numel, numqpt,
     &                       filePath, fileRoot,
     &                       XTAL_I, XTAL_E, XTAL_O, XTAL_S, XTAL_A1,
     &                       XTAL_STRESS_O, XTAL_STRAIN_O,
     &                       XTAL_EFFSS_O, XTAL_TRUESS_O,
     &                       XTAL_ITER_O, AGG_EFFSS_O,
     &                       numgrn, numslip, numvtx, kODFout,
     &                       kappa0, fCeDevVol, fCeVol,
     &                       matProp, tauSlip,
     &                       fCeDev, fCeiDev,
     &                       zBar0,
     &                       pBar0, qBar0,
     &                       pBar0Vec, qBar0Vec,
     &                       ppTBar0,
     &                       gcrot0,
     &                       sigfs,
     &                       overstress,
     &                       hardmtx,
     &                       maxIterState, maxIterNewt,
     &                       tolerState, tolerNewt,
     &                       kODF, numor, seed, angles, euler,
     &                       XTAL_TXT_OUT,
     &                       gstress,
     &                       gstress_n,
     &                       gestran,
     &                       gestran_n,
     &                       gkappa,
     &                       gkappa_n,
     &                       gstatev,
     &                       gstatev_n,
     &                       geqvalues,
     &                       ggamdot,
     &                       gcrot,
     &                       gcrot_n,
     &                       grrot,
     &                       grrot_n, gshear, ggamtar, gswitch,
     &                       tauc, gammat, grainsize, geps, geps_n, 
     &                       ggndden, ggndden_n)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetUpInitialState(
     &   numqpt, numel, noel, npt, nstatv, statev
     &   )

      implicit none
      include  'params_xtal.inc'

      integer numel, numqpt, noel, npt
      integer nstatv
      real*8  statev(nstatv)

      integer numgrn, numslip, kODF
      integer numor, seed
      real*8  gcrot0(DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  angles(DIMS, MAX_ORIEN)
      real*8  euler(DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)

      common /XtalNumGrn/  numgrn, numslip
      common /XtalODFCode/ kODF
      common /XtalOrients/ numor, seed
      common /XtalCrot0 /  gcrot0
      common /EulerAngles/ angles, euler
c
      real*8  kappa0(NKAPP)
      common  /XtalState0/ kappa0
      
      real*8  tauc(NKAPP, MAX_GRN), gammat(MAX_SLIP, MAX_GRN)
      common  /TauCrit/ tauc, gammat
c
c---------------------------------------------------------------------72
c
      call AssignCrystalODF(numel, numqpt, noel, npt,
     &                      numgrn, numslip, kODF,
     &                      numor, seed,
     &                      gcrot0,
     &                      angles,
     &                      euler, 
     &                      XTAL_O, XTAL_TXT_OUT)

      call SetUpStateVars(nstatv, statev,
     &                    numgrn, numslip,
     &                    kappa0,
     &                    gcrot0,
     &                    XTAL_O, tauc, gammat, noel)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetUpStateVars(
     &   nstatv, statev,
     &   numgrn, numslip,
     &   kappa0,
     &   gcrot0,
     &   FILE_O,
     &   tauc, gammat, noel) 

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer nstatv
      real*8  statev(nstatv)

      integer numgrn, numslip, FILE_O, noel
      real*8  kappa0(NKAPP)
      real*8  gcrot0 (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  tauc(NKAPP, MAX_GRN), gammat(MAX_SLIP, MAX_GRN)

      integer dex, ig, varsPerGrn
c
c---------------------------------------------------------------------72
c
c---- initialize state variable vector per ip
c
      varsPerGrn = 2*NVECS          ! stress, estran
c     &             + NKAPP          ! kappa
     &             + numslip        ! kappa
     &             + NSTAV          ! ee_v, p, detVe, WpNorm, ssact
     &             + NEQVA/2        ! eqps, eqstr, gam_star, gamtot 
     &             + 3*DIMS*DIMS    ! C_0, C, R
     &             + numslip        ! gamdot
     &             + 4*numslip      
     &             + 9              ! eps
     &             + 2      

      if (nstatv .ne. numgrn*varsPerGrn)
     &   call RunTimeError(FILE_O, 'nstatv .ne. numgrn*varsPerGrn')

      do ig = 1, numgrn

         dex = (ig-1)*varsPerGrn + 1           ! stress (xtal)
         call SetTensor(statev(dex), pzero, NVECS)

         dex = dex + NVECS                      ! estrain
         call SetTensor(statev(dex), pzero, NVECS)

         dex = dex + NVECS                      ! kappa
         call EqualTensors(tauc(1,noel), statev(dex), numslip)

c         dex = dex + NKAPP              ! ee_v, p, detVe, WpNorm, ssact
         dex = dex + numslip              ! ee_v, p, detVe, WpNorm, ssact
         call SetTensor(statev(dex), pzero, NSTAV)

         dex = dex + NSTAV               ! eqps, eqstr, gam_star, gamtot
         call SetTensor(statev(dex), pzero, NEQVA/2)

         dex = dex + NEQVA/2                    ! C_0
         call EqualTensors(gcrot0(1,1,ig,1,1), statev(dex), DIMS*DIMS)

         dex = dex + DIMS*DIMS                  ! C
         call EqualTensors(gcrot0(1,1,ig,1,1), statev(dex), DIMS*DIMS)

         dex = dex + DIMS*DIMS                  ! R
         call EqualTensors(Ident2nd, statev(dex), DIMS*DIMS)

         dex = dex + DIMS*DIMS                  ! gamdot
         call SetTensor(statev(dex), pzero, numslip)
         
         dex = dex + numslip                    ! shear
         call SetTensor(statev(dex), pzero, numslip)

         dex = dex + numslip                    ! gndden
         call SetTensor(statev(dex), pzero, numslip)
         
         dex = dex + numslip                    ! gamtar
         call EqualTensors(gammat(1,noel), statev(dex), numslip)
         
         dex = dex + numslip                    ! switch
         call SetTensor(statev(dex), pzero, numslip)
         
         dex = dex + numslip                    ! eps
         call SetTensor(statev(dex), pzero, 9)
         
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE EvolveMatlState(
     &   stress, statev, ddsdde, strain, dstrain, time, 
     &   dtime, temp, dtemp, ! ndi, nshr, ntens, 
     &   nstatv, drot, dfgrd0, dfgrd1, noel, npt,
     &   incr, numel, numqpt, iprint, pnewdt, coords
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer nstatv  !, ndi, nshr, ntens
      integer noel, npt, incr, numel, numqpt, iprint
      real*8  dtime, temp, dtemp, pnewdt
      real*8  stress(ntens), statev(nstatv), ddsdde(ntens, ntens)
      real*8  strain(ntens), dstrain(ntens), time(2)
      real*8  dfgrd0(3,3), dfgrd1(3,3), drot(3,3)
      real*8  coords(3)

      integer statusFlag, numIncrs
      integer numgrn, numslip
      real*8  gstress   (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran   (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa    (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev   (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geqvalues (NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot     (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot     (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot   (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps      (9, MAX_GRN, NUMQPT_T, NUMEL_T)

      real*8  gstress_n (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran_n (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa_n  (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev_n (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot_n   (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot_n   (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot0    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps_n    (9, MAX_GRN, NUMQPT_T, NUMEL_T)

      real*8  vgrad(DIMS, DIMS)
      real*8  savg_ij(DIMS, DIMS, NUMQPT_T) 
      real*8  cavg_ijkl(DIMS2, DIMS2, NUMQPT_T)
      real*8  epsxx(NUMQPT_T), epsyy(NUMQPT_T), epsth(NUMQPT_T)
      real*8  epsxy(NUMQPT_T), epsum(NUMQPT_T), spinn(NUMQPT_T)
      real*8  qptpo(NUMQPT_T), qptp(NUMQPT_T)
      real*8  sigxx(NUMQPT_T), sigyy(NUMQPT_T)
      real*8  sigth(NUMQPT_T), sigxy(NUMQPT_T)

      real*8  epsxz(NUMQPT_T),  epsyz(NUMQPT_T)
      real*8  spinxz(NUMQPT_T), spinyz(NUMQPT_T)
      real*8  sigxz(NUMQPT_T),  sigyz(NUMQPT_T)
      
      integer ielem

      common  /XtalNumGrn/ numgrn, numslip
      common  /XtalCrot0 / gcrot0
      common  /XtalVars  / gstress, gestran, gkappa, gstatev, geqvalues,
     &                     gcrot, grrot, ggamdot, geps
      common  /XtalVars_n/ gstress_n, gestran_n, gkappa_n, gstatev_n,
     &                     gcrot_n, grrot_n, geps_n
      common  /for3DProbs/ epsxz, epsyz, spinxz, spinyz, sigxz, sigyz
      
      real*8  gshear(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamtar(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gswitch(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      common  /XtalStep/ gshear, ggamtar, gswitch
      
      real*8  ggndden   (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden_n (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      common  /Dislocation/  ggndden, ggndden_n
c
c---------------------------------------------------------------------72
c
c-------- kinematic quantities - deformation path
c
      call VelocityGrad(dstrain, drot, dfgrd0, dfgrd1, vgrad, dtime, 
     &                  noel, npt, incr)
      
      call DeformPath(epsxx, epsyy, epsth, epsxy, epsum, spinn,
     &      qptpo, qptp, vgrad, temp, (temp+dtemp), numel, numqpt)
c
c-------- fetch state variables from abaqus vector array
c
      call RecoverStateVars(statev, nstatv, gstress_n, gestran_n,
     &      gkappa_n, gstatev_n, gcrot_n, grrot_n, gcrot0, ggamdot, 
     &      geqvalues, numgrn, numslip, gshear, ggamtar, gswitch,
     &      geps_n, ggndden_n)
c
c-------- compute state
c
      numIncrs = 3999
      ielem    = 1
      statusFlag = XTAL_CONVERGED
      call StressCrystal(sigxx, sigyy, sigth, sigxy, epsxx, epsyy,
     &      epsth, epsxy, epsum, spinn, qptpo, qptp, dtime, time(2),
     &      ielem, incr, numel, numqpt, iprint, savg_ij, cavg_ijkl, 
     &      statusFlag, numIncrs, noel, npt, dfgrd1, coords)
      if (statusFlag .ne. XTAL_CONVERGED) then
         write(XTAL_O, *) ' ** Umat did not converged       **'
         write(XTAL_O, *) ' ** re-scaling time step by 0.75 **'
         pnewdt = 0.75
         return
      endif
c
c-------- save state variables in abaqus vector array
c
      call SaveStateVars(statev, nstatv, gstress, gestran, gkappa, 
     &      gstatev, gcrot, grrot, gcrot0, ggamdot, geqvalues, numgrn,
     &      numslip, gshear, ggamtar, gswitch, geps, ggndden)
c
c-------- stresses and algorithmic moduli in abaqus format
c
      call SaveStressModuli(stress, ddsdde, savg_ij(1, 1, 1), 
     &                      cavg_ijkl(1, 1, 1))

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE VelocityGrad(
     &   dstran, drot_aba, dfgrd0, dfgrd1, vgrad, dtime, ielem, iqpt,
     &   incr
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

c      integer ntens, nshr
      integer ielem, iqpt, incr
      real*8  dtime
      real*8  dstran(NTENS)
      real*8  drot_aba(3,3), dfgrd0(3,3), dfgrd1(3,3), vgrad(3,3)

      integer i, j
      real*8  det
      real*8  d_aba(3,3)
      real*8  d_umat(3,3), w_umat(3,3), drot_umat(3,3), spin(3)
      real*8  matx_1(3,3), matx_2(3,3), matx_3(3,3), rdfgrd(3,3)
c
c---------------------------------------------------------------------72
c
c-------- relative deformation gradient: f=F_(n+1)*F_n^{-1}
c
      call InverseOfMat3x3(dfgrd0, matx_1, det)
      call MultAxB(dfgrd1, matx_1, rdfgrd, 3)
c
c-------- approximation to velocity gradient: L=2/dt*(f-I)*(f+I)^{-1}
c
      call AddTensors(pone, rdfgrd, -pone, Ident2nd, matx_2, 3*3)
      call AddTensors(pone, rdfgrd, pone, Ident2nd, matx_1, 3*3)
      call InverseOfMat3x3(matx_1, matx_3, det)

      call MultAxB(matx_2, matx_3, vgrad, 3)
      call SetToScaledtensor(2.0/dtime, vgrad, vgrad, 3*3)
c
c-------- rate of deformation and spin: d_umat=sym(L), w_umat=skw(L)
c
      call SetTensor(d_umat, pzero, DIMS*DIMS)
      call SetTensor(w_umat, pzero, DIMS*DIMS)

      call SymmetrizeTensor(vgrad, d_umat, 3)
      call SkewSymmetrizeTensor(vgrad, w_umat, 3)
c
c-------- incremental rotation tensor: drot_umat = exp(dt*w_umat)
c
      call Mat3x3ToVec3x1Skew(w_umat, spin, 3)
      call IncrementalRotTensor(drot_umat, spin, dtime)
c
c-------- rate of deformation from ABAQUS
c
      call SetTensor(d_aba, pzero, DIMS*DIMS)
      d_aba(1,1) = dstran(1)/dtime
      d_aba(2,2) = dstran(2)/dtime
      d_aba(3,3) = dstran(3)/dtime
      d_aba(1,2) = dstran(4)/ptwo/dtime
      d_aba(2,1) = d_aba(1,2)
      
      if (NSHR. gt. 1) then
         d_aba(1,3) = dstran(5)/ptwo/dtime
         d_aba(2,3) = dstran(6)/ptwo/dtime
         d_aba(3,1) = d_aba(1,3)
         d_aba(3,2) = d_aba(2,3)
      endif
c
c-------- print values for comparison
c
      if (ielem .eq. kPRINT_ELEM .and. iqpt .eq. kPRINT_QPT) then
         if (incr .eq. 1) write(XTAL_E, 1000) ielem, iqpt
        
         do i = 1, 3
            write(XTAL_E, 5000) (vgrad(i,j), j=1,3)
         enddo
         write(XTAL_E, *)
         do i = 1, 3
            write(XTAL_E, 5000) (d_aba(i,j), j=1,3),
     &                          (d_umat(i,j),j=1,3), 
     &                          (w_umat(i,j),j=1,3)
         enddo
         write(XTAL_E, *)
         do i = 1, 3
            write(XTAL_E, 5000) (drot_aba(i,j), j=1,3),
     &                          (drot_umat(i,j),j=1,3)
         enddo
      endif
      
1000  format(' INCR',/,18x,' L_ij',/,
     &       18x,' D_aba',32x,'D_umat',30x,'W_umat',/,
     &       18x,' Drot_aba',32x,'Drot_umat',/,
     &       38x,' (elem # ', i5, ',  qpt # ', i2, ')')
5000  format((3(2x,3(1x,e11.4))))
c
c-------- recompute velocity gradient: L = d_aba + w_umat
c
      call AddTensors(pone, d_aba, pone, w_umat, vgrad, DIMS*DIMS)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE DeformPath(
     &   epsxx, epsyy, epsth, epsxy, epsum, spinn, qptpo, qptp, vgrad, 
     &   theta_n, theta, numel, numqpt
     &   )

      implicit none
      include 'params_xtal.inc'

      integer numel, numqpt
      real*8  theta_n, theta
      real*8  vgrad(DIMS, DIMS)
      real*8  epsxx(numqpt), epsyy(numqpt), epsth(numqpt)
      real*8  epsxy(numqpt), epsum(numqpt), spinn(numqpt)
      real*8  qptpo(numqpt), qptp(numqpt)

      real*8  epsxz(NUMQPT_T),  epsyz(NUMQPT_T)
      real*8  spinxz(NUMQPT_T), spinyz(NUMQPT_T)
      real*8  sigxz(NUMQPT_T),  sigyz(NUMQPT_T)
      common /for3DProbs/ epsxz, epsyz, spinxz, spinyz, sigxz, sigyz
c
c---------------------------------------------------------------------72
c
c------- deformation rate: deviatoric and volumetric 
c
      epsxx(1) = vgrad(1, 1)
      epsyy(1) = vgrad(2, 2)
      epsth(1) = vgrad(3, 3)

      epsum(1) = epsxx(1) + epsyy(1) + epsth(1)

      epsxx(1) = epsxx(1) - epsum(1) / 3.0d0
      epsyy(1) = epsyy(1) - epsum(1) / 3.0d0
      epsth(1) = epsth(1) - epsum(1) / 3.0d0

      epsxy(1) = (vgrad(2, 1) + vgrad(1, 2))/2.d0
      epsxz(1) = (vgrad(3, 1) + vgrad(1, 3))/2.d0
      epsyz(1) = (vgrad(3, 2) + vgrad(2, 3))/2.d0
c
c------- axial vector of spin
c
      spinn(1)  = (vgrad(1, 2) - vgrad(2, 1))/2.d0
      spinxz(1) = (vgrad(1, 3) - vgrad(3, 1))/2.d0
      spinyz(1) = (vgrad(2, 3) - vgrad(3, 2))/2.d0
c
c------- temperature (K)
c
      qptpo(1) = theta_n
      qptp(1)  = theta

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE RecoverStateVars(
     &   statev, nstatv, gstress_n, gestran_n, gkappa_n, gstatev_n, 
     &   gcrot_n, grrot_n, gcrot0, ggamdot, geqvalues, numgrn, numslip,
     &   gshear, ggamtar, gswitch, geps_n, ggndden_n)

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'
      
      integer nstatv, numgrn, numslip
      real*8  statev(nstatv)

      real*8  gstress_n (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran_n (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa_n  (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev_n (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot_n   (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot_n   (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot0    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps_n    (9, MAX_GRN, NUMQPT_T, NUMEL_T)

      real*8  geqvalues (NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot   (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)

      integer varsPerGrn, ig, dex
      
      real*8  gshear(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamtar(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gswitch(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden_n (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
c
c---------------------------------------------------------------------72
c
c---- number of state variables per ip
c
      varsPerGrn = 2*NVECS          ! stress, estran
c     &             + NKAPP          ! kappa
     &             + numslip        ! kappa
     &             + NSTAV          ! ee_v, p, detVe, WpNorm, ssact
     &             + NEQVA/2        ! eqps, eqstr, gam_star, gamtot 
     &             + 3*DIMS*DIMS    ! C_0, C, R
     &             + numslip        ! gam_dot
     &             + 4*numslip      
     &             + 9              ! eps
     &             + 2      
c
c---- recover state variables from abaqus vector array
c
      do ig = 1, numgrn

         dex = (ig-1)*varsPerGrn + 1           ! stress
         call EqualTensors(statev(dex), gstress_n(1,ig,1,1), NVECS)

         dex = dex + NVECS                      ! estrain
         call EqualTensors(statev(dex), gestran_n(1,ig,1,1), NVECS)

         dex = dex + NVECS                      ! kappa
c         call EqualTensors(statev(dex), gkappa_n(1,ig,1,1), NKAPP)
         call EqualTensors(statev(dex), gkappa_n(1,ig,1,1), numslip)

c         dex = dex + NKAPP              ! ee_v, p, detVe, WpNorm, ssact
         dex = dex + numslip             ! ee_v, p, detVe, WpNorm, ssact
         call EqualTensors(statev(dex), gstatev_n(1,ig,1,1), NSTAV)

         dex = dex + NSTAV               ! eqps, eqstr, gam_star, gamtot
         call EqualTensors(statev(dex), geqvalues(1,ig,1,1), NEQVA/2)

         dex = dex + NEQVA/2                    ! C_0
         call EqualTensors(statev(dex), gcrot0(1,1,ig,1,1), DIMS*DIMS)

         dex = dex + DIMS*DIMS                  ! C
         call EqualTensors(statev(dex), gcrot_n(1,1,ig,1,1), DIMS*DIMS)

         dex = dex + DIMS*DIMS                  ! R
         call EqualTensors(statev(dex), grrot_n(1,1,ig,1,1), DIMS*DIMS)

         dex = dex + DIMS*DIMS                  ! gamdot
         call EqualTensors(statev(dex), ggamdot(1,ig,1,1), numslip)
         
         dex = dex + numslip                    ! shear
         call EqualTensors(statev(dex), gshear(1,ig,1,1), numslip)
         
         dex = dex + numslip                    ! gndden
         call EqualTensors(statev(dex), ggndden_n(1,ig,1,1), numslip)
         
         dex = dex + numslip                    ! gamtar
         call EqualTensors(statev(dex), ggamtar(1,ig,1,1), numslip)
         
         dex = dex + numslip                    ! switch
         call EqualTensors(statev(dex), gswitch(1,ig,1,1), numslip)
         
         dex = dex + numslip                    ! eps
         call EqualTensors(statev(dex), geps_n(1,ig,1,1), 9)

      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SaveStateVars(
     &   statev, nstatv, gstress, gestran, gkappa, gstatev, gcrot, 
     &   grrot, gcrot0, ggamdot, geqvalues, numgrn, numslip, gshear,
     &   ggamtar, gswitch, geps, ggndden)

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'
      
      integer nstatv, numgrn, numslip
      real*8  statev(nstatv)

      real*8  gstress   (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran   (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa    (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev   (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot     (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot     (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot0    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps      (9, MAX_GRN, NUMQPT_T, NUMEL_T)

      real*8  geqvalues (NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot   (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)

      integer varsPerGrn, ig, dex
      
      real*8  gshear(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamtar(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gswitch(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden   (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
c
c---------------------------------------------------------------------72
c
c---- number of state variables per ip
c
      varsPerGrn = 2*NVECS          ! stress, estran
c     &             + NKAPP          ! kappa
     &             + numslip        ! kappa
     &             + NSTAV          ! ee_v, p, detVe, WpNorm, ssact
     &             + NEQVA/2        ! eqps, eqstr, gam_star, gamtot 
     &             + 3*DIMS*DIMS    ! C_0, C, R
     &             + numslip        ! gam_dot
     &             + 4*numslip      ! shear, NanoGND, gamtar, switch
     &             + 9              ! geps
     &             + 2
c
c---- save state variables in abaqus vector array
c
      do ig = 1, numgrn

         dex = (ig-1)*varsPerGrn + 1           ! stress
         call EqualTensors(gstress(1,ig,1,1), statev(dex), NVECS)

         dex = dex + NVECS                      ! estrain
         call EqualTensors(gestran(1,ig,1,1), statev(dex), NVECS)

         dex = dex + NVECS                      ! kappa
c         call EqualTensors(gkappa(1,ig,1,1), statev(dex), NKAPP)
         call EqualTensors(gkappa(1,ig,1,1), statev(dex), numslip)

c         dex = dex + NKAPP              ! ee_v, p, detVe, WpNorm, ssact
         dex = dex + numslip             ! ee_v, p, detVe, WpNorm, ssact
         call EqualTensors(gstatev(1,ig,1,1), statev(dex), NSTAV)

         dex = dex + NSTAV               ! eqps, eqstr, gam_star, gamtot
         call EqualTensors(geqvalues((NEQVA/2+1),ig,1,1), statev(dex),
     &                                                        NEQVA/2)

         dex = dex + NEQVA/2                    ! C_0
c         call EqualTensors(gcrot0(1,1,ig,1,1), statev(dex), DIMS*DIMS)

         dex = dex + DIMS*DIMS                  ! C
         call EqualTensors(gcrot(1,1,ig,1,1), statev(dex), DIMS*DIMS)

         dex = dex + DIMS*DIMS                  ! R
         call EqualTensors(grrot(1,1,ig,1,1), statev(dex), DIMS*DIMS)

         dex = dex + DIMS*DIMS                  ! gamdot
         call EqualTensors(ggamdot(1,ig,1,1), statev(dex), numslip)
         
         dex = dex + numslip                    ! shear
         call EqualTensors(gshear(1,ig,1,1), statev(dex), numslip)
         
         dex = dex + numslip                    ! gndden
         call EqualTensors(ggndden(1,ig,1,1), statev(dex), numslip)
         
         dex = dex + numslip                    ! gamtar
         call EqualTensors(ggamtar(1,ig,1,1), statev(dex), numslip)
         
         dex = dex + numslip                    ! switch
         call EqualTensors(gswitch(1,ig,1,1), statev(dex), numslip)
         
         dex = dex + numslip                    ! eps
         call EqualTensors(geps(1,ig,1,1), statev(dex), 9)

      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE DeformationRate_n(
     &   dvec_n, wvec_n, d_kk_n, iqpt
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer iqpt
      real*8  d_kk_n
      real*8  dvec_n(NVECS), wvec_n(3)
c
c---------------------------------------------------------------------72
c
c------ zeros arrays
c
      d_kk_n = 0.0
      call SetTensor(dvec_n, pzero, NVECS)
      call SetTensor(wvec_n, pzero, 3)
c
c------ stop program if called during MPS runs ??
c
c      print *, 'DeformationRate_n called ... Exiting'
c      call RunTimeError(XTAL_O, 'DeformationRate_n called')

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SaveStressModuli(
     &   stress, ddsdde, savg_ij, cavg_ijkl
     &   )

      implicit none
      include 'params_xtal.inc'

      real*8  stress(NTENS), ddsdde(NTENS, NTENS)
      real*8  savg_ij(DIMS, DIMS), cavg_ijkl(DIMS2, DIMS2)

      integer i, j
c
c---------------------------------------------------------------------72
c
c---- stresses
c
      stress(1) = savg_ij(1,1)
      stress(2) = savg_ij(2,2)
      stress(3) = savg_ij(3,3)
      stress(4) = savg_ij(1,2)

      if (NSHR .gt. 1) then
         stress(5) = savg_ij(1,3)
         stress(6) = savg_ij(2,3)
      endif
c
c---- algorithmic moduli
c
      do i = 1, NTENS
         do j = 1, NTENS
            ddsdde(i,j) = cavg_ijkl(i,j)
         enddo
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE EffectiveElastStiffness(
     &   ddsdde, ntens, ndi
     &   )

      implicit none
      include 'numbers.inc'

      integer ntens, ndi
      real*8  ddsdde(ntens, ntens)

      integer k1, k2
      real*8  eg2, ebulk3, elam

      real*8  eg, ebulk
      common /ElastProps/ eg, ebulk
c
c---------------------------------------------------------------------72
c
c----- effective elastic stiffness
c
      eg2    = 2.0*eg
      ebulk3 = 3.0*ebulk
      elam   = (ebulk3-eg2)/pthree

      call SetTensor(ddsdde, pzero, ntens*ntens)

      do k1 = 1, ndi
        do k2 = 1, ndi
           ddsdde(k2,k1) = elam
        enddo
        ddsdde(k1,k1) = eg2+elam
      enddo
      do k1 = ndi+1, ntens
        ddsdde(k1,k1) = eg
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE CrystalInitialize(
     &   props, mprops, numel, numqpt,
     &   filePath, fileRoot,
     &   FILE_I, FILE_E, FILE_O, FILE_S, FILE_A1,
     &   FILE_STRESS_O, FILE_STRAIN_O,
     &   FILE_EFFSS_O, FILE_TRUESS_O,
     &   FILE_ITER_O, FILE_AGG_EFFSS_O,
     &   numgrn, numslip, numvtx, kODFout,
     &   kappa0, fCeDevVol, fCeVol,
     &   matProp, tauSlip,
     &   fCeDev, fCeiDev,
     &   zBar0,
     &   pBar0, qBar0,
     &   pBar0Vec, qBar0Vec,
     &   ppTBar0,
     &   gcrot0,
     &   sigfs,
     &   overstress,
     &   hardmtx,
     &   maxIterState, maxIterNewt,
     &   tolerState, tolerNewt,
     &   kODF, numor, seed, angles, euler,
     &   FILE_TXT_OUT,
     &   gstress,
     &   gstress_n,
     &   gestran,
     &   gestran_n,
     &   gkappa,
     &   gkappa_n,
     &   gstatev,
     &   gstatev_n,
     &   geqvalues,
     &   ggamdot,
     &   gcrot,
     &   gcrot_n,
     &   grrot,
     &   grrot_n, gshear, ggamtar, gswitch,
     &   tauc, gammat, grainsize, geps, geps_n, ggndden, ggndden_n)

      implicit none
      include 'params_xtal.inc'

      integer mprops, numel, numqpt
      real*8  props(mprops)

      character*160 filePath, fileRoot
c
      integer numgrn, numslip, numvtx, kODFout
      real*8  kappa0(NKAPP), fCeDevVol(NVECS), fCeVol
      real*8  matProp(NPROPS, MAX_SLIP), tauSlip(MAX_SLIP)
      real*8  fCeDev(NVECS, NVECS), fCeiDev(NVECS, NVECS)
      real*8  zBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0(DIMS, DIMS, MAX_SLIP), qBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0Vec(NVECS, MAX_SLIP), qBar0Vec(DIMS, MAX_SLIP)
      real*8  ppTBar0(NVECS, NVECS, MAX_SLIP)
      real*8  gcrot0(DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  sigfs(5, MAX_VTX)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)

      real*8  overstress(MAX_SLIP)

      integer maxIterState, maxIterNewt
      real*8  tolerState, tolerNewt
c
      real*8  gstress    (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstress_n  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran    (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran_n  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa     (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa_n   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev    (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev_n  (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geqvalues  (NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot    (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot      (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot_n    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot      (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot_n    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps       (9, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps_n     (9, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden    (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden_n  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
c
      integer kODF

      real*8  angles(DIMS, MAX_ORIEN)
      real*8  euler(DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)

      integer numor, seed

      integer FILE_I, FILE_E, FILE_O, FILE_S, FILE_A1,
     &        FILE_STRESS_O, FILE_STRAIN_O,
     &        FILE_EFFSS_O, FILE_TRUESS_O,
     &        FILE_ITER_O, FILE_AGG_EFFSS_O,
     &        FILE_TXT_OUT
     
      real*8  gshear(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamtar(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gswitch(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      
      real*8  tauc(NKAPP, MAX_GRN), gammat(MAX_SLIP, MAX_GRN)
      real*8  grainsize
c
c---------------------------------------------------------------------72
c
c------- open input/output files
c
      call CrystalOpenIOFiles(filePath, fileRoot, 
     &                        FILE_I, FILE_E, FILE_O, FILE_S, FILE_A1)
      call CrystalOpenIOFiles_2(filePath, fileRoot, 
     &                          FILE_STRESS_O, FILE_STRAIN_O, 
     &                          FILE_EFFSS_O, FILE_TRUESS_O, 
     &                          FILE_ITER_O, FILE_AGG_EFFSS_O)
c
c------- read/echo material data
c
      call CrystalModelData(props, mprops, numqpt, numel,
     &                      numgrn, numslip, numvtx, kODFout,
     &                      kappa0, fCeDevVol, fCeVol,
     &                      matProp, tauSlip,
     &                      fCeDev, fCeiDev,
     &                      zBar0,
     &                      pBar0, qBar0,
     &                      pBar0Vec, qBar0Vec,
     &                      ppTBar0,
     &                      gcrot0,
     &                      sigfs,
     &                      overstress,
     &                      hardmtx,
     &                      maxIterState, maxIterNewt,
     &                      tolerState, tolerNewt,
     &                      FILE_I, FILE_E, FILE_O, filePath,
     &                      kODF, numor, seed, angles, euler,
     &                      FILE_TXT_OUT, tauc, gammat, grainsize)
c
c------- initialize arrays
c
      call CrystalInitializeArrays(numqpt, numel,
     &                             numgrn, numslip,
     &                             kappa0,
     &                             gcrot0,
     &                             gstress,
     &                             gstress_n,
     &                             gestran,
     &                             gestran_n,
     &                             gkappa,
     &                             gkappa_n,
     &                             gstatev,
     &                             gstatev_n,
     &                             geqvalues,
     &                             ggamdot,
     &                             gcrot,
     &                             gcrot_n,
     &                             grrot,
     &                             grrot_n, gshear, ggamtar, gswitch,
     &                             tauc, geps, geps_n,
     &                             ggndden, ggndden_n)
c
c------- initialize useful matrices to compute algorithmic moduli Cep
c
      call CrystalInitializeMatrxCep( )
c
c------- parameters for global iterations
c
      call GlobalIterationParams( )

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE CrystalOpenIOFiles(
     &   filePath, fileRoot, 
     &   FILE_I, FILE_E, FILE_O, FILE_S, FILE_A1
     &   )

      implicit none
    
      character*160 filePath, fileRoot
      integer      FILE_I, FILE_E, FILE_O, FILE_S, FILE_A1

      integer      length1, length2
      character*160 dataFile, filename
c
c---------------------------------------------------------------------72
c
c------- root name of input/output files
c
c      write(*,*) ' Root name of crystal-model input/output files ?'
c      read 1000, dataFile
c
c------- open files
c
      length1 = index(filePath,' ') - 1
      length2 = index(fileRoot,' ') - 1

      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtali'
      open(unit=FILE_I, file=filename, status='unknown',
     &                                           access='sequential')
      rewind(FILE_I)

      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtale'
      open(unit=FILE_E, file=filename, status='unknown')
      rewind(FILE_E)

      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtalo'
      open(unit=FILE_O, file=filename, status='unknown')
      
      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtals'
      open(unit=FILE_S, file=filename, status='unknown')
      rewind(FILE_S)
      
      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtala1'
      open(unit=FILE_A1, file=filename, status='unknown')
c
c------- formats
c
1000  format(a80)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE CrystalOpenIOFiles_2(
     &   filePath, fileRoot, 
     &   FILE_STRESS_O, FILE_STRAIN_O, 
     &   FILE_EFFSS_O, FILE_TRUESS_O, 
     &   FILE_ITER_O, FILE_AGG_EFFSS_O
     &   )

      implicit none

      character*160 filePath, fileRoot
      integer FILE_STRESS_O, FILE_STRAIN_O, 
     &        FILE_EFFSS_O, FILE_TRUESS_O, 
     &        FILE_ITER_O, FILE_AGG_EFFSS_O
    
      integer      length1, length2
      character*160 dataFile, filename
c
c---------------------------------------------------------------------72
c
c------- root name of i/o files
c
c      write(*,*) ' Root name of additional output files ?'
c      read 1000, dataFile
c
c------- open files for output at specified elem/iqpt/grain
c
      length1 = index(filePath,' ') - 1
      length2 = index(fileRoot,' ') - 1

      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtal.strs'
      open(unit=FILE_STRESS_O, file=filename, status='unknown')

      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtal.strn'
      open(unit=FILE_STRAIN_O, file=filename, status='unknown')

      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtal.efss'
      open(unit=FILE_EFFSS_O, file=filename, status='unknown')

      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtal.trss'
      open(unit=FILE_TRUESS_O, file=filename, status='unknown')

      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtal.iter'
      open(unit=FILE_ITER_O, file=filename, status='unknown')
c
c------- open files for aggregate output at specified elem/iqpt
c
      filename = filePath(1:length1)//fileRoot(1:length2)//'.agg.efss'
      open(unit=FILE_AGG_EFFSS_O, file=filename, status='unknown')
c
c------- formats
c
1000  format(a80)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE CrystalModelData(
     &   props, mprops, numqpt, numel,
     &   numgrn, numslip, numvtx, kODFout,
     &   kappa0, fCeDevVol, fCeVol,
     &   matProp, tauSlip,
     &   fCeDev, fCeiDev,
     &   zBar0,
     &   pBar0, qBar0,
     &   pBar0Vec, qBar0Vec,
     &   ppTBar0,
     &   gcrot0,
     &   sigfs,
     &   overstress,
     &   hardmtx,
     &   maxIterState, maxIterNewt,
     &   tolerState, tolerNewt,
     &   FILE_I, FILE_E, FILE_O, filePath,
     &   kODF, numor, seed, angles, euler,
     &   FILE_TXT_OUT, tauc, gammat,
     &   grainsize)

      implicit none
      include 'params_xtal.inc'

      character*160 filePath
      integer mprops, numqpt, numel
      real*8  props(mprops)

      integer numgrn, numslip, numvtx, kODFout
      real*8  kappa0(NKAPP), fCeDevVol(NVECS), fCeVol
      real*8  matProp(NPROPS, MAX_SLIP), tauSlip(MAX_SLIP)
      real*8  fCeDev(NVECS, NVECS), fCeiDev(NVECS, NVECS)
      real*8  zBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0(DIMS, DIMS, MAX_SLIP), qBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0Vec(NVECS, MAX_SLIP), qBar0Vec(DIMS, MAX_SLIP)
      real*8  ppTBar0(NVECS, NVECS, MAX_SLIP)
      real*8  gcrot0(DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  sigfs(5, MAX_VTX)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)

      real*8  overstress(MAX_SLIP)

      integer maxIterState, maxIterNewt
      real*8  tolerState, tolerNewt

      integer FILE_I, FILE_E, FILE_O
      integer FILE_TXT_OUT

      integer kODF, numor, seed
      real*8  angles(DIMS, MAX_ORIEN)
      real*8  euler(DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
          
      real*8  tauc(NKAPP, MAX_GRN), gammat(MAX_SLIP, MAX_GRN)
      real*8  grainsize
c
c---------------------------------------------------------------------72
c
c------- input material data for single crystal
c
      call CrystalMaterialData(props, mprops, numqpt, numel,
     &                      numgrn, numslip, numvtx, kODFout,
     &                      kappa0, fCeDevVol, fCeVol,
     &                      matProp, tauSlip,
     &                      fCeDev, fCeiDev,
     &                      zBar0,
     &                      pBar0, qBar0,
     &                      pBar0Vec, qBar0Vec,
     &                      ppTBar0,
     &                      gcrot0,
     &                      sigfs,
     &                      overstress,
     &                      hardmtx,
     &                      FILE_I, FILE_E, FILE_O, filePath,
     &                      kODF, numor, seed, angles, euler,
     &                      FILE_TXT_OUT, tauc, gammat,
     &                      grainsize)
c
c------- input convergence data (state iterations & constitutive solver)
c
      call CrystalSolverData(props, mprops,
     &                      maxIterState, maxIterNewt,
     &                      tolerState, tolerNewt, 
     &                      FILE_I, FILE_E)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE CrystalMaterialData(
     &   props, mprops, numqpt, numel,
     &   numgrn, numslip, numvtx, kODFout,
     &   kappa0, fCeDevVol, fCeVol,
     &   matProp, tauSlip,
     &   fCeDev, fCeiDev,
     &   zBar0,
     &   pBar0, qBar0,
     &   pBar0Vec, qBar0Vec,
     &   ppTBar0,
     &   gcrot0,
     &   sigfs,
     &   overstress,
     &   hardmtx,
     &   FILE_I, FILE_E, FILE_O, filePath,
     &   kODF, numor, seed, angles, euler,
     &   FILE_TXT_OUT, tauc, gammat,
     &   grainsize)

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      character*160 filePath
      integer mprops, numqpt, numel
      real*8  props(mprops)

      integer numgrn, numslip, numvtx, kODFout
      real*8  kappa0(NKAPP), fCeDevVol(NVECS), fCeVol
      real*8  matProp(NPROPS, MAX_SLIP), tauSlip(MAX_SLIP)
      real*8  fCeDev(NVECS, NVECS), fCeiDev(NVECS, NVECS)
      real*8  zBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0(DIMS, DIMS, MAX_SLIP), qBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0Vec(NVECS, MAX_SLIP), qBar0Vec(DIMS, MAX_SLIP)
      real*8  ppTBar0(NVECS, NVECS, MAX_SLIP)
      real*8  gcrot0(DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  sigfs(5, MAX_VTX)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)

      real*8  overstress(MAX_SLIP)

      integer FILE_I, FILE_E, FILE_O
      integer FILE_TXT_OUT

      integer kODF, numor, seed
      real*8  angles(DIMS, MAX_ORIEN)
      real*8  euler(DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  xtalProp(NPROPS)

      integer is, i
      character*16 scysFile
      integer crystalID

      integer SXFILE_I
      parameter(SXFILE_I=91)

      integer    length1, length2
      character  filename*160, SingleXtalFile*160, prosa*160
      
      real*8  tauc(NKAPP, MAX_GRN), gammat(MAX_SLIP, MAX_GRN)
      real*8  grainsize
c
c---------------------------------------------------------------------72
c
c---- initialize some arrays that store material data / slip system
c
      call SetTensor(kappa0,  pzero, NKAPP)
      call SetTensor(matProp, pzero, NPROPS*MAX_SLIP)
      call SetTensor(hardmtx, pzero, MAX_SLIP*MAX_SLIP)
  
      call SetTensor(xtalProp, pzero, NPROPS)
c
c---- crystal type and number of grains per IP
c
      call SetCrystalType(crystalID, numgrn, props, mprops, 
     &                    FILE_I, FILE_E, FILE_O)
c
c---- open single crystal filename

      read(FILE_I,'(a18)') SingleXtalFile

      length1 = index(filePath,' ') - 1
      length2 = index(SingleXtalFile,' ') - 1

      filename = filePath(1:length1)//SingleXtalFile(1:length2)
      open(unit=SXFILE_I, file=filename, status='unknown',
     &                                   access='sequential')
      rewind(SXFILE_I)
c
c---- elastic stiffness matrix of single crystal
c
      call SetCrystalElasticity(crystalID, fCeDev, fCeiDev, xtalProp, 
     &                          fCeDevVol, fCeVol, props, mprops,
     &                          SXFILE_I, FILE_E, FILE_O)
c
c---- thermal expansion coefficients
c
c      call SetCrystalThermal(xtalProp, SXFILE_I, FILE_E)
       read(SXFILE_I,'(a)') prosa
c
c---- set crystal geometry: slip / twinning data
c
      call SetCrystalGeometry(crystalID, numslip, numvtx, scysFile, 
     &   zBar0, pBar0, qBar0, pBar0Vec, qBar0Vec, ppTBar0, xtalProp, 
     &   matProp, hardmtx, kappa0, tauSlip, props, mprops, SXFILE_I, 
     &   FILE_E, grainsize)
c
c---- close single crystal filename
c
      close(SXFILE_I)
c
c---- crystal orientations
c
      call SetCrystalLatticeOrient(numgrn, numqpt, numel, kODFout, 
     &                             gcrot0, props, mprops,
     &                             kODF, numor, seed, angles, euler,
     &                             FILE_I, FILE_E, FILE_O, filePath,
     &                             FILE_TXT_OUT, tauc, gammat, numslip,
     &                             grainsize)
c
c---- vertices of rate independent yield surface (single crystal)
c
      if (crystalID .eq. kFCC .or. crystalID .eq. kBCC) then
         call VerticesSCYS_Cubic(numvtx, kappa0(1), sigfs, scysFile,
     &                           FILE_E)
      else  ! crystalID=kHCP
         call VerticesSCYS_HCP(numvtx, kappa0(1), sigfs, scysFile,
     &                         FILE_E, filePath)
      endif
c
c---- For HCP crystals, the overstress will keep the same initial 
c---- difference in slip system hardness between basal/prismatic
c---- and pyramidal slip systems. 
c---- For FCC/BCC crystals, it won't have any effect.
c
      do is = 1, numslip
         overstress(is) = (tauSlip(is) - 1.0) * kappa0(1)
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE CrystalSolverData(
     &   props, mprops,
     &   maxIterState, maxIterNewt,
     &   tolerState, tolerNewt,
     &   FILE_I, FILE_E
     &   )

      implicit none

      integer mprops
      real*8  props(mprops)
     
      integer maxIterState, maxIterNewt
      real*8  tolerState, tolerNewt

      integer FILE_I, FILE_E
c
c---------------------------------------------------------------------72
c
c------- number iterations and tolerance for state iterations
c
      read(FILE_I, *) maxIterState, tolerState
c       maxIterState = nint (props(19))
c       tolerState   = props(20)
c
c------- number iterations and tolerance for newton method
c
      read(FILE_I, *) maxIterNewt, tolerNewt
c       maxIterNewt = nint (props(21))
c       tolerNewt   = props(22)
c
c------- echo input data
c
      write(FILE_E, 1000) maxIterState, maxIterNewt, 
     &                    tolerState, tolerNewt      
c
c------- format
c
1000  format(/'*-----   Local Convergence Control Parameters  -----*'/,
     &        7x,'  max iters State  = ',i5 /,
     &        7x,'  max iters Newton = ',i5 /,
     &        7x,'  tolerance State  = ',e12.5 /,
     &        7x,'  tolerance Newton = ',e12.5)
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE CrystalInitializeArrays(
     &   numqpt, numel,
     &   numgrn, numslip,
     &   kappa0,
     &   gcrot0,
     &   gstress,
     &   gstress_n,
     &   gestran,
     &   gestran_n,
     &   gkappa,
     &   gkappa_n,
     &   gstatev,
     &   gstatev_n,
     &   geqvalues,
     &   ggamdot,
     &   gcrot,
     &   gcrot_n,
     &   grrot,
     &   grrot_n, gshear, ggamtar, gswitch,
     &   tauc, geps, geps_n, ggndden, ggndden_n)

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numqpt, numel
      integer numgrn, numslip
      real*8  kappa0(NKAPP)
      real*8  gcrot0 (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)

      real*8  gstress    (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstress_n  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran    (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran_n  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa     (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa_n   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev    (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev_n  (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geqvalues  (NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot    (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot      (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot_n    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot      (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot_n    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  tauc       (NKAPP, MAX_GRN)
      real*8  geps       (9, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps_n     (9, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden    (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden_n  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)

      integer ie, ip, ig
      
      real*8  gshear(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamtar(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gswitch(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)    
c
c---------------------------------------------------------------------72
c
c------- initialize arrays used in XTAL constitutive integration method
c
      do ie = 1, numel
        do ip = 1, numqpt
          do ig = 1, numgrn

            call SetTensor(gstress  (1,ig,ip,ie), pzero, NVECS)
            call SetTensor(gestran  (1,ig,ip,ie), pzero, NVECS)
            call SetTensor(gkappa   (1,ig,ip,ie), pzero, NKAPP)
            call SetTensor(gstatev  (1,ig,ip,ie), pzero, NSTAV)
            call SetTensor(gstress_n(1,ig,ip,ie), pzero, NVECS)
            call SetTensor(gestran_n(1,ig,ip,ie), pzero, NVECS)
            call SetTensor(gkappa_n (1,ig,ip,ie), pzero, NKAPP)
            call SetTensor(gstatev_n(1,ig,ip,ie), pzero, NSTAV)
            call SetTensor(geqvalues(1,ig,ip,ie), pzero, NEQVA)
            call SetTensor(ggamdot  (1,ig,ip,ie), pzero, MAX_SLIP)
            call SetTensor(geps     (1,ig,ip,ie), pzero, 9)
            call SetTensor(geps_n   (1,ig,ip,ie), pzero, 9)
            call SetTensor(ggndden  (1,ig,ip,ie), pzero, MAX_SLIP)
            call SetTensor(ggndden_n(1,ig,ip,ie), pzero, MAX_SLIP)

            call SetTensor(gcrot  (1,1,ig,ip,ie), pzero, DIMS*DIMS)
            call SetTensor(gcrot_n(1,1,ig,ip,ie), pzero, DIMS*DIMS)
            call SetTensor(grrot  (1,1,ig,ip,ie), pzero, DIMS*DIMS)
            call SetTensor(grrot_n(1,1,ig,ip,ie), pzero, DIMS*DIMS)
            
            call SetTensor(gshear  (1,ig,ip,ie), pzero, MAX_SLIP)
            call SetTensor(ggamtar (1,ig,ip,ie), pzero, MAX_SLIP)
            call SetTensor(gswitch (1,ig,ip,ie), pzero, MAX_SLIP)  

          enddo
        enddo
      enddo
c
c------- initialize state variables and rotation tensors
c
      do ie = 1, numel
        do ip = 1, numqpt
          do ig = 1, numgrn

            call EqualTensors(tauc(1,ig), gkappa_n(1,ig,ip,ie), NKAPP)
            call EqualTensors(gcrot0(1,1,ig,ip,ie), 
     &                         gcrot_n(1,1,ig,ip,ie), DIMS*DIMS)
            call EqualTensors(Ident2nd, 
     &                         grrot_n(1,1,ig,ip,ie), DIMS*DIMS)

          enddo
        enddo
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE CrystalCloseIOFiles( )

      implicit none
      include 'params_xtal.inc'
c
c---------------------------------------------------------------------72
c
c------- close files
c
      close(XTAL_I)
      close(XTAL_O)
      close(XTAL_E)
      close(XTAL_C)
      close(XTAL_S)
c
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE CrystalCloseIOFiles_2( )

      implicit none
      include 'params_xtal.inc'
c
c---------------------------------------------------------------------72
c
c------- close files
c
      close(XTAL_STRESS_O)
      close(XTAL_STRAIN_O)
      close(XTAL_EFFSS_O)
      close(XTAL_TRUESS_O)
      close(XTAL_ITER_O)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE CrystalInitializeMatrxCep( 
     &
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer i, j

      real*8  unit4(DIMS2, DIMS2), unitDev4(DIMS2, DIMS2)
      common /unit4Tensors/ unit4, unitDev4

      real*8  fDevMat5x6(NVECS, DIMS2), fDevMat6x5(DIMS2, NVECS)
      real*8  fMatTId5x6(NVECS, DIMS2)
      common /transfMatxs/ fDevMat5x6, fDevMat6x5, fMatTId5x6
c
c---------------------------------------------------------------------72
c
c----- zero 4th order unit tensors and transformation matrices

      call SetTensor(unit4, pzero, DIMS2*DIMS2)
      call SetTensor(unitDev4, pzero, DIMS2*DIMS2)

      call SetTensor(fDevMat5x6, pzero, NVECS*DIMS2)
      call SetTensor(fDevMat6x5, pzero, DIMS2*NVECS)
      call SetTensor(fMatTId5x6, pzero, NVECS*DIMS2)
c
c----- set 4th order unit tensor (I)
c
      do i = 1, DIMS2
         unit4(i, i) = pone
      enddo
c
c----- set 4th order deviatoric unit tensor (I_dev) 
c
      do i = 1, DIMS2/2
         unitDev4(i, i) = pone
         do j = 1, DIMS2/2
            unitDev4(i, j) = unitDev4(i, j) - pthird
         enddo
      enddo
 
      do i = DIMS2/2+1, DIMS2
         unitDev4(i, i) = pone
      enddo
c
c----- set tranformation matrix [T] = [tDevMat]_5x6:
c-----    {}'_5x1 = [tDevMat]_5x6 {}'_6x1
c
      fDevMat5x6(1,1) =  pone/sqr2
      fDevMat5x6(1,2) = -pone/sqr2
      fDevMat5x6(2,3) = sqr32
      fDevMat5x6(3,4) = sqr2
      fDevMat5x6(4,5) = sqr2
      fDevMat5x6(5,6) = sqr2
c
c----- set tranformation matrix [J] = [tDevMat]_6x5:
c-----    {}'_6x1 = [tDevMat]_6x5 {}'_5x1
c
      fDevMat6x5(1,1) =  pone/sqr2
      fDevMat6x5(1,2) = -pone/(sqr2*sqr3)
      fDevMat6x5(2,1) = -pone/sqr2
      fDevMat6x5(2,2) = -pone/(sqr2*sqr3)
      fDevMat6x5(3,2) = sqr23
      fDevMat6x5(4,3) = pone/sqr2
      fDevMat6x5(5,4) = pone/sqr2
      fDevMat6x5(6,5) = pone/sqr2
c
c----- needed product: [TId]_5x6 = [T]_5x6 [I_dev]_6x6
c
      call MultAxB_G(fDevMat5x6, unitDev4, fMatTId5x6, NVECS, DIMS2, 
     &               DIMS2) 

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE GlobalIterationParams( 
     &
     &   )

      implicit none
      include 'params_xtal.inc'

      integer maxIterN
      real*8  tolerN, divTolerN, zeroTolerN
      common  /ParmsGlobalN_1/ maxIterN
      common  /ParmsGlobalN_2/ tolerN, zeroTolerN, divTolerN

      logical alwaysLineSearch
      integer searchIters
      real*8  maxStepSize, orthoToler
      common  /ParmsLS_0/ alwaysLineSearch
      common  /ParmsLS_1/ searchIters
      common  /ParmsLS_2/ maxStepSize, orthoToler

      integer maxIncrCuts, quickSolveTol, quickSeriesTol
      real*8  DEQP_ref
      common  /dtContrl_1/ maxIncrCuts, quickSolveTol, quickSeriesTol
      common  /dtContrl_2/ DEQP_ref
c
c---------------------------------------------------------------------72
c
c------- params_xtal for global newton iterations
c
      read(XTAL_I, *) maxIterN
      read(XTAL_I, *) tolerN, zeroTolerN, divTolerN
c
c------- params_xtal for global line search iterations
c
      read(XTAL_I, *) alwaysLineSearch
      read(XTAL_I, *) searchIters
      read(XTAL_I, *) maxStepSize, orthoToler
c
c------- params_xtal to control the increase/decrease of time step
c
      read(XTAL_I, *) maxIncrCuts
      read(XTAL_I, *) quickSolveTol
      read(XTAL_I, *) quickSeriesTol
      read(XTAL_I, *) DEQP_ref
c
c------- echo input data 
c
      write(XTAL_E, 1000) maxIterN, tolerN, zeroTolerN, divTolerN
      write(XTAL_E, 2000) alwaysLineSearch, searchIters, maxStepSize, 
     &                   orthoToler
      write(XTAL_E, 3000) maxIncrCuts, quickSolveTol, quickSeriesTol,
     &                   DEQP_ref
c
c------- formats
c 
1000  format(/'*-----   Global Newton Iterations   -----*'/,
     &       7x,'  maxIterN    = ',i6/
     &       7x,'  tolerN      = ',e12.5/
     &       7x,'  zeroTolerN  = ',e12.5/
     &       7x,'  divTolerN   = ',e12.5)
2000  format(/'*-----   Global Line Search Iterations   -----*'/,
     &       7x,'  alwaysLS    = ',l6/
     &       7x,'  searchIters = ',i6/
     &       7x,'  maxStepSize = ',e12.5/
     &       7x,'  orthoToler  = ',e12.5)
3000  format(/'*-----   Reset Params for controlling dtime   -----*'/,
     &       7x,'  maxIncrCuts    = ',i6/
     &       7x,'  quickSolveTol  = ',i6/
     &       7x,'  quickSeriesTol = ',i6/
     &       7x,'  DEQP_ref       = ',e12.5)
          
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetCrystalType(
     &   crystalID, numgrn, props, mprops, FILE_I, FILE_E, FILE_O
     &   )

      implicit none
      include 'params_xtal.inc'

      integer crystalID, numgrn, FILE_I, FILE_E, FILE_O
      integer mprops
      real*8  props(mprops)
c
c---------------------------------------------------------------------72
c
c------- crystal type and number of grains/IP
c
      read(FILE_I, *) crystalID, numgrn
c      crystalID = nint (props(1))
c      numgrn   = nint (props(2))
      if (crystalID .ne. kFCC .and.
     &    crystalID .ne. kBCC .and.
     &    crystalID .ne. kHCP)
     &   call RunTimeError(FILE_O, 'SetCrystalType: crystalID ?')
c
c------- echo input data
      write(FILE_E, 1000) crystalID, numgrn
c
c------- format
c
1000  format(/'*-----   Crystal Type and Number of Grns per IP -----*'/,
     &        7x,'  crystal type     = ',i5 /,
     &        7x,'  number grains/ip = ',i5)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetCrystalElasticity(
     &   crystalID, fCeDev, fCeiDev, matProp, fCeDevVol, fCeVol, 
     &   props, mprops, SXFILE_I, FILE_E, FILE_O
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer crystalID, SXFILE_I, FILE_E, FILE_O
      real*8  fCeDevVol(NVECS), fCeVol
      real*8  matProp(NPROPS)
      real*8  fCeDev(NVECS,NVECS), fCeiDev(NVECS,NVECS)

      integer mprops
      real*8  props(mprops)

      integer elastID, i, j
      real*8  eMod, eNu, gMod, eBulk
      real*8  c1, c2, c3, c4, c5, const

      common /ElastProps/ gMod, eBulk
c
c---------------------------------------------------------------------72
c
c------- type of elasticity (in crystal coordinates)
c-------    elastID = kELAS_ISO : isotropic elasticity
c-------    elastID = kELAS_ANI : anisotropic elasticity
c
      read(SXFILE_I, *) elastID
c      elastID = nint (props(3))
      if (elastID .ne. kELAS_ISO .and.
     &    elastID .ne. kELAS_ANI)
     &   call RunTimeError(FILE_O, 'SetCrystalElasticity: elastID ?')
      write(FILE_E, 1000) elastID
c
c------- initialize deviatoric elastic stiffness, its inverse and
c-------  coupled deviatoric-volumetric term
c
      call SetTensor(fCeDev, pzero, NVECS*NVECS)
      call SetTensor(fCeiDev, pzero, NVECS*NVECS)
      call SetTensor(fCeDevVol, pzero, NVECS)
c
c------- read elastic constants, i.e.,
c-------     ISO   FCC  HCP
c-------     eMod  C11  C11
c-------     eMu   C12  C12
c-------           C44  C13
c-------                C33
c-------                C44
c------- and build deviatoric elastic stiffness, as well as
c------- coupled deviatoric-volumetric and volumetric terms.
c-------  For Isotropic and Cubic Crystala: 
c-------      fCeDevVol = 0, fCeVol = Bulk Modulues
c-------  For Hexagonal Closed Packed Crystals:
c-------      si_dev = fCeDev     *ei_dev + fCeDevVol*e_kk
c-------      press  = fCeDevVol^T*ei_dev + fCeVol   *e_kk
c-------   Here: 
c-------      si_dev = {(s11-s22)/sq2,sq32*s33,sq2*s12,sq2*s13,sq2*s23}
c-------      ei_dev = {(e11-e22)/sq2,sq32*e33,sq2*e12,sq2*e13,sq2*e23}
c

      if (elastID .eq. kELAS_ISO) then
c
c---------- isotropic elastic stiffness
         read(SXFILE_I, *)  eMod, eNu 
c         eMod = props(4)
c         eNu  = props(5)
         write(FILE_E, 2000) eMod, eNu
         matProp(1) = eMod / (1. + eNu) / 2.       ! gMod
         matProp(2) = eMod / (1. - 2. * eNu) / 3.  ! eBulk
         fCeDev(1,1) = 2.0*matProp(1)
         fCeDev(2,2) = 2.0*matProp(1)
         fCeDev(3,3) = 2.0*matProp(1)
         fCeDev(4,4) = 2.0*matProp(1)
         fCeDev(5,5) = 2.0*matProp(1)
         fCeVol      = matProp(2)
      elseif ((elastID .eq. kELAS_ANI .and. crystalID .eq. kFCC) .or.
     &        (elastID .eq. kELAS_ANI .and. crystalID .eq. kBCC)) then
c
c---------- anisotropic elastic stiffness, FCC or BCC
         read(SXFILE_I, *) c1, c2, c3
c         c1 = props(4)             ! C11
c         c2 = props(5)             ! C12
c         c3 = props(6)             ! C44
         write(FILE_E, 3000) c1, c2, c3
         matProp(1) = (2.0*(c1 - c2) + 6.0*c3) / 10.0
         matProp(2) = (c1 + 2.0*c2) / 3.0
         fCeDev(1,1) = c1 - c2
         fCeDev(2,2) = c1 - c2
         fCeDev(3,3) = 2.0*c3
         fCeDev(4,4) = 2.0*c3
         fCeDev(5,5) = 2.0*c3
         fCeVol      = (c1 + 2.0*c2) / 3.0
      elseif (elastID .eq. 2 .and. crystalID .eq. kHCP) then
c
c---------- anisotropic elastic stiffness, HCP
         read(SXFILE_I, *) c1, c2, c3, c4, c5
c         c1 = props(4)             ! C11
c         c2 = props(5)             ! C12
c         c3 = props(6)             ! C13
c         c4 = props(7)             ! C33
c         c5 = props(8)             ! C44
         write(FILE_E, 4000) c1, c2, c3, c4, c5
         const = c1 + c2 - 4.0*c3 + 2.0*c4
         matProp(1) = (6.0*(c1-c2) + const + 12.0*c5) / 30.0
         matProp(2) = ((c1 + c2)*c4 - 2.0*c3*c3) / const
         fCeDev(1,1) = c1 - c2
         fCeDev(2,2) = const / 3.0
         fCeDev(3,3) = c1 - c2
         fCeDev(4,4) = 2.0*c5
         fCeDev(5,5) = 2.0*c5
         fCeDevVol(2) = - dsqrt(2.d0/27.d0) * (c1 + c2 - c3 - c4)
         fCeVol       = (2.0*c1 + 2.0*c2 + 4.0*c3 + c4) / 9.0
      endif         
c
c------- effective shear and bulk modulus (to compute effective Ce)
c
      gMod  = matProp(1)
      eBulk = matProp(2)
c
c------- inverse of deviatoric elastic stiffness
c
      do i = 1, NVECS
         fCeiDev(i,i) = 1. / fCeDev(i,i)
      enddo
c
c------- echo some computed data
c
      write(FILE_E, 5000) matProp(1), matProp(2), fCeVol,
     &                    (fCeDevVol(i),i=1,NVECS)
      write(FILE_E, 6000) ((fCeDev(i,j), j=1,5), i=1,5)
c
c------- format
c
1000  format(/'*-----   Crystal Elasticity -----*'/,
     &        7x,'  elasticity type  = ',i5)
2000  format( 7x,'  young modulus    = ',e12.5/
     &        7x,'  poisson ratio    = ',e12.5)
3000  format( 7x,'  c_1              = ',e12.5/
     &        7x,'  c_2              = ',e12.5/
     &        7x,'  c_3              = ',e12.5)
4000  format( 7x,'  c_1              = ',e12.5/
     &        7x,'  c_2              = ',e12.5/
     &        7x,'  c_3              = ',e12.5/
     &        7x,'  c_4              = ',e12.5/
     &        7x,'  c_5              = ',e12.5)
5000  format( 7x,'  shear modulus    = ',e12.5/
     &        7x,'  bulk modulus     = ',e12.5/
     &        7x,'  fCeVol           = ',e12.5/
     &        7x,'  fCeDevVol(1...5) = ',5(e12.5,1x)/
     &        7x,'  CeDev : ')
6000  format((15x,5(e12.5,2x)))

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetCrystalPlasticity(
     &   crystalID, xtalProp, crss0, props, mprops, mode, 
     &   SXFILE_I, FILE_E, grainsize
     &   )

      implicit none
      include 'params_xtal.inc'

      integer crystalID, mode, SXFILE_I, FILE_E
      real*8  crss0
      real*8  xtalProp(NPROPS)

      integer mprops
      real*8  props(mprops)

      integer i, GIndex
      common  /newAddGIndex/ GIndex
      real*8  globalgrainsize(MAX_GRN)
      common  /newAddParameters/ globalgrainsize 
      real*8  grainsize
c
c---------------------------------------------------------------------72
c
c------- slip system kinetic equation: power law (~thermal activation)
c-------   xtalProp(3) xm   : strain-rate sensitivity (m)
c-------   xtalProp(4) gam0 : reference shear rate (d_0)
c
      read(SXFILE_I, *)  xtalProp(3), xtalProp(4)
c      xtalProp(3) = props( 9)
c      xtalProp(4) = props(10)
c
c------ bounds for argument of power law
c
      call BoundForArgPowLaw(xtalProp(3))
c
c------- slip system kinetic equation: drag-controlled plastic flow
c-------   xtalProp(11) bdrag : drag coefficient (B)
c
      read(SXFILE_I, *)  xtalProp(11)
c
c------ (rss/crss)_limit to switch from (powerLaw & drag) to pure drag
c------  xtalProp(12) = rss/crss_limit
c
      call LimitRatioForDrag(xtalProp)
c
c------- slip system hardening law : Voce's type
c-------   xtalProp(5) h0    : initial work hardening rate (Theta_0)
c-------   xtalProp(6) tausi : initial slip system strength (tau_0)
c-------   xtalProp(7) taus0 : saturation threshold strength (tau_s0)
c-------   xtalProp(8) xms   : strain-rate sensitivity - state evol (m')
c-------   xtalProp(9) gamss0: ref deformation rate - state evol 
c
      read(SXFILE_I, *)  (xtalProp(i), i = 5,9)
c      xtalProp(5) = props(11)
c      xtalProp(6) = props(12)
c      xtalProp(7) = props(13)
c      xtalProp(8) = props(14)
c      xtalProp(9) = props(15)
c
c------- initial value of state variables: slip system hardness crss0
c
      read(SXFILE_I, *) crss0          ! kappa0
      read(SXFILE_I, *) grainsize      
c      crss0 = props(16)
      xtalProp(10) = crss0
	  grainsize = globalgrainsize(GIndex)
	  
c
c------- echo date
c
      write(FILE_E, 1000) mode, xtalProp(3), xtalProp(4), xtalProp(11),
     &                    xtalProp(12)
      write(FILE_E, 2000) (xtalProp(i), i=5,9), crss0, grainsize
	  !write(XTAL_O, *) 'grainsize', grainsize
	  !write(XTAL_O, *) 'globalgrainsize', globalgrainsize(GIndex)
c
c------- format
c
1000  format(/'*-----   Crystal Plasticity, mode', i3, '-----*'/,
     &        7x,'  Kinetic Equation : '/,
     &        7x,'  m                = ',e12.5/
     &        7x,'  gam0             = ',e12.5/
     &        7x,'  bdrag            = ',e12.5/
     &        7x,'  (rss/crss)_limit = ',e12.5)
2000  format( 7x,'  Hardening Law : '/,
     &        7x,'  h0               = ',e12.5/
     &        7x,'  tausi            = ',e12.5/
     &        7x,'  taus0            = ',e12.5/
     &        7x,'  xms              = ',e12.5/
     &        7x,'  gamss0           = ',e12.5/
     &        7x,'  crss0 (kappa0)   = ',e12.5/
     &        7x,'  grain size       = ',e12.5)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetCrystalSlipGeometry(
     &   crystalID, numslip, numvtx, tauSlip, scysFile, zBar0, 
     &   pBar0, qBar0, pBar0Vec, qBar0Vec, ppTBar0, FILE_E, filePath
     &   )
      
      implicit none
      include 'params_xtal.inc'

      character*160 filePath
      integer crystalID, numslip, numvtx, FILE_E
      real*8  tauSlip(MAX_SLIP)
      real*8  zBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0(DIMS, DIMS, MAX_SLIP), qBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0Vec(NVECS, MAX_SLIP), qBar0Vec(DIMS, MAX_SLIP)
      real*8  ppTBar0(NVECS, NVECS, MAX_SLIP)
      character*(*)  scysFile

      real*8  vecM(DIMS, MAX_SLIP), vecS(DIMS, MAX_SLIP)
c
c---------------------------------------------------------------------72
c
c------- Set up slip system vectors based on type of crystal
c
      if (crystalID .eq. kFCC) 
     &    call SetSlipSystemFCC(numslip, vecM, vecS, scysFile, tauSlip,
     &                          FILE_E)

      if (crystalID .eq. kBCC) 
     &    call SetSlipSystemBCC(numslip, vecM, vecS, scysFile, tauSlip,
     &                          FILE_E)
    

      if (crystalID .eq. kHCP) 
     &    call SetSlipSystemHCP(numslip, numvtx, vecM, vecS, scysFile, 
     &                          tauSlip, FILE_E, filePath)
c
c------- Set up Schmid tensors/vectors for each slip system
c
      call SetSlipSystemTensors(numslip, vecM, vecS, zBar0, pBar0, 
     &                          qBar0, pBar0Vec, qBar0Vec, ppTBar0,
     &                          FILE_E)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetSlipSystemFCC(
     &   numslip, vecM, vecS, scysFile, tauSlip, FILE_E
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip, FILE_E
      real*8  vecM(DIMS, MAX_SLIP), vecS(DIMS, MAX_SLIP)
      real*8  tauSlip(MAX_SLIP)
      character*(*)  scysFile

      integer is, i
      real*8 indexM(3,12), indexS(3,12)
      real*8 sDotm(12)
      real*8 InnerProductVec

      data indexM /1.,  1., -1.,
     &             1.,  1., -1.,
     &             1.,  1., -1.,
     &             1., -1., -1.,
     &             1., -1., -1.,
     &             1., -1., -1.,
     &             1., -1.,  1.,
     &             1., -1.,  1.,
     &             1., -1.,  1.,
     &             1.,  1.,  1.,
     &             1.,  1.,  1.,
     &             1.,  1.,  1./
      data indexS /0.,  1.,  1.,
     &             1.,  0.,  1.,
     &             1., -1.,  0.,
     &             0.,  1., -1.,
     &             1.,  0.,  1.,
     &             1.,  1.,  0.,
     &             0.,  1.,  1.,
     &             1.,  0., -1.,
     &             1.,  1.,  0.,
     &             0.,  1., -1.,
     &             1.,  0., -1.,
     &             1., -1.,  0./

c
c---------------------------------------------------------------------72
c
c------- set number of slip systems and slip system reference stress
c
      numslip = kSlipFCC
      call SetTensor(tauSlip, pone, numslip)
c
c------- slip system normals and slip directions: unit vectors
c
      do is = 1, numslip
         call UnitVector(indexS(1,is), vecS(1,is), DIMS)
         call UnitVector(indexM(1,is), vecM(1,is), DIMS)
      enddo
c
c------- file with RI vertex stresses
c
      scysFile = 'vert_fcc.028.00'     ! not used
c
c------- check normality of vecS and vecM
c
      do is = 1, numslip
         sDotm(is) = InnerProductVec(vecS(1,is), vecM(1,is), DIMS)
      enddo
c
c------- echo values
c
      write(FILE_E, 1000) numslip, scysFile
      do is = 1, numslip
         write(FILE_E, 2000) is, tauSlip(is), (vecS(i,is),i=1,3), 
     &                       (vecM(i,is),i=1,3), sDotm(is)
      enddo
c
c------- format
c
1000  format(/'*-----   Slip Systems for FCC -----*'/,
     &        7x,'  number slip syst = ',i4/
     &        7x,'  vtx stress file  = ',a15/
     &        7x,'  SS#',5x,'tau',18x,'vecS',30x,'vecM',20x,'SdotM')
2000  format( 7x, i4, 4x, f5.2, 4x, 3f10.5, 4x, 3f10.5, 4x, f10.5)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetSlipSystemBCC(
     &   numslip, vecM, vecS, scysFile, tauSlip, FILE_E
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'
      
      integer numslip, FILE_E
      real*8  vecm(DIMS, MAX_SLIP), vecs(DIMS, MAX_SLIP)
      real*8  tauSlip(MAX_SLIP)
      character*(*)  scysFile

      integer is, i
      real*8 indexM(3,12), indexS(3,12)
      real*8 sDotm(12)
      real*8 InnerProductVec

      data indexM /0.,  1.,  1.,
     &             1.,  0.,  1.,
     &             1., -1.,  0.,
     &             0.,  1., -1.,
     &             1.,  0.,  1.,
     &             1.,  1.,  0.,
     &             0.,  1.,  1.,
     &             1.,  0., -1.,
     &             1.,  1.,  0.,
     &             0.,  1., -1.,
     &             1.,  0., -1.,
     &             1., -1.,  0./
      data indexS /1.,  1., -1.,
     &             1.,  1., -1.,
     &             1.,  1., -1.,
     &             1., -1., -1.,
     &             1., -1., -1.,
     &             1., -1., -1.,
     &             1., -1.,  1.,
     &             1., -1.,  1.,
     &             1., -1.,  1.,
     &             1.,  1.,  1.,
     &             1.,  1.,  1.,
     &             1.,  1.,  1./
c
c---------------------------------------------------------------------72
c
c------- set number of slip systems and slip system reference stress
c
      numslip = kSlipBCC
      call SetTensor(tauSlip, pone, numslip)
c
c------- slip system normals and slip directions: unit vectors
c
      do is = 1, numslip
         call UnitVector(indexS(1,is), vecS(1,is), DIMS)
         call UnitVector(indexM(1,is), vecM(1,is), DIMS)
      enddo
c
c------- file with RI vertex stresses
c
      scysFile = 'vert_bcc.028.00'     ! not used
c
c------- check normality of vecS and vecM
c
      do is = 1, numslip
         sDotm(is) = InnerProductVec(vecS(1,is), vecM(1,is), DIMS)
      enddo
c
c------- echo values
c
      write(FILE_E, 1000) numslip, scysFile
      do is = 1, numslip
         write(FILE_E, 2000) is, tauSlip(is), (vecS(i,is),i=1,3), 
     &                       (vecM(i,is),i=1,3), sDotm(is)
      enddo
c
c------- format
c
1000  format(/'*-----   Slip Systems for BCC -----*'/,
     &        7x,'  number slip syst = ',i4/
     &        7x,'  vtx stress file  = ',a15/
     &        7x,'  SS#',5x,'tau',18x,'vecS',30x,'vecM',20x,'SdotM')
2000  format( 7x, i4, 4x, f5.2, 4x, 3f10.5, 4x, 3f10.5, 4x, f10.5)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetSlipSystemHCP(
     &   numslip, numvtx, vecM, vecS, scysFile, tauSlip, 
     &   FILE_E, filePath
     &   )

      implicit  none
      include   'params_xtal.inc'

      character*160 filePath
      integer   numvtx, numslip, FILE_E
      real*8    vecM(DIMS, MAX_SLIP), vecS(DIMS, MAX_SLIP)
      real*8    tauSlip(MAX_SLIP)

      integer   nmodesx, nmodes, kount, modex, nsmx, mode(5)
      integer   i, j, is, ix, io, nm
      real*8    rca, tau
      real*8    vecm4(48, 4), vecs4(48, 4), sDotm(48)
      character prosa(9)*8
      character*(*) scysFile

      real*8 InnerProductVec

      integer      lengthFile
      character*160 filename
c
c---------------------------------------------------------------------72
c
c     need input file 'slip_hcp.in'

      lengthFile = index(filePath,' ') - 1
      filename   = filePath(1:lengthFile)//'slip_hcp.in'

      io = 99
c      open(unit=io, file='slip_hcp.in', status='old')
      open(unit=io, file=filename, status='old')

      read(io,*) rca
      read(io,*) nmodesx
      read(io,*) nmodes
      read(io,*) (mode(i),i=1,nmodes)
      read(io,*) numvtx
      read(io,5) scysFile
5     format(3x,a15)
c      print *, scysFile
c      print *, '  nVertx for HCP = ', numvtx
c
      kount=1
      i=0
      do nm=1,nmodesx
         read(io,6) prosa
6        format(10a8)
         read(io,*) modex,nsmx,tau
         if(modex.ne.mode(kount)) then
            do ix=1,nsmx
               read(io,6) prosa
            enddo
         else
            kount=kount+1
            do is=1,nsmx
               i=i+1
               read(io,*) (vecm4(i,j),j=1,4),(vecs4(i,j),j=1,4)
               vecM(1,i)= vecm4(i,1)
               vecM(2,i)=(vecm4(i,1)+2.*vecm4(i,2))/sqrt(3.)
               vecM(3,i)= vecm4(i,4)/rca
               vecS(1,i)= 3./2.*vecs4(i,1)
               vecS(2,i)=(vecs4(i,1)/2.+vecs4(i,2))*sqrt(3.)
               vecS(3,i)= vecs4(i,4)*rca
               tauSlip(i) = tau
            enddo
         endif
      enddo

      close(io)
      numslip=i
c
c------- slip system normals and slip directions: unit vectors
c
      do is = 1, numslip
         call UnitVector(vecS(1,is), vecS(1,is), DIMS)
         call UnitVector(vecM(1,is), vecM(1,is), DIMS)
      enddo
c
c------- check normality of vecS and vecM
c
      do is = 1, numslip
         sDotm(is) = InnerProductVec(vecS(1,is), vecM(1,is), DIMS)
      enddo
c
c------- echo values
c
      write(FILE_E, 1000) numslip, scysFile
      do is = 1, numslip
         write(FILE_E, 2000) is, tauSlip(is), (vecS(i,is),i=1,3), 
     &                       (vecM(i,is),i=1,3), sDotm(is)
      enddo
c
c------- format
c
1000  format(/'*-----   Slip Systems for HCP -----*'/,
     &        7x,'  number slip syst = ',i4/
     &        7x,'  vtx stress file  = ',a15/
     &        7x,'  SS#',5x,'tau',18x,'vecS',30x,'vecM',20x,'SdotM')
2000  format( 7x, i4, 4x, f5.2, 4x, 3f10.5, 4x, 3f10.5, 4x, f10.5)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetSlipSystemTensors(
     &   numslip, vecM, vecS, zBar0, pBar0, qBar0, pBar0Vec, qBar0Vec,
     &   ppTBar0, FILE_E
     &   )

      implicit  none
      include  'params_xtal.inc'

      integer numslip, FILE_E
      real*8  vecM(DIMS, MAX_SLIP), vecS(DIMS, MAX_SLIP)
      real*8  zBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0(DIMS, DIMS, MAX_SLIP), pBar0Vec(NVECS, MAX_SLIP)
      real*8  qBar0(DIMS, DIMS, MAX_SLIP), qBar0Vec(DIMS, MAX_SLIP)
      real*8  ppTBar0(NVECS, NVECS, MAX_SLIP)

      integer is, j, k
c
c---------------------------------------------------------------------72
c
c------- Schmid Orientation Tensors/Vectors in Crystal Coordinates  
c
      do is = 1, numslip
c
c---------- Tensor zBar0
         call OuterProductVec(vecS(1,is), vecM(1,is), zBar0(1,1,is), 
     &                        DIMS)
c
c---------- Symmetric and Skew parts of zBar0
         call SymmetrizeTensor(zBar0(1,1,is), pBar0(1,1,is), DIMS)
         call SkewSymmetrizeTensor(zBar0(1,1,is), qBar0(1,1,is), DIMS)
c
c---------- Vector form for pBar0 and qBar0
         call Mat3x3ToVec5x1Symm(pBar0(1,1,is), pBar0Vec(1,is), DIMS)
         call Mat3x3ToVec3x1Skew(qBar0(1,1,is), qbar0Vec(1,is), DIMS)
c
c---------- Form matrix {P}{P}^T
         call OuterProductVec(pBar0Vec(1,is), pBar0Vec(1,is),
     &                        ppTBar0(1,1,is), NVECS)

      enddo
c
c------- echo Schmid orientation tensors/vectors
c
      write(FILE_E, 1000)
      do is = 1, numslip
         write(FILE_E, 1500) is
         do j = 1, DIMS
            write(FILE_E, 2000) (zBar0(j,k,is),k=1,DIMS), 
     &                          (pBar0(j,k,is),k=1,DIMS), 
     &                          (qBar0(j,k,is),k=1,DIMS)
         enddo
      enddo

      write(FILE_E, 3000)
      do is = 1, numslip
         write(FILE_E, 4000) is, (pBar0Vec(j,is),j=1,NVECS), 
     &                           (qBar0Vec(j,is),j=1,DIMS)
      enddo
c
c------- format
c
1000  format(/'*-----   Schmid Tensors -----*')
1500  format( 7x, ' SS # ', i2 / 
     &       10x, ' zBar0, pBar0, qBar0: ') 
2000  format(10x, 3(3x, 3f9.5))
3000  format(/'*-----   Schmid Vectors -----*'/
     &        7x, ' SS#, pBar0Vec, qBar0Vec: ') 
4000  format( 7x, i3, 2(3x, 5f9.5))

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetCrystalLatticeOrient(
     &   numgrn, numqpt, numel, kODFout, gcrot0, props, mprops,
     &   kODF, numor, seed, angles, euler, FILE_I, FILE_E, FILE_O, 
     &   filePath, FILE_TXT_OUT, tauc, gammat, numslip, grainsize
     &   )

      implicit none
      include 'params_xtal.inc'

      integer numgrn, numqpt, numel, kODFout, numslip
      real*8  gcrot0(DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)

      integer mprops
      real*8  props(mprops)

      character*160 filePath
      integer kODF, numor, seed, FILE_I, FILE_E, FILE_O
      integer FILE_TXT_OUT
      real*8  angles(DIMS, MAX_ORIEN)
      real*8  euler(DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)

      integer iikc, iidr, i, j
      real*8  pi, pi180, piby2, ph

      integer    length1, length2, length3
      character  textureFile*20, filename*160
      character  critFile*20
      
      real*8  tauc(NKAPP, MAX_GRN), gammat(MAX_SLIP, MAX_GRN)
      integer ig, is
      real*8  aa, bb, c, M, b, grainsize

      
c
c---------------------------------------------------------------------72
c
      pi = 4.0 * datan(1.0d+00)
      pi180 = pi/180.

      c = 1.2
      M = 50e3
      b = 0.25
c
c------------------------------------- Data read from FILE_I
c
c------- read code to assign angles at each GR/IP/EL
c
      read(FILE_I, *) kODF
c      kODF = nint (props(17))
c
c------- read multiples of increment to output texture
c
      read(FILE_I, *) kODFout
c      kODFout = nint (props(18))
c
c------- read root name for input/output texture file
c
      read(FILE_I,'(a18)') textureFile
c      textureFile = 'texture'
      critFile = 'crss'

      length1 = index(filePath,' ') - 1
      length2 = index(textureFile,' ') - 1
      length3 = index(critFile,' ') - 1

      filename = filePath(1:length1)//textureFile(1:length2)//'.txti'
      open(unit=XTAL_TXT_IN, file=filename, status='unknown',
     &                                   access='sequential')
      rewind(XTAL_TXT_IN)

      filename = filePath(1:length1)//textureFile(1:length2)//'.txto'
      open(unit=FILE_TXT_OUT, file=filename, status='unknown')

c
c------- echo data read from FILE_I
c
      write(FILE_E, 1000) kODF, kODFout, textureFile
c
c------------------------------------- Data read from XTAL_TXT_IN
c
c------- number of orientations in file
c
      read(XTAL_TXT_IN, *)  numor
      if (numor .lt. numgrn) 
     &   call RunTimeError(FILE_O, 
     &                       'SetCrystalLatticeOrient: numor < numgrn')
c
c------- read the flag for angle convention and set some constants
c-------   iikc = 0 : angles input in Kocks convention :  (psi,the,phi)
c                 1 : angles input in Canova convention : (ph,th,om)
c                 ph = 90 + phi; th = the; om = 90 - psi
c-------   iidr = 0 : angles input in degrees
c-------          1 : angles input in radians
      read(XTAL_TXT_IN, *) iikc, iidr
      piby2 = 90.
      if (iidr .eq. 1) piby2 = pi / 2.0
c
c-------- read Euler angles in file & storage them in : (Kocks, radians)
c
      do i = 1, numor
         read(XTAL_TXT_IN, *) (angles(j,i), j=1,DIMS)

         if (iikc .eq. 1) then
            ph = angles(1,i)
            angles(1,i) = piby2 - angles(3,i)
            angles(3,i) = ph - piby2
         endif

         if (iidr .eq. 0) 
     &     call SetToScaledTensor(pi180, angles(1,i), angles(1,i), DIMS)
      enddo

      if (kODF .eq. kODF_from_abq) then
         angles(1,1) = props(3)
         angles(2,1) = props(4)
         angles(3,1) = props(5)

         if (iikc .eq. 1) then
            ph = angles(1,1)
            angles(1,1) = piby2 - angles(3,1)
            angles(3,1) = ph - piby2
         endif

         if (iidr .eq. 0)
     &     call SetToScaledTensor(pi180, angles(1,1), angles(1,1), DIMS)
      endif
c
c------- seed for a random angle distribution in polycrystal (needed?)
c
      read(XTAL_TXT_IN, *) seed
c
c------- close input file with texture info
c
      close(XTAL_TXT_IN)
c
c------- read critically resolved shear stress
c  
c
c------------------------------------- Assign ODF to each GRs/IPs/ELs
c
c      call AssignCrystalODF(numel, numqpt)  ! No of args has changed!
c
c------- formats
c
1000  format(/'*-----   Lattice Orientation -----*'/
     &        7x,'  kODF             = ',i6/
     &        7x,'  kODFout          = ',i6/
     &        7x,'  root name for I/O txt file  = ',a18)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE VerticesSCYS_Cubic(
     &   numvtx, kappa0, sigfs, scys_cub, FILE_E
     &   )

      implicit none
      include  'params_xtal.inc'
      include  'numbers.inc'

      integer numvtx, FILE_E
      real*8  kappa0, sigfs(5,MAX_VTX)
      character*(*)  scys_cub

      integer nvtx, ivtx, j
      real*8  sigca(5,28), strVtx(3,3)

      data nvtx   /28/
      data sigca  /-1.73204, -1.00000,   .00000,   .00000,   .00000,
     &              1.73204, -1.00000,   .00000,   .00000,   .00000,
     &               .00000,  2.00000,   .00000,   .00000,   .00000,
     &               .00000,   .00000,  1.73204,  1.73204,  1.73204,
     &               .00000,   .00000, -1.73204,  1.73204,  1.73204,
     &               .00000,   .00000,  1.73204, -1.73204,  1.73204,
     &               .00000,   .00000,  1.73204,  1.73204, -1.73204,
     &               .00000,   .00000,  3.46410,   .00000,   .00000,
     &               .00000,   .00000,   .00000,  3.46410,   .00000,
     &               .00000,   .00000,   .00000,   .00000,  3.46410,
     &              -.86602,  -.50000,   .00000,  1.73204,  1.73204,
     &              -.86602,  -.50000,   .00000, -1.73204, -1.73204,
     &              -.86602,  -.50000,   .00000,  1.73204, -1.73204,
     &              -.86602,  -.50000,   .00000, -1.73204,  1.73204,
     &               .86602,  -.50000,  1.73204,   .00000,  1.73204,
     &               .86602,  -.50000, -1.73204,   .00000, -1.73204,
     &               .86602,  -.50000, -1.73204,   .00000,  1.73204,
     &               .86602,  -.50000,  1.73204,   .00000, -1.73204,
     &               .00000,  1.00000,  1.73204,  1.73204,   .00000,
     &               .00000,  1.00000, -1.73204, -1.73204,   .00000,
     &               .00000,  1.00000,  1.73204, -1.73204,   .00000,
     &               .00000,  1.00000, -1.73204,  1.73204,   .00000,
     &               .86602, -1.50000,  1.73204,   .00000,   .00000,
     &               .86602, -1.50000, -1.73204,   .00000,   .00000,
     &               .86602,  1.50000,   .00000,  1.73204,   .00000,
     &               .86602,  1.50000,   .00000, -1.73204,   .00000,
     &             -1.73204,   .00000,   .00000,   .00000,  1.73204,
     &             -1.73204,   .00000,   .00000,   .00000, -1.73204/
c
c---------------------------------------------------------------------72
c
c------- sigca: vertices of single crystal yield surface
c------- {sigca}={-(11-22)/sqr2,sqr32*33,sqr2*32,sqr2*31,sqr2*21}

      call SetTensor(strVtx, pzero, DIMS)

      numvtx = nvtx
      do ivtx = 1, numvtx

c         strVtx(1,1) = -(sigca(1,ivtx) + sigca(2,ivtx) / sqr3) / sqr2
c         strVtx(2,2) =  (sigca(1,ivtx) - sigca(2,ivtx) / sqr3) / sqr2
c         strVtx(3,3) = sigca(2,ivtx) * sqr2 / sqr3
c         strVtx(2,1) = sigca(5,ivtx) / sqr2
c         strVtx(3,1) = sigca(4,ivtx) / sqr2
c         strVtx(3,2) = sigca(3,ivtx) / sqr2

         sigfs(1, ivtx) = -sigca(1, ivtx)
         sigfs(2, ivtx) =  sigca(2, ivtx)
         sigfs(3, ivtx) =  sigca(5, ivtx)
         sigfs(4, ivtx) =  sigca(4, ivtx)
         sigfs(5, ivtx) =  sigca(3, ivtx)

c         call mat3x3ToVec5x1Symm(strVtx, sigfs(1,ivtx), DIMS)
         call SetToScaledTensor(kappa0, sigfs(1,ivtx), sigfs(1,ivtx), 5)

      enddo
c
c------- echo vettices
c
      write(FILE_E, 1000) numvtx
      do ivtx = 1, numvtx
         write (FILE_E, 2000) ivtx, (sigfs(j,ivtx),j=1,5)
      enddo
c
c------- formats
c
1000  format(/'*-----   Vertices of SCYS - Cubic Crystal -----*'/
     &        7x,'  No of vertices   = ',i4/
     &        7x,'  ivtx     sigfs(1...5)')
2000  format( 7x, i5, 3x, 5f10.5)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE VerticesSCYS_HCP(
     &   numvtx, kappa0, sigfs, scys_hcp, FILE_E, filePath
     &   )

      implicit  none
      include   'params_xtal.inc'

      character*160 filePath
      integer   numvtx, FILE_E
      real*8    kappa0, sigfs(5,MAX_VTX)
      character*(*)      scys_hcp

      integer   io, ivtx, j

      integer      length1, length2
      character*160 filename
c
c---------------------------------------------------------------------72
c
c------- sigca: vertices of single crystal yield surface
c-------   {sigca}={(11-22)/sqr2,sqr32*33,sqr2*21,sqr2*31,sqr2*32}
c
      length1 = index(filePath,' ') - 1
      length2 = index(scys_hcp,' ') - 1

      filename = filePath(1:length1)//scys_hcp(1:length2)

      io = 99
c      open(io, file = scys_hcp, status='old')
      open(io, file = filename, status='old')

      do ivtx = 1, numvtx
         read(io,*) (sigfs(j,ivtx),j=1,5) 
         call SetToScaledTensor(kappa0, sigfs(1,ivtx), sigfs(1,ivtx), 5)
      enddo

      close(io)
c
c------- echo vettices
c
      write(FILE_E, 1000) numvtx
      do ivtx = 1, numvtx
         write (FILE_E, 2000) ivtx, (sigfs(j,ivtx),j=1,5)
      enddo
c
c------- formats
c
1000  format(/'*-----   Vertices of SCYS - HCP Crystal -----*'/
     &        7x,'  No of vertices   = ',i4/
     &        7x,'  ivtx   sigfs(1...5)')
2000  format( 7x, i5, 3x, 5f13.5)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE AssignCrystalODF(numel, numqpt, noel, npt,
     &   numgrn, numslip, kODF,
     &   numor, seed,
     &   gcrot0,
     &   angles,
     &   euler,
     &   FILE_O, FILE_TXT_OUT)

      implicit none
      include  'params_xtal.inc'

      integer numel, numqpt, noel, npt
      integer FILE_O, FILE_TXT_OUT

      integer numgrn, numslip, kODF
      integer numor, seed
      real*8  gcrot0(DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  angles(DIMS, MAX_ORIEN)
      real*8  euler(DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)

      real    ran1_nr

      integer ie, ip, ig, i, random
      real*8  pi, pi180

      integer seed_nr
      save    seed_nr
c
c---------------------------------------------------------------------72
c
      pi = 4.0 * datan(1.0d+00)
      pi180 = pi/180.

      if (noel .eq. 1 .and. npt .eq. 1) seed_nr = -seed
      if (noel .eq. 1 .and. npt .eq. 1) write(FILE_TXT_OUT, 2000)
c
c------- same ODF in all ELs and IPs
      if (kODF .eq. kODF_same_all) then
         do ie = 1, numel
            do ip = 1, numqpt
               do ig =1, numgrn
                  call EqualTensors(angles(1,ig), 
     &                                        euler(1,ig,ip,ie), DIMS)
               enddo
            enddo
         enddo
c
c------- different ODF in all ELs; same ODF in all IPs
      elseif (kODF .eq. kODF_diff_elems) then
         do ie = 1, numel
            do ip = 1, numqpt
               random = int (ran1_nr(seed_nr) * numor)
               do ig =1, numgrn
                  call EqualTensors(angles(1,random), 
     &                                        euler(1,ig,ip,ie), DIMS)
               enddo
            enddo
         enddo
c
c------- different ODF in all ELs and IPs
      elseif (kODF .eq. kODF_diff_inpts) then
         do ie = 1, numel
            do ip = 1, numqpt
               do ig =1, numgrn
                  random = int (ran1_nr(seed_nr) * numor)
                  call EqualTensors(angles(1,random), 
     &                                        euler(1,ig,ip,ie), DIMS)
               enddo
            enddo
         enddo
c
c------- ODF read from file for all ELs, IPs, and GRs
      elseif (kODF .eq. kODF_from_file) then
c
c---------- for specific cases, e.g., BICRYSTAL  
c         if (numel .eq. numor .and. numgrn .eq. 1) then
c            do ie = 1, numel
c               do ip = 1, numqpt
c                  do ig =1, numgrn
c                     call EqualTensors(angles(1,ie), 
c     &                                        euler(1,ig,ip,ie), DIMS)
c                  enddo
c               enddo
c            enddo
cc
c---------- assign to each GR/IP/EL from file
c         elseif (numel*numqpt*numgrn .eq. numor) then
            i = 0
            do ie = 1, numel
               do ip = 1, numqpt
                  do ig =1, numgrn
                     i = i + 1
c                     call EqualTensors(angles(1,i), 
                     call EqualTensors(angles(1,noel), 
     &                                        euler(1,ig,ip,ie), DIMS)
                  enddo
               enddo
            enddo
c         else
c            call RunTimeError(FILE_O, 
c     &        'setCrystalLatticeOrient: kODF_from_file, check Inp file')
c         endif
c
c------- read ODF from abaqus input file: 1 grain / IP
      elseif (kODF .eq. kODF_from_abq) then
         do ie = 1, numel
            do ip = 1, numqpt
               do ig =1, numgrn
                  call EqualTensors(angles(1,ig),
     &                                        euler(1,ig,ip,ie), DIMS)
               enddo
            enddo
         enddo
      else
         call RunTimeError(FILE_O, 'setCrystalLatticeOrient: Bad kODF')
      endif
c
c------- build rotation matrices C0: {x}_sm = [C0] {x}_cr 
c
      do ie = 1, numel
         do ip = 1, numqpt
            do ig =1, numgrn
               call AnglesToRotMatrix(euler(1,ig,ip,ie), 
     &                                 gcrot0(1,1,ig,ip,ie), DIMS)
            enddo
         enddo
      enddo
c
c------- write initial assigned orientations
c-------  note that numel & numqpt are both dummy variables (=1)
c
c      write(*,*)
      do ie = 1, numel
         do ip = 1, numqpt
            do ig =1, numgrn
               write(FILE_TXT_OUT, 3000) 
     &                 (euler(i,ig,ip,ie)/pi180,i=1,3), ig, npt, noel
            enddo
         enddo
      enddo
c
c------- formats
c
2000  format(/'*-----   Initial Assigned Orientations -----*'/
     &        4x,'ang1',4x,'ang2',4x,'ang3',4x,' igrn',4x,'intpt',
     &        3x,'ielem ')
3000  format( 3f8.2, 3(3x,i5))

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetCrystalGeometry(
     &   crystalID, numslip, numvtx, scysFile, zBar0, pBar0, qBar0, 
     &   pBar0Vec, qBar0Vec, ppTBar0, xtalProp, matProp, hardmtx, 
     &   kappa0, tauSlip, props, mprops, SXFILE_I, FILE_E, grainsize
     &   )
      
      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      character*(*)  scysFile
      integer crystalID, numslip, numvtx, SXFILE_I, FILE_E
      real*8  kappa0(NKAPP)
      real*8  zBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0(DIMS, DIMS, MAX_SLIP), qBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0Vec(NVECS, MAX_SLIP), qBar0Vec(DIMS, MAX_SLIP)
      real*8  ppTBar0(NVECS, NVECS, MAX_SLIP)
      real*8  xtalProp(NPROPS), matProp(NPROPS, MAX_SLIP)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP), tauSlip(MAX_SLIP)

      integer mprops
      real*8  props(mprops)

      integer   nmodesx, nmodes, kount, modex, nsmx
      integer   i, j, is, js, nm, im, jm, isensex
      real*8    rca, crss0
      real*8    twshx, isectwx, thres1x, thres2x
      real*8    vecM(DIMS, MAX_SLIP), vecS(DIMS, MAX_SLIP)
      real*8    vecm4(48, 4), vecs4(48, 4), sDotm(48)
      character prosa*160
      common /SlipSystem/ vecM, vecS

      real*8 InnerProductVec

      integer   NMFILE
      parameter (NMFILE=12)
      integer   mode(NMFILE), nsm(NMFILE)
      real*8    hselfx(NMFILE), hlatex(NMFILE, NMFILE)
      real*8    grainsize
c
c---------------------------------------------------------------------72
c
c---- zero out some arrays
c 
      call SetTensor(vecM,  pzero, DIMS*MAX_SLIP)
      call SetTensor(vecS,  pzero, DIMS*MAX_SLIP)
      call SetTensor(vecm4, pzero, 48*4)
      call SetTensor(vecs4, pzero, 48*4)
c
c---- read SXFILE_I to set slip / twinning data
c
      read(SXFILE_I,*) rca
      read(SXFILE_I,*) nmodesx               ! total # modes in file
      read(SXFILE_I,*) nmodes                ! # modes used (active)
      read(SXFILE_I,*) (mode(i),i=1,nmodes)  ! labels of active modes
      read(SXFILE_I,*) numvtx
      read(SXFILE_I,'(a15)') scysFile

      kount=1                         ! counter for active modes
      i=0                             ! counter for # slip/twin systems
      do nm=1,nmodesx

         read(SXFILE_I,'(a)') prosa
         read(SXFILE_I,*) modex,nsmx,isensex
         if(modex.ne.mode(kount)) then
c........... skip data
            read(SXFILE_I,'(a)') prosa
            read(SXFILE_I,'(a)') prosa
            read(SXFILE_I,'(a)') prosa
            read(SXFILE_I,'(a)') prosa
            read(SXFILE_I,'(a)') prosa
            read(SXFILE_I,'(a)') prosa
            do is=1,nsmx
               read(SXFILE_I,*)
            enddo
         else
c........... PlastProp: 4 lines; Twindata: 1 line; LatentHard: 1 line
            call SetCrystalPlasticity(crystalID, xtalProp, crss0, 
     &                         props, mprops, modex, SXFILE_I, FILE_E,
     &                         grainsize)
            read(SXFILE_I,*) twshx, isectwx, thres1x, thres2x
            read(SXFILE_I,*) (hlatex(kount, jm), jm=1, nmodes)
            hselfx(kount) = 1.0
            nsm(kount) = nsmx

c........... indices for each slip/twin mode
            if (crystalID .eq. kFCC .or. crystalID .eq. kBCC) then
              do is=1,nsmx
                i=i+1
                kappa0(i) = crss0
                call EqualTensors(xtalProp, matProp(1,i), NPROPS)
                read(SXFILE_I,*) (vecM(j,i),j=1,3),(vecS(j,i),j=1,3)
                call UnitVector(vecS(1,i), vecS(1,i), DIMS)
                call UnitVector(vecM(1,i), vecM(1,i), DIMS)
              enddo
            endif

            if (crystalID .eq. kHCP) then
              do is=1,nsmx
                i=i+1
                kappa0(i) = crss0
                call EqualTensors(xtalProp, matProp(1,i), NPROPS)
                read(SXFILE_I,*) (vecm4(i,j),j=1,4),(vecs4(i,j),j=1,4)
                vecM(1,i)= vecm4(i,1)
                vecM(2,i)=(vecm4(i,1)+2.*vecm4(i,2))/sqrt(3.)
                vecM(3,i)= vecm4(i,4)/rca
                vecS(1,i)= 3./2.*vecs4(i,1)
                vecS(2,i)=(vecs4(i,1)/2.+vecs4(i,2))*sqrt(3.)
                vecS(3,i)= vecs4(i,4)*rca
                call UnitVector(vecS(1,i), vecS(1,i), DIMS)
                call UnitVector(vecM(1,i), vecM(1,i), DIMS)
              enddo
            endif
            kount=kount+1
         endif

      enddo

      numSlip=i
c
c------- ratio of slip system's kappa0: kappa0(is)/kappa0(1) 
c------- check normality of vecS and vecM
c
      do is = 1, numSlip
         tauSlip(is) = kappa0(is)/kappa0(1)
         sDotm(is)   = InnerProductVec(vecS(1,is), vecM(1,is), DIMS)
      enddo
c
c------- set up latent hardening matrix
c
      i=0
      do im = 1, nmodes
        do is = 1, nsm(im)
          i=i+1
          j=0
          do jm = 1, nmodes
            do js = 1, nsm(jm)
              j=j+1
              hardmtx(i,j) = hlatex(im,jm)
            enddo
          enddo
          hardmtx(i,i)=hselfx(im)
        enddo
      enddo
c
c------- echo values
c
      write(FILE_E, 1000) crystalID, numslip, scysFile
      do is = 1, numslip
         write(FILE_E, 2000) is, tauSlip(is), (vecS(i,is),i=1,3),
     &                       (vecM(i,is),i=1,3), sDotm(is)
      enddo

      write(FILE_E,'(''*-----   Latent Hardening Matrix -----*'')')
      do is = 1, numSlip
         write(FILE_E, '(8x,24F5.1)') (hardmtx(is,js), js=1,numSlip)
      enddo
c
c------- Set up Schmid tensors/vectors for each slip system
c
      call SetSlipSystemTensors(numslip, vecM, vecS, zBar0, pBar0,
     &                          qBar0, pBar0Vec, qBar0Vec, ppTBar0,
     &                          FILE_E)
c
c------- format
c
1000  format(/'*-----   Slip Systems for',i3,' (Crystal Type)-----*'/,
     &        7x,'  number slip syst = ',i4/
     &        7x,'  vtx stress file  = ',a15/
     &        7x,'  SS#',5x,'tau',18x,'vecS',30x,'vecM',20x,'SdotM')
2000  format( 7x, i4, 4x, f5.2, 4x, 3f10.5, 4x, 3f10.5, 4x, f10.5)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE stressCrystal(
     &   sigxx, sigyy, sigth, sigxy, epsxx, epsyy, epsth, epsxy, 
     &   epsum, spinn, qptpo, qptp, dtime, time, ielem, incr, numel, 
     &   numqpt, iprint, sigma, ddsdde, statusFlag, numIncrs, noel, npt,
     &   dfgrd1, coords)

      implicit none
     
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer ielem, incr, numel, numqpt, iprint, statusFlag, numIncrs
      integer noel, npt
      real*8  dtime, time
      real*8  sigxx(numqpt), sigyy(numqpt), sigth(numqpt), sigxy(numqpt)
      real*8  epsxx(numqpt), epsyy(numqpt), epsth(numqpt)
      real*8  epsxy(numqpt), epsum(numqpt), spinn(numqpt)
      real*8  qptpo(numqpt), qptp(numqpt)
      real*8  sigma(DIMS, DIMS, NUMQPT_T)
      real*8  ddsdde(DIMS2, DIMS2, NUMQPT_T), dfgrd1(3,3)
    
      integer iqpt
      real*8  thetao, theta, d_kk
      real*8  d_vec(NVECS), w_vec(3)
      real*8  coords(3)

      real*8  epsxz(NUMQPT_T),  epsyz(NUMQPT_T)
      real*8  spinxz(NUMQPT_T), spinyz(NUMQPT_T)
      real*8  sigxz(NUMQPT_T),  sigyz(NUMQPT_T)
      common /for3DProbs/ epsxz, epsyz, spinxz, spinyz, sigxz, sigyz
c
c---------------------------------------------------------------------72
c
c------ initialize vectors for symm/skew parts of velocity gradients
c
      call SetTensor(d_vec, pzero, NVECS)
      call SetTensor(w_vec, pzero, 3)
c
c------ loop over integration points
c
      do iqpt = 1, numqpt
c
c---------- recover symm/skew parts of velocity gradient at ielem/iqpt
c---------- eps_ij are deviatoric quantities
c
         d_vec(1) = (epsxx(iqpt) - epsyy(iqpt)) / sqr2 
         d_vec(2) = epsth(iqpt) * sqr32
         d_vec(3) = epsxy(iqpt) * sqr2
         if (NSHR .gt. 1) then
            d_vec(4) = epsxz(iqpt) * sqr2
            d_vec(5) = epsyz(iqpt) * sqr2
         endif
         d_kk = epsum(iqpt)

         w_vec(1) = -spinn(iqpt)
         if (NSHR .gt. 1) then
            w_vec(2) = -spinxz(iqpt)
            w_vec(3) = -spinyz(iqpt)
         endif
c
c---------- recover temperature at ielem/iqpt
c
         thetao = qptpo(iqpt)
         theta  = qptp(iqpt)
c
c---------- evolve state at ielem/iqpt
c
         call DriverCrystalEvolve(sigxx(iqpt), sigyy(iqpt), sigth(iqpt),
     &      sigxy(iqpt), sigxz(iqpt), sigyz(iqpt), d_vec, w_vec, d_kk,
     &      thetao, theta, dtime, time, iqpt, ielem, incr, numqpt, 
     &      numel, iprint, sigma(1, 1, iqpt), ddsdde(1, 1, iqpt),
     &      statusFlag, numIncrs, noel, npt, dfgrd1, coords)
         if (statusFlag .eq. kGLOBAL_FAILED) return

      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE DriverCrystalEvolve(
     &   sigxx, sigyy, sigth, sigxy, sigxz, sigyz, d_vec, w_vec, d_kk,
     &   thetao, theta, dtime, time, iqpt, ielem, incr, numqpt, numel, 
     &   iprint, savg_ij, cavg_ijkl, statusFlag, numIncrs, noel, npt,
     &   dfgrd1, coords)

      implicit  none
      include  'params_xtal.inc'
      include  'numbers.inc'
c
c---- arguments
c
      integer iqpt, ielem, incr, numqpt, numel, iprint, statusFlag
      integer numIncrs, noel, npt
      real*8  sigxx, sigyy, sigth, sigxy, sigxz, sigyz
      real*8  thetao, theta, d_kk, dtime, time
      real*8  d_vec(NVECS), w_vec(DIMS), dfgrd1(3,3)
      real*8  cavg_ijkl(DIMS2, DIMS2), savg_ij(DIMS, DIMS)
c
c---- variables passed through common blocks
c
      integer numgrn, numslip, numvtx, kODFout
      integer maxIterState, maxIterNewt
      real*8  tolerState, tolerNewt
      real*8  fCeDevVol(NVECS), fCeVol
      real*8  coords(3)

      real*8  tauSlip(MAX_SLIP)
      real*8  zBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0(DIMS, DIMS, MAX_SLIP), qBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0Vec(NVECS, MAX_SLIP), qBar0Vec(DIMS, MAX_SLIP)
      real*8  ppTBar0(NVECS, NVECS, MAX_SLIP)
      real*8  sigfs(NVECS, MAX_VTX)
      real*8  fCeDev(NVECS, NVECS), fCeiDev(NVECS, NVECS)
      real*8  matProp(NPROPS, MAX_SLIP)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)

      real*8  gstress  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev  (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstress_n(NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran_n(NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa_n (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev_n(NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geqvalues(NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot_n  (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot_n  (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot0   (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps     (9, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps_n   (9, MAX_GRN, NUMQPT_T, NUMEL_T)
c
c---- local variables
c
      integer iterCounterS, iterCounterN, ierr, igrn, nxtals_nc 
      real*8  tmp, wpnorm_avg
      real*8  s_ij(DIMS,DIMS), ee_ij(DIMS,DIMS), tau(NVECS)
      real*8  stress(NVECS), estran(NVECS), kappa(NKAPP), statev(NSTAV)
      real*8  stress_n(NVECS), estran_n(NVECS), kappa_n(NKAPP)
      real*8  statev_n(NSTAV), eqvalues(NEQVA), gamdot(MAX_SLIP)
      real*8  crot(DIMS,DIMS), rrot(DIMS,DIMS), drot(DIMS,DIMS)
      real*8  crot_n(DIMS,DIMS), rrot_n(DIMS,DIMS), crot0(DIMS,DIMS)
      real*8  c_ijkl(DIMS2,DIMS2), sdev_ij(DIMS,DIMS)
      real*8  rss(MAX_SLIP)
      real*8  eps(9), eps_n(9)
      real*8  gndden(MAX_SLIP), gndden_n(MAX_SLIP)
c
c---- common blocks
c
      common /XtalNumGrn/ numgrn, numslip
      common /XtalNumVtx/ numvtx
      common /XtalODFOut/ kODFout
      common /XtalStrVtx/ sigfs
      common /XtalProps / fCeDev, fCeiDev, matProp,
     &                    fCeDevVol, fCeVol, hardmtx
      common /XtalCrot0 / gcrot0
      common /XtalSlipG / tauSlip, zBar0, pBar0, qBar0,
     &                    pBar0Vec, qBar0Vec, ppTBar0
      common /XtalVars  / gstress, gestran, gkappa, gstatev, geqvalues,
     &                    gcrot, grrot, ggamdot, geps
      common /XtalVars_n/ gstress_n, gestran_n, gkappa_n, gstatev_n,
     &                    gcrot_n, grrot_n, geps_n

      common /nlcsolve1/  maxIterState, maxIterNewt
      common /nlcsolve2/  tolerState, tolerNewt
      
      integer numor, seed
      common  /XtalOrients/ numor, seed

      real*8 time_end
      common /timeEnd/ time_end
      
      real*8  gshear(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamtar(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gswitch(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      common  /XtalStep/ gshear, ggamtar, gswitch
      
      real*8  shear(MAX_SLIP), gamtar(MAX_SLIP), switch(MAX_SLIP)
      
      real*8  tauc(NKAPP, MAX_GRN), gammat(MAX_SLIP, MAX_GRN)
      common  /TauCrit/ tauc, gammat
      
      real*8  grainsize
      common  /XtalGrainSize/ grainsize
      
      integer activated, finished
      common  /GrainActive/ activated, finished
      data    activated, finished  /0, 0/
      
      real*8  ggndden  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden_n(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      common  /Dislocation/  ggndden, ggndden_n
      
      real*8  a(12)
      integer is
c
c---------------------------------------------------------------------72
c
c------- initialize average stress (sij) and consistent moduli (cijkl)
c
      call SetTensor(savg_ij, pzero, DIMS*DIMS)
      call SetTensor(cavg_ijkl, pzero, DIMS2*DIMS2)
c
c------- initialize average value of wp_norm
c
      wpnorm_avg = pzero
c
c------- counter for crystals that do not converge - Large Scale Appls.
      nxtals_nc = 0
c
c------- loop over all grains in aggregate at iqpt/ielem
c
      do igrn = 1, numgrn
c
c----------- zero local arrays
c
         call InitCrystalLocalArrays(tau, stress, estran, kappa, statev,
     &      stress_n, estran_n, kappa_n, statev_n, eqvalues, gamdot,
     &      rss, crot, rrot, drot, crot_n, rrot_n, crot0, c_ijkl, shear,
     &      gamtar, switch, eps, eps_n, gndden, gndden_n)
c
c---------- fetch variables at ielem/iqpt/igrn
c
         call FetchCrystalVariablesAtIP(gstress_n, gestran_n, gkappa_n,
     &      gstatev_n, geqvalues, gcrot_n, grrot_n, gcrot0, stress_n, 
     &      estran_n, kappa_n, statev_n, eqvalues, crot_n, rrot_n, 
     &      crot0, igrn, iqpt, ielem, gshear, ggamtar, gswitch, shear,
     &      gamtar, switch, geps_n, eps_n, ggndden_n, gndden_n)
c
c---------- monitor integration/convergence of XTAL constitutive eqns
c
         call MonitorCrystalIntegration(d_vec, w_vec, d_kk, s_ij, ee_ij,
     &      tau, 
     &      stress, estran, kappa, statev, eqvalues, gamdot, rss,
     &      crot, rrot,
     &      drot, stress_n, estran_n, kappa_n, statev_n, crot_n, rrot_n,
     &      crot0, fCeDev, fCeiDev, matProp, sigfs, tauSlip, zBar0, 
     &      pBar0, qBar0, pBar0Vec, qBar0Vec, ppTBar0, numslip, numvtx, 
     &      maxIterState, maxIterNewt, tolerState, tolerNewt, dtime, 
     &      fCeDevVol, fCeVol,
     &      theta, thetao, igrn, npt, noel, incr, iterCounterS, 
     &      iterCounterN, ierr, c_ijkl, hardmtx, switch, eps, eps_n,
     &      dfgrd1)
c
c---------- check convergence; if not then:
c              either reset crystal quantities or force time reduction
c
         if (ierr .ne. XTAL_CONVERGED) then
            call WriteMessage (XTAL_O, 
     &                'DriverXtalEvolve: Sub-Incrs failed!')
            call WriteMessage (XTAL_O, 
     &                'DriverXtalEvolve: nxtals_nc incremented by one')
            nxtals_nc = nxtals_nc + 1 
            write(XTAL_O, 1000) nxtals_nc, incr, igrn, npt, noel

            if (nxtals_nc .ge. NXTALS_NC_LIMIT) then      
              call WriteMessage (XTAL_O, 
     &                'DriverXtalEvolve: nxtals_nc > LIMIT')
              call WriteMessage(XTAL_O,
     &                '... will force dtime reduction by')
              call WriteMessage(XTAL_O,
     &                '... setting: status=kGLOBAL_FAILED')
              statusFlag = kGLOBAL_FAILED
              return
            endif

            call WriteMessage(XTAL_O, 'Resetting xtal quantities')
            call ResetCrystalQnts(stress, estran, kappa, statev, 
     &             eqvalues, gamdot, crot, rrot, s_ij, c_ijkl,
     &             stress_n, 
     &             estran_n, kappa_n, statev_n, crot_n, rrot_n, 
     &             matProp, pBar0Vec, fCeiDev, fCeDevVol, fCeVol, 
     &             numslip, dtime, switch, eps, eps_n, gndden, gndden_n)

         endif
c
c---------- calculate accumulated shear strain
c
         if (ierr .eq. XTAL_CONVERGED) then
c          if (noel .eq. 1000) then
              call UpdateFPandCoords()
c          endif    

      endif
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c
c---------- save computed variables at ielem/iqpt/igrn
c
         call SaveCrystalVariablesAtIP(stress, estran, kappa, statev, 
     &      gamdot, eqvalues, crot, rrot, gstress, gestran, gkappa, 
     &      gstatev, geqvalues, ggamdot, gcrot, grrot, igrn, iqpt, 
     &      ielem, shear, gamtar, switch, gshear, ggamtar, gswitch,
     &      eps, geps, gndden, ggndden)
c
c---------- output computed quantities at selected "noel,npt,igrn"
c
c         if (iprint .eq. kPRINT_MODEL .and. 
c     &        noel  .eq. kPRINT_ELEM  .and.
c     &         npt  .eq. kPRINT_QPT   .and.
c     &         igrn .eq. kPRINT_GRN        ) then

c            call OutputQnts(d_vec, s_ij, ee_ij, kappa, gamdot, rss,
c     &             statev, eqvalues, crot, rrot, c_ijkl, d_kk, dtime, 
c     &             time, numslip, incr, iterCounterS, iterCounterN, 
c     &             noel, npt, igrn, shear, gamtar, switch)

c         endif
c
c---------- average values of sij and cijkl for aggregate

         tmp = 1.d0 / numgrn
         call AddScaledTensor(tmp, s_ij, savg_ij, DIMS*DIMS)
         call AddScaledTensor(tmp, c_ijkl, cavg_ijkl, DIMS2*DIMS2)
c
c---------- average value of wp_norm
c
         wpnorm_avg = wpnorm_avg + tmp*statev(kWPNORM)

      enddo
c
c------- move aggregate deviatoric stress to main-program arrays
c
      call DeviatoricTensor(savg_ij, sdev_ij, DIMS)
      sigxx = sdev_ij(1,1)
      sigyy = sdev_ij(2,2)
      sigth = sdev_ij(3,3)
      sigxy = sdev_ij(1,2)
      if (NSHR .gt. 1) then
         sigxz = sdev_ij(1,3)
         sigyz = sdev_ij(2,3)
      endif
c
c------- output aggregate qnts at selected "ielem,iqpt"
c-------  note that ielem=1,iqpt=1 always !!
c
c      if (iprint .eq. kPRINT_MODEL .and. 
c     &     ielem .eq. kPRINT_ELEM  .and.
c     &     iqpt  .eq. kPRINT_QPT        ) then
c
c---------- write effective stress-strain curve
c
c         call WriteStressStrainCurve(sdev_ij, d_vec,
c     &                               wpnorm_avg, dtime, time,
c     &                               numgrn, incr, npt, noel)
c
c---------- write texture at specified increments
c
         if ( (incr/kODFout*kODFout) .eq. incr .or. 
     &                       incr .eq. numIncrs) then

c     &           dabs((time+dtime)-time_end) .lt. 1e-6) then
            call WriteTexture(gcrot, gkappa, d_vec, numgrn, incr, 
     &                        iqpt, ielem, npt, noel)

         endif

c      endif

1000  format(3x, 'Grains did not converged: ', i3, ':', 3x,
     &           'incr  #', i8, ';', 1x,
     &           'grain #', i5, ';', 1x,
     &           'iqpt  #', i3, ';', 1x,
     &           'elem  #', i5)

      return

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE InitCrystalLocalArrays(
     &   tau, stress, estran, kappa, statev, stress_n, estran_n, 
     &   kappa_n, statev_n, eqvalues, gamdot, rss, crot, rrot, drot, 
     &   crot_n, rrot_n, crot0, c_ijkl, shear, gamtar, switch,
     &   eps, eps_n, gndden, gndden_n)

      implicit  none
      include  'params_xtal.inc'
      include  'numbers.inc'

      real*8  tau(NVECS)
      real*8  stress(NVECS), estran(NVECS), kappa(NKAPP), statev(NSTAV)
      real*8  stress_n(NVECS), estran_n(NVECS), kappa_n(NKAPP)
      real*8  statev_n(NSTAV), eqvalues(NEQVA), gamdot(MAX_SLIP)
      real*8  crot(DIMS,DIMS), rrot(DIMS,DIMS), drot(DIMS,DIMS)
      real*8  crot_n(DIMS,DIMS), rrot_n(DIMS,DIMS), crot0(DIMS,DIMS)
      real*8  c_ijkl(DIMS2, DIMS2)
      real*8  rss(MAX_SLIP)
      real*8  shear(MAX_SLIP), gamtar(MAX_SLIP), switch(MAX_SLIP)
      real*8  eps(9), eps_n(9)
      real*8  gndden(MAX_SLIP), gndden_n(MAX_SLIP)
c
c---------------------------------------------------------------------72
c
      call SetTensor(tau,      pzero, NVECS)
      call SetTensor(stress,   pzero, NVECS)
      call SetTensor(estran,   pzero, NVECS)
      call SetTensor(kappa,    pzero, NKAPP)
      call SetTensor(statev,   pzero, NSTAV)
      call SetTensor(stress_n, pzero, NVECS)
      call SetTensor(estran_n, pzero, NVECS)
      call SetTensor(kappa_n,  pzero, NKAPP)
      call SetTensor(statev_n, pzero, NSTAV)
      call SetTensor(eqvalues, pzero, NEQVA)
      call SetTensor(gamdot,   pzero, MAX_SLIP)
      call SetTensor(crot,     pzero, DIMS*DIMS)
      call SetTensor(rrot,     pzero, DIMS*DIMS)
      call SetTensor(drot,     pzero, DIMS*DIMS)
      call SetTensor(crot_n,   pzero, DIMS*DIMS)
      call SetTensor(rrot_n,   pzero, DIMS*DIMS)
      call SetTensor(crot0,    pzero, DIMS*DIMS)
      call SetTensor(c_ijkl,   pzero, DIMS2*DIMS2)
      call SetTensor(shear,    pzero, MAX_SLIP)
      call SetTensor(gamtar,   pzero, MAX_SLIP)
      call SetTensor(switch,   pzero, MAX_SLIP)
      call SetTensor(rss,      pzero, MAX_SLIP)
      call SetTensor(eps,      pzero, 9)
      call SetTensor(eps_n,    pzero, 9)
      call SetTensor(gndden,   pzero, MAX_SLIP)
      call SetTensor(gndden_n, pzero, MAX_SLIP)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE FetchCrystalVariablesAtIP(
     &   gstress_n, gestran_n, gkappa_n, gstatev_n, geqvalues, 
     &   gcrot_n, grrot_n, gcrot0, stress_n, estran_n, kappa_n, 
     &   statev_n, eqvalues, crot_n, rrot_n, crot0, igrn, iqpt, 
     &   ielem, gshear, ggamtar, gswitch, shear, gamtar, switch,
     &   geps_n, eps_n, ggndden_n, gndden_n)

      implicit  none
      include  'params_xtal.inc'

      integer igrn, iqpt, ielem

      real*8  gstress_n (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran_n (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa_n  (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev_n (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geqvalues (NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot_n   (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot_n   (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot0    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps_n    (9, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden_n (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)

      real*8  stress_n(NVECS), estran_n(NVECS), kappa_n(NKAPP)
      real*8  statev_n(NSTAV), eqvalues(NEQVA)
      real*8  crot_n(DIMS,DIMS), rrot_n(DIMS,DIMS), crot0(DIMS,DIMS)
      real*8  eps_n(9), gndden_n(MAX_SLIP)
      
      real*8  gshear(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamtar(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gswitch(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      
      real*8  shear(MAX_SLIP), gamtar(MAX_SLIP), switch(MAX_SLIP)
c
c---------------------------------------------------------------------72
c
      call EqualTensors(gstress_n(1,igrn,iqpt,ielem), stress_n, NVECS)
      call EqualTensors(gestran_n(1,igrn,iqpt,ielem), estran_n, NVECS)
      call EqualTensors(gkappa_n (1,igrn,iqpt,ielem), kappa_n,  NKAPP)
      call EqualTensors(gstatev_n(1,igrn,iqpt,ielem), statev_n, NSTAV)
      call EqualTensors(geps_n   (1,igrn,iqpt,ielem), eps_n,    9)

      call EqualTensors(gcrot_n(1,1,igrn,iqpt,ielem), crot_n, DIMS*DIMS)
      call EqualTensors(grrot_n(1,1,igrn,iqpt,ielem), rrot_n, DIMS*DIMS)
      call EqualTensors(gcrot0(1,1,igrn,iqpt,ielem),  crot0,  DIMS*DIMS)
      
      call EqualTensors(ggndden_n(1,igrn,iqpt,ielem),gndden_n,MAX_SLIP)

      eqvalues(kEQP_n)    = geqvalues(kEQP_n,   igrn, iqpt, ielem)
      eqvalues(kMISES_n)  = geqvalues(kMISES_n, igrn, iqpt, ielem)
      eqvalues(kSHRATE_n) = geqvalues(kSHRATE_n,igrn, iqpt, ielem)
      eqvalues(kGAMTOT_n) = geqvalues(kGAMTOT_n,igrn, iqpt, ielem)
      
      call EqualTensors(gshear(1,igrn,iqpt,ielem), shear, MAX_SLIP)
      call EqualTensors(ggamtar(1,igrn,iqpt,ielem), gamtar, MAX_SLIP)
      call EqualTensors(gswitch(1,igrn,iqpt,ielem), switch, MAX_SLIP)
      
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SaveCrystalVariablesAtIP(
     &   stress, estran, kappa, statev, gamdot, eqvalues, crot, rrot, 
     &   gstress, gestran, gkappa, gstatev, geqvalues, ggamdot, gcrot, 
     &   grrot, igrn, iqpt, ielem, shear, gamtar, switch, gshear,
     &   ggamtar, gswitch, eps, geps, gndden, ggndden)

      implicit  none
      include  'params_xtal.inc'

      integer igrn, iqpt, ielem

      real*8  gstress  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev  (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geqvalues(NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps     (9, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)

      real*8  stress(NVECS), estran(NVECS), kappa(NKAPP)
      real*8  statev(NSTAV), eqvalues(NEQVA), gamdot(MAX_SLIP)
      real*8  crot(DIMS,DIMS), rrot(DIMS,DIMS)
      real*8  eps(9), gndden(MAX_SLIP)
      
      real*8  gshear(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamtar(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gswitch(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  shear(MAX_SLIP), gamtar(MAX_SLIP), switch(MAX_SLIP)
c
c---------------------------------------------------------------------72
c
      call EqualTensors(stress, gstress(1,igrn,iqpt,ielem), NVECS)
      call EqualTensors(estran, gestran(1,igrn,iqpt,ielem), NVECS)
      call EqualTensors(kappa,  gkappa (1,igrn,iqpt,ielem), NKAPP)
      call EqualTensors(statev, gstatev(1,igrn,iqpt,ielem), NSTAV)
      call EqualTensors(gamdot, ggamdot(1,igrn,iqpt,ielem), MAX_SLIP)
      call EqualTensors(crot, gcrot(1,1,igrn,iqpt,ielem), DIMS*DIMS)
      call EqualTensors(rrot, grrot(1,1,igrn,iqpt,ielem), DIMS*DIMS)
      call EqualTensors(eps, geps(1,igrn,iqpt,ielem), 9)
      call EqualTensors(gndden, ggndden(1,igrn,iqpt,ielem), MAX_SLIP)

      geqvalues(kEQP,   igrn,iqpt,ielem) = eqvalues(kEQP)
      geqvalues(kMISES, igrn,iqpt,ielem) = eqvalues(kMISES)
      geqvalues(kSHRATE,igrn,iqpt,ielem) = eqvalues(kSHRATE)
      geqvalues(kGAMTOT,igrn,iqpt,ielem) = eqvalues(kGAMTOT)
      
      call EqualTensors(shear, gshear(1,igrn,iqpt,ielem), MAX_SLIP)
      call EqualTensors(gamtar, ggamtar(1,igrn,iqpt,ielem), MAX_SLIP)
      call EqualTensors(switch, gswitch(1,igrn,iqpt,ielem), MAX_SLIP)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE UpdateCrystalHistory(
     &    numqpt, numel
     &   )

      implicit none
      include 'params_xtal.inc'
c
c---- arguments   
c
      integer numqpt, numel
c
c---- variables passed through common blocks
c
      integer numgrn, numslip
      real*8  gstress  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev  (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstress_n(NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran_n(NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa_n (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev_n(NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geqvalues(NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot_n  (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot_n  (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
c
c---- local variables
c
      integer ig, ip, ie
c
c---- common block
c
      common /XtalNumGrn/ numgrn, numslip
      common /XtalVars  / gstress, gestran, gkappa, gstatev, geqvalues,
     &                    gcrot, grrot, ggamdot
      common /XtalVars_n/ gstress_n, gestran_n, gkappa_n, gstatev_n,
     &                    gcrot_n, grrot_n
c
c---------------------------------------------------------------------72
c
c------- update variables (gamdot is not updated)
c
      do ie = 1, numel
        do ip = 1, numqpt
          do ig = 1, numgrn

            call EqualTensors(gstress  (1,ig,ip,ie), 
     &                        gstress_n(1,ig,ip,ie), NVECS)
            call EqualTensors(gestran  (1,ig,ip,ie), 
     &                        gestran_n(1,ig,ip,ie), NVECS)
            call EqualTensors(gkappa   (1,ig,ip,ie),  
     &                        gkappa_n (1,ig,ip,ie), NKAPP)
            call EqualTensors(gstatev  (1,ig,ip,ie), 
     &                        gstatev_n(1,ig,ip,ie), NSTAV)
            call EqualTensors(gcrot    (1,1,ig,ip,ie), 
     &                        gcrot_n  (1,1,ig,ip,ie), DIMS*DIMS)
            call EqualTensors(grrot    (1,1,ig,ip,ie), 
     &                        grrot_n  (1,1,ig,ip,ie), DIMS*DIMS)

            geqvalues(kEQP_n,   ig,ip,ie) = geqvalues(kEQP,   ig,ip,ie)
            geqvalues(kMISES_n, ig,ip,ie) = geqvalues(kMISES, ig,ip,ie)
            geqvalues(kSHRATE_n,ig,ip,ie) = geqvalues(kSHRATE,ig,ip,ie)
            geqvalues(kGAMTOT_n,ig,ip,ie) = geqvalues(kGAMTOT,ig,ip,ie)

          enddo
        enddo
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE ResetCrystalHistory(
     &    numqpt, numel
     &   )

      implicit none
      include 'params_xtal.inc'
c
c---- arguments
c
      integer numqpt, numel
c
c---- variables passed through common blocks
c
      integer numgrn, numslip
      real*8  gstress  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev  (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstress_n(NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran_n(NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa_n (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev_n(NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geqvalues(NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot_n  (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot_n  (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
c
c---- local variables
c
      integer ig, ip, ie
c
c---- common block
c
      common /XtalNumGrn/ numgrn, numslip
      common /XtalVars  / gstress, gestran, gkappa, gstatev, geqvalues,
     &                    gcrot, grrot, ggamdot
      common /XtalVars_n/ gstress_n, gestran_n, gkappa_n, gstatev_n,
     &                    gcrot_n, grrot_n
c
c---------------------------------------------------------------------72
c
c------- reset variables (gamdot is not reset)
c
      do ie = 1, numel
        do ip = 1, numqpt
          do ig = 1, numgrn

            call EqualTensors(gstress_n(1,ig,ip,ie), 
     &                        gstress  (1,ig,ip,ie), NVECS)
            call EqualTensors(gestran_n(1,ig,ip,ie), 
     &                        gestran  (1,ig,ip,ie), NVECS)
            call EqualTensors(gkappa_n (1,ig,ip,ie),  
     &                        gkappa   (1,ig,ip,ie), NKAPP)
            call EqualTensors(gstatev_n(1,ig,ip,ie), 
     &                        gstatev  (1,ig,ip,ie), NSTAV)
            call EqualTensors(gcrot_n  (1,1,ig,ip,ie), 
     &                        gcrot    (1,1,ig,ip,ie), DIMS*DIMS)
            call EqualTensors(grrot_n  (1,1,ig,ip,ie), 
     &                        grrot    (1,1,ig,ip,ie), DIMS*DIMS)

            geqvalues(kEQP,   ig,ip,ie) = geqvalues(kEQP_n,   ig,ip,ie)
            geqvalues(kMISES, ig,ip,ie) = geqvalues(kMISES_n, ig,ip,ie)
            geqvalues(kSHRATE,ig,ip,ie) = geqvalues(kSHRATE_n,ig,ip,ie)
            geqvalues(kGAMTOT,ig,ip,ie) = geqvalues(kGAMTOT_n,ig,ip,ie)

          enddo
        enddo
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE OutputQnts(
     &   d_vec, s_ij, ee_ij, kappa, gamdot, rss, statev, eqvalues, crot,
     &   rrot, c_ijkl, d_kk, dtime, time, numslip, incr, 
     &   iterCounterS, iterCounterN, ielem, iqpt, igrn, shear, gamtar,
     &   switch)

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip, incr, iterCounterS, iterCounterN
      integer ielem, iqpt, igrn
      real*8  d_kk, dtime, time
      real*8  d_vec(NTENS), s_ij(DIMS,DIMS), ee_ij(DIMS,DIMS)
      real*8  kappa(NKAPP), gamdot(MAX_SLIP), statev(NSTAV)
      real*8  eqvalues(NEQVA), crot(DIMS,DIMS), rrot(DIMS,DIMS)
      real*8  c_ijkl(DIMS2,DIMS2)
      real*8  rss(MAX_SLIP)

      integer i, j, is
c      character*160 message
      real*8  strain_11, stress_11, d_eff, eqstran, deqstran
      real*8  d_ij(DIMS,DIMS)

      real*8  InnerProductVec

      save    strain_11, eqstran
      data    strain_11, eqstran /0.d0, 0.d0/

      real*8  overstress(MAX_SLIP)
      common /XtalOverS/ overstress
      
      real*8  shear(MAX_SLIP), gamtar(MAX_SLIP), switch(MAX_SLIP)
      
      real*8  strain_33, stress_33
      real*8  m(MAX_SLIP), maxx
	  
      integer dex
      save    dex

      save    strain_33
      data    strain_33 /0.d0/
c
c----- NOTE: 'strain_11' makes sense only for coaxial loading
c
c---------------------------------------------------------------------72
c
c------- write headers in output files
c
      if (incr .eq. 1) then
         write(XTAL_STRESS_O,  800) ielem, iqpt, igrn
         write(XTAL_STRAIN_O,  900) ielem, iqpt, igrn
         write(XTAL_EFFSS_O,  1000) ielem, iqpt, igrn
         write(XTAL_TRUESS_O, 2000) ielem, iqpt, igrn
         write(XTAL_ITER_O,   3000) ielem, iqpt, igrn
      endif
c
c------- find the index of the maximum Schmid factor
c
         if (incr .eq. 1) then
             do is = 1, numslip
                 m(is) = rss(is) / s_ij(3,3)
             enddo
      
             maxx = dabs(m(1));
             do is = 2, numslip
                 if (dabs(m(is)) .ge. maxx) then
                     maxx = dabs(m(is))
                     dex = is
					
                 endif
             enddo
			 
         endif 		  
c
c------- equivalent total strain and strain_11
c
      d_eff = dsqrt(2./3.*InnerProductVec(d_vec, d_vec, NVECS))
      deqstran = dtime * d_eff
      eqstran  = eqstran + deqstran

      call Vec5x1ToMat3x3Symm(d_vec, d_ij, DIMS)
      strain_11 = strain_11 + dtime*(d_ij(1,1) + pthird*d_kk)
      stress_11 = s_ij(1,1)
      
      strain_33 = strain_33 + dtime*(d_ij(3,3) + pthird*d_kk)
      stress_33 = s_ij(3,3)
c
c------- write stress / consistent tangent
c
      write(XTAL_STRESS_O, 4000) incr, ((s_ij(i,j), j=i,3), i=1,3)
      write(XTAL_STRESS_O, 4500) ((c_ijkl(i,j), j=1,DIMS2), i=1,DIMS2)
c      write(XTAL_STRESS_O, 4500) (rss(i)/(kappa(1)+overstress(i)),
c     &                            i=1,numslip)
      write(XTAL_STRESS_O, 4500) (rss(i)/(kappa(i)), i=1,numslip)
c
c------- write elastic strains / crot / rotation tensors; gammadot
c
      write(XTAL_STRAIN_O, '(/i5)') incr
      do i = 1, 3
         write(XTAL_STRAIN_O, 5000) (ee_ij(i,j), j=1,3), 
     &                       (crot(i,j),j=1,3), (rrot(i,j), j=1,3)
      enddo

      write(XTAL_STRAIN_O, 5200) eqvalues(kEQP), statev(kWPNORM),
     &                           int( statev(kSSACT) )
      write(XTAL_STRAIN_O, 5500) (gamdot(i), i=1,numslip)
c
c------- write effective quantities
c
      write(XTAL_EFFSS_O, 5800) incr, int(statev(kSSACT)), dtime,
     &                time+dtime, d_eff, deqstran,
     &                statev(kWPNORM), eqvalues(kEQP),
     &                (eqvalues(kEQP)-eqvalues(kEQP_n))/dtime,
     &                eqvalues(kMISES)
c
c------- write true stress-strain curve (X-direction) (uniaxial/MPS)
c
      write(XTAL_TRUESS_O, 6000) incr, time, strain_33, shear(dex), 
     &                  abs(rss(dex)), kappa(dex), stress_33,
     &                  gamtar(dex), int(switch(dex)), dex
c
c------- write iteration counters and escalar variables
c
      write(XTAL_ITER_O, 7000) incr, iterCounterS, iterCounterN, 
     &                  statev(kEVOL), statev(kPRES), statev(kDETVe),
     &                  (kappa(i),i=1,numslip)
c     &                  (kappa(1)+overstress(i),i=1,numslip)
c
c------- formats
c
 800  format(' INCR', 15x, 'CAUCHY STRESS' ,/,
     &       20x,'CONSISTENT TANGENT ',/, 20x, 'RSS/KAPPA ',/,
     &       ' (elem # ', i5, ',  qpt # ', i2, ',  grn # ', i5, ')'/) 
 900  format(' INCR',/,18x,' EE_ij',32x,'Crot_ij',30x,'Rrot_ij',/, 
     &       19x, 'EQPS', 33x, 'WP_NORM', 30x, 'SS_ACTIV',/,
     &       48x, ' GAMMADOT(1,...,NUMSLIP)',/,38x,
     &       ' (elem # ', i5, ',  qpt # ', i2, ',  grn # ', i5, ')') 
1000  format(' INCR',3x,'SS_ACT',5x,'DTIME',8x,'TIME',8x,'D_EFF',10x,
     &      'DEQSTRAN',8x,'WP_NORM',7x,'EQP',10x,'DP_EFF',8x,'MISES',/,
     &       ' (elem # ', i5, ',  qpt # ', i2, ',  grn # ', i5, ')')
2000  format(' INCR',5x,'TIME',10x,'PE_33',5x,'SHEAR STRAIN',6x,
     &       'RSS',11x,'CRSS',10x,'S_33',9x,'Gamtar',4x,'switch',4x,
     &       'index'/,
     &       ' (elem # ', i5, ',  qpt # ', i2, ',  grn # ', i5, ')',/,
     &       ' (for uniaxial loading along Z-direction / MPS runs)' ) 
3000  format(' INCR',1x,'IT-S',1x,'IT-N',6x,'EVOL', 8x,'PRESS', 
     &       8x,'DETVe',6x,'KAPPA(1) ... KAPPA(NKAPP)'/, 
     &       ' (elem # ', i5, ',  qpt # ', i2, ',  grn # ', i5, ')') 
4000  format(i8,2x,3(1x,e11.4)/22x,2(1x,e11.4)/34x,1(1x,e11.4)/)
4500  format(6(10x,6(1x,e11.4)/))
5000  format((3(2x,3(1x,e11.4))))
5200  format(/(15x,e11.4,27x,e11.4,28x,i5))
5500  format(/(22x,6(1x,e11.4)))
5800  format(i8,2x,i3,8(2x,e12.5))
6000  format(i4,7(2x,e12.5),2x,i3,7x,i3)
7000  format(i8,1x,i3,2x,i3,1x,3(2x,e11.4),18(1x,e11.4))

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MonitorCrystalIntegration( 
     &   d_vec, w_vec, d_kk, s_ij, ee_ij, tau, stress, estran, kappa, 
     &   statev, eqvalues, gamdot, rss, crot, rrot, drot, stress_n, 
     &   estran_n,
     &   kappa_n, statev_n, crot_n, rrot_n, crot0, fCeDev, fCeiDev, 
     &   matProp, sigfs, tauSlip, zBar0, pBar0, qBar0, pBar0Vec,
     &   qBar0Vec, ppTBar0, numslip, numvtx, maxIterState, maxIterNewt, 
     &   tolerState, tolerNewt, dtime, fCeDevVol, fCeVol, 
     &   theta, thetao, igrn, iqpt, 
     &   ielem, incr, iterCounterS, iterCounterN, ierr, cepmod, hardmtx,
     &   switch, eps, eps_n, dfgrd1)

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip, numvtx, maxIterState, maxIterNewt
      integer igrn, iqpt, ielem, incr, iterCounterS, iterCounterN, ierr
      real*8  tolerState, tolerNewt, dtime, theta, thetao, d_kk
      real*8  fCeDevVol(NVECS), fCeVol
      real*8  d_vec(NVECS), w_vec(DIMS), dfgrd1(3,3)
      real*8  s_ij(DIMS,DIMS), ee_ij(DIMS,DIMS), tau(NVECS)
      real*8  stress(NVECS), estran(NVECS), kappa(NKAPP)
      real*8  statev(NSTAV), eqvalues(NEQVA), gamdot(MAX_SLIP)
      real*8  crot(DIMS,DIMS), rrot(DIMS,DIMS), drot(DIMS,DIMS)
      real*8  stress_n(NVECS), estran_n(NVECS), kappa_n(NKAPP)
      real*8  statev_n(NSTAV), crot_n(DIMS,DIMS), rrot_n(DIMS,DIMS)
      real*8  crot0(DIMS,DIMS)
      real*8  fCeDev(NVECS,NVECS), fCeiDev(NVECS,NVECS)
      real*8  matProp(NPROPS, MAX_SLIP)
      real*8  sigfs(NVECS,MAX_VTX), tauSlip(MAX_SLIP)
      real*8  zBar0(DIMS,DIMS,MAX_SLIP)
      real*8  pBar0(DIMS,DIMS,MAX_SLIP), qBar0(DIMS,DIMS,MAX_SLIP)
      real*8  pBar0Vec(NVECS,MAX_SLIP), qBar0Vec(DIMS,MAX_SLIP)
      real*8  ppTBar0(NVECS,NVECS,MAX_SLIP)
      real*8  cepmod(DIMS2,DIMS2)
      real*8  rss(MAX_SLIP)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)
      real*8  eps(9), eps_n(9)

      integer subIncr, totSubIncrs
      real*8  dtime_tr, theta_tr, tmp, d_kk_n, d_kk_tr
      real*8  d_vec_tr(NVECS), w_vec_tr(DIMS)
      real*8  d_vec_n(NVECS), w_vec_n(DIMS), tau_save(NSTAV)
      
      real*8  switch(MAX_SLIP)
c
c---------------------------------------------------------------------72
c
c------- counters for sub-increments
c
      subIncr = 1
      totSubIncrs = 1
c
c------- integrate crystal constitutive equations
c
      iterCounterS = 0
      call IntegrateCrystalEqns(d_vec, w_vec, d_kk, s_ij, ee_ij, tau, 
     &   stress, estran, kappa, statev, eqvalues, gamdot, rss,
     &   crot, rrot, 
     &   drot, stress_n, estran_n, kappa_n, statev_n, crot_n, rrot_n, 
     &   crot0, fCeDev, fCeiDev, matProp, sigfs, tauSlip, zBar0, 
     &   pBar0, qBar0, pBar0Vec, qBar0Vec, ppTBar0, numslip, numvtx, 
     &   maxIterState, maxIterNewt, tolerState, tolerNewt, dtime, 
     &   fCeDevVol, fCeVol,
     &   theta, thetao, iterCounterS, iterCounterN, incr, ierr, 
     &   cepmod, subIncr, totSubIncrs, hardmtx, switch, eps, eps_n, 
     &   dfgrd1, iqpt, ielem)
c
c------- if converged -> return, else do -> subincrementation
c------- NOTE: here, subincrementation improves the initial guess 
c-------       for the solution variables 'stress'
c
      if (ierr .eq. XTAL_CONVERGED) return

      call WriteMessage
     &     (XTAL_O, 'MonitorCrystalIntegration: using Sub-Incrs!')
c
c------- compute velocity gradient at t_n
c
      call DeformationRate_n(d_vec_n, w_vec_n, d_kk_n, iqpt)
c     MAY NEED Vec6 to Vec5 routine.........
c
c------- loop for sub-incrementation procedure
c
      do while (.true.)
c
c---------- if not converged, increase # of sub-increments 
         if (ierr .ne. XTAL_CONVERGED) then
            subIncr = 2 * subIncr - 1
            totSubIncrs = 2 * totSubIncrs
            if (totSubIncrs .gt. 128) then
               call WriteMessage(XTAL_O,
     &                  'MonitorCrystalIntegration: totSubIncrs > 128')
               return
            endif
            call SetTensor(tau, pzero, NVECS)
            if (subIncr .gt. 1) 
     &           call EqualTensors(tau_save, tau, NVECS)
c
c---------- if converged, adjust subincrements
         else if (subIncr .lt. totSubIncrs) then
            if ( (subIncr/2*2) .eq. subIncr) then
               subIncr = subIncr / 2 + 1
               totSubIncrs = totSubIncrs / 2
            else
               subIncr = subIncr + 1
            endif
c
c---------- successful return for subincrementation
         else
            call WriteMessage(XTAL_O, 'Sub-Incrs successful !!!')
            return
         endif
c
c---------- report sub-incrs status
         write(XTAL_O, 1000) igrn, iqpt, ielem, subIncr, totSubIncrs
c
c---------- trial time, velocity gradient and temperature (assumes 
c---------- linear variation of velocity gradient and temperature in dt)
         tmp = real(subIncr) / real(totSubIncrs)   
         dtime_tr = dtime * tmp
         theta_tr = (1.-tmp)*thetao + tmp*theta
         d_kk_tr  = (1.-tmp)*d_kk_n + tmp*d_kk

         call AddTensors((1.-tmp), d_vec_n, tmp, d_vec, d_vec_tr, NVECS)
         call AddTensors((1.-tmp), w_vec_n, tmp, w_vec, w_vec_tr, DIMS)
c
c---------- save current convergent solution before getting next solution
         if (subIncr .gt. 1 .and. ierr .eq. XTAL_CONVERGED)
     &                 call EqualTensors(tau, tau_save, NVECS)
c
c---------- integrate crystal constitutive equations
         iterCounterS = 0
         call IntegrateCrystalEqns(d_vec_tr, w_vec_tr, d_kk_tr, s_ij,
     &     ee_ij, tau, stress, estran, kappa, statev, eqvalues, gamdot, 
     &     rss,
     &     crot, rrot, drot, stress_n, estran_n, kappa_n, statev_n, 
     &     crot_n, rrot_n, crot0, fCeDev, fCeiDev, matProp, sigfs, 
     &     tauSlip, zBar0, pBar0, qBar0, pBar0Vec, qBar0Vec, ppTBar0, 
     &     numslip, numvtx, maxIterState, maxIterNewt, tolerNewt, 
     &     tolerState, dtime_tr, fCeDevVol, fCeVol,
     &     theta_tr, thetao, iterCounterS, iterCounterN, incr, 
     &     ierr, cepmod, subIncr, totSubIncrs, hardmtx, switch, eps,
     &     eps_n, dfgrd1, iqpt, ielem)

      enddo

1000  format(3x, 'at grain #', i5, ';', 1x, 
     &           'iqpt #', i3, ';', 1x,
     &           'elem #', i5, ';', 1x,
     &           'subInc/totSubIncrs = ', i3, '/', i3)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE WriteTexture(
     &   gcrot, gkappa, d_vec, numgrn, incr, iqpt, ielem, npt, noel
     &   )

      implicit none
      include  'params_xtal.inc'

      integer numgrn, incr, iqpt, ielem, npt, noel
      real*8  gcrot(DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa(NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  d_vec(NTENS)

      integer igrn, i
      real*8  pi, pi180, d_eff
      real*8  angle(DIMS)

      real*8  InnerProductVec
c
c---------------------------------------------------------------------72
c
      pi = 4.0 * datan(1.0d+00)
      pi180 = pi/180.
c
c-------  effective strain rate at the ip/elem
c
      d_eff = dsqrt(2./3.*InnerProductVec(d_vec, d_vec, NVECS))
c
c------- output heading
c
      if (noel .eq. 1 .and. npt .eq. 1) 
     &            write(XTAL_TXT_OUT, 1000) incr
c
c------- output orientation
c
      do igrn = 1, numgrn

         call RotMatrixToAngles(gcrot(1,1,igrn,iqpt,ielem), angle, DIMS)

         write(XTAL_TXT_OUT, 2000) (angle(i)/pi180, i=1,3), igrn,
     &                          npt, noel, gkappa(1,igrn,iqpt,ielem),
     &                          d_eff

      enddo
c
c------- formats
c
1000  format(/'*-----   Euler Angles at incr ', i8, ' -----*'/
     &        4x,'ang1',4x,'ang2',4x,'ang3',4x,' igrn',4x,'intpt',
     &        3x,'ielem ',4x,'kappa(1)',4x,'d_eff')
2000  format( 3f8.2, 3(3x,i5), 3x, f10.2, 6x, e12.2)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE WriteStressStrainCurve(
     &   sdev_ij, d_vec, wpnorm_avg, dtime, time, numgrn, 
     &   incr, npt, noel
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numgrn, incr, npt, noel
      real*8  dtime, time, wpnorm_avg
      real*8  sdev_ij(DIMS,DIMS), d_vec(NTENS)

      integer denom
      real*8  s_eff, d_eff, eqstran, deqstran
      real*8  sdev_vec(NVECS)
      real*8  d_eff_avg, s_eff_avg

      real*8  InnerProductVec

      save    eqstran
      data    eqstran /0.d0/

      save    d_eff_avg, s_eff_avg, deqstran

      integer numel_aba, numqpt_aba
      common /data_aba/ numel_aba, numqpt_aba
c
c---------------------------------------------------------------------72
c
c------- write headers in output file
c
      if (noel .eq. 1 .and. npt .eq. 1) then
         s_eff_avg = 0.0d0
         d_eff_avg = 0.0d0
         deqstran  = 0.0d0
         if (incr .eq. 1) 
     &        write(AGG_EFFSS_O, 800) numel_aba, numqpt_aba, numgrn
      endif
c
c------- effective stress
c
      call Mat3x3ToVec5x1Symm(sdev_ij, sdev_vec, DIMS)
      s_eff = dsqrt(3./2.*InnerProductVec(sdev_vec, sdev_vec, NVECS))

      s_eff_avg = s_eff_avg + s_eff
c
c-------  effective strain rate and total equivalent strain
c
      d_eff = dsqrt(2./3.*InnerProductVec(d_vec, d_vec, NVECS))
      d_eff_avg = d_eff_avg + d_eff

      deqstran = deqstran + dtime * d_eff
c
c------- write stress/strain curve
c
      if (noel .eq. numel_aba .and. npt .eq. numqpt_aba) then
c         denom = numgrn * numel_aba * numqpt_aba
         denom = numel_aba * numqpt_aba
         s_eff_avg = s_eff_avg / denom
         d_eff_avg = d_eff_avg / denom
         deqstran = deqstran / denom
         eqstran  = eqstran + deqstran
         write(AGG_EFFSS_O, 1000) incr, dtime, time+dtime, d_eff_avg, 
     &           deqstran, eqstran, wpnorm_avg, s_eff_avg
      endif
c
c------- formats
c
800   format(' INCR',9x,'DTIME',9x,'TIME',9x,'D_EFF',8x,'DEQSTRAN',
     &       8x,'EQSTRAN',8x,'WP_NORM',6x,'S_EFF'/,
     &       ' (nelem: ', i5, ',  nqpts: ', i2, ',  ngrns: ', i5, ')') 
1000  format(i8,10(2x,e12.5))

      return
      END
c
c=====================================================================72
c
      SUBROUTINE ResetCrystalQnts(
     &   stress, estran, kappa, statev, eqvalues, gamdot, crot, rrot, 
     &   s_ij, c_ijkl, stress_n, estran_n, kappa_n, statev_n, crot_n, 
     &   rrot_n, 
     &   matProp, pBar0Vec, fCeiDev, fCeDevVol, fCeVol, numslip, dtime,
     &   switch, eps, eps_n, gndden, gndden_n)

      implicit  none
      include  'params_xtal.inc'
      include  'numbers.inc'

      integer numslip
      real*8  dtime

      real*8  stress(NVECS), estran(NVECS), kappa(NKAPP)
      real*8  statev(NSTAV), eqvalues(NEQVA), gamdot(MAX_SLIP)
      real*8  crot(DIMS,DIMS), rrot(DIMS,DIMS)
      real*8  s_ij(DIMS,DIMS), c_ijkl(DIMS2,DIMS2)

      real*8  stress_n(NVECS), estran_n(NVECS), kappa_n(NKAPP)
      real*8  statev_n(NSTAV), crot_n(DIMS,DIMS), rrot_n(DIMS,DIMS)

      real*8  matProp(NPROPS, MAX_SLIP)
      real*8  pBar0Vec(NVECS,MAX_SLIP)
      real*8  fCeiDev(NVECS,NVECS), fCeDevVol(NVECS), fCeVol
      real*8  eps(9), eps_n(9)
      real*8  gndden(MAX_SLIP), gndden_n(MAX_SLIP)
c
c---- Local variables
c
      integer is
      real*8  qr5x5(NVECS,NVECS)
      real*8  pHatVec(NVECS,MAX_SLIP), ppTHat(NVECS,NVECS,MAX_SLIP)
      real*8  fCeiDevHat(NVECS,NVECS),fCeDevVolHat(NVECS)

      real*8  crss
      real*8  tau(NVECS)
      real*8  rss(MAX_SLIP), dGamdTau(MAX_SLIP)

      real*8  InnerProductVec, SSKineticEqn
c
c---- Common Blocks
c
      real*8  overstress(MAX_SLIP)
      common /XtalOverS/ overstress

      real*8  sGamPPt(NVECS,NVECS)
      common /SumGamPPt/ sGamPPt
      
      real*8  switch(MAX_SLIP)
c
c---------------------------------------------------------------------72
c
c---- Reset crystal quantities
c
      call EqualTensors(stress_n, stress, NVECS)
      call EqualTensors(estran_n, estran, NVECS)
      call EqualTensors(kappa_n , kappa , NKAPP)
      call EqualTensors(statev_n, statev, NSTAV)
      call EqualTensors(crot_n  , crot  , DIMS*DIMS)
      call EqualTensors(rrot_n  , rrot  , DIMS*DIMS)
      call EqualTensors(eps_n, eps, 9)
      call EqualTensors(gndden_n, gndden, MAX_SLIP)
                                                                    
      eqvalues(kEQP)    = eqvalues(kEQP_n)
      eqvalues(kMISES)  = eqvalues(kMISES_n)             
      eqvalues(kSHRATE) = eqvalues(kSHRATE_n)
      eqvalues(kGAMTOT) = eqvalues(kGAMTOT_n)

      call SetToScaledtensor(statev(kDETVe), stress, tau, NVECS)
c
c------- Cauchy stress (tensor form)
      call Vec5x1ToMat3x3Symm(stress, s_ij, DIMS)
      call AddScaledTensor(statev(kPRES), Ident2nd, s_ij, DIMS*DIMS)
c
c---- Compute Approximate Consistent Tangent
c
c------- Rotation tensor 
      call RotMat5x5ForSymm(crot, qr5x5, DIMS)
c
c------- Compute gamdot and sGamPPt
      call SetTensor(sGamPPt, pzero, NVECS*NVECS)
      do is = 1, numslip
c
c------- Rotate symmetric part of Schimdt tensor
         call MultAxu(qr5x5, pBar0Vec(1,is), pHatVec(1,is), NVECS)
         call OuterProductVec(pHatVec(1,is), pHatVec(1,is),
     &                        ppTHat(1,1,is), NVECS)
c
c------- Resolve shear stresses
         rss(is) = InnerProductVec(tau, pHatVec(1,is), NVECS)
c         crss = kappa(1) + overstress(is)
         crss = kappa(is)
c
c------- Shear strain rate
         gamdot(is) = SSKineticEqn(rss(is),crss,matProp(1,is),kGAMDOT,
     &                             switch(is))
c
c------- Derivative of shear strain rate
         dGamdTau(is) = SSKineticEqn(rss(is), crss, matProp(1,is),
     &                                   kdGAMdTAU, switch(is))
         call AddScaledTensor(dtime*dGamdTau(is), ppTHat(1,1,is),
     &                        sGamPPt, NVECS*NVECS)

      enddo
c
c------- Rotate elasticities
      call MultQAQT(qr5x5, fCeiDev, fCeiDevHat, NVECS)
      call MultAxu(qr5x5, fCeDevVol, fCeDevVolHat, NVECS)
c
c------- Approximate consistent tangent
      call PlasticModuli(c_ijkl, fCeiDevHat, fCeDevVolHat, fCeVol,
     &                   statev(kDETVe))

      return
      end
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE IntegrateCrystalEqns(
     &   d_vec, w_vec, d_kk, s_ij, ee_ij, tau, stress,
     &   estran, kappa, statev, eqvalues, gamdot, rss, crot, rrot, drot,
     &   stress_n, estran_n, kappa_n, statev_n, crot_n, rrot_n, 
     &   crot0, fCeDev, fCeiDev, matProp, sigfs, tauSlip, zBar0, 
     &   pBar0, qBar0, pBar0Vec, qBar0Vec, ppTBar0, numslip, numvtx, 
     &   maxIterState, maxIterNewt, tolerState, tolerNewt, dtime, 
     &   fCeDevVol, fCeVol,
     &   theta, thetao, iterCounterS, iterCounterN, globalIncr, ierr, 
     &   cepmod, subIncr, totSubIncrs, hardmtx, switch, eps, eps_n,
     &   dfgrd1, iqpt, ielem)

      implicit none
      include 'params_xtal.inc' 
      include 'numbers.inc'

      integer numslip, numvtx, maxIterState, maxIterNewt, globalIncr
      integer iterCounterS, iterCounterN, ierr, subIncr, totSubIncrs
      integer iqpt, ielem
      real*8  tolerState, tolerNewt, dtime, theta, thetao, d_kk
      real*8  fCeDevVol(NVECS), fCeVol
      real*8  d_vec(NVECS), w_vec(DIMS)
      real*8  s_ij(DIMS,DIMS), ee_ij(DIMS,DIMS), tau(NVECS)
      real*8  stress(NVECS), estran(NVECS), kappa(NKAPP)
      real*8  statev(NSTAV), eqvalues(NEQVA), gamdot(MAX_SLIP)
      real*8  crot(DIMS,DIMS), rrot(DIMS,DIMS), drot(DIMS,DIMS)
      real*8  stress_n(NVECS), estran_n(NVECS), kappa_n(NKAPP)
      real*8  statev_n(NSTAV), crot_n(DIMS,DIMS), rrot_n(DIMS,DIMS)
      real*8  crot0(DIMS,DIMS)
      real*8  fCeDev(NVECS,NVECS), fCeiDev(NVECS,NVECS)
      real*8  matProp(NPROPS,MAX_SLIP)
      real*8  sigfs(NVECS,MAX_VTX), tauSlip(MAX_SLIP)
      real*8  zBar0(DIMS,DIMS,MAX_SLIP)
      real*8  pBar0(DIMS,DIMS,MAX_SLIP), qBar0(DIMS,DIMS,MAX_SLIP)
      real*8  pBar0Vec(NVECS,MAX_SLIP), qBar0Vec(DIMS,MAX_SLIP)
      real*8  ppTBar0(NVECS,NVECS,MAX_SLIP)
      real*8  cepmod(DIMS2,DIMS2)
      real*8  rss(MAX_SLIP)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)
      real*8  eps(9), eps_n(9)
      real*8  dfgrd1(3,3)

      integer is, iterState
      logical converged
      real*8  e_kk, norm_tau0, norm_kapp0, norm_tau, norm_kapp, epsdot
      real*8  d_vec_lat(NVECS), tau_lat(NVECS)
      real*8  qr5x5(NVECS,NVECS), qr3x3(DIMS,DIMS)
      real*8  pHatVec(NVECS,MAX_SLIP), qHatVec(DIMS,MAX_SLIP)
      real*8  ppTHat(NVECS,NVECS,MAX_SLIP)
      real*8  fCeiDevHat(NVECS,NVECS),fCeDevVolHat(NVECS)

      logical ConvergeState
      real*8  InnerProductVec
      
      real*8  switch(MAX_SLIP)
c
c---------------------------------------------------------------------72
c
c------------------------- DEVIATORIC RESPONSE
c
c------- effective deformation rate 
c
      epsdot = dsqrt(2./3.*InnerProductVec(d_vec, d_vec, NVECS))
      epsdot = max(epsdot, TINY)
c
c------- rotate Scmidt tensors (vectors) to Btilde_n configuration
c
      call RotMat5x5ForSymm(crot_n, qr5x5, DIMS)
      call RotMat3x3ForSkew(crot_n, qr3x3, DIMS)
      do is = 1, numslip
         call MultAxu(qr5x5, pBar0Vec(1,is), pHatVec(1,is), NVECS)
         call MultAxu(qr3x3, qBar0Vec(1,is), qHatVec(1,is), DIMS)
      enddo     
c
c------- initialize slip hardening at time t
c
      call EqualTensors(kappa_n, kappa, NKAPP)
c
c------- initial estimate for the stress
c
      if (subIncr .eq. 1) then

         if (globalIncr .eq. 1) then
c
c------------- initial guess from viscoplastic solution
c---------------- transform D from sample to crystal coordinates
            call MultATxu(qr5x5, d_vec, d_vec_lat, NVECS)
c
c---------------- viscoplastic solution in crystal coordinates
            call StressSolveViscoPlastic(tau_lat, d_vec_lat, kappa,
     &              pBar0Vec, ppTBar0, matProp, sigfs, tauSlip, dtime, 
     &              theta, thetao, epsdot, tolerNewt, maxIterNewt, 
     &              numslip, numvtx, switch)
c
c---------------- transform computed stresses from crystal to sample axis
            call MultAxu(qr5x5, tau_lat, tau, NVECS)
         else
c
c------------- initial guess from previous solution
            call EqualTensors(stress_n, tau, NVECS)
         endif

      endif
c
c------- initial estimate for hardness and rotation tensor
c
      call HardeningRotationSolve(kappa, kappa_n, eqvalues, rrot, 
     &        rrot_n, crot, crot0, drot, tau, w_vec, gamdot, rss,
     &        pHatVec, 
     &        qHatVec, matProp, hardmtx, epsdot, dtime, numslip, 
     &        kHARD_EXPL, switch)
c
c------- initial valuies for the norm of stress and kappa
c
      norm_tau0  = dsqrt(InnerProductVec(tau, tau, NVECS))
      norm_kapp0 = dsqrt(InnerProductVec(kappa, kappa, NKAPP))
c
c------- initialized global flag to monitor Newton/State convergence
c
      ierr = XTAL_CONVERGED
c
c------- predictor volumetric elastic strain in Btilde
c
      e_kk = statev_n(kEVOL) + dtime * d_kk
c
c------- iterate for the material state
c
      iterState    = 0
      iterCounterN = 0
      converged = .false.

      do while(iterState .lt. maxIterState  .and. .not.converged)

         iterState = iterState + 1
c
c---------- rotate Scmidt tensors to Btilde configuration
c
         call RotateSlipGeometry(crot, pBar0Vec, qBar0Vec, qr5x5, qr3x3,
     &                           pHatVec, qHatVec, ppTHat, numslip)
c
c---------- rotate deviatoric & dev/volumetric elasticities to Btilde:
c
         call MultQAQT(qr5x5, fCeiDev, fCeiDevHat, NVECS)
         call MultAxu(qr5x5, fCeDevVol, fCeDevVolHat, NVECS)
c
c---------- solve for the crystal stresses
c
         call StressSolveDeviatoric(tau, d_vec, estran, estran_n, kappa,
     &           gamdot, rss, drot, fCeiDevHat, fCeDevVolHat, pHatVec,
     &           ppTHat, matProp, e_kk, dtime, tolerNewt, numslip, 
     &           maxIterNewt, iterCounterN, ierr, switch)
         if (ierr .ne. XTAL_CONVERGED) return
         norm_tau = dsqrt(InnerProductVec(tau, tau, NVECS))
c
c---------- solve for hardness and rotation tensor
c
         call HardeningRotationSolve(kappa, kappa_n, eqvalues, rrot, 
     &           rrot_n, crot, crot0, drot, tau, w_vec, gamdot, rss,
     &           pHatVec,
     &           qHatVec, matProp, hardmtx, epsdot, dtime, numslip, 
c     &           kHARD_MIDP)
     &           kHARD_ANAL, switch)
         norm_kapp = dsqrt(InnerProductVec(kappa, kappa, NKAPP))
c
c---------- check convergence
c
         converged = ConvergeState(norm_tau, norm_kapp, norm_tau0,
     &                             norm_kapp0, tolerState)

      enddo
c
c------- keep track of state iteration and check number of iterations
c
      iterCounterS = iterState
      if (iterState .ge. maxIterState) then
         call WriteWarning(XTAL_O,
     &            'StressSolveElastoViscoPlastic: iters > maxIters')
         ierr = XTAL_MAX_ITERS_HIT
         return
      endif
c
c------- continue sub-incrementation process if needed
c
      if (subIncr .lt. totSubIncrs) return
c
c------------------------ VOLUMETRIC RESPONSE
c
c------- integrate crystal volumetric  response
c
c      call RotMat5x5ForSymm(crot, qr5x5, DIMS)
c      call MultAxu(qr5x5, fCeDevVol, fCeDevVolHat, NVECS)
      call StressSolveVolumetric(statev, fCeDevVolHat, fCeVol, estran, 
     &                           e_kk, dtime)
c
c----------------------- CAUCHY STRESS - CONSISTENT TANGENT
c
c------- compute Cauchy stress s_ij (note that "tau" is deviatoric)
c-------  s_ij = stress + statev(kPRES)*I  (i.e., stress is deviatoric)
c-------  ee_ij = estran + statev(KEVOL)*I  (i.e., estran is deviatoric)
c
      call UpdateStress(s_ij, stress, tau, ee_ij, estran, statev, 
     &                  eqvalues, rrot, dfgrd1, iqpt, ielem)
c
c------- compute norm of spin, slip system activity, and eqp
c
      call UpdateDeformQnts(gamdot, eqvalues, statev, crot, zBar0,
     &                      dtime, numslip, eps, eps_n)

c      call RotMat5x5ForSymm(crot, qr5x5, DIMS)
c      call MultQAQT(qr5x5, fCeiDev, fCeiDevHat, NVECS)
c      call MultAxu(qr5x5, fCeDevVol, fCeDevVolHat, NVECS)

      call PlasticModuli(cepmod, fCeiDevHat, fCeDevVolHat, fCeVol,
     &                   statev(kDETVe))

      return
      END

c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE RotateSlipGeometry(
     &   crot, pBar0Vec, qBar0Vec, qr5x5, qr3x3, pHatVec, qHatVec,
     &   ppTHat, numslip
     &   )

      implicit none
      include 'params_xtal.inc'

      integer numslip
      real*8  crot(DIMS,DIMS)
      real*8  pBar0Vec(NVECS,MAX_SLIP), qBar0Vec(DIMS,MAX_SLIP)
      real*8  qr5x5(NVECS,NVECS), qr3x3(DIMS,DIMS)
      real*8  pHatVec(NVECS,MAX_SLIP), qHatVec(DIMS,MAX_SLIP)
      real*8  ppTHat(NVECS,NVECS,MAX_SLIP)

      integer is
c
c---------------------------------------------------------------------72
c
c------- rotation matrices for symm and skew quantities
c
      call RotMat5x5ForSymm(crot, qr5x5, DIMS)
      call RotMat3x3ForSkew(crot, qr3x3, DIMS)
c
c------- rotate Scmidt tensors (vectors) to Btilde configuration
c
      do is = 1, numslip
         call MultAxu(qr5x5, pBar0Vec(1,is), pHatVec(1,is), NVECS)
         call MultAxu(qr3x3, qBar0Vec(1,is), qHatVec(1,is), DIMS)
         call OuterProductVec(pHatVec(1,is), pHatVec(1,is), 
     &                        ppTHat(1,1,is), NVECS)
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE HardeningRotationSolve(
     &   kappa, kappa_n, eqvalues, rrot, rrot_n, crot, crot0, drot, 
     &   tau, w_vec, gamdot, rss, pHatVec, qHatVec, matProp, hardmtx, 
     &   epsdot, dtime, numslip, kInteg_Hard, switch
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip, kInteg_Hard
      real*8  epsdot, dtime
      real*8  kappa(NKAPP), kappa_n(NKAPP), eqvalues(NEQVA)
      real*8  rrot(DIMS,DIMS), rrot_n(DIMS,DIMS), crot(DIMS,DIMS)
      real*8  crot0(DIMS,DIMS), drot(DIMS,DIMS), tau(NVECS), w_vec(DIMS)
      real*8  gamdot(MAX_SLIP), matProp(NPROPS,MAX_SLIP)
      real*8  pHatVec(NVECS,MAX_SLIP), qHatVec(DIMS,MAX_SLIP)
      real*8  rss(MAX_SLIP)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)

      integer is
      real*8  InnerProductVec, SSKineticEqn

      real*8  crss
      real*8  overstress(MAX_SLIP)
      common /XtalOverS/ overstress
      
      real*8  switch(MAX_SLIP)
      integer counter
c
c---------------------------------------------------------------------72
c
c------- slip quantities: resolve shear stress and shear strain rate
c
      do is = 1, numslip
         rss(is) = InnerProductVec(tau, pHatVec(1,is), NVECS)
c         crss = kappa(1) + overstress(is)
         crss = kappa(is)
         gamdot(is) = SSKineticEqn(rss(is),crss,matProp(1,is),kGAMDOT, 
     &                             switch(is))
      enddo
c
c------- solve for hardness
c
c      call IntegrateHardening(kappa, kappa_n, eqvalues, gamdot, matProp,
c     &                     hardmtx, epsdot, dtime, numslip, kInteg_Hard)
c
c------- solve for rotation: dotR = (W - Wp)*R
c------- update slip system orientation matrix: C = R * C0
c
      call IntegrateRotation(w_vec, rrot_n, crot0, qHatVec, gamdot, 
     &                       rrot, crot, drot, dtime, numslip)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE IntegrateHardening(
     &   kappa, kappa_n, eqvalues, gamdot, matProp, hardmtx, epsdot, 
     &   dtime, numslip, kInteg_Code
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip, kInteg_Code
      real*8  epsdot, dtime
      real*8  kappa(NKAPP), kappa_n(NKAPP), eqvalues(NEQVA)
      real*8  gamdot(MAX_SLIP), matProp(NPROPS, MAX_SLIP)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)

      integer is, ik, jk
      real*8  h_0, tausi, taus0, xms, gamss0
      real*8  shr_min, shr_max, kappa_sat
      real*8  c, g_n, g_s, g
      real*8  dkappa, gamtot_n, delgam, fac
      real*8  kTHETA
      data    kTHETA /1.0d0/
c
c---------------------------------------------------------------------72
c
c------- total shear strain rate: Sum(|gamdot(i)|), i=1,numslip
c
      eqvalues(kSHRATE) = pzero
      do is = 1, numslip
         eqvalues(kSHRATE) = eqvalues(kSHRATE) + dabs(gamdot(is))
      enddo

      shr_min = 1.0d-6 * epsdot
      shr_max = 1.0d+1 * epsdot

      if (eqvalues(kSHRATE) .le. shr_min) eqvalues(kSHRATE) = shr_min
      if (eqvalues(kSHRATE) .ge. shr_max) eqvalues(kSHRATE) = shr_max
c
c------- accumulated shear strain: gamtot
c
      eqvalues(kGAMTOT) = eqvalues(kGAMTOT_n) + eqvalues(kSHRATE)*dtime
c
c------- integration of Voce's hardening law (one hardness/slip system)
c
      if (kInteg_Code .eq. kHARD_EXPL) then  
c
c---------- explicit update
c         do ik = 1, NKAPP
         do ik = 1, numslip
            h_0    = matProp(5,ik)
            tausi  = matProp(6,ik)
            taus0  = matProp(7,ik)
            xms    = matProp(8,ik)
            gamss0 = matProp(9,ik)

            kappa_sat = taus0 * ((eqvalues(kSHRATE) / gamss0)**xms)

            c = dtime*h_0
            g_s = kappa_sat - tausi
            g_n = kappa_n(ik) - tausi

            if ( (g_n/g_s) .le. 1.0 ) then
               g = g_n + c*(1.0-g_n/g_s)*eqvalues(kSHRATE)
            else
               g = g_n
            endif
            kappa(ik) = g + tausi
         enddo

      else if (kInteg_Code .eq. kHARD_MIDP) then
c
c---------- generalized mid-point rule 
c         do ik = 1, NKAPP
         do ik = 1, numslip
            h_0    = matProp(5,ik)
            tausi  = matProp(6,ik)
            taus0  = matProp(7,ik)
            xms    = matProp(8,ik)
            gamss0 = matProp(9,ik)

            kappa_sat = taus0 * ((eqvalues(kSHRATE) / gamss0)**xms)

            c = dtime*h_0
            g_s = kappa_sat - tausi
            g_n = kappa_n(ik) - tausi

            if ( (g_n/g_s) .le. 1.0 ) then
               g = g_n + c*( 
     &                (1.0-kTHETA)*(1.0-g_n/g_s)*eqvalues(kSHRATE_n)
     &                +   (kTHETA)*eqvalues(kSHRATE) 
     &                     )
               g = g / (1.0 + c*kTHETA*eqvalues(kSHRATE)/g_s)
            else
               g = g_n
            endif
            kappa(ik) = g + tausi
         enddo

      else if (kInteg_Code .eq. kHARD_ANAL) then

         gamtot_n = eqvalues(kGAMTOT_n)
         delgam   = eqvalues(kSHRATE) * dtime

         do ik = 1, numslip

            h_0    = matProp(5,ik)
            tausi  = matProp(6,ik)
            taus0  = matProp(7,ik)
            xms    = matProp(8,ik)
            gamss0 = matProp(9,ik)

            kappa_sat = taus0 * ((eqvalues(kSHRATE) / gamss0)**xms)
            g_s = kappa_sat - tausi
            fac = dabs(h_0/g_s)

            dkappa = 0.0
            do jk = 1, numslip
              dkappa = dkappa + 
     &                    hardmtx(ik,jk)*dabs(gamdot(jk))*dtime/delgam
            enddo
            dkappa = dkappa * g_s * dexp(-gamtot_n*fac) * 
     &                                        (1.0 - dexp(-delgam*fac))

            kappa(ik) = kappa_n(ik) + dkappa

         enddo
         
      else
c
c------- wrong code number
         call RunTimeError(XTAL_O, 
     &                      'IntegrateHardening: Wrong kInteg_Code!')

      endif
      write (XTAL_O, *)'running?'
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE IntegrateRotation(w_vec, rrot_n, crot0, qHatVec, 
     &   gamdot, rrot, crot, drot, dtime, numslip
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip
      real*8  dtime
      real*8  w_vec(DIMS), rrot_n(DIMS,DIMS), crot0(DIMS,DIMS)
      real*8  qHatVec(DIMS,MAX_SLIP), gamdot(MAX_SLIP)
      real*8  rrot(DIMS,DIMS), crot(DIMS,DIMS), drot(DIMS,DIMS)

      integer is
      real*8  wpHat(DIMS), spin(DIMS)
c
c---------------------------------------------------------------------72
c
c------- plastic spin from slip system activity
c
      call SetTensor(wpHat, pzero, DIMS)
      do is = 1, numslip
         call AddScaledTensor(gamdot(is), qHatVec(1,is), wpHat, DIMS)
      enddo
c
c------- net spin: (W - Wp)
c
      call SetTensor(spin, pzero, DIMS)
      call AddTensors(pone, w_vec, -pone, wpHat, spin, DIMS)
c
c------- incremental rot: dR = exp(spin*dtime)
c
      call IncrementalRotTensor(drot, spin, dtime)
c
c------- rotation tensor: R = dR*R_n
c
      call MultAxB(drot, rrot_n, rrot, DIMS)
c
c------- slip system orientation matrix crot: C = R*C_0
c
      call MultAxB(rrot, crot0, crot, DIMS)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE IncrementalRotTensor(
     &   drot, spin, dtime
     &   )

      implicit none

      real*8  dtime
      real*8  drot(3, 3), spin(3)

      real*8  th1, th2, th3, tau, taua
c
c---------------------------------------------------------------------72
c
      th1 = spin(1) * dtime
      th2 = spin(2) * dtime
      th3 = spin(3) * dtime

      tau = sqrt(th1 * th1 + th2 * th2 + th3 * th3)
      if (tau .ne. 0.0) then
         taua = tan(tau / 2.0) / tau
         th1 = taua * th1
         th2 = taua * th2
         th3 = taua * th3
         tau = taua * tau
      endif

      tau = 2.0 / (1.0 + tau * tau)

      drot(1, 1) = 1.0 - tau * (th1 * th1 + th2 * th2)
      drot(1, 2) = - tau * (th1 + th2 * th3)
      drot(1, 3) = tau * ( - th2 + th1 * th3)
      drot(2, 1) = tau * (th1 - th2 * th3)
      drot(2, 2) = 1.0 - tau * (th1 * th1 + th3 * th3)
      drot(2, 3) = - tau * (th3 + th1 * th2)
      drot(3, 1) = tau * (th2 + th1 * th3)
      drot(3, 2) = tau * (th3 - th1 * th2)
      drot(3, 3) = 1.0 - tau * (th2 * th2 + th3 * th3)
 
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      logical FUNCTION ConvergeState(
     &   norm_a, norm_b, norm_a0, norm_b0, toler
     &   )

      implicit none

      real*8  norm_a, norm_b, norm_a0, norm_b0, toler
c
c---------------------------------------------------------------------72
c
      ConvergeState = (dabs(norm_a - norm_a0) .lt. toler*norm_a0) .and. 
     &                (dabs(norm_b - norm_b0) .lt. toler*norm_b0)

      if (.not.ConvergeState) then
          norm_a0 = norm_a
          norm_b0 = norm_b
      endif

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE UpdateStress(
     &   s_ij, stress, tau, ee_ij, estran, statev, eqvalues, rrot,
     &   dfgrd1, iqpt, ielem)

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      real*8  s_ij(DIMS,DIMS), stress(NVECS), tau(NVECS)
      real*8  ee_ij(DIMS,DIMS), estran(NVECS), statev(NSTAV)
      real*8  eqvalues(NEQVA), rrot(DIMS,DIMS), dfgrd1(3,3)
      integer iqpt, ielem

      real*8  ekk3, det
      real*8  Ve_ij(DIMS, DIMS)
      real*8  matx_1(3,3), matx_2(3,3), matx_3(3,3)
      real*8  F_p(3,3)

      real*8  DetOfMat3x3, InnerProductVec
      
      real*8  plasdef1(3, 3, 1, MAX_GRN)
      common  /PlasticDefGrad1/ plasdef1
      save    /PlasticDefGrad1/
c
c---------------------------------------------------------------------72
c
c------- Elastic strain (ee) and V tensors (V = ee + I)
c
      call Vec5x1ToMat3x3Symm(estran, ee_ij, DIMS)
      ekk3 = statev(kEVOL)/3.0
      call AddScaledTensor(ekk3, Ident2nd, ee_ij, DIMS*DIMS)

      call AddTensors(pone, Ident2nd, pone, ee_ij, Ve_ij, DIMS*DIMS)
      statev(kDETVe) = DetOfMat3x3(Ve_ij)
c
c------- plastic deformation gradient
c
      call InverseOfMat3x3(rrot, matx_1, det)
      call InverseOfMat3x3(Ve_ij, matx_2, det)
      call MultAxB(matx_1, matx_2, matx_3, 3)
      call MultAxB(matx_3, dfgrd1, F_p, 3)
      call EqualTensors(F_p, plasdef1(1, 1, 1, ielem), 9)
c
c------- deviatoric (vector) and volumetric (pressure) Cauchy stress
c
      call SetToScaledtensor(1./statev(kDETVe), tau, stress, NVECS)
      statev(kPRES) = statev(kPRES)/statev(kDETVe)
c
c------- VonMises stress
c
      eqvalues(kMISES) = dsqrt(3./2.*InnerProductVec(stress, stress,
     &                                                NVECS))
c
c------- Cauchy stress (tensor form)
c
      call Vec5x1ToMat3x3Symm(stress, s_ij, DIMS)
      call AddScaledTensor(statev(kPRES), Ident2nd, s_ij, DIMS*DIMS)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE PlasticModuli(
     &   fCep, fCeiDevHat, fCeDevVolHat, fCeVol, jac_e
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      real*8  fCeVol, jac_e
      real*8  fCep(DIMS2,DIMS2)
      real*8  fCeiDevHat(NVECS,NVECS), fCeDevVolHat(NVECS)

      real*8  det
      real*8  fMat5x5(NVECS, NVECS), fVect5(NVECS)
      real*8  fOuter6x6(DIMS2, DIMS2), fVect6(DIMS2)
      real*8  fMat5x6(NVECS, DIMS2), fOuter5x6(NVECS, DIMS2)
      real*8  fCepDev(NVECS, DIMS2), fCepVol(DIMS2)

      real*8  fDevMat5x6(NVECS, DIMS2), fDevMat6x5(DIMS2, NVECS)
      real*8  fMatTId5x6(NVECS, DIMS2)
      common /transfMatxs/ fDevMat5x6, fDevMat6x5, fMatTId5x6

      real*8  lhs5x5(NVECS, NVECS)
      common /SumGamPPt/ lhs5x5
c
c
c---------------------------------------------------------------------72
c
c---- Deviatoric moduli: fCepDev_5x6
c
c-------------- ([Ce_d]^-1 + [S])^-1
      call AddTensors(pone, fCeiDevHat, pone, lhs5x5, fMat5x5, 
     &                NVECS*NVECS) 
      call MatrixInverse(fMat5x5, NVECS, NVECS, det)
c
c-------------- ([T][Pdev] + [Ce_d]^-1 {H} {1}^T)
      call MultAxu(fCeiDevHat, fCeDevVolHat, fVect5, NVECS)
      call OuterProductVec_G(fvect5, Ident, fOuter5x6, NVECS, DIMS2)
      call Addtensors(pone, fMatTId5x6, pone, fOuter5x6, fMat5x6, 
     &                NVECS*DIMS2)
c
c-------------- (([Ce_d]^-1 + [S])^-1) ([T][Pdev] + [Ce_d]^-1 {H} {1}^T)
      call MultAxB_G(fMat5x5, fMat5x6, fCepDev, NVECS, NVECS, DIMS2)
c
c---- Volumetric moduli: fCepVol_6
c
c-------------- {H}^T [T] [Pdev] + M {1}^T
      call MultATxu_G(fMatTId5x6, fCeDevVolHat, fVect6, NVECS, DIMS2)
      call AddTensors(pone, fVect6, fCeVol, Ident, fCepVol, DIMS2)
c
c-------------- {H}^T [S] [Cep_d]
      call MultAxB_G(lhs5x5, fCepDev, fMat5x6, NVECS, NVECS, DIMS2)
      call MultATxu_G(fMat5x6, fCeDevVolHat, fVect6, NVECS, DIMS2)
c
c-------------- {H}^T [T] [Pdev] + M {1}^T - {H}^T [S] [Cep_d]
      call AddScaledTensor(-pone, fVect6, fCepVol, DIMS2)
c
c---- Algorithmic Moduli: fCep_6x6
c
c-------------- [J] [Cep_d]
      call MultAxB_G(fDevMat6x5, fCepDev, fCep, DIMS2, NVECS, DIMS2)
c
c-------------- [J] [Cep_d] + {1} {Cep_v}^T -> [Cep]: {s} = [Cep] {d}
      call OuterProductVec(Ident, fCepVol, fOuter6x6, DIMS2)
      call AddScaledTensor(pone, fOuter6x6, fCep, DIMS2*DIMS2)
      call SetToScaledTensor(pone/jac_e, fCep, fCep, DIMS2*DIMS2)
c
c---- Change fCep to be consistent with abaqus representation of {d}
c
      call SetToScaledTensor(phalf, fCep(1,4), fCep(1,4), DIMS2*DIMS2/2)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE UpdateDeformQnts(
     &   gamdot, eqvalues, statev, crot, zBar0, dtime, numslip,
     &   eps, eps_n)

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip
      real*8  dtime
      real*8  statev(NSTAV), eqvalues(NEQVA), gamdot(MAX_SLIP)
      real*8  crot(DIMS,DIMS), zBar0(DIMS,DIMS,MAX_SLIP)

      integer is
      real*8  eqpdot, maxGamdot, ratio
      real*8  wp_vec(DIMS), matx_1(DIMS,DIMS)
      real*8  lp_ij(DIMS,DIMS), dp_ij(DIMS,DIMS), wp_ij(DIMS, DIMS)
      real*8  eps(9), eps_n(9)

      real*8  InnerProductVec, MaxAbsValueOfVec
c
c---------------------------------------------------------------------72
c
c------- Plastic velocity gradient lp in crystal axis
c
      call setTensor(lp_ij, pzero, DIMS*DIMS)
      do is = 1, numslip
         call AddScaledTensor(gamdot(is), zBar0(1,1,is), lp_ij, 
     &                        DIMS*DIMS)
      enddo
c
c------- Symm and Skew parts of lp: lp =  dp + wp
c
c      call EqualTensors(lp_ij, matx_1, DIMS*DIMS)
c      call MultQAQT(crot, matx_1, lp_ij, DIMS)

      call SymmetrizeTensor(lp_ij, dp_ij, DIMS)
      call SkewSymmetrizeTensor(lp_ij, wp_ij, DIMS)
c
c------- Equivalent plastic strain
c
      eps(1) = eps_n(1) + dtime * dp_ij(1,1)
      eps(2) = eps_n(2) + dtime * dp_ij(1,2)
      eps(3) = eps_n(3) + dtime * dp_ij(1,3)
      eps(4) = eps_n(4) + dtime * dp_ij(2,1)
      eps(5) = eps_n(5) + dtime * dp_ij(2,2)
      eps(6) = eps_n(6) + dtime * dp_ij(2,3)
      eps(7) = eps_n(7) + dtime * dp_ij(3,1)
      eps(8) = eps_n(8) + dtime * dp_ij(3,2)
      eps(9) = eps_n(9) + dtime * dp_ij(3,3)
      eqvalues(kEQP) = dsqrt(2./3.*InnerProductVec(eps, eps, DIMS*DIMS))
c
c------- Norm of plastic spin
c
      call Mat3x3ToVec3x1Skew(wp_ij, wp_vec, DIMS)
      statev(kWPNORM) = dsqrt(InnerProductVec(wp_vec, wp_vec, DIMS))
c
c------- Slip system activity
c
      maxGamdot = MaxAbsValueOfVec(gamdot, numslip)
      if (maxGamdot .lt. TINY) maxGamdot = TINY

      statev(kSSACT) = 0.0
      do is = 1, numslip
         ratio = dabs(gamdot(is))/maxGamdot
         if (ratio .ge. 0.05d0) statev(kSSACT) = statev(kSSACT) + pone
      enddo
        
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE StressSolveViscoPlastic(
     &   tau_lat, d_vec_lat, kappa, pBar0Vec, ppTBar0, matProp, sigfs,
     &   tauSlip, dtime, theta, thetao, epsdot, tolerNewt, maxIterNewt,
     &   numslip, numvtx, switch
     &   )

      implicit none
      include 'params_xtal.inc'

      integer maxIterNewt, numslip, numvtx
      real*8  dtime, theta, thetao, epsdot, tolerNewt
      real*8  tau_lat(NVECS), d_vec_lat(NVECS), kappa(NKAPP)
      real*8  pBar0Vec(NVECS,MAX_SLIP), ppTBar0(NVECS,NVECS,MAX_SLIP)
      real*8  sigfs(NVECS,MAX_VTX), tauSlip(MAX_SLIP)
      real*8  matProp(NPROPS,MAX_SLIP)

      integer iv, ierr
      real*8  rescale
      real*8  sdotd(MAX_VTX)
      real*8  InnerProductVec
      
      real*8  switch(MAX_SLIP)
c
c---------------------------------------------------------------------72
c
c------- solve for stresses in lattice reference frame
c
c------- use effective deformation rate to normalized (scale-down) D
c
c      epsdot = dsqrt(2./3.*InnerProductVec(d_vec_lat, d_vec_lat, NVECS))
      call SetToScaledTensor(1./epsdot, d_vec_lat, d_vec_lat, NVECS)
c
c------- compute work to select guessed values from vertex stresses
c
      do iv = 1, numvtx
         sdotd(iv) = InnerProductVec(sigfs(1,iv), d_vec_lat, NVECS) 
      enddo
c
c------- initialized local flag to monitor convergence of Newton iters
c
      ierr = XTAL_CONVERGED
c
c------- solve fot the stresses. If needed, will try all vertices
c
      do iv = 1, numvtx
c
c---------- compute and scaled initial guess
         call InitialGuessStress(tau_lat, sdotd, kappa, sigfs, pBar0Vec,
     &                           matProp, numslip, numvtx)
c
c---------- compute stresses
         call StressViscoPlastic(tau_lat, d_vec_lat, kappa, pBar0Vec, 
     &                           ppTBar0, matProp, tolerNewt,
     &                           maxIterNewt, numslip, ierr, switch)
c 
c---------- if not converged, try another initial guess
         if (ierr .ne. XTAL_CONVERGED) then
            call WriteMessage(XTAL_O, 
     &         'StressSolveViscoPlastic: did not converged !!!')
            call WriteMessage(XTAL_O, 
     &         'StressSolveViscoPlastic: will try next vertex stresses')
            ierr = XTAL_CONVERGED
c 
c---------- or else converged; then, rescale the stresses, and return
         else
            rescale = epsdot**matProp(3,1)   ! epsdot**m
            call SetToScaledTensor(rescale, tau_lat, tau_lat, NVECS)
            return
         endif

      enddo
c
c------- if gets here, then the stress solution did not converged
c-------   for all trials. 
c
      call RunTimeError(XTAL_O, 
     &       'StressSolveViscoPlastic: Newton method did not converged')

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE InitialGuessStress(
     &   s, sdotd, kappa, sigfs, pBar0Vec, matProp, numslip, numvtx
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip, numvtx
      real*8  s(NVECS), sdotd(MAX_VTX), kappa(NKAPP)
      real*8  sigfs(NVECS, MAX_VTX), pBar0Vec(NVECS, MAX_SLIP)
      real*8  matProp(NPROPS,MAX_SLIP)

      integer kmax, is
      real*8  signmax, factor
      real*8  rssk(MAX_SLIP)

      integer IndexMaxAbsValueOfVec
      real*8  MaxAbsValueOfVec, InnerProductVec, SignOf

      real*8  crss
      real*8  overstress(MAX_SLIP)
      common /XtalOverS/ overstress
c
c---------------------------------------------------------------------72
c
c------- locate sign of max value of sdotd (plastic work)
c
      kmax = IndexMaxAbsValueOfVec(sdotd, numvtx)
      signmax = SignOf(sdotd(kmax))
      sdotd(kmax) = pzero        ! avoids 'kmax' vertex being re-used
c
c------- asign correct sign to selected vertex stress (initial guess)
c
      call SetToScaledTensor(signmax, sigfs(1,kmax), s, NVECS)
c
c------ compute factor to scaled down initial guess
c------   factor = |rss|_max / kappa * gam0**m
c
      do is = 1, numslip
c         crss = kappa(1) + overstress(is)
         crss = kappa(is)
         rssk(is) = InnerProductVec(pBar0Vec(1,is), s, NVECS)/crss
      enddo
      factor = MaxAbsValueOfVec(rssk, numslip) *
     &                                      matProp(4,1)**matProp(3,1)
c
c------- scaled initial guess for stress
c
      call SetToScaledTensor(1./factor, s, s, NVECS)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE StressViscoPlastic(
     &   s, d_vec_lat, kappa, pBar0Vec, ppTBar0, matProp, toler,
     &   maxIters, numslip, ierr, switch 
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer maxIters, numslip, ierr
      real*8  toler
      real*8  s(NVECS), d_vec_lat(NVECS), kappa(NKAPP)
      real*8  pBar0Vec(NVECS,MAX_SLIP), ppTBar0(NVECS,NVECS,MAX_SLIP)
      real*8  matProp(NPROPS,MAX_SLIP)

      integer iters
      real*8  rhs_norm_0, rhs_norm, search
      real*8  delts(NVECS), s0(NVECS)
      real*8  rhs(NVECS), rss(MAX_SLIP), lhs(NVECS,NVECS)

      logical Converged
c      real*8  InnerProductVec
      real*8  NormOfVector
      
      real*8  switch(MAX_SLIP)
c
c---------------------------------------------------------------------72
c
c------- compute initial residual and its norm
c
      call FormResidual(rhs, s, rss, d_vec_lat, kappa, pBar0Vec, 
     &                  matProp, numslip, switch)
c      rhs_norm_0 = dsqrt(InnerProductVec(rhs, rhs, NVECS))
      rhs_norm_0 = NormOfVector(rhs, NVECS)
c
c------- initialize flags and start iterations
c
      iters = 0
      do while(.not.Converged(rhs, toler) .and. iters .lt. maxIters)
         iters = iters + 1
         call EqualTensors(s, s0, NVECS)
c
c---------- compute local jacobian
         call FormJacobian(lhs, rss, ppTBar0, kappa, matProp, numslip, 
     &                     switch)
c
c---------- solve for the increment of stress
c---------- use the lower half of the symmetric matrix lhs
         call lsolve(lhs, rhs, 1, ierr, NVECS)
         if (ierr .eq. XTAL_SING_JACOBIAN) return

         call EqualTensors(rhs, delts, NVECS)
c
c---------- update variables
         search = pone
         call AddTensors(pone, s0, search, delts, s, NVECS)
c
c---------- compute new residual and its norm
         call FormResidual(rhs, s, rss, d_vec_lat, kappa, pBar0Vec, 
     &                     matProp, numslip, switch)
c         rhs_norm = dsqrt(InnerProductVec(rhs, rhs, NVECS))
         rhs_norm = NormOfVector(rhs, NVECS)
c
c---------- simple line search
         do while ( rhs_norm .gt. rhs_norm_0 ) 

            search = search*0.5d0
            if (search .lt. TINY) then
               call WriteWarning(XTAL_O, 
     &                 'StressViscoPlastic: LS Failed, search < TINY')
               ierr = XTAL_LS_FAILED
               return
            endif
            call AddTensors(pone, s0, search, delts, s, NVECS)

            call FormResidual(rhs, s, rss, d_vec_lat, kappa, pBar0Vec, 
     &                        matProp, numslip, switch)
c            rhs_norm = dsqrt(InnerProductVec(rhs, rhs, NVECS))
            rhs_norm = NormOfVector(rhs, NVECS)

         enddo
c
c---------- update norm of residual
         rhs_norm_0 = rhs_norm

      enddo
c
c------- check convergence on iters < maxIters
c
      if (iters .ge. maxIters) then
        call WriteMessage(XTAL_O,'StressViscoPlastic: iters > maxIters')
        ierr = XTAL_MAX_ITERS_HIT
        return
      endif
c
c------- keep track of Newton's iterations
c
      write(XTAL_O, *) 'iterCounter (viscoplastic) = ', iters

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE FormResidual(
     &   rhs, s, rss, d_vec_lat, kappa, pBar0Vec, matProp, numslip,
     &   switch)

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip
      real*8  rhs(NVECS), s(NVECS), rss(MAX_SLIP)
      real*8  d_vec_lat(NVECS), kappa(NKAPP)
      real*8  pBar0Vec(NVECS,MAX_SLIP), matProp(NPROPS,MAX_SLIP)

      integer is
      real*8  gamdot(MAX_SLIP), d_curr(NVECS)
      real*8  InnerProductVec, SSKineticEqn

      real*8  crss
      real*8  overstress(MAX_SLIP)
      common /XtalOverS/ overstress
      
      real*8  switch(MAX_SLIP)
c
c---------------------------------------------------------------------72
c
c------- resolve shear stresses and shear strain rates
c
      do is = 1, numslip
         rss(is) = InnerProductVec(s, pBar0Vec(1,is), NVECS)
c         crss = kappa(1) + overstress(is)
         crss = kappa(is)
         gamdot(is) = SSKineticEqn(rss(is),crss,matProp(1,is),kGAMDOT, 
     &                             switch(is))
      enddo
c
c------- rate of deformation due to slip system activity
c
      call SetTensor(d_curr, pzero, NVECS)
      do is = 1, numslip
         call AddScaledTensor(gamdot(is), pBar0Vec(1,is), d_curr, NVECS)
      enddo
c
c------- residual
c
      call SetTensor(rhs, pzero, NVECS)
      call AddTensors(+pone, d_vec_lat, -pone, d_curr, rhs, NVECS)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE FormJacobian(
     &   lhs, rss, ppTBar0, kappa, matProp, numslip, switch
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip
      real*8  lhs(NVECS,NVECS), rss(MAX_SLIP), kappa(NKAPP)
      real*8  ppTBar0(NVECS,NVECS,MAX_SLIP), matProp(NPROPS,MAX_SLIP)

      integer is
      real*8  SSKineticEqn
      real*8  dGamdTau(MAX_SLIP)

      real*8  crss
      real*8  overstress(MAX_SLIP)
      common /XtalOverS/ overstress
      
      real*8  switch(MAX_SLIP)
c
c---------------------------------------------------------------------72
c
c------- d(gamdot)/d(tau)
c
      do is = 1, numslip
c         crss = kappa(1) + overstress(is)
         crss = kappa(is)
         dGamdTau(is) = SSKineticEqn(rss(is),crss,matProp(1,is),
     &                               kdGAMdTAU, switch(is))
      enddo
c
c------- jacobian
c
      call SetTensor(lhs, pzero, NVECS*NVECS)
      do is = 1, numslip
         call AddScaledTensor(dGamdTau(is), ppTBar0(1,1,is), lhs,
     &                        NVECS*NVECS)
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c

      logical FUNCTION Converged(
     &   res, toler
     &   )

      implicit none
      include 'params_xtal.inc'

      real*8  toler
      real*8  res(NVECS)

      integer i
c
c---------------------------------------------------------------------72
c     
c------- check convergence on residual
c
      Converged = ( dabs(res(1)) .lt. toler )
      do i = 2, NVECS
         Converged = ( (dabs(res(i)) .lt. toler) .and. Converged)
      enddo

      return
      END
c     
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE StressSolveDeviatoric(
     &   tau, d_vec, estran, estran_n, kappa,
     &   gamdot, rss, drot, fCeiDevHat, fCeDevVolHat, pHatVec,
     &   ppTHat, matProp, e_kk, dtime, tolerNewt, numslip,
     &   maxIterNewt, iterCounterN, ierr, switch
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip, maxIterNewt, iterCounterN, ierr
      real*8  e_kk, dtime, tolerNewt
      real*8  tau(NVECS), d_vec(NVECS), estran(NVECS), estran_n(NVECS)
      real*8  kappa(NKAPP), gamdot(MAX_SLIP), drot(DIMS,DIMS)
      real*8  fCeiDevHat(NVECS,NVECS), fCeDevVolHat(NVECS)
      real*8  pHatVec(NVECS), ppTHat(NVECS,NVECS)
      real*8  matProp(NPROPS, MAX_SLIP)
      real*8  rss(MAX_SLIP)

      integer subIncr, totSubIncrs, is
      real*8  xm_o(MAX_SLIP), tmp
      real*8  tau_o(NVECS), tau_save(NVECS)
      
      real*8  switch(MAX_SLIP)
c
c---------------------------------------------------------------------72
c
c------- counters for sub-increments / backup some variables
c
      subIncr = 1
      totSubIncrs = 1
      do is = 1, numslip
        xm_o(is) = matProp(3,is)
      enddo
      call EqualTensors(tau, tau_o, NVECS)
c
c------- integrate crystal constitutive equations
c
      call SolveDeviatoric(tau, d_vec, estran, estran_n, kappa,
     &           gamdot, rss, drot, fCeiDevHat, fCeDevVolHat, pHatVec,
     &           ppTHat, matProp, e_kk, dtime, tolerNewt, numslip,
     &           maxIterNewt, iterCounterN, ierr, switch)
c
c------- if converged -> return, else do -> continuation method in "m"
c------- NOTE: here, continuation improves the initial guess 
c-------       for the solution variables 'tau'
c
      if (ierr .eq. XTAL_CONVERGED) return

      call WriteMessage
     &     (XTAL_O, 'DriverStressSolveDev: using Contin. for m!')
c
c------- loop for sub-incrementation procedure
c
      do while (.true.)
c
c---------- if not converged, increase # of sub-increments 
         if (ierr .ne. XTAL_CONVERGED) then
            subIncr = 2 * subIncr - 1
            totSubIncrs = 2 * totSubIncrs
            if (totSubIncrs .gt. 16) then
               call WriteMessage(XTAL_O,
     &                  'DriverStressSolveDev: totSubIncrs > 16')
               do is = 1, numslip
                  matProp(3,is) = xm_o(is)
                  call BoundForArgPowLaw(matProp(3,is))
               enddo
               return
            endif
            if (subIncr .eq. 1) 
     &           call EqualTensors(tau_o, tau, NVECS)
            if (subIncr .gt. 1) 
     &           call EqualTensors(tau_save, tau, NVECS)
c
c---------- if converged, adjust subincrements
         else if (subIncr .lt. totSubIncrs) then
            if ( (subIncr/2*2) .eq. subIncr) then
               subIncr = subIncr / 2 + 1
               totSubIncrs = totSubIncrs / 2
            else
               subIncr = subIncr + 1
            endif
c
c---------- successful return for continuation method
         else
            call WriteMessage(XTAL_O, 'Contin. successful !!!')
            return
         endif
c
c---------- trial rate sensitivity coefficient m
         tmp = real(subIncr) / real(totSubIncrs)   
         do is = 1, numslip
            matProp(3,is) = xm_o(is) / tmp
            call BoundForArgPowLaw(matProp(3,is))
         enddo
c
c---------- save current convergent solution before getting next solution
         if (subIncr .gt. 1 .and. ierr .eq. XTAL_CONVERGED)
     &                 call EqualTensors(tau, tau_save, NVECS)
c
c---------- report sub-incrs status
         write(XTAL_O, 1000) subIncr, totSubIncrs, matProp(3,1)
c
c---------- integrate crystal constitutive equations
c        iterCounterN = 0
         ierr = XTAL_CONVERGED
         call SolveDeviatoric(tau, d_vec, estran, estran_n, kappa,
     &           gamdot, rss, drot, fCeiDevHat, fCeDevVolHat, pHatVec,
     &           ppTHat, matProp, e_kk, dtime, tolerNewt, numslip,
     &           maxIterNewt, iterCounterN, ierr, switch)

      enddo

1000  format(3x, 'subInc/totSubIncrs = ', i2, '/', i3, 5x, 
     &      'm(1)=', e12.5)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SolveDeviatoric(
     &   tau, d_vec, estran, estran_n, kappa, gamdot, rss, drot,
     &   fCeiDevHat, fCeDevVolHat, pHatVec, ppTHat, matProp, e_kk, 
     &   dtime, toler, numslip, maxIters, iterCounterN, ierr, switch
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip, maxIters, iterCounterN, ierr
      real*8  e_kk, dtime, toler
      real*8  tau(NVECS), d_vec(NVECS), estran(NVECS), estran_n(NVECS)
      real*8  kappa(NKAPP), gamdot(MAX_SLIP), drot(DIMS,DIMS)
      real*8  fCeiDevHat(NVECS,NVECS), fCeDevVolHat(NVECS)
      real*8  pHatVec(NVECS), ppTHat(NVECS,NVECS)
      real*8  matProp(NPROPS, MAX_SLIP)
      real*8  rss(MAX_SLIP)

      integer iters
      real*8  rhs_norm_0, rhs_norm, search
      real*8  estranHat(NVECS), estranStar(NVECS)
      real*8  deltau(NVECS), tau0(NVECS)
      real*8  rhs(NVECS), lhs(NVECS,NVECS)
      real*8  dqr5x5(NVECS,NVECS)

      logical Converged
c      real*8  InnerProductVec
      real*8  NormOfVector
      
      real*8  switch(MAX_SLIP)
c
c---------------------------------------------------------------------72
c
c------- predictor elastic strain in Btilde: dR*e_n*dR^T + D*dt 
c
      call RotMat5x5ForSymm(drot, dqr5x5, DIMS)
      call MultAxu(dqr5x5, estran_n, estranHat, NVECS)
      call AddTensors(pone, estranHat, dtime, d_vec, estranStar, NVECS)
c
c------- compute initial residual 
c
      call ComputeResidual(rhs, tau, estran, estranStar, fCeiDevHat, 
     &                     fCeDevVolHat, rss, gamdot, kappa, pHatVec, 
     &                     matProp, e_kk, dtime, numslip, switch)
c      rhs_norm_0 = dsqrt(InnerProductVec(rhs, rhs, NVECS))
      rhs_norm_0 = NormOfVector(rhs, NVECS)
c
c------- initialize flags and start Newton iterations for stress
c
      iters = 0
      do while(.not.Converged(rhs, toler) .and. iters .lt. maxIters)
         iters = iters + 1
         call EqualTensors(tau, tau0, NVECS)
c
c---------- compute local jacobian
         call ComputeJacobian(lhs, rss, kappa, ppTHat, fCeiDevHat, 
     &                        matProp, dtime, numslip, switch)
c
c---------- solve for the increment of stress
         call lsolve(lhs, rhs, 1, ierr, NVECS)
         if (ierr .eq. XTAL_SING_JACOBIAN) then
            call WriteWarning(XTAL_O,
     &             'StressSolveDeviatoric: Jacobian is singular')
            return
         endif

         call EqualTensors(rhs, deltau, NVECS)
c
c---------- update stresses
         search = pone
         call AddTensors(pone, tau0, search, deltau, tau, NVECS)
c
c---------- compute new residual
         call ComputeResidual(rhs, tau, estran, estranStar, fCeiDevHat, 
     &                        fCeDevVolHat, rss, gamdot, kappa, pHatVec,
     &                        matProp, e_kk, dtime, numslip, switch)
c         rhs_norm = dsqrt(InnerProductVec(rhs, rhs, NVECS))
         rhs_norm = NormOfVector(rhs, NVECS)
c
c---------- simple line search
         do while ( rhs_norm .gt. rhs_norm_0 )

            search = search*0.5d0
            if (search .lt. TINY) then
               call WriteWarning(XTAL_O,
     &                'StressSolveDeviatoric: LS Failed, search < TINY')
               ierr = XTAL_LS_FAILED
               return
            endif
            call AddTensors(pone, tau0, search, deltau, tau, NVECS)

            call ComputeResidual(rhs, tau, estran, estranStar, 
     &                           fCeiDevHat, fCeDevVolHat, rss, gamdot,
     &                           kappa, pHatVec, matProp, e_kk, dtime, 
     &                           numslip, switch)
c            rhs_norm = dsqrt(InnerProductVec(rhs, rhs, NVECS))
            rhs_norm = NormOfVector(rhs, NVECS)

         enddo
c
c---------- update norm of residual
         rhs_norm_0 = rhs_norm

      enddo
c
c------- keep track of max number of newton iterations
c
      iterCounterN = max(iterCounterN, iters)
c
c------- check convergence on iters < maxIters
c
      if (iters .ge. maxIters) then
         call WriteWarning(XTAL_O,
     &                'StressSolveDeviatoric: iters > maxIters')
         ierr = XTAL_MAX_ITERS_HIT
         return
      endif
         
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE ComputeResidual(
     &   rhs, tau, estran, estranStar, fCeiDevHat, fCeDevVolHat, rss, 
     &   gamdot, kappa, pHatVec, matProp, e_kk, dtime, numslip, switch
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip
      real*8  e_kk, dtime
      real*8  rhs(NVECS), tau(NVECS), estran(NVECS), estranStar(NVECS)
      real*8  fCeiDevHat(NVECS,NVECS), fCeDevVolHat(NVECS)
      real*8  rss(MAX_SLIP), gamdot(MAX_SLIP), kappa(NKAPP)
      real*8  pHatVec(NVECS,MAX_SLIP), matProp(NPROPS,MAX_SLIP)

      integer is
      real*8  tmp(NVECS)
      real*8  InnerProductVec, SSKineticEqn

      real*8  crss
      real*8  overstress(MAX_SLIP)
      common /XtalOverS/ overstress
      
      real*8  switch(MAX_SLIP)
c
c---------------------------------------------------------------------72
c
c------- elastic strain: -[Cei]({tau}-ekk*{fCeDevVol}}+{e*} = -{e}+{e*}
c
      call AddTensors(pone, tau, -e_kk, fCeDevVolHat, tmp, NVECS)
      call MultAxu(fCeiDevHat, tmp, estran, NVECS)
      call AddTensors(-pone, estran, +pone, estranStar, rhs, NVECS)
c
c------- resolve shear stresses and shear strain rates
c
      do is = 1, numslip
         rss(is) = InnerProductVec(tau, pHatVec(1,is), NVECS)
c         crss = kappa(1) + overstress(is)
         crss = kappa(is)
         gamdot(is) = SSKineticEqn(rss(is), crss, matProp(1,is), 
     &                             kGAMDOT, switch(is))
      enddo
c
c------- residual
c
      do is = 1, numslip
         call AddScaledTensor(-dtime*gamdot(is), pHatVec(1,is), rhs, 
     &                        NVECS)
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE ComputeJacobian(
     &   lhs, rss, kappa, ppTHat, fCeiDevHat, matProp, dtime, numslip, 
     &   switch)
      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip
      real*8  dtime
      real*8  lhs(NVECS,NVECS), rss(MAX_SLIP), kappa(NKAPP)
      real*8  ppTHat(NVECS,NVECS,MAX_SLIP), fCeiDevHat(NVECS,NVECS)
      real*8  matProp(NPROPS,MAX_SLIP)

      integer is
      real*8  dGamdTau(MAX_SLIP)
      real*8  SSKineticEqn

      real*8  crss
      real*8  overstress(MAX_SLIP)
      common /XtalOverS/ overstress

      real*8  sGamPPt(NVECS,NVECS)
      common /SumGamPPt/ sGamPPt
      
      real*8  switch(MAX_SLIP)
c
c---------------------------------------------------------------------72
c
c------- d(gamdot)/d(tau)
c
      do is = 1, numslip
         crss = kappa(is)
c         crss = kappa(1) + overstress(is)
         dGamdTau(is) = SSKineticEqn(rss(is), crss, matProp(1,is),
     &                               kdGAMdTAU, switch(is))
      enddo
c
c------- jacobian
c
      call SetTensor(sGamPPt, pzero, NVECS*NVECS)
      do is = 1, numslip
         call AddScaledTensor(dtime*dGamdTau(is), ppTHat(1,1,is), 
     &                        sGamPPt, NVECS*NVECS)
      enddo
      call AddTensors(pone, fCeiDevHat, pone, sGamPPt, lhs, NVECS*NVECS)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE StressSolveVolumetric(
     &   statev, fCeDevVolHat, fCeVol, estran, e_kk, dtime
     &   )

      implicit none
      include 'params_xtal.inc'

      real*8  fCeVol, e_kk, dtime
      real*8  fCeDevVolHat(NVECS), estran(NVECS), statev(NSTAV)

      real*8  InnerProductVec
c
c---------------------------------------------------------------------72
c
      statev(kEVOL) = e_kk
      statev(kPRES) = InnerProductVec(fCeDevVolHat, estran, NVECS) 
     &                + fCeVol * e_kk
      
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetTensor(
     &   tensor, value, n
     &   )
      implicit none

      integer n
      real*8  value
      real*8  tensor(n)

      integer i
c
c---------------------------------------------------------------------72
c
      do i = 1, n
         tensor(i) = value
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SymmetrizeTensor(
     &   tensor, symmTensor, n
     &   )

      implicit none
      include 'params_xtal.inc'

      integer  n
      real*8   tensor(n*n), symmTensor(n*n)
c
c---------------------------------------------------------------------72
c
      if ((n .ne. 2) .and. (n .ne. 3))
     &   call RunTimeError(XTAL_O, 'SkewSymmetrizeTensor: n =! 2 or 3')

      if (n .eq. 2) then
         symmTensor(1) = tensor(1)
         symmTensor(4) = tensor(4)
         symmTensor(3) = 0.5*(tensor(3) + tensor(2))
         symmTensor(2) = symmTensor(3)
      else   ! n = 3
         symmTensor(1) = tensor(1)                     ! 11
         symmTensor(5) = tensor(5)                     ! 22
         symmTensor(9) = tensor(9)                     ! 33
         symmTensor(4) = 0.5*(tensor(4) + tensor(2))   ! 12
         symmTensor(7) = 0.5*(tensor(7) + tensor(3))   ! 13
         symmTensor(8) = 0.5*(tensor(8) + tensor(6))   ! 23
         symmTensor(2) = symmTensor(4)                 ! 21
         symmTensor(3) = symmTensor(7)                 ! 31
         symmTensor(6) = symmTensor(8)                 ! 32
      endif

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SkewSymmetrizeTensor(
     &   tensor, skewTensor, n
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer  n
      real*8   tensor(n*n), skewTensor(n*n)
c
c---------------------------------------------------------------------72
c
      if ((n .ne. 2) .and. (n .ne. 3))
     &   call RunTimeError(XTAL_O, 'SkewSymmetrizeTensor: n =! 2 or 3')

      call SetTensor(skewTensor, pzero, n*n)

      if (n .eq. 2) then
         skewTensor(3) = 0.5*(tensor(3) - tensor(2))
         skewTensor(2) = - skewTensor(3)
      else   ! n = 3
         skewTensor(4) = 0.5*(tensor(4) - tensor(2))
         skewTensor(7) = 0.5*(tensor(7) - tensor(3))
         skewTensor(8) = 0.5*(tensor(8) - tensor(6))
         skewTensor(2) = - skewTensor(4)
         skewTensor(3) = - skewTensor(7)
         skewTensor(6) = - skewTensor(8)
      endif

      return
      END
c

c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE Mat3x3ToVec5x1Symm(
     &   matrix, vector, n
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'
     
      integer  n
      real*8   matrix(n*n), vector(5)
c
c---------------------------------------------------------------------72
c
c------- matrix is deviatoric
c
      if (n .ne. 3)
     &   call RunTimeError(XTAL_O, 'Mat3x3ToVect5x1Symm: n =! 3')

      vector(1) = (matrix(1) - matrix(5)) / sqr2
      vector(2) = matrix(9) * sqr32
      vector(3) = matrix(2) * sqr2
      vector(4) = matrix(3) * sqr2
      vector(5) = matrix(6) * sqr2
 
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE Vec5x1ToMat3x3Symm(
     &   vector, matrix, n
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'
     
      integer n
      real*8  matrix(n*n), vector(5)
c
c---------------------------------------------------------------------72
c
c------- matrix is deviatoric
c
      if (n .ne. 3)
     &   call RunTimeError(XTAL_O, 'Vect5x1ToMat3x3Symm: n =! 3')

      matrix(1) = 0.5 * (sqr2 * vector(1) - sqr23 * vector(2))
      matrix(5) =-0.5 * (sqr2 * vector(1) + sqr23 * vector(2))
      matrix(9) = vector(2) * sqr23
      matrix(2) = vector(3) / sqr2
      matrix(3) = vector(4) / sqr2
      matrix(6) = vector(5) / sqr2
      matrix(4) = matrix(2)
      matrix(7) = matrix(3)
      matrix(8) = matrix(6)
 
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE Mat3x3ToVec3x1Skew(
     &   matrix, vector, n
     &   )

      implicit none
      include 'params_xtal.inc'
     
      integer n
      real*8  matrix(n*n), vector(3)
c
c---------------------------------------------------------------------72
c
      if (n .ne. 3)
     &   call RunTimeError(XTAL_O, 'Vect5x1ToMat3x3Skew: n =! 3')
c
c------- uses tensor components below main diagonal
c
      vector(1) = matrix(2)
      vector(2) = matrix(3)
      vector(3) = matrix(6)
 
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE Vec3x1ToMat3x3Skew(
     &   vector, matrix, n
     &   )

      implicit none
      include 'params_xtal.inc'
     
      integer  n
      real*8   matrix(n*n), vector(3)
c
c---------------------------------------------------------------------72
c
      if (n .ne. 3)
     &   call RunTimeError(XTAL_O, 'Vect5x1ToMat3x3Skew: n =! 3')

      matrix(1) = 0.
      matrix(5) = 0.
      matrix(9) = 0.
      matrix(2) = vector(1)
      matrix(3) = vector(2)
      matrix(6) = vector(3)
      matrix(4) = - vector(1)
      matrix(7) = - vector(2)
      matrix(8) = - vector(3)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE RotMat5x5ForSymm(
     &   cmatrix, qr5x5, n
     &   )

      implicit 	none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer  n
      real*8   cmatrix(n*n), qr5x5(5, 5)

      real*8   c11, c21, c31, c12, c22, c32, c13, c23, c33
c
c---------------------------------------------------------------------72
c
c----- construct 5x5 rotation matrix for symm 2nd order tensors: 
c         [A_sm]=[C][A_lat][C]' <=>  {A_sm} = [qr5x5]{A_lat}
c----- with: {A}={()/sqr2,sqr32*(),sqr2*(),sqr2*(),sqr2*()}

      if (n .ne. 3)
     &   call RunTimeError(XTAL_O, 'RotMat5x5ForSymm: n =! 3')

      c11 = cmatrix(1)
      c21 = cmatrix(2)
      c31 = cmatrix(3)
      c12 = cmatrix(4)
      c22 = cmatrix(5)
      c32 = cmatrix(6)
      c13 = cmatrix(7)
      c23 = cmatrix(8)
      c33 = cmatrix(9)

      qr5x5(1, 1)  =  0.5d0 * (c11 * c11 - c12 * c12 - 
     &                                           c21 * c21 + c22 * c22)
      qr5x5(1, 2)  =  sqr3 / 2.d0 * (c13 * c13 - c23 * c23)
      qr5x5(1, 3)  =  c11 * c12 - c21 * c22
      qr5x5(1, 4)  =  c11 * c13 - c21 * c23
      qr5x5(1, 5)  =  c12 * c13 - c22 * c23
      qr5x5(2, 1)  =  sqr3 / 2.d0 * (c31 * c31 - c32 * c32)
      qr5x5(2, 2)  =  1.5d0 * c33 * c33 - 0.5d0
      qr5x5(2, 3)  =  sqr3 * c31 * c32
      qr5x5(2, 4)  =  sqr3 * c31 * c33
      qr5x5(2, 5)  =  sqr3 * c32 * c33
      qr5x5(3, 1)  =  c11 * c21 - c12 * c22
      qr5x5(3, 2)  =  sqr3 * c13 * c23
      qr5x5(3, 3)  =  c11 * c22 + c12 * c21
      qr5x5(3, 4)  =  c11 * c23 + c13 * c21
      qr5x5(3, 5)  =  c12 * c23 + c13 * c22
      qr5x5(4, 1)  =  c11 * c31 - c12 * c32
      qr5x5(4, 2)  =  sqr3 * c13 * c33
      qr5x5(4, 3)  =  c11 * c32 + c12 * c31
      qr5x5(4, 4)  =  c11 * c33 + c13 * c31
      qr5x5(4, 5)  =  c12 * c33 + c13 * c32
      qr5x5(5, 1)  =  c21 * c31 - c22 * c32
      qr5x5(5, 2)  =  sqr3 * c23 * c33
      qr5x5(5, 3)  =  c21 * c32 + c22 * c31
      qr5x5(5, 4)  =  c21 * c33 + c23 * c31
      qr5x5(5, 5)  =  c22 * c33 + c23 * c32

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE RotMat3x3ForSkew(
     &   cmatrix, qr3x3, n
     &   )

      implicit  none
      include 'params_xtal.inc'

      integer  n
      real*8   cmatrix(n*n), qr3x3(3, 3)

      real*8   c11, c21, c31, c12, c22, c32, c13, c23, c33
c
c---------------------------------------------------------------------72
c
c---- Construct 3x3 rotation matrix for skew 2nd order tensors
c        [W_sm]=[C][W_lat][C]'  <=>  {W_sm} = [qr3x3]{W_lat}

      if (n .ne. 3)
     &   call RunTimeError(XTAL_O, 'RotMat3x3ForSkew: n =! 3')

      c11 = cmatrix(1)
      c21 = cmatrix(2)
      c31 = cmatrix(3)
      c12 = cmatrix(4)
      c22 = cmatrix(5)
      c32 = cmatrix(6)
      c13 = cmatrix(7)
      c23 = cmatrix(8)
      c33 = cmatrix(9)

      qr3x3(1, 1) = c22 * c11 - c21 * c12
      qr3x3(1, 2) = c23 * c11 - c21 * c13
      qr3x3(1, 3) = c23 * c12 - c22 * c13
      qr3x3(2, 1) = c32 * c11 - c31 * c12
      qr3x3(2, 2) = c33 * c11 - c31 * c13
      qr3x3(2, 3) = c33 * c12 - c32 * c13
      qr3x3(3, 1) = c32 * c21 - c31 * c22
      qr3x3(3, 2) = c33 * c21 - c31 * c23
      qr3x3(3, 3) = c33 * c22 - c32 * c23

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MultAxu(
     &   A, u, v, n
     &   )

      implicit none
      include 'numbers.inc'

      integer n
      real*8  A(n,n), u(n), v(n)

      integer i, j
c
c---------------------------------------------------------------------72
c
      call SetTensor(v, pzero, n)

      do i = 1, n
         do j = 1, n
            v(i) = v(i) + A(i,j) * u(j)
         end do
      end do

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MultATxu(
     &   A, u, v, n
     &   )

      implicit none
      include 'numbers.inc'

      integer n
      real*8  A(n,n), u(n), v(n)

      integer i, j
c
c---------------------------------------------------------------------72
c
      call SetTensor(v, pzero, n)

      do i = 1, n
         do j = 1, n
            v(i) = v(i) + A(j,i) * u(j)
         end do
      end do

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MultAxB(
     &   A, B, C, n
     &   )

      implicit none
      include 'numbers.inc'
     
      integer n
      real*8  A(n, n), B(n, n), C(n, n)

      integer i, j, k
c
c---------------------------------------------------------------------72
c
      call SetTensor(C, pzero, n*n)

      do i = 1, n 
         do j = 1, n
            do k = 1, n
               C(i,j) = C(i,j) + A(i,k) * B(k,j)
            enddo
         enddo
      enddo
   
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MultATxB(
     &   A, B, C, n
     &   )

      implicit none
      include 'numbers.inc'
     
      integer n
      real*8  A(n, n), B(n, n), C(n, n)

      integer i, j, k
c
c---------------------------------------------------------------------72
c
      call SetTensor(C, pzero, n*n)

      do i = 1, n 
         do j = 1, n
            do k = 1, n
               C(i,j) = C(i,j) + A(k,i) * B(k,j)
            enddo
         enddo
      enddo
   
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MultAxBT(
     &   A, B, C, n
     &   )

      implicit none
      include 'numbers.inc'
     
      integer n
      real*8  A(n, n), B(n, n), C(n, n)

      integer i, j, k
c
c---------------------------------------------------------------------72
c
      call SetTensor(C, pzero, n*n)

      do i = 1, n 
         do j = 1, n
            do k = 1, n
               C(i,j) = C(i,j) + A(i,k) * B(j,k)
            enddo
         enddo
      enddo
   
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE DeviatoricTensor(
     &   tensor, tensorDev, n
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer n
      real*8  tensor(n*n), tensorDev(n*n)

      real*8  TraceOfTensor, tensHyd
c
c---------------------------------------------------------------------72
c
      if (n .ne. 3)
     &   call RunTimeError(XTAL_O, 'DeviatoricTensor: n =! 3')

      tensHyd = TraceOfTensor(tensor, n) / pthree

      tensorDev(1) = tensor(1) - tensHyd
      tensorDev(5) = tensor(5) - tensHyd
      tensorDev(9) = tensor(9) - tensHyd

      tensorDev(2) = tensor(2)
      tensorDev(3) = tensor(3)
      tensorDev(6) = tensor(6)

      tensorDev(4) = tensor(4)
      tensorDev(7) = tensor(7)
      tensorDev(8) = tensor(8)
      
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      real*8 FUNCTION TraceOfTensor(
     &   tensor, n
     &   )

      implicit none
      include 'params_xtal.inc'

      integer n
      real*8  tensor(n*n)
c
c---------------------------------------------------------------------72
c
      if (n .ne. 3)
     &   call RunTimeError(XTAL_O, 'TraceOfTensor: n =! 3')

      TraceOfTensor = tensor(1) + tensor(5) + tensor(9)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      real*8 FUNCTION ScalarProduct(
     &   tensA, tensB, n
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer n
      real*8  tensA(n*n), tensB(n*n)
c
c---------------------------------------------------------------------72
c
      if (n .ne. 3)
     &   call RunTimeError(XTAL_O, 'ScalarProduct: n =! 3')

      ScalarProduct = tensA(1)*tensB(1) + tensA(2)*tensB(2) +
     &                tensA(3)*tensB(3) + tensA(4)*tensB(4) +
     &                tensA(5)*tensB(5) + tensA(6)*tensB(6) +
     &                tensA(7)*tensB(7) + tensA(8)*tensB(8) +
     &                tensA(9)*tensB(9)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE OuterProductVec(
     &   vecU, vecV, outer, n
     &   )

      implicit none
      include 'numbers.inc'

      integer n
      real*8  vecU(n), vecV(n)
      real*8  outer(n, n)
     
      integer i, j
c
c---------------------------------------------------------------------72
c
      do i = 1, n
         do j = 1, n
            outer(i,j) = vecU(i) * vecV(j)
         enddo
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      real*8 FUNCTION InnerProductVec(
     &   vecU, vecV, n
     &   )

      implicit none
      include 'numbers.inc'

      integer n
      real*8  vecU(n), vecV(n)

      integer i
c
c---------------------------------------------------------------------72
c
c---- compute inner product of two vectors: product = {u}^T*{v}
c
      InnerProductVec = pzero
      do i = 1, n
         InnerproductVec = InnerProductVec + vecU(i) * vecV(i)
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MultQAQT(
     &   Q, A, B, n
     &   )

      implicit none

      integer n
      real*8  Q(n,n), A(n,n), B(n,n)

      real*8  tmp(n,n)
c
c---------------------------------------------------------------------72
c
      call MultAxB(Q, A, tmp, n)       ! tmp = Q*A
      call MultAxBT(tmp, Q, B, n)      ! B = tmp*QT

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MultQTAQ(
     &   Q, A, B, n
     &   )

      implicit none

      integer n
      real*8  Q(n,n), A(n,n), B(n,n)

      real*8  tmp(n,n)
c
c---------------------------------------------------------------------72
c
      call MultATxB(Q, A, tmp, n)      ! tmp = QT*A
      call MultAxB(tmp, Q, B, n)       ! B = tmp*Q

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE EqualTensors(
     &   tensA, tensB, n
     &   )

      implicit none

      integer n
      real*8  tensA(n), tensB(n)

      integer i
c
c---------------------------------------------------------------------72
c
      do i = 1, n
         tensB(i) = tensA(i)
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE AddTensors(
     &   coef_A, tensA, coef_B, tensB, tensC, n
     &   )

      implicit none

      integer n
      real*8  coef_A, coef_B
      real*8  tensA(n), tensB(n), tensC(n)

      integer i
c
c---------------------------------------------------------------------72
c
       do i = 1, n
          tensC(i) = coef_A * tensA(i) + coef_B * tensB(i) 
       enddo

       return
       END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetToScaledTensor(
     &   fac_A, tensA, tensB, n
     &   )

      implicit none

      integer n
      real*8  fac_A
      real*8  tensA(n), tensB(n)

      integer i
c
c---------------------------------------------------------------------72
c
      do i = 1, n
          tensB(i) = fac_A * tensA(i)
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE AddScaledTensor(
     &   fac_A, tensA, tensB, n
     &   )

      implicit none

      integer n
      real*8  fac_A
      real*8  tensA(n), tensB(n)

      integer i
c
c---------------------------------------------------------------------72
c
      do i = 1, n
          tensB(i) = tensB(i) + fac_A * tensA(i)
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      real*8 FUNCTION DetOfMat3x3(
     &   tensor
     &   )

      implicit none
      
      real*8  tensor(3, 3)

      real*8  det11, det12, det13, det21, det22, det23
c
c---------------------------------------------------------------------72
c
c  Determinant of a second order tensor (3x3 matrix)
 
      det11 = tensor(1, 1) * tensor(2, 2) * tensor(3, 3)
      det12 = tensor(1, 2) * tensor(2, 3) * tensor(3, 1)
      det13 = tensor(1, 3) * tensor(2, 1) * tensor(3, 2)
      det21 = tensor(1, 1) * tensor(2, 3) * tensor(3, 2)
      det22 = tensor(1, 2) * tensor(2, 1) * tensor(3, 3)
      det23 = tensor(1, 3) * tensor(2, 2) * tensor(3, 1)
      
      DetOfMat3x3 = det11 + det12 + det13 - det21 - det22 - det23

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE InverseOfMat3x3(
     &   tensor, tensor_inv, tensor_det
     &   )

      implicit 	none

      real*8  tensor_det
      real*8  tensor(3, 3), tensor_inv(3, 3)

      real*8  tinv11, tinv12, tinv13
      real*8  tinv21, tinv22, tinv23
      real*8  tinv31, tinv32, tinv33
     
      real*8  DetOfMat3x3
c
c---------------------------------------------------------------------72
c
      tensor_det = DetOfMat3x3(tensor)

      tinv11 = tensor(2, 2) * tensor(3, 3) - tensor(2, 3) * tensor(3, 2)
      tinv12 = tensor(1, 3) * tensor(3, 2) - tensor(1, 2) * tensor(3, 3)
      tinv13 = tensor(1, 2) * tensor(2, 3) - tensor(1, 3) * tensor(2, 2)
      tinv21 = tensor(2, 3) * tensor(3, 1) - tensor(2, 1) * tensor(3, 3)
      tinv22 = tensor(1, 1) * tensor(3, 3) - tensor(1, 3) * tensor(3, 1)
      tinv23 = tensor(1, 3) * tensor(2, 1) - tensor(1, 1) * tensor(2, 3)
      tinv31 = tensor(2, 1) * tensor(3, 2) - tensor(2, 2) * tensor(3, 1)
      tinv32 = tensor(3, 1) * tensor(1, 2) - tensor(1, 1) * tensor(3, 2)
      tinv33 = tensor(1, 1) * tensor(2, 2) - tensor(1, 2) * tensor(2, 1)

      tensor_inv(1, 1) = tinv11 / tensor_det
      tensor_inv(1, 2) = tinv12 / tensor_det
      tensor_inv(1, 3) = tinv13 / tensor_det
      tensor_inv(2, 1) = tinv21 / tensor_det
      tensor_inv(2, 2) = tinv22 / tensor_det
      tensor_inv(2, 3) = tinv23 / tensor_det
      tensor_inv(3, 1) = tinv31 / tensor_det
      tensor_inv(3, 2) = tinv32 / tensor_det
      tensor_inv(3, 3) = tinv33 / tensor_det

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      real*8 FUNCTION SignOf(
     &   value
     &   )

      implicit none

      real*8  value
c
c---------------------------------------------------------------------72
c
c------- compute the sign of a number
c
      SignOf = 1.0
      if (value .lt. 0.0) SignOf = -1.0

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      real*8 FUNCTION Power(
     &   x, y
     &   )

      implicit none

      real*8  x, y
c
c---------------------------------------------------------------------72
c
c---- evaluates  x^y
c
      if (x .eq. 0.0) then
         if (y .gt. 0.0) then
            Power = 0.d0
         elseif (y .lt. 0.0) then
            Power = 1.d+300
         else
            Power = 1.d0
         endif
      else
         Power = y * log10(dabs(x))
         if (Power .gt. 300.0) then
            Power = 1.d+300
         else
            Power = 10.d0 ** Power
         endif
         if (x .lt. 0.0) Power = -Power
      endif

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE UnitVector(
     &   vector, unitV, n
     &   )
      implicit none
      include 'params_xtal.inc'

      integer n
      real*8  vector(n), unitV(n)

      real*8  magnitude, InnerProductVec
c
c---------------------------------------------------------------------72
c
      magnitude = dsqrt(InnerProductVec(vector, vector, n))
      if (magnitude .le. 0.d0) 
     &        call RunTimeError(XTAL_O, 'UnitVector: magnitude <= 0.0')
      call SetToScaledTensor(1./magnitude, vector, unitV, n)

      return
      END
c
c=====================================================================72
c
c=====================================================================72
c
      integer FUNCTION IndexMaxAbsValueOfVec(
     &   vector, n
     &   )

      implicit none

      integer n
      real*8  vector(n)

      integer indexVal, i
      real*8  maxValue
c
c---------------------------------------------------------------------72
c
      indexVal = 1
      maxValue = dabs(vector(1))
      do i = 2, n
         if (dabs(vector(i)) .gt. maxValue) then
            indexVal = i
            maxValue = dabs(vector(i))
         endif
      enddo

      IndexMaxAbsValueOfVec = indexVal

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      real*8 FUNCTION MaxAbsValueOfVec(
     &   vector, n
     &   )

      implicit none

      integer n
      real*8  vector(n)

      integer i
      real*8  maxValue
c
c---------------------------------------------------------------------72
c
      maxValue = dabs(vector(1))
      do i = 2, n
         if (dabs(vector(i)) .gt. maxValue) maxValue = dabs(vector(i))
      enddo

      MaxAbsValueOfVec = maxValue

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MatrixInverse(
     &   a, n, np, det
     &   )
c
c---- routine borrowed from Hammid Youssef
c---- input:  a   : nonsym matrix to invert
c----         n   : order of matrix to invert
c----         np  : dimension of matrix in calling routine
c---- output: a   : matrix inverse
c----         det : determinant of matrix
c---- local:  k   : VECTEUR DE TRAVAIL ENTIER 
c
      IMPLICIT DOUBLE  PRECISION (A-H,O-Z)
      DIMENSION A(NP,NP),K(6)
      DATA ZERO,UN,EPS/0.D0,1.D0,1.D-13/
c
c---------------------------------------------------------------------72
c
c-------  Initialization
c
      DET=UN
      DO I=1,N
         K(I)=I
      END DO
c
c-------  inversion de la matrix A
c
      DO II=1,N
c
c-------  search for non-zero pivot in column II
c
         DO I=II,N
            XPIV=A(I,II)
            IF(DABS(XPIV).GT.EPS) GO TO 10
         END DO
         DET=ZERO
         GOTO 1000
c
c-------  exchange rows II and I
c
 10      DET=DET*XPIV
         IF(I.EQ.II) GO TO 20
         I1=K(II)
         K(II)=K(I)
         K(I)=I1
         DO J=1,N
            C=A(I,J)
            A(I,J)=A(II,J)
            A(II,J)=C
         END DO
         DET=-DET
c
c-------  normalize the row of the pivot
c
 20      C=UN/XPIV
         A(II,II)=UN
         DO J=1,N
            A(II,J)=A(II,J)*C
         END DO
c
c-------  Elimination
c
         DO I=1,N
            IF(I.EQ.II) GO TO 30
            C=A(I,II)
            A(I,II)=ZERO
            DO J=1,N
               A(I,J)=A(I,J)-C*A(II,J)
            END DO
 30      END DO
      END DO
c
c-------  re-order columns of inverse
c
      DO J=1,N
c-------  search J1 until K(J1)=J
         DO J1=J,N
            JJ=K(J1)
            IF(JJ.EQ.J) GO TO 100
         END DO
100      IF(J.EQ.J1) GO TO 110
c-------  exchange columns J and J1
         K(J1)=K(J)
         DO I=1,N
            C=A(I,J)
            A(I,J)=A(I,J1)
            A(I,J1)=C
         END DO
110   END DO

1000  CONTINUE

      return
      END
c
c=====================================================================72
c
c     
c=====================================================================72
c     
      SUBROUTINE lsolve(a, x, nrhs, mystatus, n)
c     
c---- routine borrowed from Paul Dawson (Cornell University)
c---- Solve symmetric positive definite 5x5 linear systems. By calling 
c----  with the identity as right hand sides, you can use this routine
c----  to invert the matrix.
c
c     a    -- "n x n" matrix                                   (in/out)     
c     x    -- "n x nrhs" matrix, the array of right hand sides (in/out)
c     nrhs -- number of right hand sides                       (in)
c     mystatus -- return status                                (out)
c     n    -- size of system (must be equal to 5)
c     
      implicit none
      include 'params_xtal.inc'

      integer nrhs, mystatus, n
      real*8  a(n,n), x(n,*)
     
      integer i
      real*8  v1, v2, v3, v4, sum, not_zero, amax
c     
c---------------------------------------------------------------------72
c     
c------- check size of matrix
c
      if (n .ne. 5)
     &   call RunTimeError(XTAL_O, 'lsolve : n .ne. 5')
c     
c------- set the scale according to the largest entry of a.
c     
      amax = 0.0d0
      amax = max(a(1,1), a(2,1), a(3,1), a(4,1), a(5,1))
      amax = max(a(2,2), a(3,2), a(4,2), a(5,2), amax)
      amax = max(a(3,3), a(4,3), a(5,3), amax)
      amax = max(a(4,4), a(5,4), amax)
      amax = max(a(5,5), amax)
     
      not_zero = amax * TINY
c     
c------- A = LDL'.
c
c------- At each step, check the size of each pivot
c---------- j = 1.
      if (a(1,1) .LT. not_zero) then
        mystatus = XTAL_SING_JACOBIAN
        return
      endif
     
      sum = 1.0/a(1,1)

      a(2,1) = a(2,1) * sum
      a(3,1) = a(3,1) * sum
      a(4,1) = a(4,1) * sum
      a(5,1) = a(5,1) * sum
c
c---------- j = 2.
      v1 = a(2,1)*a(1,1)
      a(2,2) = a(2,2) - a(2,1)*v1
      
      a(3,2) = a(3,2)-a(3,1)*v1
      a(4,2) = a(4,2)-a(4,1)*v1
      a(5,2) = a(5,2)-a(5,1)*v1

      if (a(2,2) .LT. not_zero) then
        mystatus = XTAL_SING_JACOBIAN
        return
      endif

      sum = 1.0/a(2,2)

      a(3,2) = a(3,2) * sum
      a(4,2) = a(4,2) * sum
      a(5,2) = a(5,2) * sum
c
c---------- j = 3.
      v1 = a(3,1)*a(1,1)
      v2 = a(3,2)*a(2,2)
      sum = a(3,1)*v1 + a(3,2)*v2

      a(3,3) = a(3,3) - sum
      
      if (a(3,3) .LT. not_zero) then
        mystatus = XTAL_SING_JACOBIAN
        return
      endif

      sum = 1.0/a(3,3)

      a(4,3) = a(4,3)-a(4,1)*v1-a(4,2)*v2
      a(5,3) = a(5,3)-a(5,1)*v1-a(5,2)*v2
      a(4,3) = a(4,3) * sum
      a(5,3) = a(5,3) * sum
c
c---------- j = 4.
      v1 = a(4,1)*a(1,1)
      v2 = a(4,2)*a(2,2)
      v3 = a(4,3)*a(3,3)
      sum = a(4,1)*v1+a(4,2)*v2+a(4,3)*v3

      a(4,4) = a(4,4) - sum

      if (a(4,4) .LT. not_zero) then
        mystatus = XTAL_SING_JACOBIAN
        return
      endif

      a(5,4) = a(5,4)-a(5,1)*v1-a(5,2)*v2-a(5,3)*v3
      a(5,4) = a(5,4)/a(4,4)
c
c---------- j = 5.
      v1 = a(5,1)*a(1,1)
      v2 = a(5,2)*a(2,2)
      v3 = a(5,3)*a(3,3)
      v4 = a(5,4)*a(4,4)
      sum = a(5,1)*v1+a(5,2)*v2+a(5,3)*v3+a(5,4)*v4

      a(5,5) = a(5,5) - sum

      if (a(5,5) .LT. not_zero) then
        mystatus = XTAL_SING_JACOBIAN
        return
      endif
c     
c------- solve for RHS
c
      do i=1, nrhs
c
c---------- Ly=b. 
        x(2, i) = x(2, i)
     &     - a(2,1)*x(1, i)
        x(3, i) = x(3, i)
     &     - a(3,1)*x(1, i) - a(3,2)*x(2, i)
        x(4, i) = x(4, i)
     &     - a(4,1)*x(1, i) - a(4,2)*x(2, i) - a(4,3)*x(3, i)
        x(5, i) = x(5, i)
     &     - a(5,1)*x(1, i) - a(5,2)*x(2, i) - a(5,3)*x(3, i)
     &     - a(5,4)*x(4, i)
c
c---------- Dz=y.
        x(1, i) = x(1, i)/a(1,1)
        x(2, i) = x(2, i)/a(2,2)
        x(3, i) = x(3, i)/a(3,3)
        x(4, i) = x(4, i)/a(4,4)
        x(5, i) = x(5, i)/a(5,5)
c
c---------- L'x=z.
        x(4, i) = x(4, i)
     &     - a(5,4)*x(5, i)
        x(3, i) = x(3, i)
     &     - a(4,3)*x(4, i) - a(5,3)*x(5, i)
        x(2, i) = x(2, i)
     &     - a(3,2)*x(3, i) - a(4,2)*x(4, i) - a(5,2)*x(5, i)
        x(1, i) = x(1, i)
     &     - a(2,1)*x(2, i) - a(3,1)*x(3, i) - a(4,1)*x(4, i)
     &     - a(5,1)*x(5, i)

      enddo

      return
      END
c     
c=====================================================================72
c     
c
c=====================================================================72
c
      SUBROUTINE Vec5x1ToVec6x1(
     &   vec5x1, vec6x1
     &   )

      implicit none
      include 'numbers.inc'

      real*8  vec5x1(5), vec6x1(6)
c
c---------------------------------------------------------------------72
c
      vec6x1(1) = 0.5 * (sqr2*vec5x1(1) - sqr23*vec5x1(2)) 
      vec6x1(2) =-0.5 * (sqr2*vec5x1(1) + sqr23*vec5x1(2)) 
      vec6x1(3) = vec5x1(2) * sqr23
      vec6x1(4) = vec5x1(3) / sqr2
      vec6x1(5) = vec5x1(4) / sqr2
      vec6x1(6) = vec5x1(5) / sqr2

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MultAxB_G(
     &   A, B, C, m1, m2, m3
     &   )

      implicit none
      include 'numbers.inc'
     
      integer m1, m2, m3
      real*8  A(m1, m2), B(m2, m3), C(m1, m3)

      integer i, j, k
c
c---------------------------------------------------------------------72
c
c---- Dimensions: A is m1 x m2; B is m2 x m3; C is m1 x m3 
c
      call SetTensor(C, pzero, m1*m3)

      do i = 1, m1 
         do j = 1, m3
            do k = 1, m2
               C(i,j) = C(i,j) + A(i,k) * B(k,j)
            enddo
         enddo
      enddo
   
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MultAxu_G(
     &   A, u, v, m1, m2
     &   )

      implicit none
      include 'numbers.inc'

      integer m1, m2
      real*8  A(m1,m2), u(m2), v(m1)

      integer i, j
c
c---------------------------------------------------------------------72
c
c---- Dimensions: A is m1 x m2; u is m2; v is m1
c
      call SetTensor(v, pzero, m1)

      do i = 1, m1
         do j = 1, m2
            v(i) = v(i) + A(i,j) * u(j)
         end do
      end do

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MultATxu_G(
     &   A, u, v, m1, m2
     &   )

      implicit none
      include 'numbers.inc'

      integer m1, m2
      real*8  A(m1,m2), u(m1), v(m2)

      integer i, j
c
c---------------------------------------------------------------------72
c
c---- Dimensions: A is m2 x m1; u is m1; v is m2
c
      call SetTensor(v, pzero, m2)

      do i = 1, m2
         do j = 1, m1
            v(i) = v(i) + A(j,i) * u(j)
         end do
      end do

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE OuterProductVec_G(
     &   vecU, vecV, outer, m1, m2
     &   )

      implicit none
      include 'numbers.inc'

      integer m1, m2
      real*8  vecU(m1), vecV(m2)
      real*8  outer(m1, m2)
     
      integer i, j
c
c---------------------------------------------------------------------72
c
      do i = 1, m1
         do j = 1, m2
            outer(i,j) = vecU(i) * vecV(j)
         enddo
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      real*8 FUNCTION NormOfVector(
     &   vec,  n
     &   )

      implicit none
      include 'numbers.inc'

      integer n
      real*8  vec(n)

      real*8  max_vec
      real*8  MaxAbsValueOfVec, InnerProductVec
c
c---------------------------------------------------------------------72
c
      max_vec = MaxAbsValueOfVec(vec, n)
      if (max_vec .gt. pzero)
     &    call SetToScaledTensor(pone/max_vec, vec, vec, n)

      NormOfVector = dsqrt(InnerProductVec(vec, vec, n))

      if (max_vec .gt. pzero) then
         NormOfVector = max_vec*NormOfVector
         call SetToScaledTensor(max_vec, vec, vec, n)
      endif

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
c---- random number generator from numerical recipes
      FUNCTION ran1_nr(idum)
      IMPLICIT NONE
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1_nr,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     &    NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy/0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do j=NTAB+8,1,-1
           k=idum/IQ
           idum=IA*(idum-k*IQ)-IR*k
           if (idum.lt.0) idum=idum+IM
           if (j.le.NTAB) iv(j)=idum
        enddo
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1_nr=min(AM*iy,RNMX)
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE AnglesToRotMatrix(
     &   angle, crot, n
     &   )

      implicit none

      integer n
      real*8  angle(n)
      real*8  crot(n,n)

      real*8  sps, cps, sth, cth, sph, cph
c
c---------------------------------------------------------------------72
c
c------- Construct [C] matrix from euler angles (in radians)
c-------     {a}_sm = [C] {a}_xtal
c
      sps = dsin(angle(1))
      cps = dcos(angle(1))
      sth = dsin(angle(2))
      cth = dcos(angle(2))
      sph = dsin(angle(3))
      cph = dcos(angle(3))

      crot(1,1) = -sps * sph - cps * cph * cth
      crot(2,1) =  cps * sph - sps * cph * cth
      crot(3,1) =  cph * sth
      crot(1,2) =  cph * sps - sph * cps * cth
      crot(2,2) = -cps * cph - sps * sph * cth
      crot(3,2) =  sph * sth
      crot(1,3) =  cps * sth
      crot(2,3) =  sps * sth
      crot(3,3) =  cth

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE RotMatrixToAngles(
     &   crot, angle, n
     &   )

      implicit none

      integer n
      real*8  angle(n)
      real*8  crot(n,n)
      
      real*8  sth
c
c---------------------------------------------------------------------72
c
c------- compute euler angles from [C] (in radians)
c
      angle(2) = acos(crot(3,3))
      if (dabs(crot(3,3)) .ne. 1.0) then
         sth = dsin(angle(2))
         angle(1) = datan2(crot(2,3)/sth, crot(1,3)/sth)
         angle(3) = datan2(crot(3,2)/sth, crot(3,1)/sth)
      else
         angle(1) = 0.
         angle(3) = datan2(-crot(1,2), -crot(1,1))
      endif

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      real*8 FUNCTION SSKineticEqn(
     &   rss, crss, matProp, kflag, rswitch
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer kflag
      real*8  crss
      real*8  matProp(NPROPS)

      real*8  xm, gam0, arg, pow, rss
      real*8  CheckArgPowLaw, SignOf, Power
      
      real*8  rswitch, C
c
c---------------------------------------------------------------------72
c
c------- recover power law material constants
c
      xm   = matProp(3) 
      gam0 = matProp(4)
      C = (1./crss)**(1./xm)
c
c------- check argument of power law
c
      arg = rss/crss
      if (dabs(arg) .ne. 0) arg = CheckArgPowLaw(arg, SignOf(arg))
c
c------- evaluate slip system kinetic equation 
c
      pow = Power(dabs(arg), 1./xm-1.)
      if (kflag .eq. kGAMDOT) then
         if (rswitch .eq. pone) then
             SSKineticEqn = gam0 * arg * pow
         else if (rswitch .eq. pzero) then
             SSKineticEqn = gam0 * C * arg * pow
         else
             call RunTimeError(XTAL_O, 'Unknown switch status')
         endif
      else if (kflag .eq. kdGAMdTAU) then
         if (rswitch .eq. pone) then
             SSKineticEqn = gam0 / (xm*crss) * pow
         else if (rswitch .eq. pzero) then
             SSKineticEqn = gam0 / (xm*crss) * C * pow
         else
             call RunTimeError(XTAL_O, 'Unknown switch status')
         endif
      else if (kflag .eq. kdGAMdKAPP) then
         if (rswitch .eq. pone) then
             SSKineticEqn = -gam0 / (xm*crss) * arg * pow
         else if (rswitch .eq. pzero) then
             SSKineticEqn = -gam0 / (xm*crss) * C * arg * pow
         else
             call RunTimeError(XTAL_O, 'Unknown switch status')
         endif
      else
         call RunTimeError(XTAL_O, 'SSKineticEqn: unknown kflag')
      endif

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      real*8 FUNCTION SSKineticEqnDrag(
     &   rss, crss, matProp, kflag
     &   )

      implicit none
      include 'params_xtal.inc'

      integer kflag
      real*8  rss, crss
      real*8  matProp(NPROPS)

      real*8  xm, gam0, arg, pow, bdrag, ratio_drag, tmp1, tmp2
      real*8  gamdot_1, gamdot_2, gamdot_e, dgamdot_1, dgamdot_2
      real*8  CheckArgPowLaw, SignOf, Power
c
c---------------------------------------------------------------------72
c
c------- recover power law material constants
c
      xm    = matProp(3) 
      gam0  = matProp(4)
      bdrag = matProp(11)
      ratio_drag = matProp(12)

      arg = rss/crss
      if (dabs(arg) .ge. ratio_drag) then
c
c------- flow rule based on dislocation drag kinetics
c
c---------- evaluate slip system kinetic equation 
         if (kflag .eq. kGAMDOT) then 
            SSKineticEqnDrag = arg / bdrag
         else if (kflag .eq. kdGAMdTAU) then
            SSKineticEqnDrag = 1.0 / (bdrag*crss)
         else if (kflag .eq. kdGAMdKAPP) then
            SSKineticEqnDrag = -1.0 / (bdrag*crss) * arg
         else
            call RunTimeError(XTAL_O, 'SSKineticEqnDrag: unkown kflag')
         endif

      else
c
c------- flow rule based on power law & dislocation drag kinetics
c
c---------- check argument of power law
         if (dabs(arg) .ne. 0) arg = CheckArgPowLaw(arg, SignOf(arg))
c
c---------- evaluate slip system kinetic equation 
         pow = Power(dabs(arg), 1./xm-1.)
         gamdot_1 = gam0 * arg * pow
         gamdot_2 = arg / bdrag
         gamdot_e =  gamdot_1*gamdot_2 / (gamdot_1 + gamdot_2)
         if (kflag .eq. kGAMDOT) then 
            SSKineticEqnDrag = gamdot_e
         else if (kflag .eq. kdGAMdTAU) then
            dgamdot_1 = gam0 / (xm*crss) * pow
            dgamdot_2 = 1.0 / (bdrag*crss)
            tmp1 = dgamdot_1/gamdot_1*gamdot_e/gamdot_1
            tmp2 = dgamdot_2/gamdot_2*gamdot_e/gamdot_2
            SSKineticEqnDrag = gamdot_e * (tmp1 + tmp2)
         else if (kflag .eq. kdGAMdKAPP) then
            dgamdot_1 = -gam0 / (xm*crss) * arg * pow
            dgamdot_2 = -1.0 / (bdrag*crss) * arg
            tmp1 = dgamdot_1/gamdot_1*gamdot_e/gamdot_1
            tmp2 = dgamdot_2/gamdot_2*gamdot_e/gamdot_2
            SSKineticEqnDrag = gamdot_e * (tmp1 + tmp2)
         else
            call RunTimeError(XTAL_O, 'SSKineticEqnDrag: unkown kflag')
         endif
      endif

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE BoundForArgPowLaw(
     &   xm
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      real*8  xm

      real*8  xmm

      real*8  argMin, argMax
      common /POWERLAWARG/ argMin, argMax

      real*8  EXPON
      data EXPON /280.d0/
c
c---------------------------------------------------------------------72
c
c---- limits (bounds) on the value of the argument for power law
c---- note: In general: xm <= 1.0 (1/xm >= 1.0)
c
      xmm = pone/xm - pone
      if (dabs(xm-pone) .lt. TINY) xmm = pone

      argMin = dexp(-EXPON/(xmm)*dlog(10.d0))
      argMax = dexp( EXPON/(xmm)*dlog(10.d0))

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      real*8 FUNCTION CheckArgPowLaw(
     &   arg, sgn
     &   )
 
      implicit none
      
      real*8  arg, sgn

      real*8  argMin, argMax                 ! values assigned in
      common /POWERLAWARG/ argMin, argMax    ! BoundForArgPowerLaw
c
c---------------------------------------------------------------------72
c
c---- check range of argument for power law
c
      if (dabs(arg) .lt. argMin) then
         CheckArgPowLaw = argMin * sgn
      else if (dabs(arg) .ge. argMax) then
         CheckArgPowLaw = argMax * sgn
      else
         CheckArgPowLaw = arg
      endif

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE LimitRatioForDrag(
     &   matProp
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      real*8  matProp(NPROPS)

      real*8  xm, gam0, bdrag, arg, argMax

      real*8  FACTOR
      data FACTOR /1.0d+10/

      real*8  EXPON
      data EXPON /280.d0/
c
c---------------------------------------------------------------------72
c
c------ (rss/crss)_limit to switch from (powerLaw & drag) to pure drag
c
      xm    = matProp(3)
      gam0  = matProp(4)
      bdrag = matProp(11)

      if (dabs(xm-pone) .lt. TINY) then
         matProp(12)=1.0d+300
         return
      endif

      arg = FACTOR/(bdrag*gam0)

      if (xm .gt. 0.700) then
         argMax = dexp(EXPON*(pone/xm-pone)*dlog(10.0d0))
         if (arg .gt. argMax) then
            matProp(12)=1.0d+300
            return
         endif
      endif

      matProp(12) = dexp(xm/(pone-xm)*dlog(arg))

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE RunTimeError(
     &   io, message
     &   )

      implicit none
      
      character*(*) message
      integer io
c
c---------------------------------------------------------------------72
c
      write(io, 1000) message
      write(*, 1000) message
c      call closfl( )
      call CrystalCloseIOFiles( )
      call CrystalCloseIOFiles_2( )
      stop

1000  format(/,'***ERROR Message: '/, 3x, a)

      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE WriteWarning(
     &   io, message
     &   )

      implicit none
      
      character*(*) message
      integer io
c
c---------------------------------------------------------------------72
c
      write(io, 1000) message

1000  format(/,'***WARNING Message: '/, 3x, a)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE WriteMessage(
     &   io, message
     &   )

      implicit none
      
      character*(*) message
      integer io
c
c---------------------------------------------------------------------72
c
      write(io, 1000) message

1000  format('***Message: ', a)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetUpCrystalProps2(
     &   props, mprops, numqpt, numel
     &   )
     
      implicit none
      include 'params_xtal.inc'

      integer numel, numqpt, mprops
      real*8  props(mprops)

      character*160 filePath, fileRoot
      common /XtalWorkDir/ filePath  /XtalRootFile/ fileRoot
c
      integer numgrn, numslip, numvtx, kODFout
      real*8  kappa0(NKAPP), fCeDevVol(NVECS), fCeVol
      real*8  matProp(NPROPS, MAX_SLIP), tauSlip(MAX_SLIP)
      real*8  fCeDev(NVECS, NVECS), fCeiDev(NVECS, NVECS)
      real*8  zBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0(DIMS, DIMS, MAX_SLIP), qBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0Vec(NVECS, MAX_SLIP), qBar0Vec(DIMS, MAX_SLIP)
      real*8  ppTBar0(NVECS, NVECS, MAX_SLIP)
      real*8  gcrot0(DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  sigfs(5, MAX_VTX)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)

      common /XtalNumGrn/ numgrn, numslip
      common /XtalNumVtx/ numvtx
      common /XtalODFOut/ kODFout
      common /XtalState0/ kappa0
      common /XtalProps / fCeDev, fCeiDev, matProp,
     &                    fCeDevVol, fCeVol, hardmtx
      common /XtalSlipG / tauSlip, zBar0, pBar0, qBar0,
     &                    pBar0Vec, qBar0Vec, ppTBar0
      common /XtalCrot0 / gcrot0
      common /XtalStrVtx/ sigfs

      real*8  overstress(MAX_SLIP)
      common /XtalOverS/ overstress

      integer maxIterState, maxIterNewt
      real*8  tolerState, tolerNewt

      common /nlcsolve1/ maxIterState, maxIterNewt
      common /nlcsolve2/ tolerState, tolerNewt
c
      real*8  gstress    (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstress_n  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran    (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran_n  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa     (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa_n   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev    (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev_n  (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geqvalues  (NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot    (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot      (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot_n    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot      (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot_n    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps       (9, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps_n     (9, MAX_GRN, NUMQPT_T, NUMEL_T)

      common  /XtalVars  / gstress, gestran, gkappa, gstatev, geqvalues,
     &                     gcrot, grrot, ggamdot, geps
      common  /XtalVars_n/ gstress_n, gestran_n, gkappa_n, gstatev_n,
     &                     gcrot_n, grrot_n, geps_n
      integer kODF
      common /XtalODFCode/ kODF

      real*8  angles(DIMS, MAX_ORIEN)
      real*8  euler(DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      common /EulerAngles/ angles, euler

      integer numor, seed
      common  /XtalOrients/ numor, seed
      
      real*8  gforden    (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gforden_n  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden    (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden_n  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gsub       (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gsub_n     (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gdensity (1, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gdensity_n (1, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gfdensity (1, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gfdensity_n (1, MAX_GRN, NUMQPT_T, NUMEL_T)
      
      common  /Dislocation/  gforden, gforden_n, ggndden, ggndden_n,
     &                       gsub, gsub_n,
     &                       gdensity, gdensity_n, gfdensity, 
     &                       gfdensity_n
c
c---------------------------------------------------------------------72
c
      call CrystalInitialize2(props, mprops, numel, numqpt,
     &                       filePath, fileRoot,
     &                       XTAL_I, XTAL_E, XTAL_O, XTAL_G1, XTAL_G2,
     &                       XTAL_Z1, XTAL_Z2, XTAL_STRESS_O,
     &                       XTAL_STRAIN_O, XTAL_EFFSS_O, XTAL_TRUESS_O,
     &                       XTAL_ITER_O, AGG_EFFSS_O,
     &                       numgrn, numslip, numvtx, kODFout,
     &                       kappa0, fCeDevVol, fCeVol,
     &                       matProp, tauSlip,
     &                       fCeDev, fCeiDev,
     &                       zBar0,
     &                       pBar0, qBar0,
     &                       pBar0Vec, qBar0Vec,
     &                       ppTBar0,
     &                       gcrot0,
     &                       sigfs,
     &                       overstress,
     &                       hardmtx,
     &                       maxIterState, maxIterNewt,
     &                       tolerState, tolerNewt,
     &                       kODF, numor, seed, angles, euler,
     &                       XTAL_TXT_OUT,
     &                       gstress,
     &                       gstress_n,
     &                       gestran,
     &                       gestran_n,
     &                       gkappa,
     &                       gkappa_n,
     &                       gstatev,
     &                       gstatev_n,
     &                       geqvalues,
     &                       ggamdot,
     &                       gcrot,
     &                       gcrot_n,
     &                       grrot,
     &                       grrot_n,
     &                       gforden, gforden_n, geps, geps_n, ggndden,
     &                       ggndden_n, gsub, gsub_n, gdensity,
     &                       gdensity_n, gfdensity, gfdensity_n)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetUpInitialState2(
     &   numqpt, numel, noel, npt, nstatv, statev
     &   )

      implicit none
      include  'params_xtal.inc'

      integer numel, numqpt, noel, npt
      integer nstatv
      real*8  statev(nstatv)

      integer numgrn, numslip, kODF
      integer numor, seed
      real*8  gcrot0(DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  angles(DIMS, MAX_ORIEN)
      real*8  euler(DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)

      common /XtalNumGrn/  numgrn, numslip
      common /XtalODFCode/ kODF
      common /XtalOrients/ numor, seed
      common /XtalCrot0 /  gcrot0
      common /EulerAngles/ angles, euler
c
      real*8  kappa0(NKAPP)
      common  /XtalState0/ kappa0
c
c---------------------------------------------------------------------72
c
      call AssignCrystalODF2(numel, numqpt, noel, npt,
     &                      numgrn, numslip, kODF,
     &                      numor, seed,
     &                      gcrot0,
     &                      angles,
     &                      euler, 
     &                      XTAL_O, XTAL_TXT_OUT)

      call SetUpStateVars2(nstatv, statev,
     &                    numgrn, numslip,
     &                    kappa0,
     &                    gcrot0,
     &                    XTAL_O)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetUpStateVars2(
     &   nstatv, statev,
     &   numgrn, numslip,
     &   kappa0,
     &   gcrot0,
     &   FILE_O
     &   ) 

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer nstatv
      real*8  statev(nstatv)

      integer numgrn, numslip, FILE_O
      real*8  kappa0(NKAPP)
      real*8  gcrot0 (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)

      integer dex, ig, varsPerGrn
c
c---------------------------------------------------------------------72
c
c---- initialize state variable vector per ip
c
      varsPerGrn = 2*NVECS          ! stress, estran
c     &             + NKAPP          ! kappa
     &             + numslip        ! kappa
     &             + NSTAV          ! ee_v, p, detVe, WpNorm, ssact
     &             + NEQVA/2        ! eqps, eqstr, gam_star, gamtot 
     &             + 3*DIMS*DIMS    ! C_0, C, R
     &             + numslip        ! gamdot
     &             + 4*numslip    
     &             + 9
     &             + 1            
     &             + 1 

      if (nstatv .ne. numgrn*varsPerGrn)
     &   call RunTimeError2(FILE_O, 'nstatv .ne. numgrn*varsPerGrn')

      do ig = 1, numgrn

         dex = (ig-1)*varsPerGrn + 1           ! stress (xtal)
         call SetTensor2(statev(dex), pzero, NVECS)

         dex = dex + NVECS                      ! estrain
         call SetTensor2(statev(dex), pzero, NVECS)

         dex = dex + NVECS                      ! kappa
c         call SetTensor(statev(dex), kappa0, NKAPP)
         call EqualTensors2(kappa0, statev(dex), numslip)

c         dex = dex + NKAPP              ! ee_v, p, detVe, WpNorm, ssact
         dex = dex + numslip              ! ee_v, p, detVe, WpNorm, ssact
         call SetTensor2(statev(dex), pzero, NSTAV)

         dex = dex + NSTAV               ! eqps, eqstr, gam_star, gamtot
         call SetTensor2(statev(dex), pzero, NEQVA/2)

         dex = dex + NEQVA/2                    ! C_0
         call EqualTensors2(gcrot0(1,1,ig,1,1), statev(dex), DIMS*DIMS)

         dex = dex + DIMS*DIMS                  ! C
         call EqualTensors2(gcrot0(1,1,ig,1,1), statev(dex), DIMS*DIMS)

         dex = dex + DIMS*DIMS                  ! R
         call EqualTensors2(Ident2nd, statev(dex), DIMS*DIMS)

         dex = dex + DIMS*DIMS                  ! gamdot
         call SetTensor2(statev(dex), pzero, numslip)
         
         dex = dex + numslip                    ! forden   SDV 71
         call SetTensor2(statev(dex), 1.d12, numslip)           ! changed the value from 5.d12
         
         dex = dex + numslip                    ! gndden   SDV 83
         call SetTensor2(statev(dex), pzero, numslip)
         
         dex = dex + numslip                    ! sub
         call SetTensor2(statev(dex), 5.d1, numslip)         
         
         dex = dex + 2 * numslip                    ! eps
         call SetTensor2(statev(dex), pzero, 9)
         
         dex = dex + 9                    ! Average forest density   SDV 128 means average for 1 element
         call SetTensor2(statev(dex), pzero, 1)
         
         dex = dex + 1                    ! Average GND density   SDV 129
         call SetTensor2(statev(dex), pzero, 1)      
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE EvolveMatlState2(
     &   stress, statev, ddsdde, strain, dstrain, time, 
     &   dtime, temp, dtemp, ! ndi, nshr, ntens, 
     &   nstatv, drot, dfgrd0, dfgrd1, noel, npt,
     &   incr, numel, numqpt, iprint, pnewdt, coords
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer nstatv  !, ndi, nshr, ntens
      integer noel, npt, incr, numel, numqpt, iprint
      real*8  dtime, temp, dtemp, pnewdt
      real*8  stress(ntens), statev(nstatv), ddsdde(ntens, ntens)
      real*8  strain(ntens), dstrain(ntens), time(2)
      real*8  dfgrd0(3,3), dfgrd1(3,3), drot(3,3), coords(3)

      integer statusFlag, numIncrs
      integer numgrn, numslip
      real*8  gstress   (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran   (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa    (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev   (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geqvalues (NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot     (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot     (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot   (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)

      real*8  gstress_n (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran_n (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa_n  (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev_n (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot_n   (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot_n   (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot0    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps      (9, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps_n    (9, MAX_GRN, NUMQPT_T, NUMEL_T)

      real*8  vgrad(DIMS, DIMS)
      real*8  savg_ij(DIMS, DIMS, NUMQPT_T) 
      real*8  cavg_ijkl(DIMS2, DIMS2, NUMQPT_T)
      real*8  epsxx(NUMQPT_T), epsyy(NUMQPT_T), epsth(NUMQPT_T)
      real*8  epsxy(NUMQPT_T), epsum(NUMQPT_T), spinn(NUMQPT_T)
      real*8  qptpo(NUMQPT_T), qptp(NUMQPT_T)
      real*8  sigxx(NUMQPT_T), sigyy(NUMQPT_T)
      real*8  sigth(NUMQPT_T), sigxy(NUMQPT_T)

      real*8  epsxz(NUMQPT_T),  epsyz(NUMQPT_T)
      real*8  spinxz(NUMQPT_T), spinyz(NUMQPT_T)
      real*8  sigxz(NUMQPT_T),  sigyz(NUMQPT_T)

      integer ielem

      common  /XtalNumGrn/ numgrn, numslip
      common  /XtalCrot0 / gcrot0
      common  /XtalVars  / gstress, gestran, gkappa, gstatev, geqvalues,
     &                     gcrot, grrot, ggamdot, geps
      common  /XtalVars_n/ gstress_n, gestran_n, gkappa_n, gstatev_n,
     &                     gcrot_n, grrot_n, geps_n
      common  /for3DProbs/ epsxz, epsyz, spinxz, spinyz, sigxz, sigyz
      
      real*8  gforden   (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gforden_n (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden   (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden_n (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gsub      (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gsub_n    (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gdensity (1, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gdensity_n (1, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gfdensity (1, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gfdensity_n (1, MAX_GRN, NUMQPT_T, NUMEL_T)
      common  /Dislocation/  gforden, gforden_n, ggndden, ggndden_n,
     &                       gsub, gsub_n,
     &                       gdensity, gdensity_n, gfdensity, 
     &                       gfdensity_n 
c
c---------------------------------------------------------------------72
c
c-------- kinematic quantities - deformation path
c
      call VelocityGrad2(dstrain, drot, dfgrd0, dfgrd1, vgrad, dtime, 
     &                  noel, npt, incr)
      
      call DeformPath2(epsxx, epsyy, epsth, epsxy, epsum, spinn,
     &      qptpo, qptp, vgrad, temp, (temp+dtemp), numel, numqpt)
c
c-------- fetch state variables from abaqus vector array
c
      call RecoverStateVars2(statev, nstatv, gstress_n, gestran_n,
     &      gkappa_n, gstatev_n, gcrot_n, grrot_n, gcrot0, ggamdot, 
     &      geqvalues, numgrn, numslip, gforden_n, geps_n, ggndden_n,
     &      gsub_n, gdensity_n, gfdensity_n)
c
c-------- compute state
c
      numIncrs = 1000
      ielem    = 1
      statusFlag = XTAL_CONVERGED
      call StressCrystal2(sigxx, sigyy, sigth, sigxy, epsxx, epsyy,
     &      epsth, epsxy, epsum, spinn, qptpo, qptp, dtime, time(2),
     &      ielem, incr, numel, numqpt, iprint, savg_ij, cavg_ijkl, 
     &      statusFlag, numIncrs, noel, npt, dfgrd1, coords)
      if (statusFlag .ne. XTAL_CONVERGED) then
         write(XTAL_O, *) ' ** Umat did not converged       **'
         write(XTAL_O, *) ' ** re-scaling time step by 0.75 **'
         pnewdt = 0.75
         return
      endif
c
c-------- save state variables in abaqus vector array
c
      call SaveStateVars2(statev, nstatv, gstress, gestran, gkappa, 
     &      gstatev, gcrot, grrot, gcrot0, ggamdot, geqvalues, numgrn,
     &      numslip, gforden, geps, ggndden, gsub, gdensity, gfdensity)
c
c-------- stresses and algorithmic moduli in abaqus format
c
      call SaveStressModuli2(stress, ddsdde, savg_ij(1, 1, 1), 
     &                      cavg_ijkl(1, 1, 1))

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE VelocityGrad2(
     &   dstran, drot_aba, dfgrd0, dfgrd1, vgrad, dtime, ielem, iqpt,
     &   incr
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

c      integer ntens, nshr
      integer ielem, iqpt, incr
      real*8  dtime
      real*8  dstran(NTENS)
      real*8  drot_aba(3,3), dfgrd0(3,3), dfgrd1(3,3), vgrad(3,3)

      integer i, j
      real*8  det
      real*8  d_aba(3,3)
      real*8  d_umat(3,3), w_umat(3,3), drot_umat(3,3), spin(3)
      real*8  matx_1(3,3), matx_2(3,3), matx_3(3,3), rdfgrd(3,3)
      
      common  /MacroStrainRate/ d_aba
c
c---------------------------------------------------------------------72
c
c-------- relative deformation gradient: f=F_(n+1)*F_n^{-1}
c
      call InverseOfMat3x32(dfgrd0, matx_1, det)
      call MultAxB2(dfgrd1, matx_1, rdfgrd, 3)
c
c-------- approximation to velocity gradient: L=2/dt*(f-I)*(f+I)^{-1}
c
      call AddTensors2(pone, rdfgrd, -pone, Ident2nd, matx_2, 3*3)
      call AddTensors2(pone, rdfgrd, pone, Ident2nd, matx_1, 3*3)
      call InverseOfMat3x32(matx_1, matx_3, det)

      call MultAxB2(matx_2, matx_3, vgrad, 3)
      call SetToScaledtensor2(2.0/dtime, vgrad, vgrad, 3*3)
c
c-------- rate of deformation and spin: d_umat=sym(L), w_umat=skw(L)
c
      call SetTensor2(d_umat, pzero, DIMS*DIMS)
      call SetTensor2(w_umat, pzero, DIMS*DIMS)

      call SymmetrizeTensor2(vgrad, d_umat, 3)
      call SkewSymmetrizeTensor2(vgrad, w_umat, 3)
c
c-------- incremental rotation tensor: drot_umat = exp(dt*w_umat)
c
      call Mat3x3ToVec3x1Skew2(w_umat, spin, 3)
      call IncrementalRotTensor2(drot_umat, spin, dtime)
c
c-------- rate of deformation from ABAQUS
c
      call SetTensor2(d_aba, pzero, DIMS*DIMS)
      d_aba(1,1) = dstran(1)/dtime
      d_aba(2,2) = dstran(2)/dtime
      d_aba(3,3) = dstran(3)/dtime
      d_aba(1,2) = dstran(4)/ptwo/dtime
      d_aba(2,1) = d_aba(1,2)

      if (NSHR. gt. 1) then
         d_aba(1,3) = dstran(5)/ptwo/dtime
         d_aba(2,3) = dstran(6)/ptwo/dtime
         d_aba(3,1) = d_aba(1,3)
         d_aba(3,2) = d_aba(2,3)
      endif
c
c-------- print values for comparison
c
      if (ielem .eq. kPRINT_ELEM .and. iqpt .eq. kPRINT_QPT) then
         if (incr .eq. 1) write(XTAL_E, 1000) ielem, iqpt
         write(XTAL_E, '(/i5)') incr
         do i = 1, 3
            write(XTAL_E, 5000) (vgrad(i,j), j=1,3)
         enddo
         write(XTAL_E, *)
         do i = 1, 3
            write(XTAL_E, 5000) (d_aba(i,j), j=1,3),
     &                          (d_umat(i,j),j=1,3), 
     &                          (w_umat(i,j),j=1,3)
         enddo
         write(XTAL_E, *)
         do i = 1, 3
            write(XTAL_E, 5000) (drot_aba(i,j), j=1,3),
     &                          (drot_umat(i,j),j=1,3)
         enddo
      endif

1000  format(' INCR',/,18x,' L_ij',/,
     &       18x,' D_aba',32x,'D_umat',30x,'W_umat',/,
     &       18x,' Drot_aba',32x,'Drot_umat',/,
     &       38x,' (elem # ', i5, ',  qpt # ', i2, ')')
5000  format((3(2x,3(1x,e11.4))))
c
c-------- recompute velocity gradient: L = d_aba + w_umat
c
      call AddTensors2(pone, d_aba, pone, w_umat, vgrad, DIMS*DIMS)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE DeformPath2(
     &   epsxx, epsyy, epsth, epsxy, epsum, spinn, qptpo, qptp, vgrad, 
     &   theta_n, theta, numel, numqpt
     &   )

      implicit none
      include 'params_xtal.inc'

      integer numel, numqpt
      real*8  theta_n, theta
      real*8  vgrad(DIMS, DIMS)
      real*8  epsxx(numqpt), epsyy(numqpt), epsth(numqpt)
      real*8  epsxy(numqpt), epsum(numqpt), spinn(numqpt)
      real*8  qptpo(numqpt), qptp(numqpt)

      real*8  epsxz(NUMQPT_T),  epsyz(NUMQPT_T)
      real*8  spinxz(NUMQPT_T), spinyz(NUMQPT_T)
      real*8  sigxz(NUMQPT_T),  sigyz(NUMQPT_T)
      common /for3DProbs/ epsxz, epsyz, spinxz, spinyz, sigxz, sigyz
c
c---------------------------------------------------------------------72
c
c------- deformation rate: deviatoric and volumetric 
c
      epsxx(1) = vgrad(1, 1)
      epsyy(1) = vgrad(2, 2)
      epsth(1) = vgrad(3, 3)

      epsum(1) = epsxx(1) + epsyy(1) + epsth(1)

      epsxx(1) = epsxx(1) - epsum(1) / 3.0d0
      epsyy(1) = epsyy(1) - epsum(1) / 3.0d0
      epsth(1) = epsth(1) - epsum(1) / 3.0d0

      epsxy(1) = (vgrad(2, 1) + vgrad(1, 2))/2.d0
      epsxz(1) = (vgrad(3, 1) + vgrad(1, 3))/2.d0
      epsyz(1) = (vgrad(3, 2) + vgrad(2, 3))/2.d0
c
c------- axial vector of spin
c
      spinn(1)  = (vgrad(1, 2) - vgrad(2, 1))/2.d0
      spinxz(1) = (vgrad(1, 3) - vgrad(3, 1))/2.d0
      spinyz(1) = (vgrad(2, 3) - vgrad(3, 2))/2.d0
c
c------- temperature (K)
c
      qptpo(1) = theta_n
      qptp(1)  = theta

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE RecoverStateVars2(
     &   statev, nstatv, gstress_n, gestran_n, gkappa_n, gstatev_n, 
     &   gcrot_n, grrot_n, gcrot0, ggamdot, geqvalues, numgrn, numslip,
     &   gforden_n, geps_n, ggndden_n, gsub_n, gdensity_n, gfdensity_n)

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'
      
      integer nstatv, numgrn, numslip
	  real*8  maxx(1)
      common  /newAddParameters/ maxx
      real*8  statev(nstatv)

      real*8  gstress_n (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran_n (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa_n  (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev_n (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot_n   (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot_n   (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot0    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps_n    (9, MAX_GRN, NUMQPT_T, NUMEL_T)

      real*8  geqvalues (NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot   (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gforden_n (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden_n (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gsub_n    (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gdensity_n (1, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gfdensity_n (1, MAX_GRN, NUMQPT_T, NUMEL_T)

      integer varsPerGrn, ig, dex
c
c---------------------------------------------------------------------72
c
c---- number of state variables per ip
c
      varsPerGrn = 2*NVECS          ! stress, estran
c     &             + NKAPP          ! kappa
     &             + numslip        ! kappa
     &             + NSTAV          ! ee_v, p, detVe, WpNorm, ssact
     &             + NEQVA/2        ! eqps, eqstr, gam_star, gamtot 
     &             + 3*DIMS*DIMS    ! C_0, C, R
     &             + numslip        ! gam_dot
     &             + 4*numslip     
     &             + 9
     &             + 1 
     &             + 1 
c
c---- recover state variables from abaqus vector array
c
      do ig = 1, numgrn

         dex = (ig-1)*varsPerGrn + 1           ! stress
         call EqualTensors2(statev(dex), gstress_n(1,ig,1,1), NVECS)

         dex = dex + NVECS                      ! estrain
         call EqualTensors2(statev(dex), gestran_n(1,ig,1,1), NVECS)

         dex = dex + NVECS                      ! kappa
c         call EqualTensors(statev(dex), gkappa_n(1,ig,1,1), NKAPP)
         call EqualTensors2(statev(dex), gkappa_n(1,ig,1,1), numslip)

c         dex = dex + NKAPP              ! ee_v, p, detVe, WpNorm, ssact
         dex = dex + numslip             ! ee_v, p, detVe, WpNorm, ssact
         call EqualTensors2(statev(dex), gstatev_n(1,ig,1,1), NSTAV)

         dex = dex + NSTAV               ! eqps, eqstr, gam_star, gamtot
         call EqualTensors2(statev(dex), geqvalues(1,ig,1,1), NEQVA/2)

         dex = dex + NEQVA/2                    ! C_0
         call EqualTensors2(statev(dex), gcrot0(1,1,ig,1,1), DIMS*DIMS)

         dex = dex + DIMS*DIMS                  ! C
         call EqualTensors2(statev(dex), gcrot_n(1,1,ig,1,1), DIMS*DIMS)

         dex = dex + DIMS*DIMS                  ! R
         call EqualTensors2(statev(dex), grrot_n(1,1,ig,1,1), DIMS*DIMS)

         dex = dex + DIMS*DIMS                  ! gamdot
         call EqualTensors2(statev(dex), ggamdot(1,ig,1,1), numslip)
         
         dex = dex + numslip                    ! forden
         call EqualTensors2(statev(dex), gforden_n(1,ig,1,1), numslip)
         
         dex = dex + numslip                    ! gndden
         call EqualTensors2(statev(dex), ggndden_n(1,ig,1,1), numslip)
         
         dex = dex + numslip                    ! sub
         call EqualTensors2(statev(dex), gsub_n(1,ig,1,1), numslip)
         
         dex = dex + 2 * numslip                    ! eps
         call EqualTensors2(statev(dex), geps_n(1,ig,1,1), 9)
         
         dex = dex + 9                    ! Average forest density
         call EqualTensors2(statev(dex), gdensity_n(1,ig,1,1), 1)
         
         dex = dex + 1                    ! Average GND density
         call EqualTensors2(statev(dex), gfdensity_n(1,ig,1,1), 1)
		 
      enddo
      
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SaveStateVars2(
     &   statev, nstatv, gstress, gestran, gkappa, gstatev, gcrot, 
     &   grrot, gcrot0, ggamdot, geqvalues, numgrn, numslip, 
     &   gforden, geps, ggndden, gsub, gdensity, gfdensity)

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'
      
      integer nstatv, numgrn, numslip 
      real*8  statev(nstatv)

      real*8  gstress   (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran   (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa    (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev   (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot     (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot     (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot0    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps      (9, MAX_GRN, NUMQPT_T, NUMEL_T)

      real*8  geqvalues (NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot   (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gforden   (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden   (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gsub      (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gdensity (1, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gfdensity (1, MAX_GRN, NUMQPT_T, NUMEL_T)

	  

      integer varsPerGrn, ig, dex
c
c---------------------------------------------------------------------72
c
c---- number of state variables per ip
c
      varsPerGrn = 2*NVECS          ! stress, estran
c     &             + NKAPP          ! kappa
     &             + numslip        ! kappa
     &             + NSTAV          ! ee_v, p, detVe, WpNorm, ssact
     &             + NEQVA/2        ! eqps, eqstr, gam_star, gamtot 
     &             + 3*DIMS*DIMS    ! C_0, C, R
     &             + numslip        ! gam_dot
     &             + 4*numslip     
     &             + 9
     &             + 1
     &             + 1 
c
c---- save state variables in abaqus vector array
c
      do ig = 1, numgrn

         dex = (ig-1)*varsPerGrn + 1           ! stress
         call EqualTensors2(gstress(1,ig,1,1), statev(dex), NVECS)

         dex = dex + NVECS                      ! estrain
         call EqualTensors2(gestran(1,ig,1,1), statev(dex), NVECS)

         dex = dex + NVECS                      ! kappa
c         call EqualTensors(gkappa(1,ig,1,1), statev(dex), NKAPP)
         call EqualTensors2(gkappa(1,ig,1,1), statev(dex), numslip)

c         dex = dex + NKAPP              ! ee_v, p, detVe, WpNorm, ssact
         dex = dex + numslip             ! ee_v, p, detVe, WpNorm, ssact
         call EqualTensors2(gstatev(1,ig,1,1), statev(dex), NSTAV)

         dex = dex + NSTAV               ! eqps, eqstr, gam_star, gamtot
         call EqualTensors2(geqvalues((NEQVA/2+1),ig,1,1), statev(dex),
     &                                                        NEQVA/2)

         dex = dex + NEQVA/2                    ! C_0
c         call EqualTensors(gcrot0(1,1,ig,1,1), statev(dex), DIMS*DIMS)

         dex = dex + DIMS*DIMS                  ! C
         call EqualTensors2(gcrot(1,1,ig,1,1), statev(dex), DIMS*DIMS)

         dex = dex + DIMS*DIMS                  ! R
         call EqualTensors2(grrot(1,1,ig,1,1), statev(dex), DIMS*DIMS)

         dex = dex + DIMS*DIMS                  ! gamdot
         call EqualTensors2(ggamdot(1,ig,1,1), statev(dex), numslip)
         
         dex = dex + numslip                    ! forden
         call EqualTensors2(gforden(1,ig,1,1), statev(dex), numslip)
         
         dex = dex + numslip                    ! gndden
         call EqualTensors2(ggndden(1,ig,1,1), statev(dex), numslip)
         
         dex = dex + numslip                    ! sub
         call EqualTensors2(gsub(1,ig,1,1), statev(dex), numslip)
         
         dex = dex + 2 * numslip                    ! eps
         call EqualTensors2(geps(1,ig,1,1), statev(dex), 9)
         
         dex = dex + 9                    ! Average forest density
         call EqualTensors2(gdensity(1,ig,1,1), statev(dex), 1)
         
         dex = dex + 1                    ! Average GND density
         call EqualTensors2(gfdensity(1,ig,1,1), statev(dex), 1)

      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE DeformationRate_n2(
     &   dvec_n, wvec_n, d_kk_n, iqpt
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer iqpt
      real*8  d_kk_n
      real*8  dvec_n(NVECS), wvec_n(3)
c
c---------------------------------------------------------------------72
c
c------ zeros arrays
c
      d_kk_n = 0.0
      call SetTensor2(dvec_n, pzero, NVECS)
      call SetTensor2(wvec_n, pzero, 3)
c
c------ stop program if called during MPS runs ??
c
c      print *, 'DeformationRate_n called ... Exiting'
c      call RunTimeError(XTAL_O, 'DeformationRate_n called')

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SaveStressModuli2(
     &   stress, ddsdde, savg_ij, cavg_ijkl
     &   )

      implicit none
      include 'params_xtal.inc'

      real*8  stress(NTENS), ddsdde(NTENS, NTENS)
      real*8  savg_ij(DIMS, DIMS), cavg_ijkl(DIMS2, DIMS2)

      integer i, j
c
c---------------------------------------------------------------------72
c
c---- stresses
c
      stress(1) = savg_ij(1,1)
      stress(2) = savg_ij(2,2)
      stress(3) = savg_ij(3,3)
      stress(4) = savg_ij(1,2)

      if (NSHR .gt. 1) then
         stress(5) = savg_ij(1,3)
         stress(6) = savg_ij(2,3)
      endif
c
c---- algorithmic moduli
c
      do i = 1, NTENS
         do j = 1, NTENS
            ddsdde(i,j) = cavg_ijkl(i,j)
         enddo
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE EffectiveElastStiffness2(
     &   ddsdde, ntens, ndi
     &   )

      implicit none
      include 'numbers.inc'

      integer ntens, ndi
      real*8  ddsdde(ntens, ntens)

      integer k1, k2
      real*8  eg2, ebulk3, elam

      real*8  eg, ebulk
      common /ElastProps/ eg, ebulk
c
c---------------------------------------------------------------------72
c
c----- effective elastic stiffness
c
      eg2    = 2.0*eg
      ebulk3 = 3.0*ebulk
      elam   = (ebulk3-eg2)/pthree

      call SetTensor2(ddsdde, pzero, ntens*ntens)

      do k1 = 1, ndi
        do k2 = 1, ndi
           ddsdde(k2,k1) = elam
        enddo
        ddsdde(k1,k1) = eg2+elam
      enddo
      do k1 = ndi+1, ntens
        ddsdde(k1,k1) = eg
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE CrystalInitialize2(
     &   props, mprops, numel, numqpt,
     &   filePath, fileRoot,
     &   FILE_I, FILE_E, FILE_O, FILE_G1, FILE_G2,
     &   FILE_Z1, FILE_Z2, FILE_STRESS_O, FILE_STRAIN_O,
     &   FILE_EFFSS_O, FILE_TRUESS_O,
     &   FILE_ITER_O, FILE_AGG_EFFSS_O,
     &   numgrn, numslip, numvtx, kODFout,
     &   kappa0, fCeDevVol, fCeVol,
     &   matProp, tauSlip,
     &   fCeDev, fCeiDev,
     &   zBar0,
     &   pBar0, qBar0,
     &   pBar0Vec, qBar0Vec,
     &   ppTBar0,
     &   gcrot0,
     &   sigfs,
     &   overstress,
     &   hardmtx,
     &   maxIterState, maxIterNewt,
     &   tolerState, tolerNewt,
     &   kODF, numor, seed, angles, euler,
     &   FILE_TXT_OUT,
     &   gstress,
     &   gstress_n,
     &   gestran,
     &   gestran_n,
     &   gkappa,
     &   gkappa_n,
     &   gstatev,
     &   gstatev_n,
     &   geqvalues,
     &   ggamdot,
     &   gcrot,
     &   gcrot_n,
     &   grrot,
     &   grrot_n, 
     &   gforden, gforden_n, geps, geps_n,
     &   ggndden, ggndden_n, gsub, gsub_n, gdensity, gdensity_n,
     &   gfdensity, gfdensity_n)

      implicit none
      include 'params_xtal.inc'

      integer mprops, numel, numqpt
      real*8  props(mprops)

      character*160 filePath, fileRoot
c
      integer numgrn, numslip, numvtx, kODFout
      real*8  kappa0(NKAPP), fCeDevVol(NVECS), fCeVol
      real*8  matProp(NPROPS, MAX_SLIP), tauSlip(MAX_SLIP)
      real*8  fCeDev(NVECS, NVECS), fCeiDev(NVECS, NVECS)
      real*8  zBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0(DIMS, DIMS, MAX_SLIP), qBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0Vec(NVECS, MAX_SLIP), qBar0Vec(DIMS, MAX_SLIP)
      real*8  ppTBar0(NVECS, NVECS, MAX_SLIP)
      real*8  gcrot0(DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  sigfs(5, MAX_VTX)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)

      real*8  overstress(MAX_SLIP)

      integer maxIterState, maxIterNewt
      real*8  tolerState, tolerNewt
c
      real*8  gstress    (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstress_n  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran    (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran_n  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa     (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa_n   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev    (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev_n  (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geqvalues  (NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot    (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot      (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot_n    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot      (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot_n    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gforden    (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gforden_n  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps       (9, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps_n     (9, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden    (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden_n  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gsub       (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gsub_n     (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gdensity   (1, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gdensity_n (1, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gfdensity   (1, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gfdensity_n (1, MAX_GRN, NUMQPT_T, NUMEL_T)
c
      integer kODF

      real*8  angles(DIMS, MAX_ORIEN)
      real*8  euler(DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)

      integer numor, seed

      integer FILE_I, FILE_E, FILE_O, FILE_G1, FILE_G2,
     &        FILE_Z1, FILE_Z2, FILE_STRESS_O, FILE_STRAIN_O,
     &        FILE_EFFSS_O, FILE_TRUESS_O,
     &        FILE_ITER_O, FILE_AGG_EFFSS_O,
     &        FILE_TXT_OUT
c
c---------------------------------------------------------------------72
c
c------- open input/output files
c
      call CrystalOpenIOFiles2(filePath, fileRoot, 
     &                        FILE_I, FILE_E, FILE_O, FILE_G1, FILE_G2,
     &                        FILE_Z1, FILE_Z2)
      call CrystalOpenIOFiles_22(filePath, fileRoot, 
     &                          FILE_STRESS_O, FILE_STRAIN_O, 
     &                          FILE_EFFSS_O, FILE_TRUESS_O, 
     &                          FILE_ITER_O, FILE_AGG_EFFSS_O)
c
c------- read/echo material data
c
      call CrystalModelData2(props, mprops, numqpt, numel,
     &                      numgrn, numslip, numvtx, kODFout,
     &                      kappa0, fCeDevVol, fCeVol,
     &                      matProp, tauSlip,
     &                      fCeDev, fCeiDev,
     &                      zBar0,
     &                      pBar0, qBar0,
     &                      pBar0Vec, qBar0Vec,
     &                      ppTBar0,
     &                      gcrot0,
     &                      sigfs,
     &                      overstress,
     &                      hardmtx,
     &                      maxIterState, maxIterNewt,
     &                      tolerState, tolerNewt,
     &                      FILE_I, FILE_E, FILE_O, filePath,
     &                      kODF, numor, seed, angles, euler,
     &                      FILE_TXT_OUT)
c
c------- initialize arrays
c
      call CrystalInitializeArrays2(numqpt, numel,
     &                             numgrn, numslip,
     &                             kappa0,
     &                             gcrot0,
     &                             gstress,
     &                             gstress_n,
     &                             gestran,
     &                             gestran_n,
     &                             gkappa,
     &                             gkappa_n,
     &                             gstatev,
     &                             gstatev_n,
     &                             geqvalues,
     &                             ggamdot,
     &                             gcrot,
     &                             gcrot_n,
     &                             grrot,
     &                             grrot_n,
     &                             gforden, gforden_n, geps, geps_n,
     &                             ggndden, ggndden_n, gsub, gsub_n,
     &                             gdensity,gdensity_n, gfdensity,
     &                             gfdensity_n)
c
c------- initialize useful matrices to compute algorithmic moduli Cep
c
      call CrystalInitializeMatrxCep2( )
c
c------- parameters for global iterations
c
      call GlobalIterationParams2( )

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE CrystalOpenIOFiles2(
     &   filePath, fileRoot, 
     &   FILE_I, FILE_E, FILE_O, FILE_G1, FILE_G2,
     &   FILE_Z1, FILE_Z2)

      implicit none
    
      character*160 filePath, fileRoot
      integer      FILE_I, FILE_E, FILE_O, FILE_G1, FILE_G2
      integer      FILE_Z1, FILE_Z2

      integer      length1, length2
      character*160 dataFile, filename
c
c---------------------------------------------------------------------72
c
c------- root name of input/output files
c
c      write(*,*) ' Root name of crystal-model input/output files ?'
c      read 1000, dataFile
c
c------- open files
c
      length1 = index(filePath,' ') - 1
      length2 = index(fileRoot,' ') - 1

      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtali'
      open(unit=FILE_I, file=filename, status='unknown',
     &                                           access='sequential')
      rewind(FILE_I)

      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtale'
      open(unit=FILE_E, file=filename, status='unknown')
      rewind(FILE_E)

      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtalo'
      open(unit=FILE_O, file=filename, status='unknown')
      
      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtalg1'
      open(unit=FILE_G1, file=filename, status='unknown')
      
      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtalg2'
      open(unit=FILE_G2, file=filename, status='unknown')
      
      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtalz1'
      open(unit=FILE_Z1, file=filename, status='unknown')
      
      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtalz2'
      open(unit=FILE_Z2, file=filename, status='unknown')
c
c------- formats
c
1000  format(a80)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE CrystalOpenIOFiles_22(
     &   filePath, fileRoot, 
     &   FILE_STRESS_O, FILE_STRAIN_O, 
     &   FILE_EFFSS_O, FILE_TRUESS_O, 
     &   FILE_ITER_O, FILE_AGG_EFFSS_O
     &   )

      implicit none

      character*160 filePath, fileRoot
      integer FILE_STRESS_O, FILE_STRAIN_O, 
     &        FILE_EFFSS_O, FILE_TRUESS_O, 
     &        FILE_ITER_O, FILE_AGG_EFFSS_O
    
      integer      length1, length2
      character*160 dataFile, filename
c
c---------------------------------------------------------------------72
c
c------- root name of i/o files
c
c      write(*,*) ' Root name of additional output files ?'
c      read 1000, dataFile
c
c------- open files for output at specified elem/iqpt/grain
c
      length1 = index(filePath,' ') - 1
      length2 = index(fileRoot,' ') - 1

      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtal.strs'
      open(unit=FILE_STRESS_O, file=filename, status='unknown')

      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtal.strn'
      open(unit=FILE_STRAIN_O, file=filename, status='unknown')

      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtal.efss'
      open(unit=FILE_EFFSS_O, file=filename, status='unknown')

      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtal.trss'
      open(unit=FILE_TRUESS_O, file=filename, status='unknown')

      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtal.iter'
      open(unit=FILE_ITER_O, file=filename, status='unknown')
c
c------- open files for aggregate output at specified elem/iqpt
c
      filename = filePath(1:length1)//fileRoot(1:length2)//'.agg.efss'
      open(unit=FILE_AGG_EFFSS_O, file=filename, status='unknown')
c
c------- formats
c
1000  format(a80)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE CrystalModelData2(
     &   props, mprops, numqpt, numel,
     &   numgrn, numslip, numvtx, kODFout,
     &   kappa0, fCeDevVol, fCeVol,
     &   matProp, tauSlip,
     &   fCeDev, fCeiDev,
     &   zBar0,
     &   pBar0, qBar0,
     &   pBar0Vec, qBar0Vec,
     &   ppTBar0,
     &   gcrot0,
     &   sigfs,
     &   overstress,
     &   hardmtx,
     &   maxIterState, maxIterNewt,
     &   tolerState, tolerNewt,
     &   FILE_I, FILE_E, FILE_O, filePath,
     &   kODF, numor, seed, angles, euler,
     &   FILE_TXT_OUT
     &   )

      implicit none
      include 'params_xtal.inc'

      character*160 filePath
      integer mprops, numqpt, numel
      real*8  props(mprops)

      integer numgrn, numslip, numvtx, kODFout
      real*8  kappa0(NKAPP), fCeDevVol(NVECS), fCeVol
      real*8  matProp(NPROPS, MAX_SLIP), tauSlip(MAX_SLIP)
      real*8  fCeDev(NVECS, NVECS), fCeiDev(NVECS, NVECS)
      real*8  zBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0(DIMS, DIMS, MAX_SLIP), qBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0Vec(NVECS, MAX_SLIP), qBar0Vec(DIMS, MAX_SLIP)
      real*8  ppTBar0(NVECS, NVECS, MAX_SLIP)
      real*8  gcrot0(DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  sigfs(5, MAX_VTX)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)

      real*8  overstress(MAX_SLIP)

      integer maxIterState, maxIterNewt
      real*8  tolerState, tolerNewt

      integer FILE_I, FILE_E, FILE_O
      integer FILE_TXT_OUT

      integer kODF, numor, seed
      real*8  angles(DIMS, MAX_ORIEN)
      real*8  euler(DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
c
c---------------------------------------------------------------------72
c
c------- input material data for single crystal
c
      call CrystalMaterialData2(props, mprops, numqpt, numel,
     &                      numgrn, numslip, numvtx, kODFout,
     &                      kappa0, fCeDevVol, fCeVol,
     &                      matProp, tauSlip,
     &                      fCeDev, fCeiDev,
     &                      zBar0,
     &                      pBar0, qBar0,
     &                      pBar0Vec, qBar0Vec,
     &                      ppTBar0,
     &                      gcrot0,
     &                      sigfs,
     &                      overstress,
     &                      hardmtx,
     &                      FILE_I, FILE_E, FILE_O, filePath,
     &                      kODF, numor, seed, angles, euler,
     &                      FILE_TXT_OUT)
c
c------- input convergence data (state iterations & constitutive solver)
c
      call CrystalSolverData2(props, mprops,
     &                      maxIterState, maxIterNewt,
     &                      tolerState, tolerNewt, 
     &                      FILE_I, FILE_E)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE CrystalMaterialData2(
     &   props, mprops, numqpt, numel,
     &   numgrn, numslip, numvtx, kODFout,
     &   kappa0, fCeDevVol, fCeVol,
     &   matProp, tauSlip,
     &   fCeDev, fCeiDev,
     &   zBar0,
     &   pBar0, qBar0,
     &   pBar0Vec, qBar0Vec,
     &   ppTBar0,
     &   gcrot0,
     &   sigfs,
     &   overstress,
     &   hardmtx,
     &   FILE_I, FILE_E, FILE_O, filePath,
     &   kODF, numor, seed, angles, euler,
     &   FILE_TXT_OUT
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      character*160 filePath
      integer mprops, numqpt, numel
      real*8  props(mprops)

      integer numgrn, numslip, numvtx, kODFout
      real*8  kappa0(NKAPP), fCeDevVol(NVECS), fCeVol
      real*8  matProp(NPROPS, MAX_SLIP), tauSlip(MAX_SLIP)
      real*8  fCeDev(NVECS, NVECS), fCeiDev(NVECS, NVECS)
      real*8  zBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0(DIMS, DIMS, MAX_SLIP), qBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0Vec(NVECS, MAX_SLIP), qBar0Vec(DIMS, MAX_SLIP)
      real*8  ppTBar0(NVECS, NVECS, MAX_SLIP)
      real*8  gcrot0(DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  sigfs(5, MAX_VTX)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)

      real*8  overstress(MAX_SLIP)

      integer FILE_I, FILE_E, FILE_O
      integer FILE_TXT_OUT

      integer kODF, numor, seed
      real*8  angles(DIMS, MAX_ORIEN)
      real*8  euler(DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  xtalProp(NPROPS)

      integer is, i
      character*16 scysFile
      integer crystalID

      integer SXFILE_I
      parameter(SXFILE_I=91)

      integer    length1, length2
      character  filename*160, SingleXtalFile*160, prosa*160
c
c---------------------------------------------------------------------72
c
c---- initialize some arrays that store material data / slip system
c
      call SetTensor2(kappa0,  pzero, NKAPP)
      call SetTensor2(matProp, pzero, NPROPS*MAX_SLIP)
      call SetTensor2(hardmtx, pzero, MAX_SLIP*MAX_SLIP)
  
      call SetTensor2(xtalProp, pzero, NPROPS)
c
c---- crystal type and number of grains per IP
c
      call SetCrystalType2(crystalID, numgrn, props, mprops, 
     &                    FILE_I, FILE_E, FILE_O)
c
c---- open single crystal filename

      read(FILE_I,'(a18)') SingleXtalFile

      length1 = index(filePath,' ') - 1
      length2 = index(SingleXtalFile,' ') - 1

      filename = filePath(1:length1)//SingleXtalFile(1:length2)
      open(unit=SXFILE_I, file=filename, status='unknown',
     &                                   access='sequential')
      rewind(SXFILE_I)
c
c---- elastic stiffness matrix of single crystal
c
      call SetCrystalElasticity2(crystalID, fCeDev, fCeiDev, xtalProp, 
     &                          fCeDevVol, fCeVol, props, mprops,
     &                          SXFILE_I, FILE_E, FILE_O)
c
c---- thermal expansion coefficients
c
c      call SetCrystalThermal(xtalProp, SXFILE_I, FILE_E)
       read(SXFILE_I,'(a)') prosa
c
c---- set crystal geometry: slip / twinning data
c
      call SetCrystalGeometry2(crystalID, numslip, numvtx, scysFile, 
     &   zBar0, pBar0, qBar0, pBar0Vec, qBar0Vec, ppTBar0, xtalProp, 
     &   matProp, hardmtx, kappa0, tauSlip, props, mprops, SXFILE_I, 
     &   FILE_E)
c
c---- close single crystal filename
c
      close(SXFILE_I)
c
c---- crystal orientations
c
      call SetCrystalLatticeOrient2(numgrn, numqpt, numel, kODFout, 
     &                             gcrot0, props, mprops,
     &                             kODF, numor, seed, angles, euler,
     &                             FILE_I, FILE_E, FILE_O, filePath,
     &                             FILE_TXT_OUT)
c
c---- vertices of rate independent yield surface (single crystal)
c
      if (crystalID .eq. kFCC .or. crystalID .eq. kBCC) then
         call VerticesSCYS_Cubic2(numvtx, kappa0(1), sigfs, scysFile,
     &                           FILE_E)
      else  ! crystalID=kHCP
         call VerticesSCYS_HCP2(numvtx, kappa0(1), sigfs, scysFile,
     &                         FILE_E, filePath)
      endif
c
c---- For HCP crystals, the overstress will keep the same initial 
c---- difference in slip system hardness between basal/prismatic
c---- and pyramidal slip systems. 
c---- For FCC/BCC crystals, it won't have any effect.
c
      do is = 1, numslip
         overstress(is) = (tauSlip(is) - 1.0) * kappa0(1)
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE CrystalSolverData2(
     &   props, mprops,
     &   maxIterState, maxIterNewt,
     &   tolerState, tolerNewt,
     &   FILE_I, FILE_E
     &   )

      implicit none

      integer mprops
      real*8  props(mprops)
     
      integer maxIterState, maxIterNewt
      real*8  tolerState, tolerNewt

      integer FILE_I, FILE_E
c
c---------------------------------------------------------------------72
c
c------- number iterations and tolerance for state iterations
c
      read(FILE_I, *) maxIterState, tolerState
c       maxIterState = nint (props(19))
c       tolerState   = props(20)
c
c------- number iterations and tolerance for newton method
c
      read(FILE_I, *) maxIterNewt, tolerNewt
c       maxIterNewt = nint (props(21))
c       tolerNewt   = props(22)
c
c------- echo input data
c
      write(FILE_E, 1000) maxIterState, maxIterNewt, 
     &                    tolerState, tolerNewt      
c
c------- format
c
1000  format(/'*-----   Local Convergence Control Parameters  -----*'/,
     &        7x,'  max iters State  = ',i5 /,
     &        7x,'  max iters Newton = ',i5 /,
     &        7x,'  tolerance State  = ',e12.5 /,
     &        7x,'  tolerance Newton = ',e12.5)
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE CrystalInitializeArrays2(
     &   numqpt, numel,
     &   numgrn, numslip,
     &   kappa0,
     &   gcrot0,
     &   gstress,
     &   gstress_n,
     &   gestran,
     &   gestran_n,
     &   gkappa,
     &   gkappa_n,
     &   gstatev,
     &   gstatev_n,
     &   geqvalues,
     &   ggamdot,
     &   gcrot,
     &   gcrot_n,
     &   grrot,
     &   grrot_n,
     &   gforden, gforden_n, geps, geps_n,
     &   ggndden, ggndden_n, gsub, gsub_n,
     &   gdensity, gdensity_n,
     &   gfdensity, gfdensity_n)

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numqpt, numel
      integer numgrn, numslip
      real*8  kappa0(NKAPP)
      real*8  gcrot0 (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)

      real*8  gstress    (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstress_n  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran    (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran_n  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa     (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa_n   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev    (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev_n  (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geqvalues  (NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot    (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot      (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot_n    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot      (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot_n    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gforden    (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gforden_n  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps       (9, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps_n     (9, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden    (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden_n  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gsub       (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gsub_n     (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gdensity   (1, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gdensity_n (1, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gfdensity   (1, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gfdensity_n (1, MAX_GRN, NUMQPT_T, NUMEL_T)

      integer ie, ip, ig
c
c---------------------------------------------------------------------72
c
c------- initialize arrays used in XTAL constitutive integration method
c
      do ie = 1, numel
        do ip = 1, numqpt
          do ig = 1, numgrn

            call SetTensor2(gstress  (1,ig,ip,ie), pzero, NVECS)
            call SetTensor2(gestran  (1,ig,ip,ie), pzero, NVECS)
            call SetTensor2(gkappa   (1,ig,ip,ie), pzero, NKAPP)
            call SetTensor2(gstatev  (1,ig,ip,ie), pzero, NSTAV)
            call SetTensor2(gstress_n(1,ig,ip,ie), pzero, NVECS)
            call SetTensor2(gestran_n(1,ig,ip,ie), pzero, NVECS)
            call SetTensor2(gkappa_n (1,ig,ip,ie), pzero, NKAPP)
            call SetTensor2(gstatev_n(1,ig,ip,ie), pzero, NSTAV)
            call SetTensor2(geqvalues(1,ig,ip,ie), pzero, NEQVA)
            call SetTensor2(ggamdot  (1,ig,ip,ie), pzero, MAX_SLIP)
            call SetTensor2(gforden  (1,ig,ip,ie), pzero, MAX_SLIP)
            call SetTensor2(gforden_n(1,ig,ip,ie), pzero, MAX_SLIP)
            call SetTensor2(geps     (1,ig,ip,ie), pzero, 9)
            call SetTensor2(geps_n   (1,ig,ip,ie), pzero, 9)
            call SetTensor2(ggndden  (1,ig,ip,ie), pzero, MAX_SLIP)
            call SetTensor2(ggndden_n(1,ig,ip,ie), pzero, MAX_SLIP)
            call SetTensor2(gsub     (1,ig,ip,ie), pzero, MAX_SLIP)
            call SetTensor2(gsub_n   (1,ig,ip,ie), pzero, MAX_SLIP)
            call SetTensor2(gdensity  (1,ig,ip,ie), pzero, 1)
            call SetTensor2(gdensity_n(1,ig,ip,ie), pzero, 1)
            call SetTensor2(gfdensity  (1,ig,ip,ie), pzero, 1)
            call SetTensor2(gfdensity_n(1,ig,ip,ie), pzero, 1)

            call SetTensor2(gcrot  (1,1,ig,ip,ie), pzero, DIMS*DIMS)
            call SetTensor2(gcrot_n(1,1,ig,ip,ie), pzero, DIMS*DIMS)
            call SetTensor2(grrot  (1,1,ig,ip,ie), pzero, DIMS*DIMS)
            call SetTensor2(grrot_n(1,1,ig,ip,ie), pzero, DIMS*DIMS)

          enddo
        enddo
      enddo
c
c------- initialize state variables and rotation tensors
c
      do ie = 1, numel
        do ip = 1, numqpt
          do ig = 1, numgrn

            call EqualTensors2(kappa0, gkappa_n(1,ig,ip,ie), NKAPP)
            call EqualTensors2(gcrot0(1,1,ig,ip,ie), 
     &                         gcrot_n(1,1,ig,ip,ie), DIMS*DIMS)
            call EqualTensors2(Ident2nd, 
     &                         grrot_n(1,1,ig,ip,ie), DIMS*DIMS)

          enddo
        enddo
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE CrystalCloseIOFiles2( )

      implicit none
      include 'params_xtal.inc'
c
c---------------------------------------------------------------------72
c
c------- close files
c
      close(XTAL_I)
      close(XTAL_O)
      close(XTAL_E)
      close(XTAL_G1)
      close(XTAL_G2)
c
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE CrystalCloseIOFiles_22( )

      implicit none
      include 'params_xtal.inc'
c
c---------------------------------------------------------------------72
c
c------- close files
c
      close(XTAL_STRESS_O)
      close(XTAL_STRAIN_O)
      close(XTAL_EFFSS_O)
      close(XTAL_TRUESS_O)
      close(XTAL_ITER_O)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE CrystalInitializeMatrxCep2( 
     &
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer i, j

      real*8  unit4(DIMS2, DIMS2), unitDev4(DIMS2, DIMS2)
      common /unit4Tensors/ unit4, unitDev4

      real*8  fDevMat5x6(NVECS, DIMS2), fDevMat6x5(DIMS2, NVECS)
      real*8  fMatTId5x6(NVECS, DIMS2)
      common /transfMatxs/ fDevMat5x6, fDevMat6x5, fMatTId5x6
c
c---------------------------------------------------------------------72
c
c----- zero 4th order unit tensors and transformation matrices

      call SetTensor2(unit4, pzero, DIMS2*DIMS2)
      call SetTensor2(unitDev4, pzero, DIMS2*DIMS2)

      call SetTensor2(fDevMat5x6, pzero, NVECS*DIMS2)
      call SetTensor2(fDevMat6x5, pzero, DIMS2*NVECS)
      call SetTensor2(fMatTId5x6, pzero, NVECS*DIMS2)
c
c----- set 4th order unit tensor (I)
c
      do i = 1, DIMS2
         unit4(i, i) = pone
      enddo
c
c----- set 4th order deviatoric unit tensor (I_dev) 
c
      do i = 1, DIMS2/2
         unitDev4(i, i) = pone
         do j = 1, DIMS2/2
            unitDev4(i, j) = unitDev4(i, j) - pthird
         enddo
      enddo
 
      do i = DIMS2/2+1, DIMS2
         unitDev4(i, i) = pone
      enddo
c
c----- set tranformation matrix [T] = [tDevMat]_5x6:
c-----    {}'_5x1 = [tDevMat]_5x6 {}'_6x1
c
      fDevMat5x6(1,1) =  pone/sqr2
      fDevMat5x6(1,2) = -pone/sqr2
      fDevMat5x6(2,3) = sqr32
      fDevMat5x6(3,4) = sqr2
      fDevMat5x6(4,5) = sqr2
      fDevMat5x6(5,6) = sqr2
c
c----- set tranformation matrix [J] = [tDevMat]_6x5:
c-----    {}'_6x1 = [tDevMat]_6x5 {}'_5x1
c
      fDevMat6x5(1,1) =  pone/sqr2
      fDevMat6x5(1,2) = -pone/(sqr2*sqr3)
      fDevMat6x5(2,1) = -pone/sqr2
      fDevMat6x5(2,2) = -pone/(sqr2*sqr3)
      fDevMat6x5(3,2) = sqr23
      fDevMat6x5(4,3) = pone/sqr2
      fDevMat6x5(5,4) = pone/sqr2
      fDevMat6x5(6,5) = pone/sqr2
c
c----- needed product: [TId]_5x6 = [T]_5x6 [I_dev]_6x6
c
      call MultAxB_G2(fDevMat5x6, unitDev4, fMatTId5x6, NVECS, DIMS2, 
     &               DIMS2) 

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE GlobalIterationParams2( 
     &
     &   )

      implicit none
      include 'params_xtal.inc'

      integer maxIterN
      real*8  tolerN, divTolerN, zeroTolerN
      common  /ParmsGlobalN_1/ maxIterN
      common  /ParmsGlobalN_2/ tolerN, zeroTolerN, divTolerN

      logical alwaysLineSearch
      integer searchIters
      real*8  maxStepSize, orthoToler
      common  /ParmsLS_0/ alwaysLineSearch
      common  /ParmsLS_1/ searchIters
      common  /ParmsLS_2/ maxStepSize, orthoToler

      integer maxIncrCuts, quickSolveTol, quickSeriesTol
      real*8  DEQP_ref
      common  /dtContrl_1/ maxIncrCuts, quickSolveTol, quickSeriesTol
      common  /dtContrl_2/ DEQP_ref
c
c---------------------------------------------------------------------72
c
c------- params_xtal for global newton iterations
c
      read(XTAL_I, *) maxIterN
      read(XTAL_I, *) tolerN, zeroTolerN, divTolerN
c
c------- params_xtal for global line search iterations
c
      read(XTAL_I, *) alwaysLineSearch
      read(XTAL_I, *) searchIters
      read(XTAL_I, *) maxStepSize, orthoToler
c
c------- params_xtal to control the increase/decrease of time step
c
      read(XTAL_I, *) maxIncrCuts
      read(XTAL_I, *) quickSolveTol
      read(XTAL_I, *) quickSeriesTol
      read(XTAL_I, *) DEQP_ref
c
c------- echo input data 
c
      write(XTAL_E, 1000) maxIterN, tolerN, zeroTolerN, divTolerN
      write(XTAL_E, 2000) alwaysLineSearch, searchIters, maxStepSize, 
     &                   orthoToler
      write(XTAL_E, 3000) maxIncrCuts, quickSolveTol, quickSeriesTol,
     &                   DEQP_ref
c
c------- formats
c 
1000  format(/'*-----   Global Newton Iterations   -----*'/,
     &       7x,'  maxIterN    = ',i6/
     &       7x,'  tolerN      = ',e12.5/
     &       7x,'  zeroTolerN  = ',e12.5/
     &       7x,'  divTolerN   = ',e12.5)
2000  format(/'*-----   Global Line Search Iterations   -----*'/,
     &       7x,'  alwaysLS    = ',l6/
     &       7x,'  searchIters = ',i6/
     &       7x,'  maxStepSize = ',e12.5/
     &       7x,'  orthoToler  = ',e12.5)
3000  format(/'*-----   Reset Params for controlling dtime   -----*'/,
     &       7x,'  maxIncrCuts    = ',i6/
     &       7x,'  quickSolveTol  = ',i6/
     &       7x,'  quickSeriesTol = ',i6/
     &       7x,'  DEQP_ref       = ',e12.5)
          
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetCrystalType2(
     &   crystalID, numgrn, props, mprops, FILE_I, FILE_E, FILE_O
     &   )

      implicit none
      include 'params_xtal.inc'

      integer crystalID, numgrn, FILE_I, FILE_E, FILE_O
      integer mprops
      real*8  props(mprops)
c
c---------------------------------------------------------------------72
c
c------- crystal type and number of grains/IP
c
      read(FILE_I, *) crystalID, numgrn
c      crystalID = nint (props(1))
c      numgrn   = nint (props(2))
      if (crystalID .ne. kFCC .and.
     &    crystalID .ne. kBCC .and.
     &    crystalID .ne. kHCP)
     &   call RunTimeError2(FILE_O, 'SetCrystalType: crystalID ?')
c
c------- echo input data
      write(FILE_E, 1000) crystalID, numgrn
c
c------- format
c
1000  format(/'*-----   Crystal Type and Number of Grns per IP -----*'/,
     &        7x,'  crystal type     = ',i5 /,
     &        7x,'  number grains/ip = ',i5)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetCrystalElasticity2(
     &   crystalID, fCeDev, fCeiDev, matProp, fCeDevVol, fCeVol, 
     &   props, mprops, SXFILE_I, FILE_E, FILE_O
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer crystalID, SXFILE_I, FILE_E, FILE_O
      real*8  fCeDevVol(NVECS), fCeVol
      real*8  matProp(NPROPS)
      real*8  fCeDev(NVECS,NVECS), fCeiDev(NVECS,NVECS)

      integer mprops
      real*8  props(mprops)

      integer elastID, i, j
      real*8  eMod, eNu, gMod, eBulk
      real*8  c1, c2, c3, c4, c5, const

      common /ElastProps/ gMod, eBulk
c
c---------------------------------------------------------------------72
c
c------- type of elasticity (in crystal coordinates)
c-------    elastID = kELAS_ISO : isotropic elasticity
c-------    elastID = kELAS_ANI : anisotropic elasticity
c
      read(SXFILE_I, *) elastID
c      elastID = nint (props(3))
      if (elastID .ne. kELAS_ISO .and.
     &    elastID .ne. kELAS_ANI)
     &   call RunTimeError2(FILE_O, 'SetCrystalElasticity: elastID ?')
      write(FILE_E, 1000) elastID
c
c------- initialize deviatoric elastic stiffness, its inverse and
c-------  coupled deviatoric-volumetric term
c
      call SetTensor2(fCeDev, pzero, NVECS*NVECS)
      call SetTensor2(fCeiDev, pzero, NVECS*NVECS)
      call SetTensor2(fCeDevVol, pzero, NVECS)
c
c------- read elastic constants, i.e.,
c-------     ISO   FCC  HCP
c-------     eMod  C11  C11
c-------     eMu   C12  C12
c-------           C44  C13
c-------                C33
c-------                C44
c------- and build deviatoric elastic stiffness, as well as
c------- coupled deviatoric-volumetric and volumetric terms.
c-------  For Isotropic and Cubic Crystala: 
c-------      fCeDevVol = 0, fCeVol = Bulk Modulues
c-------  For Hexagonal Closed Packed Crystals:
c-------      si_dev = fCeDev     *ei_dev + fCeDevVol*e_kk
c-------      press  = fCeDevVol^T*ei_dev + fCeVol   *e_kk
c-------   Here: 
c-------      si_dev = {(s11-s22)/sq2,sq32*s33,sq2*s12,sq2*s13,sq2*s23}
c-------      ei_dev = {(e11-e22)/sq2,sq32*e33,sq2*e12,sq2*e13,sq2*e23}
c

      if (elastID .eq. kELAS_ISO) then
c
c---------- isotropic elastic stiffness
         read(SXFILE_I, *)  eMod, eNu 
c         eMod = props(4)
c         eNu  = props(5)
         write(FILE_E, 2000) eMod, eNu
         matProp(1) = eMod / (1. + eNu) / 2.       ! gMod
         matProp(2) = eMod / (1. - 2. * eNu) / 3.  ! eBulk
         fCeDev(1,1) = 2.0*matProp(1)
         fCeDev(2,2) = 2.0*matProp(1)
         fCeDev(3,3) = 2.0*matProp(1)
         fCeDev(4,4) = 2.0*matProp(1)
         fCeDev(5,5) = 2.0*matProp(1)
         fCeVol      = matProp(2)
      elseif ((elastID .eq. kELAS_ANI .and. crystalID .eq. kFCC) .or.
     &        (elastID .eq. kELAS_ANI .and. crystalID .eq. kBCC)) then
c
c---------- anisotropic elastic stiffness, FCC or BCC
         read(SXFILE_I, *) c1, c2, c3
c         c1 = props(4)             ! C11
c         c2 = props(5)             ! C12
c         c3 = props(6)             ! C44
         write(FILE_E, 3000) c1, c2, c3
         matProp(1) = (2.0*(c1 - c2) + 6.0*c3) / 10.0
         matProp(2) = (c1 + 2.0*c2) / 3.0
         fCeDev(1,1) = c1 - c2
         fCeDev(2,2) = c1 - c2
         fCeDev(3,3) = 2.0*c3
         fCeDev(4,4) = 2.0*c3
         fCeDev(5,5) = 2.0*c3
         fCeVol      = (c1 + 2.0*c2) / 3.0
      elseif (elastID .eq. 2 .and. crystalID .eq. kHCP) then
c
c---------- anisotropic elastic stiffness, HCP
         read(SXFILE_I, *) c1, c2, c3, c4, c5
c         c1 = props(4)             ! C11
c         c2 = props(5)             ! C12
c         c3 = props(6)             ! C13
c         c4 = props(7)             ! C33
c         c5 = props(8)             ! C44
         write(FILE_E, 4000) c1, c2, c3, c4, c5
         const = c1 + c2 - 4.0*c3 + 2.0*c4
         matProp(1) = (6.0*(c1-c2) + const + 12.0*c5) / 30.0
         matProp(2) = ((c1 + c2)*c4 - 2.0*c3*c3) / const
         fCeDev(1,1) = c1 - c2
         fCeDev(2,2) = const / 3.0
         fCeDev(3,3) = c1 - c2
         fCeDev(4,4) = 2.0*c5
         fCeDev(5,5) = 2.0*c5
         fCeDevVol(2) = - dsqrt(2.d0/27.d0) * (c1 + c2 - c3 - c4)
         fCeVol       = (2.0*c1 + 2.0*c2 + 4.0*c3 + c4) / 9.0
      endif         
c
c------- effective shear and bulk modulus (to compute effective Ce)
c
      gMod  = matProp(1)
      eBulk = matProp(2)
c
c------- inverse of deviatoric elastic stiffness
c
      do i = 1, NVECS
         fCeiDev(i,i) = 1. / fCeDev(i,i)
      enddo
c
c------- echo some computed data
c
      write(FILE_E, 5000) matProp(1), matProp(2), fCeVol,
     &                    (fCeDevVol(i),i=1,NVECS)
      write(FILE_E, 6000) ((fCeDev(i,j), j=1,5), i=1,5)
c
c------- format
c
1000  format(/'*-----   Crystal Elasticity -----*'/,
     &        7x,'  elasticity type  = ',i5)
2000  format( 7x,'  young modulus    = ',e12.5/
     &        7x,'  poisson ratio    = ',e12.5)
3000  format( 7x,'  c_1              = ',e12.5/
     &        7x,'  c_2              = ',e12.5/
     &        7x,'  c_3              = ',e12.5)
4000  format( 7x,'  c_1              = ',e12.5/
     &        7x,'  c_2              = ',e12.5/
     &        7x,'  c_3              = ',e12.5/
     &        7x,'  c_4              = ',e12.5/
     &        7x,'  c_5              = ',e12.5)
5000  format( 7x,'  shear modulus    = ',e12.5/
     &        7x,'  bulk modulus     = ',e12.5/
     &        7x,'  fCeVol           = ',e12.5/
     &        7x,'  fCeDevVol(1...5) = ',5(e12.5,1x)/
     &        7x,'  CeDev : ')
6000  format((15x,5(e12.5,2x)))

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetCrystalPlasticity2(
     &   crystalID, xtalProp, crss0, props, mprops, mode, 
     &   SXFILE_I, FILE_E
     &   )

      implicit none
      include 'params_xtal.inc'

      integer crystalID, mode, SXFILE_I, FILE_E
      real*8  crss0
      real*8  xtalProp(NPROPS)

      integer mprops
      real*8  props(mprops)

      integer i
c
c---------------------------------------------------------------------72
c
c------- slip system kinetic equation: power law (~thermal activation)
c-------   xtalProp(3) xm   : strain-rate sensitivity (m)
c-------   xtalProp(4) gam0 : reference shear rate (d_0)
c
      read(SXFILE_I, *)  xtalProp(3), xtalProp(4)
c      xtalProp(3) = props( 9)
c      xtalProp(4) = props(10)
c
c------ bounds for argument of power law
c
      call BoundForArgPowLaw2(xtalProp(3))
c
c------- slip system kinetic equation: drag-controlled plastic flow
c-------   xtalProp(11) bdrag : drag coefficient (B)
c
      read(SXFILE_I, *)  xtalProp(11)
c
c------ (rss/crss)_limit to switch from (powerLaw & drag) to pure drag
c------  xtalProp(12) = rss/crss_limit
c
      call LimitRatioForDrag2(xtalProp)
c
c------- slip system hardening law : Voce's type
c-------   xtalProp(5) h0    : initial work hardening rate (Theta_0)
c-------   xtalProp(6) tausi : initial slip system strength (tau_0)
c-------   xtalProp(7) taus0 : saturation threshold strength (tau_s0)
c-------   xtalProp(8) xms   : strain-rate sensitivity - state evol (m')
c-------   xtalProp(9) gamss0: ref deformation rate - state evol 
c
      read(SXFILE_I, *)  (xtalProp(i), i = 5,9)
c      xtalProp(5) = props(11)
c      xtalProp(6) = props(12)
c      xtalProp(7) = props(13)
c      xtalProp(8) = props(14)
c      xtalProp(9) = props(15)
c
c------- initial value of state variables: slip system hardness crss0
c
      read(SXFILE_I, *) crss0          ! kappa0
      read(SXFILE_I, *)         
c      crss0 = props(16)
      xtalProp(10) = crss0
c
c------- echo date
c
      write(FILE_E, 1000) mode, xtalProp(3), xtalProp(4), xtalProp(11),
     &                    xtalProp(12)
      write(FILE_E, 2000) (xtalProp(i), i=5,9), crss0
c
c------- format
c
1000  format(/'*-----   Crystal Plasticity, mode', i3, '-----*'/,
     &        7x,'  Kinetic Equation : '/,
     &        7x,'  m                = ',e12.5/
     &        7x,'  gam0             = ',e12.5/
     &        7x,'  bdrag            = ',e12.5/
     &        7x,'  (rss/crss)_limit = ',e12.5)
2000  format( 7x,'  Hardening Law : '/,
     &        7x,'  h0               = ',e12.5/
     &        7x,'  tausi            = ',e12.5/
     &        7x,'  taus0            = ',e12.5/
     &        7x,'  xms              = ',e12.5/
     &        7x,'  gamss0           = ',e12.5/
     &        7x,'  crss0 (kappa0)   = ',e12.5)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetCrystalSlipGeometry2(
     &   crystalID, numslip, numvtx, tauSlip, scysFile, zBar0, 
     &   pBar0, qBar0, pBar0Vec, qBar0Vec, ppTBar0, FILE_E, filePath
     &   )
      
      implicit none
      include 'params_xtal.inc'

      character*160 filePath
      integer crystalID, numslip, numvtx, FILE_E
      real*8  tauSlip(MAX_SLIP)
      real*8  zBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0(DIMS, DIMS, MAX_SLIP), qBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0Vec(NVECS, MAX_SLIP), qBar0Vec(DIMS, MAX_SLIP)
      real*8  ppTBar0(NVECS, NVECS, MAX_SLIP)
      character*(*)  scysFile

      real*8  vecM(DIMS, MAX_SLIP), vecS(DIMS, MAX_SLIP)
c
c---------------------------------------------------------------------72
c
c------- Set up slip system vectors based on type of crystal
c
      if (crystalID .eq. kFCC) 
     &    call SetSlipSystemFCC2(numslip, vecM, vecS, scysFile, tauSlip,
     &                          FILE_E)

      if (crystalID .eq. kBCC) 
     &    call SetSlipSystemBCC2(numslip, vecM, vecS, scysFile, tauSlip,
     &                          FILE_E)
    

      if (crystalID .eq. kHCP) 
     &    call SetSlipSystemHCP2(numslip, numvtx, vecM, vecS, scysFile, 
     &                          tauSlip, FILE_E, filePath)
c
c------- Set up Schmid tensors/vectors for each slip system
c
      call SetSlipSystemTensors2(numslip, vecM, vecS, zBar0, pBar0, 
     &                          qBar0, pBar0Vec, qBar0Vec, ppTBar0,
     &                          FILE_E)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetSlipSystemFCC2(
     &   numslip, vecM, vecS, scysFile, tauSlip, FILE_E
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip, FILE_E
      real*8  vecM(DIMS, MAX_SLIP), vecS(DIMS, MAX_SLIP)
      real*8  tauSlip(MAX_SLIP)
      character*(*)  scysFile

      integer is, i
      real*8 indexM(3,12), indexS(3,12)
      real*8 sDotm(12)
      real*8 InnerProductVec2

      data indexM /1.,  1., -1.,
     &             1.,  1., -1.,
     &             1.,  1., -1.,
     &             1., -1., -1.,
     &             1., -1., -1.,
     &             1., -1., -1.,
     &             1., -1.,  1.,
     &             1., -1.,  1.,
     &             1., -1.,  1.,
     &             1.,  1.,  1.,
     &             1.,  1.,  1.,
     &             1.,  1.,  1./
      data indexS /0.,  1.,  1.,
     &             1.,  0.,  1.,
     &             1., -1.,  0.,
     &             0.,  1., -1.,
     &             1.,  0.,  1.,
     &             1.,  1.,  0.,
     &             0.,  1.,  1.,
     &             1.,  0., -1.,
     &             1.,  1.,  0.,
     &             0.,  1., -1.,
     &             1.,  0., -1.,
     &             1., -1.,  0./

c
c---------------------------------------------------------------------72
c
c------- set number of slip systems and slip system reference stress
c
      numslip = kSlipFCC
      call SetTensor2(tauSlip, pone, numslip)
c
c------- slip system normals and slip directions: unit vectors
c
      do is = 1, numslip
         call UnitVector2(indexS(1,is), vecS(1,is), DIMS)
         call UnitVector2(indexM(1,is), vecM(1,is), DIMS)
      enddo
c
c------- file with RI vertex stresses
c
      scysFile = 'vert_fcc.028.00'     ! not used
c
c------- check normality of vecS and vecM
c
      do is = 1, numslip
         sDotm(is) = InnerProductVec2(vecS(1,is), vecM(1,is), DIMS)
      enddo
c
c------- echo values
c
      write(FILE_E, 1000) numslip, scysFile
      do is = 1, numslip
         write(FILE_E, 2000) is, tauSlip(is), (vecS(i,is),i=1,3), 
     &                       (vecM(i,is),i=1,3), sDotm(is)
      enddo
c
c------- format
c
1000  format(/'*-----   Slip Systems for FCC -----*'/,
     &        7x,'  number slip syst = ',i4/
     &        7x,'  vtx stress file  = ',a15/
     &        7x,'  SS#',5x,'tau',18x,'vecS',30x,'vecM',20x,'SdotM')
2000  format( 7x, i4, 4x, f5.2, 4x, 3f10.5, 4x, 3f10.5, 4x, f10.5)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetSlipSystemBCC2(
     &   numslip, vecM, vecS, scysFile, tauSlip, FILE_E
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'
      
      integer numslip, FILE_E
      real*8  vecm(DIMS, MAX_SLIP), vecs(DIMS, MAX_SLIP)
      real*8  tauSlip(MAX_SLIP)
      character*(*)  scysFile

      integer is, i
      real*8 indexM(3,12), indexS(3,12)
      real*8 sDotm(12)
      real*8 InnerProductVec2

      data indexM /0.,  1.,  1.,
     &             1.,  0.,  1.,
     &             1., -1.,  0.,
     &             0.,  1., -1.,
     &             1.,  0.,  1.,
     &             1.,  1.,  0.,
     &             0.,  1.,  1.,
     &             1.,  0., -1.,
     &             1.,  1.,  0.,
     &             0.,  1., -1.,
     &             1.,  0., -1.,
     &             1., -1.,  0./
      data indexS /1.,  1., -1.,
     &             1.,  1., -1.,
     &             1.,  1., -1.,
     &             1., -1., -1.,
     &             1., -1., -1.,
     &             1., -1., -1.,
     &             1., -1.,  1.,
     &             1., -1.,  1.,
     &             1., -1.,  1.,
     &             1.,  1.,  1.,
     &             1.,  1.,  1.,
     &             1.,  1.,  1./
c
c---------------------------------------------------------------------72
c
c------- set number of slip systems and slip system reference stress
c
      numslip = kSlipBCC
      call SetTensor2(tauSlip, pone, numslip)
c
c------- slip system normals and slip directions: unit vectors
c
      do is = 1, numslip
         call UnitVector2(indexS(1,is), vecS(1,is), DIMS)
         call UnitVector2(indexM(1,is), vecM(1,is), DIMS)
      enddo
c
c------- file with RI vertex stresses
c
      scysFile = 'vert_bcc.028.00'     ! not used
c
c------- check normality of vecS and vecM
c
      do is = 1, numslip
         sDotm(is) = InnerProductVec2(vecS(1,is), vecM(1,is), DIMS)
      enddo
c
c------- echo values
c
      write(FILE_E, 1000) numslip, scysFile
      do is = 1, numslip
         write(FILE_E, 2000) is, tauSlip(is), (vecS(i,is),i=1,3), 
     &                       (vecM(i,is),i=1,3), sDotm(is)
      enddo
c
c------- format
c
1000  format(/'*-----   Slip Systems for BCC -----*'/,
     &        7x,'  number slip syst = ',i4/
     &        7x,'  vtx stress file  = ',a15/
     &        7x,'  SS#',5x,'tau',18x,'vecS',30x,'vecM',20x,'SdotM')
2000  format( 7x, i4, 4x, f5.2, 4x, 3f10.5, 4x, 3f10.5, 4x, f10.5)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetSlipSystemHCP2(
     &   numslip, numvtx, vecM, vecS, scysFile, tauSlip, 
     &   FILE_E, filePath
     &   )

      implicit  none
      include   'params_xtal.inc'

      character*160 filePath
      integer   numvtx, numslip, FILE_E
      real*8    vecM(DIMS, MAX_SLIP), vecS(DIMS, MAX_SLIP)
      real*8    tauSlip(MAX_SLIP)

      integer   nmodesx, nmodes, kount, modex, nsmx, mode(5)
      integer   i, j, is, ix, io, nm
      real*8    rca, tau
      real*8    vecm4(48, 4), vecs4(48, 4), sDotm(48)
      character prosa(9)*8
      character*(*) scysFile

      real*8 InnerProductVec2

      integer      lengthFile
      character*160 filename
c
c---------------------------------------------------------------------72
c
c     need input file 'slip_hcp.in'

      lengthFile = index(filePath,' ') - 1
      filename   = filePath(1:lengthFile)//'slip_hcp.in'

      io = 99
c      open(unit=io, file='slip_hcp.in', status='old')
      open(unit=io, file=filename, status='old')

      read(io,*) rca
      read(io,*) nmodesx
      read(io,*) nmodes
      read(io,*) (mode(i),i=1,nmodes)
      read(io,*) numvtx
      read(io,5) scysFile
5     format(3x,a15)
c      print *, scysFile
c      print *, '  nVertx for HCP = ', numvtx
c
      kount=1
      i=0
      do nm=1,nmodesx
         read(io,6) prosa
6        format(10a8)
         read(io,*) modex,nsmx,tau
         if(modex.ne.mode(kount)) then
            do ix=1,nsmx
               read(io,6) prosa
            enddo
         else
            kount=kount+1
            do is=1,nsmx
               i=i+1
               read(io,*) (vecm4(i,j),j=1,4),(vecs4(i,j),j=1,4)
               vecM(1,i)= vecm4(i,1)
               vecM(2,i)=(vecm4(i,1)+2.*vecm4(i,2))/sqrt(3.)
               vecM(3,i)= vecm4(i,4)/rca
               vecS(1,i)= 3./2.*vecs4(i,1)
               vecS(2,i)=(vecs4(i,1)/2.+vecs4(i,2))*sqrt(3.)
               vecS(3,i)= vecs4(i,4)*rca
               tauSlip(i) = tau
            enddo
         endif
      enddo

      close(io)
      numslip=i
c
c------- slip system normals and slip directions: unit vectors
c
      do is = 1, numslip
         call UnitVector2(vecS(1,is), vecS(1,is), DIMS)
         call UnitVector2(vecM(1,is), vecM(1,is), DIMS)
      enddo
c
c------- check normality of vecS and vecM
c
      do is = 1, numslip
         sDotm(is) = InnerProductVec2(vecS(1,is), vecM(1,is), DIMS)
      enddo
c
c------- echo values
c
      write(FILE_E, 1000) numslip, scysFile
      do is = 1, numslip
         write(FILE_E, 2000) is, tauSlip(is), (vecS(i,is),i=1,3), 
     &                       (vecM(i,is),i=1,3), sDotm(is)
      enddo
c
c------- format
c
1000  format(/'*-----   Slip Systems for HCP -----*'/,
     &        7x,'  number slip syst = ',i4/
     &        7x,'  vtx stress file  = ',a15/
     &        7x,'  SS#',5x,'tau',18x,'vecS',30x,'vecM',20x,'SdotM')
2000  format( 7x, i4, 4x, f5.2, 4x, 3f10.5, 4x, 3f10.5, 4x, f10.5)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetSlipSystemTensors2(
     &   numslip, vecM, vecS, zBar0, pBar0, qBar0, pBar0Vec, qBar0Vec,
     &   ppTBar0, FILE_E
     &   )

      implicit  none
      include  'params_xtal.inc'

      integer numslip, FILE_E
      real*8  vecM(DIMS, MAX_SLIP), vecS(DIMS, MAX_SLIP)
      real*8  zBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0(DIMS, DIMS, MAX_SLIP), pBar0Vec(NVECS, MAX_SLIP)
      real*8  qBar0(DIMS, DIMS, MAX_SLIP), qBar0Vec(DIMS, MAX_SLIP)
      real*8  ppTBar0(NVECS, NVECS, MAX_SLIP)

      integer is, j, k
c
c---------------------------------------------------------------------72
c
c------- Schmid Orientation Tensors/Vectors in Crystal Coordinates  
c
      do is = 1, numslip
c
c---------- Tensor zBar0
         call OuterProductVec2(vecS(1,is), vecM(1,is), zBar0(1,1,is), 
     &                        DIMS)
c
c---------- Symmetric and Skew parts of zBar0
         call SymmetrizeTensor2(zBar0(1,1,is), pBar0(1,1,is), DIMS)
         call SkewSymmetrizeTensor2(zBar0(1,1,is), qBar0(1,1,is), DIMS)
c
c---------- Vector form for pBar0 and qBar0
         call Mat3x3ToVec5x1Symm2(pBar0(1,1,is), pBar0Vec(1,is), DIMS)
         call Mat3x3ToVec3x1Skew2(qBar0(1,1,is), qbar0Vec(1,is), DIMS)
c
c---------- Form matrix {P}{P}^T
         call OuterProductVec2(pBar0Vec(1,is), pBar0Vec(1,is),
     &                        ppTBar0(1,1,is), NVECS)

      enddo
c
c------- echo Schmid orientation tensors/vectors
c
      write(FILE_E, 1000)
      do is = 1, numslip
         write(FILE_E, 1500) is
         do j = 1, DIMS
            write(FILE_E, 2000) (zBar0(j,k,is),k=1,DIMS), 
     &                          (pBar0(j,k,is),k=1,DIMS), 
     &                          (qBar0(j,k,is),k=1,DIMS)
         enddo
      enddo

      write(FILE_E, 3000)
      do is = 1, numslip
         write(FILE_E, 4000) is, (pBar0Vec(j,is),j=1,NVECS), 
     &                           (qBar0Vec(j,is),j=1,DIMS)
      enddo
c
c------- format
c
1000  format(/'*-----   Schmid Tensors -----*')
1500  format( 7x, ' SS # ', i2 / 
     &       10x, ' zBar0, pBar0, qBar0: ') 
2000  format(10x, 3(3x, 3f9.5))
3000  format(/'*-----   Schmid Vectors -----*'/
     &        7x, ' SS#, pBar0Vec, qBar0Vec: ') 
4000  format( 7x, i3, 2(3x, 5f9.5))

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetCrystalLatticeOrient2(
     &   numgrn, numqpt, numel, kODFout, gcrot0, props, mprops,
     &   kODF, numor, seed, angles, euler, FILE_I, FILE_E, FILE_O, 
     &   filePath, FILE_TXT_OUT
     &   )

      implicit none
      include 'params_xtal.inc'

      integer numgrn, numqpt, numel, kODFout
      real*8  gcrot0(DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)

      integer mprops
      real*8  props(mprops)

      character*160 filePath
      integer kODF, numor, seed, FILE_I, FILE_E, FILE_O
      integer FILE_TXT_OUT
      real*8  angles(DIMS, MAX_ORIEN)
      real*8  euler(DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)

      integer iikc, iidr, i, j
      real*8  pi, pi180, piby2, ph

      integer    length1, length2
      character  textureFile*20, filename*160
c
c---------------------------------------------------------------------72
c
      pi = 4.0 * datan(1.0d+00)
      pi180 = pi/180.
c
c------------------------------------- Data read from FILE_I
c
c------- read code to assign angles at each GR/IP/EL
c
      read(FILE_I, *) kODF
c      kODF = nint (props(17))
c
c------- read multiples of increment to output texture
c
      read(FILE_I, *) kODFout
c      kODFout = nint (props(18))
c
c------- read root name for input/output texture file
c
      read(FILE_I,'(a18)') textureFile
c      textureFile = 'texture'

      length1 = index(filePath,' ') - 1
      length2 = index(textureFile,' ') - 1

      filename = filePath(1:length1)//textureFile(1:length2)//'.txti'
      open(unit=XTAL_TXT_IN, file=filename, status='unknown',
     &                                   access='sequential')
      rewind(XTAL_TXT_IN)

      filename = filePath(1:length1)//textureFile(1:length2)//'.txto'
      open(unit=FILE_TXT_OUT, file=filename, status='unknown')
c
c------- echo data read from FILE_I
c
      write(FILE_E, 1000) kODF, kODFout, textureFile
c
c------------------------------------- Data read from XTAL_TXT_IN
c
c------- number of orientations in file
c
      read(XTAL_TXT_IN, *)  numor
      if (numor .lt. numgrn) 
     &   call RunTimeError2(FILE_O, 
     &                       'SetCrystalLatticeOrient: numor < numgrn')
c
c------- read the flag for angle convention and set some constants
c-------   iikc = 0 : angles input in Kocks convention :  (psi,the,phi)
c                 1 : angles input in Canova convention : (ph,th,om)
c                 ph = 90 + phi; th = the; om = 90 - psi
c-------   iidr = 0 : angles input in degrees
c-------          1 : angles input in radians
      read(XTAL_TXT_IN, *) iikc, iidr
      piby2 = 90.
      if (iidr .eq. 1) piby2 = pi / 2.0
c
c-------- read Euler angles in file & storage them in : (Kocks, radians)
c
      do i = 1, numor
         read(XTAL_TXT_IN, *) (angles(j,i), j=1,DIMS)

         if (iikc .eq. 1) then
            ph = angles(1,i)
            angles(1,i) = piby2 - angles(3,i)
            angles(3,i) = ph - piby2
         endif

         if (iidr .eq. 0) 
     &     call SetToScaledTensor2(pi180, angles(1,i), angles(1,i),DIMS)
      enddo

      if (kODF .eq. kODF_from_abq) then
         angles(1,1) = props(3)
         angles(2,1) = props(4)
         angles(3,1) = props(5)

         if (iikc .eq. 1) then
            ph = angles(1,1)
            angles(1,1) = piby2 - angles(3,1)
            angles(3,1) = ph - piby2
         endif

         if (iidr .eq. 0)
     &     call SetToScaledTensor2(pi180, angles(1,1), angles(1,1),DIMS)
      endif
c
c------- seed for a random angle distribution in polycrystal (needed?)
c
      read(XTAL_TXT_IN, *) seed
c
c------- close input file with texture info
c
      close(XTAL_TXT_IN)
c
c------------------------------------- Assign ODF to each GRs/IPs/ELs
c
c      call AssignCrystalODF(numel, numqpt)  ! No of args has changed!
c
c------- formats
c
1000  format(/'*-----   Lattice Orientation -----*'/
     &        7x,'  kODF             = ',i6/
     &        7x,'  kODFout          = ',i6/
     &        7x,'  root name for I/O txt file  = ',a18)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE VerticesSCYS_Cubic2(
     &   numvtx, kappa0, sigfs, scys_cub, FILE_E
     &   )

      implicit none
      include  'params_xtal.inc'
      include  'numbers.inc'

      integer numvtx, FILE_E
      real*8  kappa0, sigfs(5,MAX_VTX)
      character*(*)  scys_cub

      integer nvtx, ivtx, j
      real*8  sigca(5,28), strVtx(3,3)

      data nvtx   /28/
      data sigca  /-1.73204, -1.00000,   .00000,   .00000,   .00000,
     &              1.73204, -1.00000,   .00000,   .00000,   .00000,
     &               .00000,  2.00000,   .00000,   .00000,   .00000,
     &               .00000,   .00000,  1.73204,  1.73204,  1.73204,
     &               .00000,   .00000, -1.73204,  1.73204,  1.73204,
     &               .00000,   .00000,  1.73204, -1.73204,  1.73204,
     &               .00000,   .00000,  1.73204,  1.73204, -1.73204,
     &               .00000,   .00000,  3.46410,   .00000,   .00000,
     &               .00000,   .00000,   .00000,  3.46410,   .00000,
     &               .00000,   .00000,   .00000,   .00000,  3.46410,
     &              -.86602,  -.50000,   .00000,  1.73204,  1.73204,
     &              -.86602,  -.50000,   .00000, -1.73204, -1.73204,
     &              -.86602,  -.50000,   .00000,  1.73204, -1.73204,
     &              -.86602,  -.50000,   .00000, -1.73204,  1.73204,
     &               .86602,  -.50000,  1.73204,   .00000,  1.73204,
     &               .86602,  -.50000, -1.73204,   .00000, -1.73204,
     &               .86602,  -.50000, -1.73204,   .00000,  1.73204,
     &               .86602,  -.50000,  1.73204,   .00000, -1.73204,
     &               .00000,  1.00000,  1.73204,  1.73204,   .00000,
     &               .00000,  1.00000, -1.73204, -1.73204,   .00000,
     &               .00000,  1.00000,  1.73204, -1.73204,   .00000,
     &               .00000,  1.00000, -1.73204,  1.73204,   .00000,
     &               .86602, -1.50000,  1.73204,   .00000,   .00000,
     &               .86602, -1.50000, -1.73204,   .00000,   .00000,
     &               .86602,  1.50000,   .00000,  1.73204,   .00000,
     &               .86602,  1.50000,   .00000, -1.73204,   .00000,
     &             -1.73204,   .00000,   .00000,   .00000,  1.73204,
     &             -1.73204,   .00000,   .00000,   .00000, -1.73204/
c
c---------------------------------------------------------------------72
c
c------- sigca: vertices of single crystal yield surface
c------- {sigca}={-(11-22)/sqr2,sqr32*33,sqr2*32,sqr2*31,sqr2*21}

      call SetTensor2(strVtx, pzero, DIMS)

      numvtx = nvtx
      do ivtx = 1, numvtx

c         strVtx(1,1) = -(sigca(1,ivtx) + sigca(2,ivtx) / sqr3) / sqr2
c         strVtx(2,2) =  (sigca(1,ivtx) - sigca(2,ivtx) / sqr3) / sqr2
c         strVtx(3,3) = sigca(2,ivtx) * sqr2 / sqr3
c         strVtx(2,1) = sigca(5,ivtx) / sqr2
c         strVtx(3,1) = sigca(4,ivtx) / sqr2
c         strVtx(3,2) = sigca(3,ivtx) / sqr2

         sigfs(1, ivtx) = -sigca(1, ivtx)
         sigfs(2, ivtx) =  sigca(2, ivtx)
         sigfs(3, ivtx) =  sigca(5, ivtx)
         sigfs(4, ivtx) =  sigca(4, ivtx)
         sigfs(5, ivtx) =  sigca(3, ivtx)

c         call mat3x3ToVec5x1Symm(strVtx, sigfs(1,ivtx), DIMS)
         call SetToScaledTensor2(kappa0, sigfs(1,ivtx), sigfs(1,ivtx),5)

      enddo
c
c------- echo vettices
c
      write(FILE_E, 1000) numvtx
      do ivtx = 1, numvtx
         write (FILE_E, 2000) ivtx, (sigfs(j,ivtx),j=1,5)
      enddo
c
c------- formats
c
1000  format(/'*-----   Vertices of SCYS - Cubic Crystal -----*'/
     &        7x,'  No of vertices   = ',i4/
     &        7x,'  ivtx     sigfs(1...5)')
2000  format( 7x, i5, 3x, 5f10.5)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE VerticesSCYS_HCP2(
     &   numvtx, kappa0, sigfs, scys_hcp, FILE_E, filePath
     &   )

      implicit  none
      include   'params_xtal.inc'

      character*160 filePath
      integer   numvtx, FILE_E
      real*8    kappa0, sigfs(5,MAX_VTX)
      character*(*)      scys_hcp

      integer   io, ivtx, j

      integer      length1, length2
      character*160 filename
c
c---------------------------------------------------------------------72
c
c------- sigca: vertices of single crystal yield surface
c-------   {sigca}={(11-22)/sqr2,sqr32*33,sqr2*21,sqr2*31,sqr2*32}
c
      length1 = index(filePath,' ') - 1
      length2 = index(scys_hcp,' ') - 1

      filename = filePath(1:length1)//scys_hcp(1:length2)

      io = 99
c      open(io, file = scys_hcp, status='old')
      open(io, file = filename, status='old')

      do ivtx = 1, numvtx
         read(io,*) (sigfs(j,ivtx),j=1,5) 
         call SetToScaledTensor2(kappa0, sigfs(1,ivtx), sigfs(1,ivtx),5)
      enddo

      close(io)
c
c------- echo vettices
c
      write(FILE_E, 1000) numvtx
      do ivtx = 1, numvtx
         write (FILE_E, 2000) ivtx, (sigfs(j,ivtx),j=1,5)
      enddo
c
c------- formats
c
1000  format(/'*-----   Vertices of SCYS - HCP Crystal -----*'/
     &        7x,'  No of vertices   = ',i4/
     &        7x,'  ivtx   sigfs(1...5)')
2000  format( 7x, i5, 3x, 5f13.5)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE AssignCrystalODF2(numel, numqpt, noel, npt,
     &   numgrn, numslip, kODF,
     &   numor, seed,
     &   gcrot0,
     &   angles,
     &   euler,
     &   FILE_O, FILE_TXT_OUT
     &   )

      implicit none
      include  'params_xtal.inc'

      integer numel, numqpt, noel, npt
      integer FILE_O, FILE_TXT_OUT

      integer numgrn, numslip, kODF
      integer numor, seed
      real*8  gcrot0(DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  angles(DIMS, MAX_ORIEN)
      real*8  euler(DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)

      real    ran1_nr2

      integer ie, ip, ig, i, random
      real*8  pi, pi180

      integer seed_nr
      save    seed_nr
c
c---------------------------------------------------------------------72
c
      pi = 4.0 * datan(1.0d+00)
      pi180 = pi/180.

      if (noel .eq. 1 .and. npt .eq. 1) seed_nr = -seed
      if (noel .eq. 1 .and. npt .eq. 1) write(FILE_TXT_OUT, 2000)
c
c------- same ODF in all ELs and IPs
      if (kODF .eq. kODF_same_all) then
         do ie = 1, numel
            do ip = 1, numqpt
               do ig =1, numgrn
                  call EqualTensors2(angles(1,ig), 
     &                                        euler(1,ig,ip,ie), DIMS)
               enddo
            enddo
         enddo
c
c------- different ODF in all ELs; same ODF in all IPs
      elseif (kODF .eq. kODF_diff_elems) then
         do ie = 1, numel
            do ip = 1, numqpt
               random = int (ran1_nr2(seed_nr) * numor)
               do ig =1, numgrn
                  call EqualTensors2(angles(1,random), 
     &                                        euler(1,ig,ip,ie), DIMS)
               enddo
            enddo
         enddo
c
c------- different ODF in all ELs and IPs
      elseif (kODF .eq. kODF_diff_inpts) then
         do ie = 1, numel
            do ip = 1, numqpt
               do ig =1, numgrn
                  random = int (ran1_nr2(seed_nr) * numor)
                  call EqualTensors2(angles(1,random), 
     &                                        euler(1,ig,ip,ie), DIMS)
               enddo
            enddo
         enddo
c
c------- ODF read from file for all ELs, IPs, and GRs
      elseif (kODF .eq. kODF_from_file) then
c
c---------- for specific cases, e.g., BICRYSTAL  
c         if (numel .eq. numor .and. numgrn .eq. 1) then
c            do ie = 1, numel
c               do ip = 1, numqpt
c                  do ig =1, numgrn
c                     call EqualTensors(angles(1,ie), 
c     &                                        euler(1,ig,ip,ie), DIMS)
c                  enddo
c               enddo
c            enddo
cc
c---------- assign to each GR/IP/EL from file
c         elseif (numel*numqpt*numgrn .eq. numor) then
            i = 0
            do ie = 1, numel
               do ip = 1, numqpt
                  do ig =1, numgrn
                     i = i + 1
c                     call EqualTensors(angles(1,i), 
                     call EqualTensors2(angles(1,noel), 
     &                                        euler(1,ig,ip,ie), DIMS)
                  enddo
               enddo
            enddo
c         else
c            call RunTimeError(FILE_O, 
c     &        'setCrystalLatticeOrient: kODF_from_file, check Inp file')
c         endif
c
c------- read ODF from abaqus input file: 1 grain / IP
      elseif (kODF .eq. kODF_from_abq) then
         do ie = 1, numel
            do ip = 1, numqpt
               do ig =1, numgrn
                  call EqualTensors2(angles(1,ig),
     &                                        euler(1,ig,ip,ie), DIMS)
               enddo
            enddo
         enddo
      else
         call RunTimeError2(FILE_O, 'setCrystalLatticeOrient: Bad kODF')
      endif
c
c------- build rotation matrices C0: {x}_sm = [C0] {x}_cr 
c
      do ie = 1, numel
         do ip = 1, numqpt
            do ig =1, numgrn
               call AnglesToRotMatrix2(euler(1,ig,ip,ie), 
     &                                 gcrot0(1,1,ig,ip,ie), DIMS)
            enddo
         enddo
      enddo
c
c------- write initial assigned orientations
c-------  note that numel & numqpt are both dummy variables (=1)
c
c      write(*,*)
      do ie = 1, numel
         do ip = 1, numqpt
            do ig =1, numgrn
               write(FILE_TXT_OUT, 3000) 
     &                 (euler(i,ig,ip,ie)/pi180,i=1,3), ig, npt, noel
            enddo
         enddo
      enddo
c
c------- formats
c
2000  format(/'*-----   Initial Assigned Orientations -----*'/
     &        4x,'ang1',4x,'ang2',4x,'ang3',4x,' igrn',4x,'intpt',
     &        3x,'ielem ')
3000  format( 3f8.2, 3(3x,i5))

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetCrystalGeometry2(
     &   crystalID, numslip, numvtx, scysFile, zBar0, pBar0, qBar0, 
     &   pBar0Vec, qBar0Vec, ppTBar0, xtalProp, matProp, hardmtx, 
     &   kappa0, tauSlip, props, mprops, SXFILE_I, FILE_E
     &   )
      
      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      character*(*)  scysFile
      integer crystalID, numslip, numvtx, SXFILE_I, FILE_E
      real*8  kappa0(NKAPP)
      real*8  zBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0(DIMS, DIMS, MAX_SLIP), qBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0Vec(NVECS, MAX_SLIP), qBar0Vec(DIMS, MAX_SLIP)
      real*8  ppTBar0(NVECS, NVECS, MAX_SLIP)
      real*8  xtalProp(NPROPS), matProp(NPROPS, MAX_SLIP)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP), tauSlip(MAX_SLIP)

      integer mprops
      real*8  props(mprops)

      integer   nmodesx, nmodes, kount, modex, nsmx
      integer   i, j, is, js, nm, im, jm, isensex
      real*8    rca, crss0
      real*8    twshx, isectwx, thres1x, thres2x
      real*8    vecM(DIMS, MAX_SLIP), vecS(DIMS, MAX_SLIP)
      real*8    vecm4(48, 4), vecs4(48, 4), sDotm(48)
      character prosa*160
      common /SlipSystem/ vecM, vecS

      real*8 InnerProductVec2

      integer   NMFILE
      parameter (NMFILE=12)
      integer   mode(NMFILE), nsm(NMFILE)
      real*8    hselfx(NMFILE), hlatex(NMFILE, NMFILE)
c
c---------------------------------------------------------------------72
c
c---- zero out some arrays
c 
      call SetTensor2(vecM,  pzero, DIMS*MAX_SLIP)
      call SetTensor2(vecS,  pzero, DIMS*MAX_SLIP)
      call SetTensor2(vecm4, pzero, 48*4)
      call SetTensor2(vecs4, pzero, 48*4)
c
c---- read SXFILE_I to set slip / twinning data
c
      read(SXFILE_I,*) rca
      read(SXFILE_I,*) nmodesx               ! total # modes in file
      read(SXFILE_I,*) nmodes                ! # modes used (active)
      read(SXFILE_I,*) (mode(i),i=1,nmodes)  ! labels of active modes
      read(SXFILE_I,*) numvtx
      read(SXFILE_I,'(a15)') scysFile

      kount=1                         ! counter for active modes
      i=0                             ! counter for # slip/twin systems
      do nm=1,nmodesx

         read(SXFILE_I,'(a)') prosa
         read(SXFILE_I,*) modex,nsmx,isensex
         if(modex.ne.mode(kount)) then
c........... skip data
            read(SXFILE_I,'(a)') prosa
            read(SXFILE_I,'(a)') prosa
            read(SXFILE_I,'(a)') prosa
            read(SXFILE_I,'(a)') prosa
            read(SXFILE_I,'(a)') prosa
            read(SXFILE_I,'(a)') prosa
            do is=1,nsmx
               read(SXFILE_I,*)
            enddo
         else
c........... PlastProp: 4 lines; Twindata: 1 line; LatentHard: 1 line
            call SetCrystalPlasticity2(crystalID, xtalProp, crss0, 
     &                         props, mprops, modex, SXFILE_I, FILE_E)
            read(SXFILE_I,*) twshx, isectwx, thres1x, thres2x
            read(SXFILE_I,*) (hlatex(kount, jm), jm=1, nmodes)
            hselfx(kount) = 1.0
            nsm(kount) = nsmx

c........... indices for each slip/twin mode
            if (crystalID .eq. kFCC .or. crystalID .eq. kBCC) then
              do is=1,nsmx
                i=i+1
                kappa0(i) = crss0
                call EqualTensors2(xtalProp, matProp(1,i), NPROPS)
                read(SXFILE_I,*) (vecM(j,i),j=1,3),(vecS(j,i),j=1,3)
                call UnitVector2(vecS(1,i), vecS(1,i), DIMS)
                call UnitVector2(vecM(1,i), vecM(1,i), DIMS)
              enddo
            endif

            if (crystalID .eq. kHCP) then
              do is=1,nsmx
                i=i+1
                kappa0(i) = crss0
                call EqualTensors2(xtalProp, matProp(1,i), NPROPS)
                read(SXFILE_I,*) (vecm4(i,j),j=1,4),(vecs4(i,j),j=1,4)
                vecM(1,i)= vecm4(i,1)
                vecM(2,i)=(vecm4(i,1)+2.*vecm4(i,2))/sqrt(3.)
                vecM(3,i)= vecm4(i,4)/rca
                vecS(1,i)= 3./2.*vecs4(i,1)
                vecS(2,i)=(vecs4(i,1)/2.+vecs4(i,2))*sqrt(3.)
                vecS(3,i)= vecs4(i,4)*rca
                call UnitVector2(vecS(1,i), vecS(1,i), DIMS)
                call UnitVector2(vecM(1,i), vecM(1,i), DIMS)
              enddo
            endif
            kount=kount+1
         endif

      enddo

      numSlip=i
c
c------- ratio of slip system's kappa0: kappa0(is)/kappa0(1) 
c------- check normality of vecS and vecM
c
      do is = 1, numSlip
         tauSlip(is) = kappa0(is)/kappa0(1)
         sDotm(is)   = InnerProductVec2(vecS(1,is), vecM(1,is), DIMS)
      enddo
c
c------- set up latent hardening matrix
c
      i=0
      do im = 1, nmodes
        do is = 1, nsm(im)
          i=i+1
          j=0
          do jm = 1, nmodes
            do js = 1, nsm(jm)
              j=j+1
              hardmtx(i,j) = hlatex(im,jm)
            enddo
          enddo
          hardmtx(i,i)=hselfx(im)
        enddo
      enddo
c
c------- echo values
c
      write(FILE_E, 1000) crystalID, numslip, scysFile
      do is = 1, numslip
         write(FILE_E, 2000) is, tauSlip(is), (vecS(i,is),i=1,3),
     &                       (vecM(i,is),i=1,3), sDotm(is)
      enddo

      write(FILE_E,'(''*-----   Latent Hardening Matrix -----*'')')
      do is = 1, numSlip
         write(FILE_E, '(8x,24F5.1)') (hardmtx(is,js), js=1,numSlip)
      enddo
c
c------- Set up Schmid tensors/vectors for each slip system
c
      call SetSlipSystemTensors2(numslip, vecM, vecS, zBar0, pBar0,
     &                          qBar0, pBar0Vec, qBar0Vec, ppTBar0,
     &                          FILE_E)
c
c------- format
c
1000  format(/'*-----   Slip Systems for',i3,' (Crystal Type)-----*'/,
     &        7x,'  number slip syst = ',i4/
     &        7x,'  vtx stress file  = ',a15/
     &        7x,'  SS#',5x,'tau',18x,'vecS',30x,'vecM',20x,'SdotM')
2000  format( 7x, i4, 4x, f5.2, 4x, 3f10.5, 4x, 3f10.5, 4x, f10.5)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE stressCrystal2(
     &   sigxx, sigyy, sigth, sigxy, epsxx, epsyy, epsth, epsxy, 
     &   epsum, spinn, qptpo, qptp, dtime, time, ielem, incr, numel, 
     &   numqpt, iprint, sigma, ddsdde, statusFlag, numIncrs, noel, npt,
     &   dfgrd1, coords)

      implicit none
     
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer ielem, incr, numel, numqpt, iprint, statusFlag, numIncrs
      integer noel, npt
      real*8  dtime, time
      real*8  sigxx(numqpt), sigyy(numqpt), sigth(numqpt), sigxy(numqpt)
      real*8  epsxx(numqpt), epsyy(numqpt), epsth(numqpt)
      real*8  epsxy(numqpt), epsum(numqpt), spinn(numqpt)
      real*8  qptpo(numqpt), qptp(numqpt)
      real*8  sigma(DIMS, DIMS, NUMQPT_T)
      real*8  ddsdde(DIMS2, DIMS2, NUMQPT_T), dfgrd1(3,3)
      real*8  coords(3)
    
      integer iqpt
      real*8  thetao, theta, d_kk
      real*8  d_vec(NVECS), w_vec(3)

      real*8  epsxz(NUMQPT_T),  epsyz(NUMQPT_T)
      real*8  spinxz(NUMQPT_T), spinyz(NUMQPT_T)
      real*8  sigxz(NUMQPT_T),  sigyz(NUMQPT_T)
      common /for3DProbs/ epsxz, epsyz, spinxz, spinyz, sigxz, sigyz
c
c---------------------------------------------------------------------72
c
c------ initialize vectors for symm/skew parts of velocity gradients
c
      call SetTensor2(d_vec, pzero, NVECS)
      call SetTensor2(w_vec, pzero, 3)
c
c------ loop over integration points
c
      do iqpt = 1, numqpt
c
c---------- recover symm/skew parts of velocity gradient at ielem/iqpt
c---------- eps_ij are deviatoric quantities
c
         d_vec(1) = (epsxx(iqpt) - epsyy(iqpt)) / sqr2 
         d_vec(2) = epsth(iqpt) * sqr32
         d_vec(3) = epsxy(iqpt) * sqr2
         if (NSHR .gt. 1) then
            d_vec(4) = epsxz(iqpt) * sqr2
            d_vec(5) = epsyz(iqpt) * sqr2
         endif
         d_kk = epsum(iqpt)

         w_vec(1) = -spinn(iqpt)
         if (NSHR .gt. 1) then
            w_vec(2) = -spinxz(iqpt)
            w_vec(3) = -spinyz(iqpt)
         endif
c
c---------- recover temperature at ielem/iqpt
c
         thetao = qptpo(iqpt)
         theta  = qptp(iqpt)
c
c---------- evolve state at ielem/iqpt
c
         call DriverCrystalEvolve2(sigxx(iqpt), sigyy(iqpt),sigth(iqpt),
     &      sigxy(iqpt), sigxz(iqpt), sigyz(iqpt), d_vec, w_vec, d_kk,
     &      thetao, theta, dtime, time, iqpt, ielem, incr, numqpt, 
     &      numel, iprint, sigma(1, 1, iqpt), ddsdde(1, 1, iqpt),
     &      statusFlag, numIncrs, noel, npt, dfgrd1, coords)
         if (statusFlag .eq. kGLOBAL_FAILED) return

      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE DriverCrystalEvolve2(
     &   sigxx, sigyy, sigth, sigxy, sigxz, sigyz, d_vec, w_vec, d_kk,
     &   thetao, theta, dtime, time, iqpt, ielem, incr, numqpt, numel, 
     &   iprint, savg_ij, cavg_ijkl, statusFlag, numIncrs, noel, npt,
     &   dfgrd1, coords)

      implicit  none
      include  'params_xtal.inc'
      include  'numbers.inc'
c
c---- arguments
c
      integer iqpt, ielem, incr, numqpt, numel, iprint, statusFlag
      integer numIncrs, noel, npt
      real*8  sigxx, sigyy, sigth, sigxy, sigxz, sigyz
      real*8  thetao, theta, d_kk, dtime, time
      real*8  d_vec(NVECS), w_vec(DIMS), dfgrd1(3,3)
      real*8  cavg_ijkl(DIMS2, DIMS2), savg_ij(DIMS, DIMS)
      real*8  coords(3)
c
c---- variables passed through common blocks
c
      integer numgrn, numslip, numvtx, kODFout
      integer maxIterState, maxIterNewt
      real*8  tolerState, tolerNewt
      real*8  fCeDevVol(NVECS), fCeVol

      real*8  tauSlip(MAX_SLIP)
      real*8  zBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0(DIMS, DIMS, MAX_SLIP), qBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0Vec(NVECS, MAX_SLIP), qBar0Vec(DIMS, MAX_SLIP)
      real*8  ppTBar0(NVECS, NVECS, MAX_SLIP)
      real*8  sigfs(NVECS, MAX_VTX)
      real*8  fCeDev(NVECS, NVECS), fCeiDev(NVECS, NVECS)
      real*8  matProp(NPROPS, MAX_SLIP)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)

      real*8  gstress  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev  (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstress_n(NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran_n(NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa_n (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev_n(NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geqvalues(NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot_n  (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot_n  (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot0   (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps     (9, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps_n   (9, MAX_GRN, NUMQPT_T, NUMEL_T)
c
c---- local variables
c
      integer iterCounterS, iterCounterN, ierr, igrn, nxtals_nc 
      real*8  tmp, wpnorm_avg
      real*8  s_ij(DIMS,DIMS), ee_ij(DIMS,DIMS), tau(NVECS)
      real*8  stress(NVECS), estran(NVECS), kappa(NKAPP), statev(NSTAV)
      real*8  stress_n(NVECS), estran_n(NVECS), kappa_n(NKAPP)
      real*8  statev_n(NSTAV), eqvalues(NEQVA), gamdot(MAX_SLIP)
      real*8  crot(DIMS,DIMS), rrot(DIMS,DIMS), drot(DIMS,DIMS)
      real*8  crot_n(DIMS,DIMS), rrot_n(DIMS,DIMS), crot0(DIMS,DIMS)
      real*8  c_ijkl(DIMS2,DIMS2), sdev_ij(DIMS,DIMS)
      real*8  rss(MAX_SLIP)
      real*8  forden(MAX_SLIP), forden_n(MAX_SLIP)
      real*8  eps(9), eps_n(9)
      real*8  gndden(MAX_SLIP), gndden_n(MAX_SLIP)
      real*8  sub(MAX_SLIP), sub_n(MAX_SLIP)
      real*8  density, density_n, fdensity, fdensity_n
c
c---- common blocks
c
      common /XtalNumGrn/ numgrn, numslip
      common /XtalNumVtx/ numvtx
      common /XtalODFOut/ kODFout
      common /XtalStrVtx/ sigfs
      common /XtalProps / fCeDev, fCeiDev, matProp,
     &                    fCeDevVol, fCeVol, hardmtx
      common /XtalCrot0 / gcrot0
      common /XtalSlipG / tauSlip, zBar0, pBar0, qBar0,
     &                    pBar0Vec, qBar0Vec, ppTBar0
      common /XtalVars  / gstress, gestran, gkappa, gstatev, geqvalues,
     &                    gcrot, grrot, ggamdot, geps
      common /XtalVars_n/ gstress_n, gestran_n, gkappa_n, gstatev_n,
     &                    gcrot_n, grrot_n, geps_n

      common /nlcsolve1/  maxIterState, maxIterNewt
      common /nlcsolve2/  tolerState, tolerNewt

      real*8 time_end
      common /timeEnd/ time_end
      
      real*8  gforden  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gforden_n(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden_n(MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gsub     (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gsub_n   (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gdensity   (1, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gdensity_n (1, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gfdensity   (1, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gfdensity_n (1, MAX_GRN, NUMQPT_T, NUMEL_T)
      common  /Dislocation/  gforden, gforden_n, ggndden, ggndden_n,
     &                       gsub, gsub_n,
     &                    gdensity, gdensity_n, gfdensity, gfdensity_n
      
      integer numel_aba, numqpt_aba
      common /data_aba/ numel_aba, numqpt_aba
c
c---------------------------------------------------------------------72
c
c------- initialize average stress (sij) and consistent moduli (cijkl)
c
      call SetTensor2(savg_ij, pzero, DIMS*DIMS)
      call SetTensor2(cavg_ijkl, pzero, DIMS2*DIMS2)
c
c------- initialize average value of wp_norm
c
      wpnorm_avg = pzero
c
c------- counter for crystals that do not converge - Large Scale Appls.
      nxtals_nc = 0
c
c------- loop over all grains in aggregate at iqpt/ielem
c
      do igrn = 1, numgrn
c
c----------- zero local arrays
c
         call InitCrystalLocalArrays2(tau, stress, estran, kappa,statev,
     &      stress_n, estran_n, kappa_n, statev_n, eqvalues, gamdot,
     &      rss, crot, rrot, drot, crot_n, rrot_n, crot0, c_ijkl,
     &      forden, forden_n, eps, eps_n, gndden, gndden_n, sub, sub_n,
     &      density, density_n, fdensity, fdensity_n)
c
c---------- fetch variables at ielem/iqpt/igrn
c
         call FetchCrystalVariablesAtIP2(gstress_n, gestran_n, gkappa_n,
     &      gstatev_n, geqvalues, gcrot_n, grrot_n, gcrot0, stress_n, 
     &      estran_n, kappa_n, statev_n, eqvalues, crot_n, rrot_n, 
     &      crot0, igrn, iqpt, ielem, gforden_n, forden_n, geps_n,
     &      eps_n, ggndden_n, gndden_n, gsub_n, sub_n, gdensity_n,
     &      density_n, gfdensity_n, fdensity_n)
c
c---------- monitor integration/convergence of XTAL constitutive eqns
c
         call MonitorCrystalIntegration2(d_vec, w_vec, d_kk, s_ij,ee_ij,
     &      tau, 
     &      stress, estran, kappa, statev, eqvalues, gamdot, rss,
     &      crot, rrot,
     &      drot, stress_n, estran_n, kappa_n, statev_n, crot_n, rrot_n,
     &      crot0, fCeDev, fCeiDev, matProp, sigfs, tauSlip, zBar0, 
     &      pBar0, qBar0, pBar0Vec, qBar0Vec, ppTBar0, numslip, numvtx, 
     &      maxIterState, maxIterNewt, tolerState, tolerNewt, dtime, 
     &      fCeDevVol, fCeVol,
     &      theta, thetao, igrn, npt, noel, incr, iterCounterS, 
     &      iterCounterN, ierr, c_ijkl, hardmtx,
     &      eps, eps_n, dfgrd1, density, density_n,
     &      fdensity, fdensity_n)
c
c---------- check convergence; if not then:
c              either reset crystal quantities or force time reduction
c
         if (ierr .ne. XTAL_CONVERGED) then
            call WriteMessage2 (XTAL_O, 
     &                'DriverXtalEvolve: Sub-Incrs failed!')
            call WriteMessage2 (XTAL_O, 
     &                'DriverXtalEvolve: nxtals_nc incremented by one')
            nxtals_nc = nxtals_nc + 1 
            write(XTAL_O, 1000) nxtals_nc, incr, igrn, npt, noel

            if (nxtals_nc .ge. NXTALS_NC_LIMIT) then      
              call WriteMessage2 (XTAL_O, 
     &                'DriverXtalEvolve: nxtals_nc > LIMIT')
              call WriteMessage2(XTAL_O,
     &                '... will force dtime reduction by')
              call WriteMessage2(XTAL_O,
     &                '... setting: status=kGLOBAL_FAILED')
              statusFlag = kGLOBAL_FAILED
              return
            endif

            call WriteMessage2(XTAL_O, 'Resetting xtal quantities')
            call ResetCrystalQnts2(stress, estran, kappa, statev, 
     &             eqvalues, gamdot, crot, rrot, s_ij, c_ijkl,
     &             stress_n, 
     &             estran_n, kappa_n, statev_n, crot_n, rrot_n, 
     &             matProp, pBar0Vec, fCeiDev, fCeDevVol, fCeVol, 
     &             numslip, dtime, eps, eps_n, forden, forden_n,
     &             gndden, gndden_n, sub, sub_n, density, density_n,
     &             fdensity, fdensity_n)

         endif      
c
c--------- Output
c--------------------------------------------------------------
c          call Initialelementcoordinate2(coords, noel)
c          if(noel .eq. 1000) stop
c--------------------------------------------------------------
         if (ierr .eq. XTAL_CONVERGED) then
c             call UpdateCoordinates(coords, npt, noel) 
             
             call CalculateSSD(gamdot, dtime, numslip, forden,
     &                         forden_n, incr, sub, sub_n)
             
         if ( incr .gt. 1) then
c				call CalculateGND2(gamdot, numslip, npt, noel, gndden, 
c     &                        gndden_n, incr, dtime)
         endif
                 
      call UpdateKappa(kappa, forden, gndden, numslip, sub, incr, numgrn)
      
      call CalculateForestDensity(density, numslip, forden, incr, noel,
     &                            npt)
      
      call CalculateGNDDensity(fdensity, numslip, gndden, incr, noel,
     &                          npt)
             
             if (noel .eq. numel_aba .and. npt .eq. numqpt_aba) then
                 call UpdateFPandCoords()
c                 if (incr .eq. 1) call FindNeighbors()
             endif
         endif        
c
c---------- save computed variables at ielem/iqpt/igrn
c
         call SaveCrystalVariablesAtIP2(stress, estran, kappa, statev, 
     &      gamdot, eqvalues, crot, rrot, gstress, gestran, gkappa, 
     &      gstatev, geqvalues, ggamdot, gcrot, grrot, igrn, iqpt, 
     &      ielem, forden, gforden, eps, geps, gndden, ggndden, sub,
     &      gsub, density, gdensity, fdensity, gfdensity)
c
c---------- output computed quantities at selected "noel,npt,igrn"
c
c         if (iprint .eq. kPRINT_MODEL .and. 
c     &        noel  .eq. kPRINT_ELEM  .and.
c     &         npt  .eq. kPRINT_QPT   .and.
c     &         igrn .eq. kPRINT_GRN        ) then

c            call OutputQnts2(d_vec, s_ij, ee_ij, kappa, gamdot, rss,
c     &             statev, eqvalues, crot, rrot, c_ijkl, d_kk, dtime, 
c     &             time, numslip, incr, iterCounterS, iterCounterN, 
c     &             noel, npt, igrn)

c         endif
c
c---------- average values of sij and cijkl for aggregate

         tmp = 1.d0 / numgrn
         call AddScaledTensor2(tmp, s_ij, savg_ij, DIMS*DIMS)
         call AddScaledTensor2(tmp, c_ijkl, cavg_ijkl, DIMS2*DIMS2)
c
c---------- average value of wp_norm
c
         wpnorm_avg = wpnorm_avg + tmp*statev(kWPNORM)

      enddo
c
c------- move aggregate deviatoric stress to main-program arrays
c
      call DeviatoricTensor2(savg_ij, sdev_ij, DIMS)
      sigxx = sdev_ij(1,1)
      sigyy = sdev_ij(2,2)
      sigth = sdev_ij(3,3)
      sigxy = sdev_ij(1,2)
      if (NSHR .gt. 1) then
         sigxz = sdev_ij(1,3)
         sigyz = sdev_ij(2,3)
      endif
c
c------- output aggregate qnts at selected "ielem,iqpt"
c-------  note that ielem=1,iqpt=1 always !!
c
c      if (iprint .eq. kPRINT_MODEL .and. 
c     &     ielem .eq. kPRINT_ELEM  .and.
c     &     iqpt  .eq. kPRINT_QPT        ) then
c
c---------- write effective stress-strain curve
c
c         call WriteStressStrainCurve2(sdev_ij, d_vec,
c     &                               wpnorm_avg, dtime, time,
c     &                               numgrn, incr, npt, noel)
c
c---------- write texture at specified increments
c
         if ( (incr/kODFout*kODFout) .eq. incr .or. 
     &                       incr .eq. numIncrs) then
c     &           dabs((time+dtime)-time_end) .lt. 1e-6) then
            call WriteTexture2(gcrot, gkappa, d_vec, numgrn, incr, 
     &                        iqpt, ielem, npt, noel)

         endif

c      endif

1000  format(3x, 'Grains did not converged: ', i3, ':', 3x,
     &           'incr  #', i8, ';', 1x,
     &           'grain #', i5, ';', 1x,
     &           'iqpt  #', i3, ';', 1x,
     &           'elem  #', i5)

      return

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE InitCrystalLocalArrays2(
     &   tau, stress, estran, kappa, statev, stress_n, estran_n, 
     &   kappa_n, statev_n, eqvalues, gamdot, rss, crot, rrot, drot, 
     &   crot_n, rrot_n, crot0, c_ijkl, forden, forden_n, eps, eps_n,
     &   gndden, gndden_n, sub, sub_n, density, density_n, fdensity,
     &   fdensity_n)

      implicit  none
      include  'params_xtal.inc'
      include  'numbers.inc'

      real*8  tau(NVECS)
      real*8  stress(NVECS), estran(NVECS), kappa(NKAPP), statev(NSTAV)
      real*8  stress_n(NVECS), estran_n(NVECS), kappa_n(NKAPP)
      real*8  statev_n(NSTAV), eqvalues(NEQVA), gamdot(MAX_SLIP)
      real*8  crot(DIMS,DIMS), rrot(DIMS,DIMS), drot(DIMS,DIMS)
      real*8  crot_n(DIMS,DIMS), rrot_n(DIMS,DIMS), crot0(DIMS,DIMS)
      real*8  c_ijkl(DIMS2, DIMS2)
      real*8  rss(MAX_SLIP)
      real*8  forden(MAX_SLIP), forden_n(MAX_SLIP)
      real*8  eps(9), eps_n(9)
      real*8  gndden(MAX_SLIP), gndden_n(MAX_SLIP)
      real*8  sub(MAX_SLIP), sub_n(MAX_SLIP)
      real*8  density, density_n
      real*8  fdensity, fdensity_n
c
c---------------------------------------------------------------------72
c
      call SetTensor2(tau,      pzero, NVECS)
      call SetTensor2(stress,   pzero, NVECS)
      call SetTensor2(estran,   pzero, NVECS)
      call SetTensor2(kappa,    pzero, NKAPP)
      call SetTensor2(statev,   pzero, NSTAV)
      call SetTensor2(stress_n, pzero, NVECS)
      call SetTensor2(estran_n, pzero, NVECS)
      call SetTensor2(kappa_n,  pzero, NKAPP)
      call SetTensor2(statev_n, pzero, NSTAV)
      call SetTensor2(eqvalues, pzero, NEQVA)
      call SetTensor2(gamdot,   pzero, MAX_SLIP)
      call SetTensor2(crot,     pzero, DIMS*DIMS)
      call SetTensor2(rrot,     pzero, DIMS*DIMS)
      call SetTensor2(drot,     pzero, DIMS*DIMS)
      call SetTensor2(crot_n,   pzero, DIMS*DIMS)
      call SetTensor2(rrot_n,   pzero, DIMS*DIMS)
      call SetTensor2(crot0,    pzero, DIMS*DIMS)
      call SetTensor2(c_ijkl,   pzero, DIMS2*DIMS2)
      call SetTensor2(forden,   pzero, MAX_SLIP)
      call SetTensor2(forden_n, pzero, MAX_SLIP)
      call SetTensor2(eps,      pzero, 9)
      call SetTensor2(eps_n,    pzero, 9)
      call SetTensor2(gndden,   pzero, MAX_SLIP)
      call SetTensor2(gndden_n, pzero, MAX_SLIP)
      call SetTensor2(sub,   pzero, MAX_SLIP)
      call SetTensor2(sub_n, pzero, MAX_SLIP)
      call SetTensor2(density,   pzero, 1)
      call SetTensor2(density_n, pzero, 1)
      call SetTensor2(fdensity,   pzero, 1)
      call SetTensor2(fdensity_n, pzero, 1)

      call SetTensor2(rss,      pzero, MAX_SLIP)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE FetchCrystalVariablesAtIP2(
     &   gstress_n, gestran_n, gkappa_n, gstatev_n, geqvalues, 
     &   gcrot_n, grrot_n, gcrot0, stress_n, estran_n, kappa_n, 
     &   statev_n, eqvalues, crot_n, rrot_n, crot0, igrn, iqpt, 
     &   ielem, gforden_n, forden_n, geps_n, eps_n, ggndden_n, gndden_n,
     &   gsub_n, sub_n, gdensity_n, density_n, gfdensity_n, fdensity_n)

      implicit  none
      include  'params_xtal.inc'

      integer igrn, iqpt, ielem

      real*8  gstress_n (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran_n (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa_n  (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev_n (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geqvalues (NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot_n   (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot_n   (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot0    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gforden_n (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps_n    (9, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden_n (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gsub_n    (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gdensity_n (1, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gfdensity_n (1, MAX_GRN, NUMQPT_T, NUMEL_T)

      real*8  stress_n(NVECS), estran_n(NVECS), kappa_n(NKAPP)
      real*8  statev_n(NSTAV), eqvalues(NEQVA)
      real*8  crot_n(DIMS,DIMS), rrot_n(DIMS,DIMS), crot0(DIMS,DIMS)
      real*8  forden_n(MAX_SLIP), gndden_n(MAX_SLIP), sub_n(MAX_SLIP)
      real*8  eps_n(9), density_n, fdensity_n
c
c---------------------------------------------------------------------72
c
      call EqualTensors2(gstress_n(1,igrn,iqpt,ielem), stress_n, NVECS)
      call EqualTensors2(gestran_n(1,igrn,iqpt,ielem), estran_n, NVECS)
      call EqualTensors2(gkappa_n (1,igrn,iqpt,ielem), kappa_n,  NKAPP)
      call EqualTensors2(gstatev_n(1,igrn,iqpt,ielem), statev_n, NSTAV)
      call EqualTensors2(geps_n   (1,igrn,iqpt,ielem), eps_n, 9)

      call EqualTensors2(gcrot_n(1,1,igrn,iqpt,ielem), crot_n,DIMS*DIMS)
      call EqualTensors2(grrot_n(1,1,igrn,iqpt,ielem), rrot_n,DIMS*DIMS)
      call EqualTensors2(gcrot0(1,1,igrn,iqpt,ielem),  crot0, DIMS*DIMS)
      
      call EqualTensors2(gforden_n(1,igrn,iqpt,ielem),forden_n,MAX_SLIP)
      call EqualTensors2(ggndden_n(1,igrn,iqpt,ielem),gndden_n,MAX_SLIP)
      call EqualTensors2(gsub_n(1,igrn,iqpt,ielem),sub_n,MAX_SLIP)
      call EqualTensors2(gdensity_n(1,igrn,iqpt,ielem),density_n,1)
      call EqualTensors2(gfdensity_n(1,igrn,iqpt,ielem),fdensity_n,1)

      eqvalues(kEQP_n)    = geqvalues(kEQP_n,   igrn, iqpt, ielem)
      eqvalues(kMISES_n)  = geqvalues(kMISES_n, igrn, iqpt, ielem)
      eqvalues(kSHRATE_n) = geqvalues(kSHRATE_n,igrn, iqpt, ielem)
      eqvalues(kGAMTOT_n) = geqvalues(kGAMTOT_n,igrn, iqpt, ielem)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SaveCrystalVariablesAtIP2(
     &   stress, estran, kappa, statev, gamdot, eqvalues, crot, rrot, 
     &   gstress, gestran, gkappa, gstatev, geqvalues, ggamdot, gcrot, 
     &   grrot, igrn, iqpt, ielem, forden, gforden, eps, geps,
     &   gndden, ggndden, sub, gsub, density, gdensity, fdensity,
     &   gfdensity)

      implicit  none
      include  'params_xtal.inc'

      integer igrn, iqpt, ielem

      real*8  gstress  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev  (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geqvalues(NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gforden  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geps     (9, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggndden  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gsub     (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gdensity (1, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gfdensity (1, MAX_GRN, NUMQPT_T, NUMEL_T)

      real*8  stress(NVECS), estran(NVECS), kappa(NKAPP)
      real*8  statev(NSTAV), eqvalues(NEQVA), gamdot(MAX_SLIP)
      real*8  crot(DIMS,DIMS), rrot(DIMS,DIMS)
      real*8  forden(MAX_SLIP), gndden(MAX_SLIP), sub(MAX_SLIP)
      real*8  eps(9), density, fdensity
c
c---------------------------------------------------------------------72
c
      call EqualTensors2(stress, gstress(1,igrn,iqpt,ielem), NVECS)
      call EqualTensors2(estran, gestran(1,igrn,iqpt,ielem), NVECS)
      call EqualTensors2(kappa,  gkappa (1,igrn,iqpt,ielem), NKAPP)
      call EqualTensors2(statev, gstatev(1,igrn,iqpt,ielem), NSTAV)
      call EqualTensors2(gamdot, ggamdot(1,igrn,iqpt,ielem), MAX_SLIP)
      call EqualTensors2(crot, gcrot(1,1,igrn,iqpt,ielem), DIMS*DIMS)
      call EqualTensors2(rrot, grrot(1,1,igrn,iqpt,ielem), DIMS*DIMS)
      call EqualTensors2(forden, gforden(1,igrn,iqpt,ielem), MAX_SLIP)
      call EqualTensors2(eps, geps(1,igrn,iqpt,ielem), 9)
      call EqualTensors2(gndden, ggndden(1,igrn,iqpt,ielem), MAX_SLIP)
      call EqualTensors2(sub, gsub(1,igrn,iqpt,ielem), MAX_SLIP)
      call EqualTensors2(density, gdensity(1,igrn,iqpt,ielem), 1)
      call EqualTensors2(fdensity, gfdensity(1,igrn,iqpt,ielem), 1)

      geqvalues(kEQP,   igrn,iqpt,ielem) = eqvalues(kEQP)
      geqvalues(kMISES, igrn,iqpt,ielem) = eqvalues(kMISES)
      geqvalues(kSHRATE,igrn,iqpt,ielem) = eqvalues(kSHRATE)
      geqvalues(kGAMTOT,igrn,iqpt,ielem) = eqvalues(kGAMTOT)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE UpdateCrystalHistory2(
     &    numqpt, numel
     &   )

      implicit none
      include 'params_xtal.inc'
c
c---- arguments   
c
      integer numqpt, numel
c
c---- variables passed through common blocks
c
      integer numgrn, numslip
      real*8  gstress  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev  (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstress_n(NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran_n(NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa_n (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev_n(NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geqvalues(NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot_n  (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot_n  (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
c
c---- local variables
c
      integer ig, ip, ie
c
c---- common block
c
      common /XtalNumGrn/ numgrn, numslip
      common /XtalVars  / gstress, gestran, gkappa, gstatev, geqvalues,
     &                    gcrot, grrot, ggamdot
      common /XtalVars_n/ gstress_n, gestran_n, gkappa_n, gstatev_n,
     &                    gcrot_n, grrot_n
c
c---------------------------------------------------------------------72
c
c------- update variables (gamdot is not updated)
c
      do ie = 1, numel
        do ip = 1, numqpt
          do ig = 1, numgrn

            call EqualTensors2(gstress  (1,ig,ip,ie), 
     &                        gstress_n(1,ig,ip,ie), NVECS)
            call EqualTensors2(gestran  (1,ig,ip,ie), 
     &                        gestran_n(1,ig,ip,ie), NVECS)
            call EqualTensors2(gkappa   (1,ig,ip,ie),  
     &                        gkappa_n (1,ig,ip,ie), NKAPP)
            call EqualTensors2(gstatev  (1,ig,ip,ie), 
     &                        gstatev_n(1,ig,ip,ie), NSTAV)
            call EqualTensors2(gcrot    (1,1,ig,ip,ie), 
     &                        gcrot_n  (1,1,ig,ip,ie), DIMS*DIMS)
            call EqualTensors2(grrot    (1,1,ig,ip,ie), 
     &                        grrot_n  (1,1,ig,ip,ie), DIMS*DIMS)

            geqvalues(kEQP_n,   ig,ip,ie) = geqvalues(kEQP,   ig,ip,ie)
            geqvalues(kMISES_n, ig,ip,ie) = geqvalues(kMISES, ig,ip,ie)
            geqvalues(kSHRATE_n,ig,ip,ie) = geqvalues(kSHRATE,ig,ip,ie)
            geqvalues(kGAMTOT_n,ig,ip,ie) = geqvalues(kGAMTOT,ig,ip,ie)

          enddo
        enddo
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE ResetCrystalHistory2(
     &    numqpt, numel
     &   )

      implicit none
      include 'params_xtal.inc'
c
c---- arguments
c
      integer numqpt, numel
c
c---- variables passed through common blocks
c
      integer numgrn, numslip
      real*8  gstress  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev  (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstress_n(NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran_n(NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa_n (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev_n(NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geqvalues(NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot_n  (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot_n  (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
c
c---- local variables
c
      integer ig, ip, ie
c
c---- common block
c
      common /XtalNumGrn/ numgrn, numslip
      common /XtalVars  / gstress, gestran, gkappa, gstatev, geqvalues,
     &                    gcrot, grrot, ggamdot
      common /XtalVars_n/ gstress_n, gestran_n, gkappa_n, gstatev_n,
     &                    gcrot_n, grrot_n
c
c---------------------------------------------------------------------72
c
c------- reset variables (gamdot is not reset)
c
      do ie = 1, numel
        do ip = 1, numqpt
          do ig = 1, numgrn

            call EqualTensors2(gstress_n(1,ig,ip,ie), 
     &                        gstress  (1,ig,ip,ie), NVECS)
            call EqualTensors2(gestran_n(1,ig,ip,ie), 
     &                        gestran  (1,ig,ip,ie), NVECS)
            call EqualTensors2(gkappa_n (1,ig,ip,ie),  
     &                        gkappa   (1,ig,ip,ie), NKAPP)
            call EqualTensors2(gstatev_n(1,ig,ip,ie), 
     &                        gstatev  (1,ig,ip,ie), NSTAV)
            call EqualTensors2(gcrot_n  (1,1,ig,ip,ie), 
     &                        gcrot    (1,1,ig,ip,ie), DIMS*DIMS)
            call EqualTensors2(grrot_n  (1,1,ig,ip,ie), 
     &                        grrot    (1,1,ig,ip,ie), DIMS*DIMS)

            geqvalues(kEQP,   ig,ip,ie) = geqvalues(kEQP_n,   ig,ip,ie)
            geqvalues(kMISES, ig,ip,ie) = geqvalues(kMISES_n, ig,ip,ie)
            geqvalues(kSHRATE,ig,ip,ie) = geqvalues(kSHRATE_n,ig,ip,ie)
            geqvalues(kGAMTOT,ig,ip,ie) = geqvalues(kGAMTOT_n,ig,ip,ie)

          enddo
        enddo
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE OutputQnts2(
     &   d_vec, s_ij, ee_ij, kappa, gamdot, rss, statev, eqvalues, crot,
     &   rrot, c_ijkl, d_kk, dtime, time, numslip, incr, 
     &   iterCounterS, iterCounterN, ielem, iqpt, igrn
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip, incr, iterCounterS, iterCounterN
      integer ielem, iqpt, igrn
      real*8  d_kk, dtime, time
      real*8  d_vec(NTENS), s_ij(DIMS,DIMS), ee_ij(DIMS,DIMS)
      real*8  kappa(NKAPP), gamdot(MAX_SLIP), statev(NSTAV)
      real*8  eqvalues(NEQVA), crot(DIMS,DIMS), rrot(DIMS,DIMS)
      real*8  c_ijkl(DIMS2,DIMS2)
      real*8  rss(MAX_SLIP)
	  
      integer i, j
c      character*160 message
      real*8  strain_11, stress_11, d_eff, eqstran, deqstran
      real*8  d_ij(DIMS,DIMS)

      real*8  InnerProductVec2

      save    strain_11, eqstran
      data    strain_11, eqstran /0.d0, 0.d0/

      real*8  overstress(MAX_SLIP)
      common /XtalOverS/ overstress
	  
c
c----- NOTE: 'strain_11' makes sense only for coaxial loading
c
c---------------------------------------------------------------------72
c
c------- write headers in output files
c
      if (incr .eq. 1) then
         write(XTAL_STRESS_O,  800) ielem, iqpt, igrn
         write(XTAL_STRAIN_O,  900) ielem, iqpt, igrn
         write(XTAL_EFFSS_O,  1000) ielem, iqpt, igrn
         write(XTAL_TRUESS_O, 2000) ielem, iqpt, igrn
         write(XTAL_ITER_O,   3000) ielem, iqpt, igrn
      endif

c
c------- equivalent total strain and strain_11
c
      d_eff = dsqrt(2./3.*InnerProductVec2(d_vec, d_vec, NVECS))
      deqstran = dtime * d_eff
      eqstran  = eqstran + deqstran

      call Vec5x1ToMat3x3Symm2(d_vec, d_ij, DIMS)
      strain_11 = strain_11 + dtime*(d_ij(1,1) + pthird*d_kk)
      stress_11 = s_ij(1,1)
c
c------- write stress / consistent tangent
c
      write(XTAL_STRESS_O, 4000) incr, ((s_ij(i,j), j=i,3), i=1,3)
      write(XTAL_STRESS_O, 4500) ((c_ijkl(i,j), j=1,DIMS2), i=1,DIMS2)
c      write(XTAL_STRESS_O, 4500) (rss(i)/(kappa(1)+overstress(i)),
c     &                            i=1,numslip)
      write(XTAL_STRESS_O, 4500) (rss(i)/(kappa(i)), i=1,numslip)
c
c------- write elastic strains / crot / rotation tensors; gammadot
c
      write(XTAL_STRAIN_O, '(/i5)') incr
      do i = 1, 3
         write(XTAL_STRAIN_O, 5000) (ee_ij(i,j), j=1,3), 
     &                       (crot(i,j),j=1,3), (rrot(i,j), j=1,3)
      enddo

      write(XTAL_STRAIN_O, 5200) eqvalues(kEQP), statev(kWPNORM),
     &                           int( statev(kSSACT) )
      write(XTAL_STRAIN_O, 5500) (gamdot(i), i=1,numslip)
c
c------- write effective quantities
c
      write(XTAL_EFFSS_O, 5800) incr, int(statev(kSSACT)), dtime,
     &                time+dtime, d_eff, deqstran,
     &                statev(kWPNORM), eqvalues(kEQP),
     &                (eqvalues(kEQP)-eqvalues(kEQP_n))/dtime,
     &                eqvalues(kMISES)
c
c------- write true stress-strain curve (X-direction) (uniaxial/MPS)
c
      write(XTAL_TRUESS_O, 6000) incr, time, strain_11, stress_11
c
c------- write iteration counters and escalar variables
c
      write(XTAL_ITER_O, 7000) incr, iterCounterS, iterCounterN, 
     &                  statev(kEVOL), statev(kPRES), statev(kDETVe),
     &                  (kappa(i),i=1,numslip)
c     &                  (kappa(1)+overstress(i),i=1,numslip)
c
c------- formats
c
 800  format(' INCR', 15x, 'CAUCHY STRESS' ,/,
     &       20x,'CONSISTENT TANGENT ',/, 20x, 'RSS/KAPPA ',/,
     &       ' (elem # ', i5, ',  qpt # ', i2, ',  grn # ', i5, ')'/) 
 900  format(' INCR',/,18x,' EE_ij',32x,'Crot_ij',30x,'Rrot_ij',/, 
     &       19x, 'EQPS', 33x, 'WP_NORM', 30x, 'SS_ACTIV',/,
     &       48x, ' GAMMADOT(1,...,NUMSLIP)',/,38x,
     &       ' (elem # ', i5, ',  qpt # ', i2, ',  grn # ', i5, ')') 
1000  format(' INCR',3x,'SS_ACT',5x,'DTIME',8x,'TIME',8x,'D_EFF',10x,
     &      'DEQSTRAN',8x,'WP_NORM',7x,'EQP',10x,'DP_EFF',8x,'MISES',/,
     &       ' (elem # ', i5, ',  qpt # ', i2, ',  grn # ', i5, ')')
2000  format(' INCR',7x,'TIME',9x,'E_11',9x,'SIGMA_11',/,
     &       ' (elem # ', i5, ',  qpt # ', i2, ',  grn # ', i5, ')',/,
     &       ' (for uniaxial loading along X-direction / MPS runs)' ) 
3000  format(' INCR',1x,'IT-S',1x,'IT-N',6x,'EVOL', 8x,'PRESS', 
     &       8x,'DETVe',6x,'KAPPA(1) ... KAPPA(NKAPP)'/, 
     &       ' (elem # ', i5, ',  qpt # ', i2, ',  grn # ', i5, ')') 
4000  format(i8,2x,3(1x,e11.4)/22x,2(1x,e11.4)/34x,1(1x,e11.4)/)
4500  format(6(10x,6(1x,e11.4)/))
5000  format((3(2x,3(1x,e11.4))))
5200  format(/(15x,e11.4,27x,e11.4,28x,i5))
5500  format(/(22x,6(1x,e11.4)))
5800  format(i8,2x,i3,8(2x,e12.5))
6000  format(i8,8(2x,e12.5))
7000  format(i8,1x,i3,2x,i3,1x,3(2x,e11.4),18(1x,e11.4))

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MonitorCrystalIntegration2( 
     &   d_vec, w_vec, d_kk, s_ij, ee_ij, tau, stress, estran, kappa, 
     &   statev, eqvalues, gamdot, rss, crot, rrot, drot, stress_n, 
     &   estran_n,
     &   kappa_n, statev_n, crot_n, rrot_n, crot0, fCeDev, fCeiDev, 
     &   matProp, sigfs, tauSlip, zBar0, pBar0, qBar0, pBar0Vec,
     &   qBar0Vec, ppTBar0, numslip, numvtx, maxIterState, maxIterNewt, 
     &   tolerState, tolerNewt, dtime, fCeDevVol, fCeVol, 
     &   theta, thetao, igrn, iqpt, 
     &   ielem, incr, iterCounterS, iterCounterN, ierr, cepmod, hardmtx,
     &   eps, eps_n, dfgrd1, density, density_n, fdensity, fdensity_n)

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip, numvtx, maxIterState, maxIterNewt
      integer igrn, iqpt, ielem, incr, iterCounterS, iterCounterN, ierr
      real*8  tolerState, tolerNewt, dtime, theta, thetao, d_kk
      real*8  fCeDevVol(NVECS), fCeVol
      real*8  d_vec(NVECS), w_vec(DIMS), dfgrd1(3,3)
      real*8  s_ij(DIMS,DIMS), ee_ij(DIMS,DIMS), tau(NVECS)
      real*8  stress(NVECS), estran(NVECS), kappa(NKAPP)
      real*8  statev(NSTAV), eqvalues(NEQVA), gamdot(MAX_SLIP)
      real*8  crot(DIMS,DIMS), rrot(DIMS,DIMS), drot(DIMS,DIMS)
      real*8  stress_n(NVECS), estran_n(NVECS), kappa_n(NKAPP)
      real*8  statev_n(NSTAV), crot_n(DIMS,DIMS), rrot_n(DIMS,DIMS)
      real*8  crot0(DIMS,DIMS)
      real*8  fCeDev(NVECS,NVECS), fCeiDev(NVECS,NVECS)
      real*8  matProp(NPROPS, MAX_SLIP)
      real*8  sigfs(NVECS,MAX_VTX), tauSlip(MAX_SLIP)
      real*8  zBar0(DIMS,DIMS,MAX_SLIP)
      real*8  pBar0(DIMS,DIMS,MAX_SLIP), qBar0(DIMS,DIMS,MAX_SLIP)
      real*8  pBar0Vec(NVECS,MAX_SLIP), qBar0Vec(DIMS,MAX_SLIP)
      real*8  ppTBar0(NVECS,NVECS,MAX_SLIP)
      real*8  cepmod(DIMS2,DIMS2)
      real*8  rss(MAX_SLIP)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)
      real*8  eps(9), eps_n(9)
      real*8  density, density_n
      real*8  fdensity, fdensity_n

      integer subIncr, totSubIncrs
      real*8  dtime_tr, theta_tr, tmp, d_kk_n, d_kk_tr
      real*8  d_vec_tr(NVECS), w_vec_tr(DIMS)
      real*8  d_vec_n(NVECS), w_vec_n(DIMS), tau_save(NSTAV)
c
c---------------------------------------------------------------------72
c
c------- counters for sub-increments
c
      subIncr = 1
      totSubIncrs = 1
c
c------- integrate crystal constitutive equations
c
      iterCounterS = 0
      call IntegrateCrystalEqns2(d_vec, w_vec, d_kk, s_ij, ee_ij, tau, 
     &   stress, estran, kappa, statev, eqvalues, gamdot, rss,
     &   crot, rrot, 
     &   drot, stress_n, estran_n, kappa_n, statev_n, crot_n, rrot_n, 
     &   crot0, fCeDev, fCeiDev, matProp, sigfs, tauSlip, zBar0, 
     &   pBar0, qBar0, pBar0Vec, qBar0Vec, ppTBar0, numslip, numvtx, 
     &   maxIterState, maxIterNewt, tolerState, tolerNewt, dtime, 
     &   fCeDevVol, fCeVol,
     &   theta, thetao, iterCounterS, iterCounterN, incr, ierr, 
     &   cepmod, subIncr, totSubIncrs, hardmtx,
     &   eps, eps_n, dfgrd1, iqpt, ielem, density, density_n,
     &   fdensity, fdensity_n)
c
c------- if converged -> return, else do -> subincrementation
c------- NOTE: here, subincrementation improves the initial guess 
c-------       for the solution variables 'stress'
c
      if (ierr .eq. XTAL_CONVERGED) return

      call WriteMessage2
     &     (XTAL_O, 'MonitorCrystalIntegration: using Sub-Incrs!')
c
c------- compute velocity gradient at t_n
c
      call DeformationRate_n2(d_vec_n, w_vec_n, d_kk_n, iqpt)
c     MAY NEED Vec6 to Vec5 routine.........
c
c------- loop for sub-incrementation procedure
c
      do while (.true.)
c
c---------- if not converged, increase # of sub-increments 
         if (ierr .ne. XTAL_CONVERGED) then
            subIncr = 2 * subIncr - 1
            totSubIncrs = 2 * totSubIncrs
            if (totSubIncrs .gt. 128) then
               call WriteMessage2(XTAL_O,
     &                  'MonitorCrystalIntegration: totSubIncrs > 128')
               return
            endif
            call SetTensor2(tau, pzero, NVECS)
            if (subIncr .gt. 1) 
     &           call EqualTensors2(tau_save, tau, NVECS)
c
c---------- if converged, adjust subincrements
         else if (subIncr .lt. totSubIncrs) then
            if ( (subIncr/2*2) .eq. subIncr) then
               subIncr = subIncr / 2 + 1
               totSubIncrs = totSubIncrs / 2
            else
               subIncr = subIncr + 1
            endif
c
c---------- successful return for subincrementation
         else
            call WriteMessage2(XTAL_O, 'Sub-Incrs successful !!!')
            return
         endif
c
c---------- report sub-incrs status
         write(XTAL_O, 1000) igrn, iqpt, ielem, subIncr, totSubIncrs
c
c---------- trial time, velocity gradient and temperature (assumes 
c---------- linear variation of velocity gradient and temperature in dt)
         tmp = real(subIncr) / real(totSubIncrs)   
         dtime_tr = dtime * tmp
         theta_tr = (1.-tmp)*thetao + tmp*theta
         d_kk_tr  = (1.-tmp)*d_kk_n + tmp*d_kk

         call AddTensors2((1.-tmp), d_vec_n, tmp, d_vec, d_vec_tr,NVECS)
         call AddTensors2((1.-tmp), w_vec_n, tmp, w_vec, w_vec_tr, DIMS)
c
c---------- save current convergent solution before getting next solution
         if (subIncr .gt. 1 .and. ierr .eq. XTAL_CONVERGED)
     &                 call EqualTensors2(tau, tau_save, NVECS)
c
c---------- integrate crystal constitutive equations
         iterCounterS = 0
         call IntegrateCrystalEqns2(d_vec_tr, w_vec_tr, d_kk_tr, s_ij,
     &     ee_ij, tau, stress, estran, kappa, statev, eqvalues, gamdot, 
     &     rss,
     &     crot, rrot, drot, stress_n, estran_n, kappa_n, statev_n, 
     &     crot_n, rrot_n, crot0, fCeDev, fCeiDev, matProp, sigfs, 
     &     tauSlip, zBar0, pBar0, qBar0, pBar0Vec, qBar0Vec, ppTBar0, 
     &     numslip, numvtx, maxIterState, maxIterNewt, tolerNewt, 
     &     tolerState, dtime_tr, fCeDevVol, fCeVol,
     &     theta_tr, thetao, iterCounterS, iterCounterN, incr, 
     &     ierr, cepmod, subIncr, totSubIncrs, hardmtx,
     &     eps, eps_n, dfgrd1, iqpt, ielem, density, density_n,
     &     fdensity, fdensity_n)

      enddo

1000  format(3x, 'at grain #', i5, ';', 1x, 
     &           'iqpt #', i3, ';', 1x,
     &           'elem #', i5, ';', 1x,
     &           'subInc/totSubIncrs = ', i3, '/', i3)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE WriteTexture2(
     &   gcrot, gkappa, d_vec, numgrn, incr, iqpt, ielem, npt, noel
     &   )

      implicit none
      include  'params_xtal.inc'

      integer numgrn, incr, iqpt, ielem, npt, noel
      real*8  gcrot(DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gkappa(NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  d_vec(NTENS)

      integer igrn, i
      real*8  pi, pi180, d_eff
      real*8  angle(DIMS)

      real*8  InnerProductVec2
c
c---------------------------------------------------------------------72
c
      pi = 4.0 * datan(1.0d+00)
      pi180 = pi/180.
c
c-------  effective strain rate at the ip/elem
c
      d_eff = dsqrt(2./3.*InnerProductVec2(d_vec, d_vec, NVECS))
c
c------- output heading
c
      if (noel .eq. 1 .and. npt .eq. 1) 
     &            write(XTAL_TXT_OUT, 1000) incr
c
c------- output orientation
c
      do igrn = 1, numgrn

         call RotMatrixToAngles2(gcrot(1,1,igrn,iqpt,ielem), angle,DIMS)

         write(XTAL_TXT_OUT, 2000) (angle(i)/pi180, i=1,3), igrn,
     &                          npt, noel, gkappa(1,igrn,iqpt,ielem),
     &                          d_eff

      enddo
c
c------- formats
c
1000  format(/'*-----   Euler Angles at incr ', i8, ' -----*'/
     &        4x,'ang1',4x,'ang2',4x,'ang3',4x,' igrn',4x,'intpt',
     &        3x,'ielem ',4x,'kappa(1)',4x,'d_eff')
2000  format( 3f8.2, 3(3x,i5), 3x, f10.2, 6x, e12.2)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE WriteStressStrainCurve2(
     &   sdev_ij, d_vec, wpnorm_avg, dtime, time, numgrn, 
     &   incr, npt, noel
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numgrn, incr, npt, noel
      real*8  dtime, time, wpnorm_avg
      real*8  sdev_ij(DIMS,DIMS), d_vec(NTENS)

      integer denom
      real*8  s_eff, d_eff, eqstran, deqstran
      real*8  sdev_vec(NVECS)
      real*8  d_eff_avg, s_eff_avg

      real*8  InnerProductVec2

      save    eqstran
      data    eqstran /0.d0/

      save    d_eff_avg, s_eff_avg, deqstran

      integer numel_aba, numqpt_aba
      common /data_aba/ numel_aba, numqpt_aba
c
c---------------------------------------------------------------------72
c
c------- write headers in output file
c
      if (noel .eq. 1 .and. npt .eq. 1) then
         s_eff_avg = 0.0d0
         d_eff_avg = 0.0d0
         deqstran  = 0.0d0
         if (incr .eq. 1) 
     &        write(AGG_EFFSS_O, 800) numel_aba, numqpt_aba, numgrn
      endif
c
c------- effective stress
c
      call Mat3x3ToVec5x1Symm2(sdev_ij, sdev_vec, DIMS)
      s_eff = dsqrt(3./2.*InnerProductVec2(sdev_vec, sdev_vec, NVECS))

      s_eff_avg = s_eff_avg + s_eff
c
c-------  effective strain rate and total equivalent strain
c
      d_eff = dsqrt(2./3.*InnerProductVec2(d_vec, d_vec, NVECS))
      d_eff_avg = d_eff_avg + d_eff

      deqstran = deqstran + dtime * d_eff
c
c------- write stress/strain curve
c
      if (noel .eq. numel_aba .and. npt .eq. numqpt_aba) then
c         denom = numgrn * numel_aba * numqpt_aba
         denom = numel_aba * numqpt_aba
         s_eff_avg = s_eff_avg / denom
         d_eff_avg = d_eff_avg / denom
         deqstran = deqstran / denom
         eqstran  = eqstran + deqstran
         write(AGG_EFFSS_O, 1000) incr, dtime, time+dtime, d_eff_avg, 
     &           deqstran, eqstran, wpnorm_avg, s_eff_avg
      endif
c
c------- formats
c
800   format(' INCR',9x,'DTIME',9x,'TIME',9x,'D_EFF',8x,'DEQSTRAN',
     &       8x,'EQSTRAN',8x,'WP_NORM',6x,'S_EFF'/,
     &       ' (nelem: ', i5, ',  nqpts: ', i2, ',  ngrns: ', i5, ')') 
1000  format(i8,10(2x,e12.5))

      return
      END
c
c=====================================================================72
c
      SUBROUTINE ResetCrystalQnts2(
     &   stress, estran, kappa, statev, eqvalues, gamdot, crot, rrot, 
     &   s_ij, c_ijkl, stress_n, estran_n, kappa_n, statev_n, crot_n, 
     &   rrot_n, 
     &   matProp, pBar0Vec, fCeiDev, fCeDevVol, fCeVol, numslip, dtime,
     &   eps, eps_n, forden, forden_n, gndden, gndden_n, sub, sub_n,
     &   density, density_n, fdensity, fdensity_n)

      implicit  none
      include  'params_xtal.inc'
      include  'numbers.inc'

      integer numslip
      real*8  dtime

      real*8  stress(NVECS), estran(NVECS), kappa(NKAPP)
      real*8  statev(NSTAV), eqvalues(NEQVA), gamdot(MAX_SLIP)
      real*8  crot(DIMS,DIMS), rrot(DIMS,DIMS)
      real*8  s_ij(DIMS,DIMS), c_ijkl(DIMS2,DIMS2)

      real*8  stress_n(NVECS), estran_n(NVECS), kappa_n(NKAPP)
      real*8  statev_n(NSTAV), crot_n(DIMS,DIMS), rrot_n(DIMS,DIMS)

      real*8  matProp(NPROPS, MAX_SLIP)
      real*8  pBar0Vec(NVECS,MAX_SLIP)
      real*8  fCeiDev(NVECS,NVECS), fCeDevVol(NVECS), fCeVol
      real*8  eps(9), eps_n(9)
      real*8  forden(MAX_SLIP), forden_n(MAX_SLIP)
      real*8  gndden(MAX_SLIP), gndden_n(MAX_SLIP)
      real*8  sub(MAX_SLIP), sub_n(MAX_SLIP)
      real*8  density, density_n
      real*8  fdensity, fdensity_n
c
c---- Local variables
c
      integer is
      real*8  qr5x5(NVECS,NVECS)
      real*8  pHatVec(NVECS,MAX_SLIP), ppTHat(NVECS,NVECS,MAX_SLIP)
      real*8  fCeiDevHat(NVECS,NVECS),fCeDevVolHat(NVECS)

      real*8  crss
      real*8  tau(NVECS)
      real*8  rss(MAX_SLIP), dGamdTau(MAX_SLIP)

      real*8  InnerProductVec2, SSKineticEqn2
c
c---- Common Blocks
c
      real*8  overstress(MAX_SLIP)
      common /XtalOverS/ overstress

      real*8  sGamPPt(NVECS,NVECS)
      common /SumGamPPt/ sGamPPt
c
c---------------------------------------------------------------------72
c
c---- Reset crystal quantities
c
      call EqualTensors2(stress_n, stress, NVECS)
      call EqualTensors2(estran_n, estran, NVECS)
      call EqualTensors2(kappa_n , kappa , NKAPP)
      call EqualTensors2(statev_n, statev, NSTAV)
      call EqualTensors2(crot_n  , crot  , DIMS*DIMS)
      call EqualTensors2(rrot_n  , rrot  , DIMS*DIMS)
      call EqualTensors2(eps_n, eps, 9)
      call EqualTensors2(forden_n, forden, MAX_SLIP)
      call EqualTensors2(gndden_n, gndden, MAX_SLIP)
      call EqualTensors2(sub_n, sub, MAX_SLIP)
      call EqualTensors2(density_n, density, 1)
      call EqualTensors2(fdensity_n, fdensity, 1)
                                                                        
      eqvalues(kEQP)    = eqvalues(kEQP_n)
      eqvalues(kMISES)  = eqvalues(kMISES_n)             
      eqvalues(kSHRATE) = eqvalues(kSHRATE_n)
      eqvalues(kGAMTOT) = eqvalues(kGAMTOT_n)

      call SetToScaledtensor2(statev(kDETVe), stress, tau, NVECS)
c
c------- Cauchy stress (tensor form)
      call Vec5x1ToMat3x3Symm2(stress, s_ij, DIMS)
      call AddScaledTensor2(statev(kPRES), Ident2nd, s_ij, DIMS*DIMS)
c
c---- Compute Approximate Consistent Tangent
c
c------- Rotation tensor 
      call RotMat5x5ForSymm2(crot, qr5x5, DIMS)
c
c------- Compute gamdot and sGamPPt
      call SetTensor2(sGamPPt, pzero, NVECS*NVECS)
      do is = 1, numslip
c
c------- Rotate symmetric part of Schimdt tensor
         call MultAxu2(qr5x5, pBar0Vec(1,is), pHatVec(1,is), NVECS)
         call OuterProductVec2(pHatVec(1,is), pHatVec(1,is),
     &                        ppTHat(1,1,is), NVECS)
c
c------- Resolve shear stresses
         rss(is) = InnerProductVec2(tau, pHatVec(1,is), NVECS)
c         crss = kappa(1) + overstress(is)
         crss = kappa(is)
c
c------- Shear strain rate
         gamdot(is) = SSKineticEqn2(rss(is),crss,matProp(1,is),kGAMDOT)
c
c------- Derivative of shear strain rate
         dGamdTau(is) = SSKineticEqn2(rss(is), crss, matProp(1,is),
     &                                   kdGAMdTAU)
         call AddScaledTensor2(dtime*dGamdTau(is), ppTHat(1,1,is),
     &                        sGamPPt, NVECS*NVECS)

      enddo
c
c------- Rotate elasticities
      call MultQAQT2(qr5x5, fCeiDev, fCeiDevHat, NVECS)
      call MultAxu2(qr5x5, fCeDevVol, fCeDevVolHat, NVECS)
c
c------- Approximate consistent tangent
      call PlasticModuli2(c_ijkl, fCeiDevHat, fCeDevVolHat, fCeVol,
     &                   statev(kDETVe))

      return
      end
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE IntegrateCrystalEqns2(
     &   d_vec, w_vec, d_kk, s_ij, ee_ij, tau, stress,
     &   estran, kappa, statev, eqvalues, gamdot, rss, crot, rrot, drot,
     &   stress_n, estran_n, kappa_n, statev_n, crot_n, rrot_n, 
     &   crot0, fCeDev, fCeiDev, matProp, sigfs, tauSlip, zBar0, 
     &   pBar0, qBar0, pBar0Vec, qBar0Vec, ppTBar0, numslip, numvtx, 
     &   maxIterState, maxIterNewt, tolerState, tolerNewt, dtime, 
     &   fCeDevVol, fCeVol,
     &   theta, thetao, iterCounterS, iterCounterN, globalIncr, ierr, 
     &   cepmod, subIncr, totSubIncrs, hardmtx,
     &   eps, eps_n, dfgrd1, iqpt, ielem, density, density_n,
     &   fdensity, fdensity_n)

      implicit none
      include 'params_xtal.inc' 
      include 'numbers.inc'

      integer numslip, numvtx, maxIterState, maxIterNewt, globalIncr
      integer iterCounterS, iterCounterN, ierr, subIncr, totSubIncrs
      integer iqpt, ielem
      real*8  tolerState, tolerNewt, dtime, theta, thetao, d_kk
      real*8  fCeDevVol(NVECS), fCeVol
      real*8  d_vec(NVECS), w_vec(DIMS)
      real*8  s_ij(DIMS,DIMS), ee_ij(DIMS,DIMS), tau(NVECS)
      real*8  stress(NVECS), estran(NVECS), kappa(NKAPP)
      real*8  statev(NSTAV), eqvalues(NEQVA), gamdot(MAX_SLIP)
      real*8  crot(DIMS,DIMS), rrot(DIMS,DIMS), drot(DIMS,DIMS)
      real*8  stress_n(NVECS), estran_n(NVECS), kappa_n(NKAPP)
      real*8  statev_n(NSTAV), crot_n(DIMS,DIMS), rrot_n(DIMS,DIMS)
      real*8  crot0(DIMS,DIMS)
      real*8  fCeDev(NVECS,NVECS), fCeiDev(NVECS,NVECS)
      real*8  matProp(NPROPS,MAX_SLIP)
      real*8  sigfs(NVECS,MAX_VTX), tauSlip(MAX_SLIP)
      real*8  zBar0(DIMS,DIMS,MAX_SLIP)
      real*8  pBar0(DIMS,DIMS,MAX_SLIP), qBar0(DIMS,DIMS,MAX_SLIP)
      real*8  pBar0Vec(NVECS,MAX_SLIP), qBar0Vec(DIMS,MAX_SLIP)
      real*8  ppTBar0(NVECS,NVECS,MAX_SLIP)
      real*8  cepmod(DIMS2,DIMS2)
      real*8  rss(MAX_SLIP)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)
      real*8  eps(9), eps_n(9)
      real*8  dfgrd1(3,3)
      real*8  density, density_n
      real*8  fdensity, fdensity_n

      integer is, iterState
      logical converged2
      real*8  e_kk, norm_tau0, norm_kapp0, norm_tau, norm_kapp, epsdot
      real*8  d_vec_lat(NVECS), tau_lat(NVECS)
      real*8  qr5x5(NVECS,NVECS), qr3x3(DIMS,DIMS)
      real*8  pHatVec(NVECS,MAX_SLIP), qHatVec(DIMS,MAX_SLIP)
      real*8  ppTHat(NVECS,NVECS,MAX_SLIP)
      real*8  fCeiDevHat(NVECS,NVECS),fCeDevVolHat(NVECS)

      logical ConvergeState2
      real*8  InnerProductVec2
c
c---------------------------------------------------------------------72
c
c------------------------- DEVIATORIC RESPONSE
c
c------- effective deformation rate 
c
      epsdot = dsqrt(2./3.*InnerProductVec2(d_vec, d_vec, NVECS))
      epsdot = max(epsdot, TINY)
c
c------- rotate Scmidt tensors (vectors) to Btilde_n configuration
c
      call RotMat5x5ForSymm2(crot_n, qr5x5, DIMS)
      call RotMat3x3ForSkew2(crot_n, qr3x3, DIMS)
      do is = 1, numslip
         call MultAxu2(qr5x5, pBar0Vec(1,is), pHatVec(1,is), NVECS)
         call MultAxu2(qr3x3, qBar0Vec(1,is), qHatVec(1,is), DIMS)
      enddo
c
c------- initialze slip hardening at time t
c
      call EqualTensors2(kappa_n, kappa, NKAPP)
c
c------- initial estimate for the stress
c
      if (subIncr .eq. 1) then

         if (globalIncr .eq. 1) then
c
c------------- initial guess from viscoplastic solution
c---------------- transform D from sample to crystal coordinates
            call MultATxu2(qr5x5, d_vec, d_vec_lat, NVECS)
c
c---------------- viscoplastic solution in crystal coordinates
            call StressSolveViscoPlastic2(tau_lat, d_vec_lat, kappa,
     &              pBar0Vec, ppTBar0, matProp, sigfs, tauSlip, dtime, 
     &              theta, thetao, epsdot, tolerNewt, maxIterNewt, 
     &              numslip, numvtx)
c
c---------------- transform computed stresses from crystal to sample axis
            call MultAxu2(qr5x5, tau_lat, tau, NVECS)
         else
c
c------------- initial guess from previous solution
            call EqualTensors2(stress_n, tau, NVECS)
         endif

      endif
c
c------- initial estimate for hardness and rotation tensor
c
      call HardeningRotationSolve2(kappa, kappa_n, eqvalues, rrot, 
     &        rrot_n, crot, crot0, drot, tau, w_vec, gamdot, rss,
     &        pHatVec, 
     &        qHatVec, matProp, hardmtx, epsdot, dtime, numslip, 
     &        kHARD_EXPL)
c
c------- initial valuies for the norm of stress and kappa
c
      norm_tau0  = dsqrt(InnerProductVec2(tau, tau, NVECS))
      norm_kapp0 = dsqrt(InnerProductVec2(kappa, kappa, NKAPP))
c
c------- initialized global flag to monitor Newton/State convergence
c
      ierr = XTAL_CONVERGED
c
c------- predictor volumetric elastic strain in Btilde
c
      e_kk = statev_n(kEVOL) + dtime * d_kk
c
c------- iterate for the material state
c
      iterState    = 0
      iterCounterN = 0
      converged2 = .false.

      do while(iterState .lt. maxIterState  .and. .not.converged2)

         iterState = iterState + 1
c
c---------- rotate Scmidt tensors to Btilde configuration
c
         call RotateSlipGeometry2(crot, pBar0Vec, qBar0Vec, qr5x5,qr3x3,
     &                           pHatVec, qHatVec, ppTHat, numslip)
c
c---------- rotate deviatoric & dev/volumetric elasticities to Btilde:
c
         call MultQAQT2(qr5x5, fCeiDev, fCeiDevHat, NVECS)
         call MultAxu2(qr5x5, fCeDevVol, fCeDevVolHat, NVECS)
c
c---------- solve for the crystal stresses
c
         call StressSolveDeviatoric2(tau, d_vec, estran, estran_n,kappa,
     &           gamdot, rss, drot, fCeiDevHat, fCeDevVolHat, pHatVec,
     &           ppTHat, matProp, e_kk, dtime, tolerNewt, numslip, 
     &           maxIterNewt, iterCounterN, ierr)
         if (ierr .ne. XTAL_CONVERGED) return
         norm_tau = dsqrt(InnerProductVec2(tau, tau, NVECS))
c
c---------- solve for hardness and rotation tensor
c
         call HardeningRotationSolve2(kappa, kappa_n, eqvalues, rrot, 
     &           rrot_n, crot, crot0, drot, tau, w_vec, gamdot, rss,
     &           pHatVec,
     &           qHatVec, matProp, hardmtx, epsdot, dtime, numslip, 
c     &           kHARD_MIDP)
     &           kHARD_ANAL)
         norm_kapp = dsqrt(InnerProductVec2(kappa, kappa, NKAPP))
c
c---------- check convergence
c
         converged2 = ConvergeState2(norm_tau, norm_kapp, norm_tau0,
     &                             norm_kapp0, tolerState)

      enddo
c
c------- keep track of state iteration and check number of iterations
c
      iterCounterS = iterState
      if (iterState .ge. maxIterState) then
         call WriteWarning2(XTAL_O,
     &            'StressSolveElastoViscoPlastic: iters > maxIters')
         ierr = XTAL_MAX_ITERS_HIT
         return
      endif
c
c------- continue sub-incrementation process if needed
c
      if (subIncr .lt. totSubIncrs) return
c
c------------------------ VOLUMETRIC RESPONSE
c
c------- integrate crystal volumetric  response
c
c      call RotMat5x5ForSymm(crot, qr5x5, DIMS)
c      call MultAxu(qr5x5, fCeDevVol, fCeDevVolHat, NVECS)
      call StressSolveVolumetric2(statev, fCeDevVolHat, fCeVol, estran, 
     &                           e_kk, dtime)
c
c----------------------- CAUCHY STRESS - CONSISTENT TANGENT
c
c------- compute Cauchy stress s_ij (note that "tau" is deviatoric)
c-------  s_ij = stress + statev(kPRES)*I  (i.e., stress is deviatoric)
c-------  ee_ij = estran + statev(KEVOL)*I  (i.e., estran is deviatoric)
c
      call UpdateStress2(s_ij, stress, tau, ee_ij, estran, statev, 
     &                  eqvalues, rrot, dfgrd1, iqpt, ielem)
c
c------- compute norm of spin, slip system activity, and eqp
c
      call UpdateDeformQnts2(gamdot, eqvalues, statev, crot, zBar0,
     &                      dtime, numslip, eps, eps_n)

c      call RotMat5x5ForSymm(crot, qr5x5, DIMS)
c      call MultQAQT(qr5x5, fCeiDev, fCeiDevHat, NVECS)
c      call MultAxu(qr5x5, fCeDevVol, fCeDevVolHat, NVECS)

      call PlasticModuli2(cepmod, fCeiDevHat, fCeDevVolHat, fCeVol,
     &                   statev(kDETVe))

      return
      END

c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE RotateSlipGeometry2(
     &   crot, pBar0Vec, qBar0Vec, qr5x5, qr3x3, pHatVec, qHatVec,
     &   ppTHat, numslip
     &   )

      implicit none
      include 'params_xtal.inc'

      integer numslip
      real*8  crot(DIMS,DIMS)
      real*8  pBar0Vec(NVECS,MAX_SLIP), qBar0Vec(DIMS,MAX_SLIP)
      real*8  qr5x5(NVECS,NVECS), qr3x3(DIMS,DIMS)
      real*8  pHatVec(NVECS,MAX_SLIP), qHatVec(DIMS,MAX_SLIP)
      real*8  ppTHat(NVECS,NVECS,MAX_SLIP)

      integer is
c
c---------------------------------------------------------------------72
c
c------- rotation matrices for symm and skew quantities
c
      call RotMat5x5ForSymm2(crot, qr5x5, DIMS)
      call RotMat3x3ForSkew2(crot, qr3x3, DIMS)
c
c------- rotate Scmidt tensors (vectors) to Btilde configuration
c
      do is = 1, numslip
         call MultAxu2(qr5x5, pBar0Vec(1,is), pHatVec(1,is), NVECS)
         call MultAxu2(qr3x3, qBar0Vec(1,is), qHatVec(1,is), DIMS)
         call OuterProductVec2(pHatVec(1,is), pHatVec(1,is), 
     &                        ppTHat(1,1,is), NVECS)
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE HardeningRotationSolve2(
     &   kappa, kappa_n, eqvalues, rrot, rrot_n, crot, crot0, drot, 
     &   tau, w_vec, gamdot, rss, pHatVec, qHatVec, matProp, hardmtx, 
     &   epsdot, dtime, numslip, kInteg_Hard)

      implicit none
      include 'params_xtal.inc'

      integer numslip, kInteg_Hard
      real*8  epsdot, dtime
      real*8  kappa(NKAPP), kappa_n(NKAPP), eqvalues(NEQVA)
      real*8  rrot(DIMS,DIMS), rrot_n(DIMS,DIMS), crot(DIMS,DIMS)
      real*8  crot0(DIMS,DIMS), drot(DIMS,DIMS), tau(NVECS), w_vec(DIMS)
      real*8  gamdot(MAX_SLIP), matProp(NPROPS,MAX_SLIP)
      real*8  pHatVec(NVECS,MAX_SLIP), qHatVec(DIMS,MAX_SLIP)
      real*8  rss(MAX_SLIP)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)

      integer is
      real*8  InnerProductVec2, SSKineticEqn2

      real*8  crss
      real*8  overstress(MAX_SLIP)
      common /XtalOverS/ overstress
c
c---------------------------------------------------------------------72
c
c------- slip quantities: resolve shear stress and shear strain rate
c
      do is = 1, numslip
         rss(is) = InnerProductVec2(tau, pHatVec(1,is), NVECS)
c         crss = kappa(1) + overstress(is)
         crss = kappa(is)
         gamdot(is) = SSKineticEqn2(rss(is),crss,matProp(1,is),kGAMDOT)
      enddo
c
c------- solve for hardness
c
c      call IntegrateHardening(kappa, kappa_n, eqvalues, gamdot, matProp,
c     &                     hardmtx, epsdot, dtime, numslip, kInteg_Hard)
c
c------- solve for rotation: dotR = (W - Wp)*R
c------- update slip system orientation matrix: C = R * C0
c
      call IntegrateRotation2(w_vec, rrot_n, crot0, qHatVec, gamdot, 
     &                       rrot, crot, drot, dtime, numslip)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE IntegrateHardening2(
     &   kappa, kappa_n, eqvalues, gamdot, matProp, hardmtx, epsdot, 
     &   dtime, numslip, kInteg_Code
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip, kInteg_Code
      real*8  epsdot, dtime
      real*8  kappa(NKAPP), kappa_n(NKAPP), eqvalues(NEQVA)
      real*8  gamdot(MAX_SLIP), matProp(NPROPS, MAX_SLIP)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)

      integer is, ik, jk
      real*8  h_0, tausi, taus0, xms, gamss0
      real*8  shr_min, shr_max, kappa_sat
      real*8  c, g_n, g_s, g
      real*8  dkappa, gamtot_n, delgam, fac
      real*8  kTHETA
      data    kTHETA /1.0d0/
c
c---------------------------------------------------------------------72
c
c------- total shear strain rate: Sum(|gamdot(i)|), i=1,numslip
c
      eqvalues(kSHRATE) = pzero
      do is = 1, numslip
         eqvalues(kSHRATE) = eqvalues(kSHRATE) + dabs(gamdot(is))
      enddo

      shr_min = 1.0d-6 * epsdot
      shr_max = 1.0d+1 * epsdot

      if (eqvalues(kSHRATE) .le. shr_min) eqvalues(kSHRATE) = shr_min
      if (eqvalues(kSHRATE) .ge. shr_max) eqvalues(kSHRATE) = shr_max
c
c------- accumulated shear strain: gamtot
c
      eqvalues(kGAMTOT) = eqvalues(kGAMTOT_n) + eqvalues(kSHRATE)*dtime
c
c------- integration of Voce's hardening law (one hardness/slip system)
c
      if (kInteg_Code .eq. kHARD_EXPL) then  
c
c---------- explicit update
c         do ik = 1, NKAPP
         do ik = 1, numslip
            h_0    = matProp(5,ik)
            tausi  = matProp(6,ik)
            taus0  = matProp(7,ik)
            xms    = matProp(8,ik)
            gamss0 = matProp(9,ik)

            kappa_sat = taus0 * ((eqvalues(kSHRATE) / gamss0)**xms)

            c = dtime*h_0
            g_s = kappa_sat - tausi
            g_n = kappa_n(ik) - tausi

            if ( (g_n/g_s) .le. 1.0 ) then
               g = g_n + c*(1.0-g_n/g_s)*eqvalues(kSHRATE)
            else
               g = g_n
            endif
            kappa(ik) = g + tausi
         enddo

      else if (kInteg_Code .eq. kHARD_MIDP) then
c
c---------- generalized mid-point rule 
c         do ik = 1, NKAPP
         do ik = 1, numslip
            h_0    = matProp(5,ik)
            tausi  = matProp(6,ik)
            taus0  = matProp(7,ik)
            xms    = matProp(8,ik)
            gamss0 = matProp(9,ik)

            kappa_sat = taus0 * ((eqvalues(kSHRATE) / gamss0)**xms)

            c = dtime*h_0
            g_s = kappa_sat - tausi
            g_n = kappa_n(ik) - tausi

            if ( (g_n/g_s) .le. 1.0 ) then
               g = g_n + c*( 
     &                (1.0-kTHETA)*(1.0-g_n/g_s)*eqvalues(kSHRATE_n)
     &                +   (kTHETA)*eqvalues(kSHRATE) 
     &                     )
               g = g / (1.0 + c*kTHETA*eqvalues(kSHRATE)/g_s)
            else
               g = g_n
            endif
            kappa(ik) = g + tausi
         enddo

      else if (kInteg_Code .eq. kHARD_ANAL) then

         gamtot_n = eqvalues(kGAMTOT_n)
         delgam   = eqvalues(kSHRATE) * dtime

         do ik = 1, numslip

            h_0    = matProp(5,ik)
            tausi  = matProp(6,ik)
            taus0  = matProp(7,ik)
            xms    = matProp(8,ik)
            gamss0 = matProp(9,ik)

            kappa_sat = taus0 * ((eqvalues(kSHRATE) / gamss0)**xms)
            g_s = kappa_sat - tausi
            fac = dabs(h_0/g_s)

            dkappa = 0.0
            do jk = 1, numslip
              dkappa = dkappa + 
     &                    hardmtx(ik,jk)*dabs(gamdot(jk))*dtime/delgam
            enddo
            dkappa = dkappa * g_s * dexp(-gamtot_n*fac) * 
     &                                        (1.0 - dexp(-delgam*fac))

            kappa(ik) = kappa_n(ik) + dkappa

         enddo
         
      else
c
c------- wrong code number
         call RunTimeError2(XTAL_O, 
     &                      'IntegrateHardening: Wrong kInteg_Code!')

      endif

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE IntegrateRotation2(w_vec, rrot_n, crot0, qHatVec, 
     &   gamdot, rrot, crot, drot, dtime, numslip
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip
      real*8  dtime
      real*8  w_vec(DIMS), rrot_n(DIMS,DIMS), crot0(DIMS,DIMS)
      real*8  qHatVec(DIMS,MAX_SLIP), gamdot(MAX_SLIP)
      real*8  rrot(DIMS,DIMS), crot(DIMS,DIMS), drot(DIMS,DIMS)

      integer is
      real*8  wpHat(DIMS), spin(DIMS)
c
c---------------------------------------------------------------------72
c
c------- plastic spin from slip system activity
c
      call SetTensor2(wpHat, pzero, DIMS)
      do is = 1, numslip
         call AddScaledTensor2(gamdot(is), qHatVec(1,is), wpHat, DIMS)
      enddo
c
c------- net spin: (W - Wp)
c
      call SetTensor2(spin, pzero, DIMS)
      call AddTensors2(pone, w_vec, -pone, wpHat, spin, DIMS)
c
c------- incremental rot: dR = exp(spin*dtime)
c
      call IncrementalRotTensor2(drot, spin, dtime)
c
c------- rotation tensor: R = dR*R_n
c
      call MultAxB2(drot, rrot_n, rrot, DIMS)
c
c------- slip system orientation matrix crot: C = R*C_0
c
      call MultAxB2(rrot, crot0, crot, DIMS)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE IncrementalRotTensor2(
     &   drot, spin, dtime
     &   )

      implicit none

      real*8  dtime
      real*8  drot(3, 3), spin(3)

      real*8  th1, th2, th3, tau, taua
c
c---------------------------------------------------------------------72
c
      th1 = spin(1) * dtime
      th2 = spin(2) * dtime
      th3 = spin(3) * dtime

      tau = sqrt(th1 * th1 + th2 * th2 + th3 * th3)
      if (tau .ne. 0.0) then
         taua = tan(tau / 2.0) / tau
         th1 = taua * th1
         th2 = taua * th2
         th3 = taua * th3
         tau = taua * tau
      endif

      tau = 2.0 / (1.0 + tau * tau)

      drot(1, 1) = 1.0 - tau * (th1 * th1 + th2 * th2)
      drot(1, 2) = - tau * (th1 + th2 * th3)
      drot(1, 3) = tau * ( - th2 + th1 * th3)
      drot(2, 1) = tau * (th1 - th2 * th3)
      drot(2, 2) = 1.0 - tau * (th1 * th1 + th3 * th3)
      drot(2, 3) = - tau * (th3 + th1 * th2)
      drot(3, 1) = tau * (th2 + th1 * th3)
      drot(3, 2) = tau * (th3 - th1 * th2)
      drot(3, 3) = 1.0 - tau * (th2 * th2 + th3 * th3)
 
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      logical FUNCTION ConvergeState2(
     &   norm_a, norm_b, norm_a0, norm_b0, toler
     &   )

      implicit none

      real*8  norm_a, norm_b, norm_a0, norm_b0, toler
c
c---------------------------------------------------------------------72
c
      ConvergeState2 = (dabs(norm_a - norm_a0) .lt. toler*norm_a0) .and. 
     &                (dabs(norm_b - norm_b0) .lt. toler*norm_b0)

      if (.not.ConvergeState2) then
          norm_a0 = norm_a
          norm_b0 = norm_b
      endif

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE UpdateStress2(
     &   s_ij, stress, tau, ee_ij, estran, statev, eqvalues, rrot, 
     &   dfgrd1, iqpt, ielem)

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      real*8  s_ij(DIMS,DIMS), stress(NVECS), tau(NVECS)
      real*8  ee_ij(DIMS,DIMS), estran(NVECS), statev(NSTAV)
      real*8  eqvalues(NEQVA), rrot(DIMS,DIMS), dfgrd1(3,3)
      integer iqpt, ielem

      real*8  ekk3, det
      real*8  Ve_ij(DIMS, DIMS)
      real*8  matx_1(3,3), matx_2(3,3), matx_3(3,3)
      real*8  F_p(3,3)

      real*8  DetOfMat3x32, InnerProductVec2
      
      real*8  plasdef1(3, 3, 1, MAX_GRN)
      common  /PlasticDefGrad1/ plasdef1
      save    /PlasticDefGrad1/
c
c---------------------------------------------------------------------72
c
c------- Elastic strain (ee) and V tensors (V = ee + I)
c
      call Vec5x1ToMat3x3Symm2(estran, ee_ij, DIMS)
      ekk3 = statev(kEVOL)/3.0
      call AddScaledTensor2(ekk3, Ident2nd, ee_ij, DIMS*DIMS)

      call AddTensors2(pone, Ident2nd, pone, ee_ij, Ve_ij, DIMS*DIMS)
      statev(kDETVe) = DetOfMat3x32(Ve_ij)
c
c------- plastic deformation gradient
c
      call InverseOfMat3x32(rrot, matx_1, det)
      call InverseOfMat3x32(Ve_ij, matx_2, det)
      call MultAxB2(matx_1, matx_2, matx_3, 3)
      call MultAxB2(matx_3, dfgrd1, F_p, 3)
      call EqualTensors2(F_p, plasdef1(1, 1, 1, ielem), 9)
c
c------- deviatoric (vector) and volumetric (pressure) Cauchy stress
c
      call SetToScaledtensor2(1./statev(kDETVe), tau, stress, NVECS)
      statev(kPRES) = statev(kPRES)/statev(kDETVe)
c
c------- VonMises stress
c
      eqvalues(kMISES) = dsqrt(3./2.*InnerProductVec2(stress, stress,
     &                                                NVECS))
c
c------- Cauchy stress (tensor form)
c
      call Vec5x1ToMat3x3Symm2(stress, s_ij, DIMS)
      call AddScaledTensor2(statev(kPRES), Ident2nd, s_ij, DIMS*DIMS)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE PlasticModuli2(
     &   fCep, fCeiDevHat, fCeDevVolHat, fCeVol, jac_e
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      real*8  fCeVol, jac_e
      real*8  fCep(DIMS2,DIMS2)
      real*8  fCeiDevHat(NVECS,NVECS), fCeDevVolHat(NVECS)

      real*8  det
      real*8  fMat5x5(NVECS, NVECS), fVect5(NVECS)
      real*8  fOuter6x6(DIMS2, DIMS2), fVect6(DIMS2)
      real*8  fMat5x6(NVECS, DIMS2), fOuter5x6(NVECS, DIMS2)
      real*8  fCepDev(NVECS, DIMS2), fCepVol(DIMS2)

      real*8  fDevMat5x6(NVECS, DIMS2), fDevMat6x5(DIMS2, NVECS)
      real*8  fMatTId5x6(NVECS, DIMS2)
      common /transfMatxs/ fDevMat5x6, fDevMat6x5, fMatTId5x6

      real*8  lhs5x5(NVECS, NVECS)
      common /SumGamPPt/ lhs5x5
c
c
c---------------------------------------------------------------------72
c
c---- Deviatoric moduli: fCepDev_5x6
c
c-------------- ([Ce_d]^-1 + [S])^-1
      call AddTensors2(pone, fCeiDevHat, pone, lhs5x5, fMat5x5, 
     &                NVECS*NVECS) 
      call MatrixInverse2(fMat5x5, NVECS, NVECS, det)
c
c-------------- ([T][Pdev] + [Ce_d]^-1 {H} {1}^T)
      call MultAxu2(fCeiDevHat, fCeDevVolHat, fVect5, NVECS)
      call OuterProductVec_G2(fvect5, Ident, fOuter5x6, NVECS, DIMS2)
      call Addtensors2(pone, fMatTId5x6, pone, fOuter5x6, fMat5x6, 
     &                NVECS*DIMS2)
c
c-------------- (([Ce_d]^-1 + [S])^-1) ([T][Pdev] + [Ce_d]^-1 {H} {1}^T)
      call MultAxB_G2(fMat5x5, fMat5x6, fCepDev, NVECS, NVECS, DIMS2)
c
c---- Volumetric moduli: fCepVol_6
c
c-------------- {H}^T [T] [Pdev] + M {1}^T
      call MultATxu_G2(fMatTId5x6, fCeDevVolHat, fVect6, NVECS, DIMS2)
      call AddTensors2(pone, fVect6, fCeVol, Ident, fCepVol, DIMS2)
c
c-------------- {H}^T [S] [Cep_d]
      call MultAxB_G2(lhs5x5, fCepDev, fMat5x6, NVECS, NVECS, DIMS2)
      call MultATxu_G2(fMat5x6, fCeDevVolHat, fVect6, NVECS, DIMS2)
c
c-------------- {H}^T [T] [Pdev] + M {1}^T - {H}^T [S] [Cep_d]
      call AddScaledTensor2(-pone, fVect6, fCepVol, DIMS2)
c
c---- Algorithmic Moduli: fCep_6x6
c
c-------------- [J] [Cep_d]
      call MultAxB_G2(fDevMat6x5, fCepDev, fCep, DIMS2, NVECS, DIMS2)
c
c-------------- [J] [Cep_d] + {1} {Cep_v}^T -> [Cep]: {s} = [Cep] {d}
      call OuterProductVec2(Ident, fCepVol, fOuter6x6, DIMS2)
      call AddScaledTensor2(pone, fOuter6x6, fCep, DIMS2*DIMS2)
      call SetToScaledTensor2(pone/jac_e, fCep, fCep, DIMS2*DIMS2)
c
c---- Change fCep to be consistent with abaqus representation of {d}
c
      call SetToScaledTensor2(phalf, fCep(1,4), fCep(1,4),DIMS2*DIMS2/2)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE UpdateDeformQnts2(
     &   gamdot, eqvalues, statev, crot, zBar0, dtime, numslip,
     &   eps, eps_n)

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip
      real*8  dtime
      real*8  statev(NSTAV), eqvalues(NEQVA), gamdot(MAX_SLIP)
      real*8  crot(DIMS,DIMS), zBar0(DIMS,DIMS,MAX_SLIP)
      real*8  eps(9), eps_n(9)

      integer is
      real*8  eqpdot, maxGamdot, ratio
      real*8  wp_vec(DIMS), matx_1(DIMS,DIMS)
      real*8  lp_ij(DIMS,DIMS), dp_ij(DIMS,DIMS), wp_ij(DIMS, DIMS)

      real*8  InnerProductVec2, MaxAbsValueOfVec2
c
c---------------------------------------------------------------------72
c
c------- Plastic velocity gradient lp in crystal axis
c
      call setTensor2(lp_ij, pzero, DIMS*DIMS)
      do is = 1, numslip
         call AddScaledTensor2(gamdot(is), zBar0(1,1,is), lp_ij, 
     &                        DIMS*DIMS)
      enddo
c
c------- Symm and Skew parts of lp: lp =  dp + wp
c
c      call EqualTensors(lp_ij, matx_1, DIMS*DIMS)
c      call MultQAQT(crot, matx_1, lp_ij, DIMS)

      call SymmetrizeTensor2(lp_ij, dp_ij, DIMS)
      call SkewSymmetrizeTensor2(lp_ij, wp_ij, DIMS)
c
c------- Equivalent plastic strain
c
      eps(1) = eps_n(1) + dtime * dp_ij(1,1)
      eps(2) = eps_n(2) + dtime * dp_ij(1,2)
      eps(3) = eps_n(3) + dtime * dp_ij(1,3)
      eps(4) = eps_n(4) + dtime * dp_ij(2,1)
      eps(5) = eps_n(5) + dtime * dp_ij(2,2)
      eps(6) = eps_n(6) + dtime * dp_ij(2,3)
      eps(7) = eps_n(7) + dtime * dp_ij(3,1)
      eps(8) = eps_n(8) + dtime * dp_ij(3,2)
      eps(9) = eps_n(9) + dtime * dp_ij(3,3)
      eqvalues(kEQP) = dsqrt(2./3.*InnerProductVec2(eps, eps, 
     &                 DIMS*DIMS))
c      eqpdot = dsqrt(2./3.*InnerProductVec2(dp_ij, dp_ij, DIMS*DIMS))
c      eqvalues(kEQP) = eqvalues(kEQP_n) + dtime*eqpdot
c
c------- Norm of plastic spin
c
c      call Mat3x3ToVec3x1Skew2(wp_ij, wp_vec, DIMS)
      statev(kWPNORM) = dsqrt(InnerProductVec2(wp_vec, wp_vec, DIMS))
c
c------- Slip system activity
c
      maxGamdot = MaxAbsValueOfVec2(gamdot, numslip)
      if (maxGamdot .lt. TINY) maxGamdot = TINY

      statev(kSSACT) = 0.0
      do is = 1, numslip
         ratio = dabs(gamdot(is))/maxGamdot
         if (ratio .ge. 0.05d0) statev(kSSACT) = statev(kSSACT) + pone
      enddo
        
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE StressSolveViscoPlastic2(
     &   tau_lat, d_vec_lat, kappa, pBar0Vec, ppTBar0, matProp, sigfs,
     &   tauSlip, dtime, theta, thetao, epsdot, tolerNewt, maxIterNewt,
     &   numslip, numvtx
     &   )

      implicit none
      include 'params_xtal.inc'

      integer maxIterNewt, numslip, numvtx
      real*8  dtime, theta, thetao, epsdot, tolerNewt
      real*8  tau_lat(NVECS), d_vec_lat(NVECS), kappa(NKAPP)
      real*8  pBar0Vec(NVECS,MAX_SLIP), ppTBar0(NVECS,NVECS,MAX_SLIP)
      real*8  sigfs(NVECS,MAX_VTX), tauSlip(MAX_SLIP)
      real*8  matProp(NPROPS,MAX_SLIP)

      integer iv, ierr
      real*8  rescale
      real*8  sdotd(MAX_VTX)
      real*8  InnerProductVec2
c
c---------------------------------------------------------------------72
c
c------- solve for stresses in lattice reference frame
c
c------- use effective deformation rate to normalized (scale-down) D
c
c      epsdot = dsqrt(2./3.*InnerProductVec(d_vec_lat, d_vec_lat, NVECS))
      call SetToScaledTensor2(1./epsdot, d_vec_lat, d_vec_lat, NVECS)
c
c------- compute work to select guessed values from vertex stresses
c
      do iv = 1, numvtx
         sdotd(iv) = InnerProductVec2(sigfs(1,iv), d_vec_lat, NVECS) 
      enddo
c
c------- initialized local flag to monitor convergence of Newton iters
c
      ierr = XTAL_CONVERGED
c
c------- solve fot the stresses. If needed, will try all vertices
c
      do iv = 1, numvtx
c
c---------- compute and scaled initial guess
         call InitialGuessStress2(tau_lat, sdotd, kappa, sigfs,pBar0Vec,
     &                           matProp, numslip, numvtx)
c
c---------- compute stresses
         call StressViscoPlastic2(tau_lat, d_vec_lat, kappa, pBar0Vec, 
     &                           ppTBar0, matProp, tolerNewt,
     &                           maxIterNewt, numslip, ierr)
c 
c---------- if not converged, try another initial guess
         if (ierr .ne. XTAL_CONVERGED) then
            call WriteMessage2(XTAL_O, 
     &         'StressSolveViscoPlastic: did not converged !!!')
            call WriteMessage2(XTAL_O, 
     &         'StressSolveViscoPlastic: will try next vertex stresses')
            ierr = XTAL_CONVERGED
c 
c---------- or else converged; then, rescale the stresses, and return
         else
            rescale = epsdot**matProp(3,1)   ! epsdot**m
            call SetToScaledTensor2(rescale, tau_lat, tau_lat, NVECS)
            return
         endif

      enddo
c
c------- if gets here, then the stress solution did not converged
c-------   for all trials. 
c
      call RunTimeError2(XTAL_O, 
     &       'StressSolveViscoPlastic: Newton method did not converged')

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE InitialGuessStress2(
     &   s, sdotd, kappa, sigfs, pBar0Vec, matProp, numslip, numvtx
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip, numvtx
      real*8  s(NVECS), sdotd(MAX_VTX), kappa(NKAPP)
      real*8  sigfs(NVECS, MAX_VTX), pBar0Vec(NVECS, MAX_SLIP)
      real*8  matProp(NPROPS,MAX_SLIP)

      integer kmax, is
      real*8  signmax, factor
      real*8  rssk(MAX_SLIP)

      integer IndexMaxAbsValueOfVec2
      real*8  MaxAbsValueOfVec2, InnerProductVec2, SignOf2

      real*8  crss
      real*8  overstress(MAX_SLIP)
      common /XtalOverS/ overstress
c
c---------------------------------------------------------------------72
c
c------- locate sign of max value of sdotd (plastic work)
c
      kmax = IndexMaxAbsValueOfVec2(sdotd, numvtx)
      signmax = SignOf2(sdotd(kmax))
      sdotd(kmax) = pzero        ! avoids 'kmax' vertex being re-used
c
c------- asign correct sign to selected vertex stress (initial guess)
c
      call SetToScaledTensor2(signmax, sigfs(1,kmax), s, NVECS)
c
c------ compute factor to scaled down initial guess
c------   factor = |rss|_max / kappa * gam0**m
c
      do is = 1, numslip
c         crss = kappa(1) + overstress(is)
         crss = kappa(is)
         rssk(is) = InnerProductVec2(pBar0Vec(1,is), s, NVECS)/crss
      enddo
      factor = MaxAbsValueOfVec2(rssk, numslip) *
     &                                      matProp(4,1)**matProp(3,1)
c
c------- scaled initial guess for stress
c
      call SetToScaledTensor2(1./factor, s, s, NVECS)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE StressViscoPlastic2(
     &   s, d_vec_lat, kappa, pBar0Vec, ppTBar0, matProp, toler,
     &   maxIters, numslip, ierr 
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer maxIters, numslip, ierr
      real*8  toler
      real*8  s(NVECS), d_vec_lat(NVECS), kappa(NKAPP)
      real*8  pBar0Vec(NVECS,MAX_SLIP), ppTBar0(NVECS,NVECS,MAX_SLIP)
      real*8  matProp(NPROPS,MAX_SLIP)

      integer iters
      real*8  rhs_norm_0, rhs_norm, search
      real*8  delts(NVECS), s0(NVECS)
      real*8  rhs(NVECS), rss(MAX_SLIP), lhs(NVECS,NVECS)

      logical Converged2
c      real*8  InnerProductVec
      real*8  NormOfVector2
c
c---------------------------------------------------------------------72
c
c------- compute initial residual and its norm
c
      call FormResidual2(rhs, s, rss, d_vec_lat, kappa, pBar0Vec, 
     &                  matProp, numslip)
c      rhs_norm_0 = dsqrt(InnerProductVec(rhs, rhs, NVECS))
      rhs_norm_0 = NormOfVector2(rhs, NVECS)
c
c------- initialize flags and start iterations
c
      iters = 0
      do while(.not.Converged2(rhs, toler) .and. iters .lt. maxIters)
         iters = iters + 1
         call EqualTensors2(s, s0, NVECS)
c
c---------- compute local jacobian
         call FormJacobian2(lhs, rss, ppTBar0, kappa, matProp, numslip)
c
c---------- solve for the increment of stress
c---------- use the lower half of the symmetric matrix lhs
         call lsolve2(lhs, rhs, 1, ierr, NVECS)
         if (ierr .eq. XTAL_SING_JACOBIAN) return

         call EqualTensors2(rhs, delts, NVECS)
c
c---------- update variables
         search = pone
         call AddTensors2(pone, s0, search, delts, s, NVECS)
c
c---------- compute new residual and its norm
         call FormResidual2(rhs, s, rss, d_vec_lat, kappa, pBar0Vec, 
     &                     matProp, numslip)
c         rhs_norm = dsqrt(InnerProductVec(rhs, rhs, NVECS))
         rhs_norm = NormOfVector2(rhs, NVECS)
c
c---------- simple line search
         do while ( rhs_norm .gt. rhs_norm_0 ) 

            search = search*0.5d0
            if (search .lt. TINY) then
               call WriteWarning2(XTAL_O, 
     &                 'StressViscoPlastic: LS Failed, search < TINY')
               ierr = XTAL_LS_FAILED
               return
            endif
            call AddTensors2(pone, s0, search, delts, s, NVECS)

            call FormResidual2(rhs, s, rss, d_vec_lat, kappa, pBar0Vec, 
     &                        matProp, numslip)
c            rhs_norm = dsqrt(InnerProductVec(rhs, rhs, NVECS))
            rhs_norm = NormOfVector2(rhs, NVECS)

         enddo
c
c---------- update norm of residual
         rhs_norm_0 = rhs_norm

      enddo
c
c------- check convergence on iters < maxIters
c
      if (iters .ge. maxIters) then
        call WriteMessage2(XTAL_O,'StressViscoPlastic: iters>maxIters')
        ierr = XTAL_MAX_ITERS_HIT
        return
      endif
c
c------- keep track of Newton's iterations
c
      write(XTAL_O, *) 'iterCounter (viscoplastic) = ', iters

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE FormResidual2(
     &   rhs, s, rss, d_vec_lat, kappa, pBar0Vec, matProp, numslip
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip
      real*8  rhs(NVECS), s(NVECS), rss(MAX_SLIP)
      real*8  d_vec_lat(NVECS), kappa(NKAPP)
      real*8  pBar0Vec(NVECS,MAX_SLIP), matProp(NPROPS,MAX_SLIP)

      integer is
      real*8  gamdot(MAX_SLIP), d_curr(NVECS)
      real*8  InnerProductVec2, SSKineticEqn2

      real*8  crss
      real*8  overstress(MAX_SLIP)
      common /XtalOverS/ overstress
c
c---------------------------------------------------------------------72
c
c------- resolve shear stresses and shear strain rates
c
      do is = 1, numslip
         rss(is) = InnerProductVec2(s, pBar0Vec(1,is), NVECS)
c         crss = kappa(1) + overstress(is)
         crss = kappa(is)
         gamdot(is) = SSKineticEqn2(rss(is),crss,matProp(1,is),kGAMDOT)
      enddo
c
c------- rate of deformation due to slip system activity
c
      call SetTensor2(d_curr, pzero, NVECS)
      do is = 1, numslip
         call AddScaledTensor2(gamdot(is), pBar0Vec(1,is), d_curr,NVECS)
      enddo
c
c------- residual
c
      call SetTensor2(rhs, pzero, NVECS)
      call AddTensors2(+pone, d_vec_lat, -pone, d_curr, rhs, NVECS)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE FormJacobian2(
     &   lhs, rss, ppTBar0, kappa, matProp, numslip
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip
      real*8  lhs(NVECS,NVECS), rss(MAX_SLIP), kappa(NKAPP)
      real*8  ppTBar0(NVECS,NVECS,MAX_SLIP), matProp(NPROPS,MAX_SLIP)

      integer is
      real*8  SSKineticEqn2
      real*8  dGamdTau(MAX_SLIP)

      real*8  crss
      real*8  overstress(MAX_SLIP)
      common /XtalOverS/ overstress
c
c---------------------------------------------------------------------72
c
c------- d(gamdot)/d(tau)
c
      do is = 1, numslip
c         crss = kappa(1) + overstress(is)
         crss = kappa(is)
         dGamdTau(is) = SSKineticEqn2(rss(is),crss,matProp(1,is),
     &                               kdGAMdTAU)
      enddo
c
c------- jacobian
c
      call SetTensor2(lhs, pzero, NVECS*NVECS)
      do is = 1, numslip
         call AddScaledTensor2(dGamdTau(is), ppTBar0(1,1,is), lhs,
     &                        NVECS*NVECS)
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      logical FUNCTION Converged2(
     &   res, toler
     &   )

      implicit none
      include 'params_xtal.inc'

      real*8  toler
      real*8  res(NVECS)

      integer i
c
c---------------------------------------------------------------------72
c     
c------- check convergence on residual
c
      Converged2 = ( dabs(res(1)) .lt. toler )
      do i = 2, NVECS
         Converged2 = ( (dabs(res(i)) .lt. toler) .and. Converged2)
      enddo

      return
      END
c     
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE StressSolveDeviatoric2(
     &   tau, d_vec, estran, estran_n, kappa,
     &   gamdot, rss, drot, fCeiDevHat, fCeDevVolHat, pHatVec,
     &   ppTHat, matProp, e_kk, dtime, tolerNewt, numslip,
     &   maxIterNewt, iterCounterN, ierr
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip, maxIterNewt, iterCounterN, ierr
      real*8  e_kk, dtime, tolerNewt
      real*8  tau(NVECS), d_vec(NVECS), estran(NVECS), estran_n(NVECS)
      real*8  kappa(NKAPP), gamdot(MAX_SLIP), drot(DIMS,DIMS)
      real*8  fCeiDevHat(NVECS,NVECS), fCeDevVolHat(NVECS)
      real*8  pHatVec(NVECS), ppTHat(NVECS,NVECS)
      real*8  matProp(NPROPS, MAX_SLIP)
      real*8  rss(MAX_SLIP)

      integer subIncr, totSubIncrs, is
      real*8  xm_o(MAX_SLIP), tmp
      real*8  tau_o(NVECS), tau_save(NVECS)
c
c---------------------------------------------------------------------72
c
c------- counters for sub-increments / backup some variables
c
      subIncr = 1
      totSubIncrs = 1
      do is = 1, numslip
        xm_o(is) = matProp(3,is)
      enddo
      call EqualTensors2(tau, tau_o, NVECS)
c
c------- integrate crystal constitutive equations
c
      call SolveDeviatoric2(tau, d_vec, estran, estran_n, kappa,
     &           gamdot, rss, drot, fCeiDevHat, fCeDevVolHat, pHatVec,
     &           ppTHat, matProp, e_kk, dtime, tolerNewt, numslip,
     &           maxIterNewt, iterCounterN, ierr)
c
c------- if converged -> return, else do -> continuation method in "m"
c------- NOTE: here, continuation improves the initial guess 
c-------       for the solution variables 'tau'
c
      if (ierr .eq. XTAL_CONVERGED) return

      call WriteMessage2
     &     (XTAL_O, 'DriverStressSolveDev: using Contin. for m!')
c
c------- loop for sub-incrementation procedure
c
      do while (.true.)
c
c---------- if not converged, increase # of sub-increments 
         if (ierr .ne. XTAL_CONVERGED) then
            subIncr = 2 * subIncr - 1
            totSubIncrs = 2 * totSubIncrs
            if (totSubIncrs .gt. 16) then
               call WriteMessage2(XTAL_O,
     &                  'DriverStressSolveDev: totSubIncrs > 16')
               do is = 1, numslip
                  matProp(3,is) = xm_o(is)
                  call BoundForArgPowLaw2(matProp(3,is))
               enddo
               return
            endif
            if (subIncr .eq. 1) 
     &           call EqualTensors2(tau_o, tau, NVECS)
            if (subIncr .gt. 1) 
     &           call EqualTensors2(tau_save, tau, NVECS)
c
c---------- if converged, adjust subincrements
         else if (subIncr .lt. totSubIncrs) then
            if ( (subIncr/2*2) .eq. subIncr) then
               subIncr = subIncr / 2 + 1
               totSubIncrs = totSubIncrs / 2
            else
               subIncr = subIncr + 1
            endif
c
c---------- successful return for continuation method
         else
            call WriteMessage2(XTAL_O, 'Contin. successful !!!')
            return
         endif
c
c---------- trial rate sensitivity coefficient m
         tmp = real(subIncr) / real(totSubIncrs)   
         do is = 1, numslip
            matProp(3,is) = xm_o(is) / tmp
            call BoundForArgPowLaw2(matProp(3,is))
         enddo
c
c---------- save current convergent solution before getting next solution
         if (subIncr .gt. 1 .and. ierr .eq. XTAL_CONVERGED)
     &                 call EqualTensors2(tau, tau_save, NVECS)
c
c---------- report sub-incrs status
         write(XTAL_O, 1000) subIncr, totSubIncrs, matProp(3,1)
c
c---------- integrate crystal constitutive equations
c        iterCounterN = 0
         ierr = XTAL_CONVERGED
         call SolveDeviatoric2(tau, d_vec, estran, estran_n, kappa,
     &           gamdot, rss, drot, fCeiDevHat, fCeDevVolHat, pHatVec,
     &           ppTHat, matProp, e_kk, dtime, tolerNewt, numslip,
     &           maxIterNewt, iterCounterN, ierr)

      enddo

1000  format(3x, 'subInc/totSubIncrs = ', i2, '/', i3, 5x, 
     &      'm(1)=', e12.5)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SolveDeviatoric2(
     &   tau, d_vec, estran, estran_n, kappa, gamdot, rss, drot,
     &   fCeiDevHat, fCeDevVolHat, pHatVec, ppTHat, matProp, e_kk, 
     &   dtime, toler, numslip, maxIters, iterCounterN, ierr
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip, maxIters, iterCounterN, ierr
      real*8  e_kk, dtime, toler
      real*8  tau(NVECS), d_vec(NVECS), estran(NVECS), estran_n(NVECS)
      real*8  kappa(NKAPP), gamdot(MAX_SLIP), drot(DIMS,DIMS)
      real*8  fCeiDevHat(NVECS,NVECS), fCeDevVolHat(NVECS)
      real*8  pHatVec(NVECS), ppTHat(NVECS,NVECS)
      real*8  matProp(NPROPS, MAX_SLIP)
      real*8  rss(MAX_SLIP)

      integer iters
      real*8  rhs_norm_0, rhs_norm, search
      real*8  estranHat(NVECS), estranStar(NVECS)
      real*8  deltau(NVECS), tau0(NVECS)
      real*8  rhs(NVECS), lhs(NVECS,NVECS)
      real*8  dqr5x5(NVECS,NVECS)

      logical Converged2
c      real*8  InnerProductVec
      real*8  NormOfVector2
c
c---------------------------------------------------------------------72
c
c------- predictor elastic strain in Btilde: dR*e_n*dR^T + D*dt 
c
      call RotMat5x5ForSymm2(drot, dqr5x5, DIMS)
      call MultAxu2(dqr5x5, estran_n, estranHat, NVECS)
      call AddTensors2(pone, estranHat, dtime, d_vec, estranStar, NVECS)
c
c------- compute initial residual 
c
      call ComputeResidual2(rhs, tau, estran, estranStar, fCeiDevHat, 
     &                     fCeDevVolHat, rss, gamdot, kappa, pHatVec, 
     &                     matProp, e_kk, dtime, numslip)
c      rhs_norm_0 = dsqrt(InnerProductVec(rhs, rhs, NVECS))
      rhs_norm_0 = NormOfVector2(rhs, NVECS)
c
c------- initialize flags and start Newton iterations for stress
c
      iters = 0
      do while(.not.Converged2(rhs, toler) .and. iters .lt. maxIters)
         iters = iters + 1
         call EqualTensors2(tau, tau0, NVECS)
c
c---------- compute local jacobian
         call ComputeJacobian2(lhs, rss, kappa, ppTHat, fCeiDevHat, 
     &                        matProp, dtime, numslip)
c
c---------- solve for the increment of stress
         call lsolve2(lhs, rhs, 1, ierr, NVECS)
         if (ierr .eq. XTAL_SING_JACOBIAN) then
            call WriteWarning2(XTAL_O,
     &             'StressSolveDeviatoric: Jacobian is singular')
            return
         endif

         call EqualTensors2(rhs, deltau, NVECS)
c
c---------- update stresses
         search = pone
         call AddTensors2(pone, tau0, search, deltau, tau, NVECS)
c
c---------- compute new residual
         call ComputeResidual2(rhs, tau, estran, estranStar, fCeiDevHat, 
     &                        fCeDevVolHat, rss, gamdot, kappa, pHatVec,
     &                        matProp, e_kk, dtime, numslip)
c         rhs_norm = dsqrt(InnerProductVec(rhs, rhs, NVECS))
         rhs_norm = NormOfVector2(rhs, NVECS)
c
c---------- simple line search
         do while ( rhs_norm .gt. rhs_norm_0 )

            search = search*0.5d0
            if (search .lt. TINY) then
               call WriteWarning2(XTAL_O,
     &                'StressSolveDeviatoric: LS Failed, search < TINY')
               ierr = XTAL_LS_FAILED
               return
            endif
            call AddTensors2(pone, tau0, search, deltau, tau, NVECS)

            call ComputeResidual2(rhs, tau, estran, estranStar, 
     &                           fCeiDevHat, fCeDevVolHat, rss, gamdot,
     &                           kappa, pHatVec, matProp, e_kk, dtime, 
     &                           numslip)
c            rhs_norm = dsqrt(InnerProductVec(rhs, rhs, NVECS))
            rhs_norm = NormOfVector2(rhs, NVECS)

         enddo
c
c---------- update norm of residual
         rhs_norm_0 = rhs_norm

      enddo
c
c------- keep track of max number of newton iterations
c
      iterCounterN = max(iterCounterN, iters)
c
c------- check convergence on iters < maxIters
c
      if (iters .ge. maxIters) then
         call WriteWarning2(XTAL_O,
     &                'StressSolveDeviatoric: iters > maxIters')
         ierr = XTAL_MAX_ITERS_HIT
         return
      endif
         
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE ComputeResidual2(
     &   rhs, tau, estran, estranStar, fCeiDevHat, fCeDevVolHat, rss, 
     &   gamdot, kappa, pHatVec, matProp, e_kk, dtime, numslip
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip
      real*8  e_kk, dtime
      real*8  rhs(NVECS), tau(NVECS), estran(NVECS), estranStar(NVECS)
      real*8  fCeiDevHat(NVECS,NVECS), fCeDevVolHat(NVECS)
      real*8  rss(MAX_SLIP), gamdot(MAX_SLIP), kappa(NKAPP)
      real*8  pHatVec(NVECS,MAX_SLIP), matProp(NPROPS,MAX_SLIP)

      integer is
      real*8  tmp(NVECS)
      real*8  InnerProductVec2, SSKineticEqn2

      real*8  crss
      real*8  overstress(MAX_SLIP)
      common /XtalOverS/ overstress
c
c---------------------------------------------------------------------72
c
c------- elastic strain: -[Cei]({tau}-ekk*{fCeDevVol}}+{e*} = -{e}+{e*}
c
      call AddTensors2(pone, tau, -e_kk, fCeDevVolHat, tmp, NVECS)
      call MultAxu2(fCeiDevHat, tmp, estran, NVECS)
      call AddTensors2(-pone, estran, +pone, estranStar, rhs, NVECS)
c
c------- resolve shear stresses and shear strain rates
c
      do is = 1, numslip
         rss(is) = InnerProductVec2(tau, pHatVec(1,is), NVECS)
c         crss = kappa(1) + overstress(is)
         crss = kappa(is)
         gamdot(is) = SSKineticEqn2(rss(is), crss, matProp(1,is), 
     &                             kGAMDOT)
      enddo
c
c------- residual
c
      do is = 1, numslip
         call AddScaledTensor2(-dtime*gamdot(is), pHatVec(1,is), rhs, 
     &                        NVECS)
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE ComputeJacobian2(
     &   lhs, rss, kappa, ppTHat, fCeiDevHat, matProp, dtime, numslip 
     &   )
      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip
      real*8  dtime
      real*8  lhs(NVECS,NVECS), rss(MAX_SLIP), kappa(NKAPP)
      real*8  ppTHat(NVECS,NVECS,MAX_SLIP), fCeiDevHat(NVECS,NVECS)
      real*8  matProp(NPROPS,MAX_SLIP)

      integer is
      real*8  dGamdTau(MAX_SLIP)
      real*8  SSKineticEqn2

      real*8  crss
      real*8  overstress(MAX_SLIP)
      common /XtalOverS/ overstress

      real*8  sGamPPt(NVECS,NVECS)
      common /SumGamPPt/ sGamPPt
c
c---------------------------------------------------------------------72
c
c------- d(gamdot)/d(tau)
c
      do is = 1, numslip
         crss = kappa(is)
c         crss = kappa(1) + overstress(is)
         dGamdTau(is) = SSKineticEqn2(rss(is), crss, matProp(1,is),
     &                               kdGAMdTAU)
      enddo
c
c------- jacobian
c
      call SetTensor2(sGamPPt, pzero, NVECS*NVECS)
      do is = 1, numslip
         call AddScaledTensor2(dtime*dGamdTau(is), ppTHat(1,1,is), 
     &                        sGamPPt, NVECS*NVECS)
      enddo
      call AddTensors2(pone, fCeiDevHat, pone, sGamPPt, lhs,NVECS*NVECS)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE StressSolveVolumetric2(
     &   statev, fCeDevVolHat, fCeVol, estran, e_kk, dtime
     &   )

      implicit none
      include 'params_xtal.inc'

      real*8  fCeVol, e_kk, dtime
      real*8  fCeDevVolHat(NVECS), estran(NVECS), statev(NSTAV)

      real*8  InnerProductVec2
c
c---------------------------------------------------------------------72
c
      statev(kEVOL) = e_kk
      statev(kPRES) = InnerProductVec2(fCeDevVolHat, estran, NVECS) 
     &                + fCeVol * e_kk
      
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetTensor2(
     &   tensor, value, n
     &   )
      implicit none

      integer n
      real*8  value
      real*8  tensor(n)

      integer i
c
c---------------------------------------------------------------------72
c
      do i = 1, n
         tensor(i) = value
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SymmetrizeTensor2(
     &   tensor, symmTensor, n
     &   )

      implicit none
      include 'params_xtal.inc'

      integer  n
      real*8   tensor(n*n), symmTensor(n*n)
c
c---------------------------------------------------------------------72
c
      if ((n .ne. 2) .and. (n .ne. 3))
     &   call RunTimeError2(XTAL_O,'SkewSymmetrizeTensor: n =! 2 or 3')

      if (n .eq. 2) then
         symmTensor(1) = tensor(1)
         symmTensor(4) = tensor(4)
         symmTensor(3) = 0.5*(tensor(3) + tensor(2))
         symmTensor(2) = symmTensor(3)
      else   ! n = 3
         symmTensor(1) = tensor(1)                     ! 11
         symmTensor(5) = tensor(5)                     ! 22
         symmTensor(9) = tensor(9)                     ! 33
         symmTensor(4) = 0.5*(tensor(4) + tensor(2))   ! 12
         symmTensor(7) = 0.5*(tensor(7) + tensor(3))   ! 13
         symmTensor(8) = 0.5*(tensor(8) + tensor(6))   ! 23
         symmTensor(2) = symmTensor(4)                 ! 21
         symmTensor(3) = symmTensor(7)                 ! 31
         symmTensor(6) = symmTensor(8)                 ! 32
      endif

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SkewSymmetrizeTensor2(
     &   tensor, skewTensor, n
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer  n
      real*8   tensor(n*n), skewTensor(n*n)
c
c---------------------------------------------------------------------72
c
      if ((n .ne. 2) .and. (n .ne. 3))
     &   call RunTimeError2(XTAL_O,'SkewSymmetrizeTensor: n =! 2 or 3')

      call SetTensor2(skewTensor, pzero, n*n)

      if (n .eq. 2) then
         skewTensor(3) = 0.5*(tensor(3) - tensor(2))
         skewTensor(2) = - skewTensor(3)
      else   ! n = 3
         skewTensor(4) = 0.5*(tensor(4) - tensor(2))
         skewTensor(7) = 0.5*(tensor(7) - tensor(3))
         skewTensor(8) = 0.5*(tensor(8) - tensor(6))
         skewTensor(2) = - skewTensor(4)
         skewTensor(3) = - skewTensor(7)
         skewTensor(6) = - skewTensor(8)
      endif

      return
      END
c

c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE Mat3x3ToVec5x1Symm2(
     &   matrix, vector, n
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'
     
      integer  n
      real*8   matrix(n*n), vector(5)
c
c---------------------------------------------------------------------72
c
c------- matrix is deviatoric
c
      if (n .ne. 3)
     &   call RunTimeError2(XTAL_O, 'Mat3x3ToVect5x1Symm: n =! 3')

      vector(1) = (matrix(1) - matrix(5)) / sqr2
      vector(2) = matrix(9) * sqr32
      vector(3) = matrix(2) * sqr2
      vector(4) = matrix(3) * sqr2
      vector(5) = matrix(6) * sqr2
 
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE Vec5x1ToMat3x3Symm2(
     &   vector, matrix, n
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'
     
      integer n
      real*8  matrix(n*n), vector(5)
c
c---------------------------------------------------------------------72
c
c------- matrix is deviatoric
c
      if (n .ne. 3)
     &   call RunTimeError2(XTAL_O, 'Vect5x1ToMat3x3Symm: n =! 3')

      matrix(1) = 0.5 * (sqr2 * vector(1) - sqr23 * vector(2))
      matrix(5) =-0.5 * (sqr2 * vector(1) + sqr23 * vector(2))
      matrix(9) = vector(2) * sqr23
      matrix(2) = vector(3) / sqr2
      matrix(3) = vector(4) / sqr2
      matrix(6) = vector(5) / sqr2
      matrix(4) = matrix(2)
      matrix(7) = matrix(3)
      matrix(8) = matrix(6)
 
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE Mat3x3ToVec3x1Skew2(
     &   matrix, vector, n
     &   )

      implicit none
      include 'params_xtal.inc'
     
      integer n
      real*8  matrix(n*n), vector(3)
c
c---------------------------------------------------------------------72
c
      if (n .ne. 3)
     &   call RunTimeError2(XTAL_O, 'Vect5x1ToMat3x3Skew: n =! 3')
c
c------- uses tensor components below main diagonal
c
      vector(1) = matrix(2)
      vector(2) = matrix(3)
      vector(3) = matrix(6)
 
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE Vec3x1ToMat3x3Skew2(
     &   vector, matrix, n
     &   )

      implicit none
      include 'params_xtal.inc'
     
      integer  n
      real*8   matrix(n*n), vector(3)
c
c---------------------------------------------------------------------72
c
      if (n .ne. 3)
     &   call RunTimeError2(XTAL_O, 'Vect5x1ToMat3x3Skew: n =! 3')

      matrix(1) = 0.
      matrix(5) = 0.
      matrix(9) = 0.
      matrix(2) = vector(1)
      matrix(3) = vector(2)
      matrix(6) = vector(3)
      matrix(4) = - vector(1)
      matrix(7) = - vector(2)
      matrix(8) = - vector(3)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE RotMat5x5ForSymm2(
     &   cmatrix, qr5x5, n
     &   )

      implicit 	none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer  n
      real*8   cmatrix(n*n), qr5x5(5, 5)

      real*8   c11, c21, c31, c12, c22, c32, c13, c23, c33
c
c---------------------------------------------------------------------72
c
c----- construct 5x5 rotation matrix for symm 2nd order tensors: 
c         [A_sm]=[C][A_lat][C]' <=>  {A_sm} = [qr5x5]{A_lat}
c----- with: {A}={()/sqr2,sqr32*(),sqr2*(),sqr2*(),sqr2*()}

      if (n .ne. 3)
     &   call RunTimeError2(XTAL_O, 'RotMat5x5ForSymm: n =! 3')

      c11 = cmatrix(1)
      c21 = cmatrix(2)
      c31 = cmatrix(3)
      c12 = cmatrix(4)
      c22 = cmatrix(5)
      c32 = cmatrix(6)
      c13 = cmatrix(7)
      c23 = cmatrix(8)
      c33 = cmatrix(9)

      qr5x5(1, 1)  =  0.5d0 * (c11 * c11 - c12 * c12 - 
     &                                           c21 * c21 + c22 * c22)
      qr5x5(1, 2)  =  sqr3 / 2.d0 * (c13 * c13 - c23 * c23)
      qr5x5(1, 3)  =  c11 * c12 - c21 * c22
      qr5x5(1, 4)  =  c11 * c13 - c21 * c23
      qr5x5(1, 5)  =  c12 * c13 - c22 * c23
      qr5x5(2, 1)  =  sqr3 / 2.d0 * (c31 * c31 - c32 * c32)
      qr5x5(2, 2)  =  1.5d0 * c33 * c33 - 0.5d0
      qr5x5(2, 3)  =  sqr3 * c31 * c32
      qr5x5(2, 4)  =  sqr3 * c31 * c33
      qr5x5(2, 5)  =  sqr3 * c32 * c33
      qr5x5(3, 1)  =  c11 * c21 - c12 * c22
      qr5x5(3, 2)  =  sqr3 * c13 * c23
      qr5x5(3, 3)  =  c11 * c22 + c12 * c21
      qr5x5(3, 4)  =  c11 * c23 + c13 * c21
      qr5x5(3, 5)  =  c12 * c23 + c13 * c22
      qr5x5(4, 1)  =  c11 * c31 - c12 * c32
      qr5x5(4, 2)  =  sqr3 * c13 * c33
      qr5x5(4, 3)  =  c11 * c32 + c12 * c31
      qr5x5(4, 4)  =  c11 * c33 + c13 * c31
      qr5x5(4, 5)  =  c12 * c33 + c13 * c32
      qr5x5(5, 1)  =  c21 * c31 - c22 * c32
      qr5x5(5, 2)  =  sqr3 * c23 * c33
      qr5x5(5, 3)  =  c21 * c32 + c22 * c31
      qr5x5(5, 4)  =  c21 * c33 + c23 * c31
      qr5x5(5, 5)  =  c22 * c33 + c23 * c32

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE RotMat3x3ForSkew2(
     &   cmatrix, qr3x3, n
     &   )

      implicit  none
      include 'params_xtal.inc'

      integer  n
      real*8   cmatrix(n*n), qr3x3(3, 3)

      real*8   c11, c21, c31, c12, c22, c32, c13, c23, c33
c
c---------------------------------------------------------------------72
c
c---- Construct 3x3 rotation matrix for skew 2nd order tensors
c        [W_sm]=[C][W_lat][C]'  <=>  {W_sm} = [qr3x3]{W_lat}

      if (n .ne. 3)
     &   call RunTimeError2(XTAL_O, 'RotMat3x3ForSkew: n =! 3')

      c11 = cmatrix(1)
      c21 = cmatrix(2)
      c31 = cmatrix(3)
      c12 = cmatrix(4)
      c22 = cmatrix(5)
      c32 = cmatrix(6)
      c13 = cmatrix(7)
      c23 = cmatrix(8)
      c33 = cmatrix(9)

      qr3x3(1, 1) = c22 * c11 - c21 * c12
      qr3x3(1, 2) = c23 * c11 - c21 * c13
      qr3x3(1, 3) = c23 * c12 - c22 * c13
      qr3x3(2, 1) = c32 * c11 - c31 * c12
      qr3x3(2, 2) = c33 * c11 - c31 * c13
      qr3x3(2, 3) = c33 * c12 - c32 * c13
      qr3x3(3, 1) = c32 * c21 - c31 * c22
      qr3x3(3, 2) = c33 * c21 - c31 * c23
      qr3x3(3, 3) = c33 * c22 - c32 * c23

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MultAxu2(
     &   A, u, v, n
     &   )

      implicit none
      include 'numbers.inc'

      integer n
      real*8  A(n,n), u(n), v(n)

      integer i, j
c
c---------------------------------------------------------------------72
c
      call SetTensor2(v, pzero, n)

      do i = 1, n
         do j = 1, n
            v(i) = v(i) + A(i,j) * u(j)
         end do
      end do

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MultATxu2(
     &   A, u, v, n
     &   )

      implicit none
      include 'numbers.inc'

      integer n
      real*8  A(n,n), u(n), v(n)

      integer i, j
c
c---------------------------------------------------------------------72
c
      call SetTensor2(v, pzero, n)

      do i = 1, n
         do j = 1, n
            v(i) = v(i) + A(j,i) * u(j)
         end do
      end do

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MultAxB2(
     &   A, B, C, n
     &   )

      implicit none
      include 'numbers.inc'
     
      integer n
      real*8  A(n, n), B(n, n), C(n, n)

      integer i, j, k
c
c---------------------------------------------------------------------72
c
      call SetTensor2(C, pzero, n*n)

      do i = 1, n 
         do j = 1, n
            do k = 1, n
               C(i,j) = C(i,j) + A(i,k) * B(k,j)
            enddo
         enddo
      enddo
   
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MultATxB2(
     &   A, B, C, n
     &   )

      implicit none
      include 'numbers.inc'
     
      integer n
      real*8  A(n, n), B(n, n), C(n, n)

      integer i, j, k
c
c---------------------------------------------------------------------72
c
      call SetTensor2(C, pzero, n*n)

      do i = 1, n 
         do j = 1, n
            do k = 1, n
               C(i,j) = C(i,j) + A(k,i) * B(k,j)
            enddo
         enddo
      enddo
   
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MultAxBT2(
     &   A, B, C, n
     &   )

      implicit none
      include 'numbers.inc'
     
      integer n
      real*8  A(n, n), B(n, n), C(n, n)

      integer i, j, k
c
c---------------------------------------------------------------------72
c
      call SetTensor2(C, pzero, n*n)

      do i = 1, n 
         do j = 1, n
            do k = 1, n
               C(i,j) = C(i,j) + A(i,k) * B(j,k)
            enddo
         enddo
      enddo
   
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE DeviatoricTensor2(
     &   tensor, tensorDev, n
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer n
      real*8  tensor(n*n), tensorDev(n*n)

      real*8  TraceOfTensor2, tensHyd
c
c---------------------------------------------------------------------72
c
      if (n .ne. 3)
     &   call RunTimeError2(XTAL_O, 'DeviatoricTensor: n =! 3')

      tensHyd = TraceOfTensor2(tensor, n) / pthree

      tensorDev(1) = tensor(1) - tensHyd
      tensorDev(5) = tensor(5) - tensHyd
      tensorDev(9) = tensor(9) - tensHyd

      tensorDev(2) = tensor(2)
      tensorDev(3) = tensor(3)
      tensorDev(6) = tensor(6)

      tensorDev(4) = tensor(4)
      tensorDev(7) = tensor(7)
      tensorDev(8) = tensor(8)
      
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      real*8 FUNCTION TraceOfTensor2(
     &   tensor, n
     &   )

      implicit none
      include 'params_xtal.inc'

      integer n
      real*8  tensor(n*n)
c
c---------------------------------------------------------------------72
c
      if (n .ne. 3)
     &   call RunTimeError2(XTAL_O, 'TraceOfTensor: n =! 3')

      TraceOfTensor2 = tensor(1) + tensor(5) + tensor(9)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      real*8 FUNCTION ScalarProduct2(
     &   tensA, tensB, n
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer n
      real*8  tensA(n*n), tensB(n*n)
c
c---------------------------------------------------------------------72
c
      if (n .ne. 3)
     &   call RunTimeError2(XTAL_O, 'ScalarProduct: n =! 3')

      ScalarProduct2 = tensA(1)*tensB(1) + tensA(2)*tensB(2) +
     &                tensA(3)*tensB(3) + tensA(4)*tensB(4) +
     &                tensA(5)*tensB(5) + tensA(6)*tensB(6) +
     &                tensA(7)*tensB(7) + tensA(8)*tensB(8) +
     &                tensA(9)*tensB(9)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE OuterProductVec2(
     &   vecU, vecV, outer, n
     &   )

      implicit none
      include 'numbers.inc'

      integer n
      real*8  vecU(n), vecV(n)
      real*8  outer(n, n)
     
      integer i, j
c
c---------------------------------------------------------------------72
c
      do i = 1, n
         do j = 1, n
            outer(i,j) = vecU(i) * vecV(j)
         enddo
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      real*8 FUNCTION InnerProductVec2(
     &   vecU, vecV, n
     &   )

      implicit none
      include 'numbers.inc'

      integer n
      real*8  vecU(n), vecV(n)

      integer i
c
c---------------------------------------------------------------------72
c
c---- compute inner product of two vectors: product = {u}^T*{v}
c
      InnerProductVec2 = pzero
      do i = 1, n
         InnerproductVec2 = InnerProductVec2 + vecU(i) * vecV(i)
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MultQAQT2(
     &   Q, A, B, n
     &   )

      implicit none

      integer n
      real*8  Q(n,n), A(n,n), B(n,n)

      real*8  tmp(n,n)
c
c---------------------------------------------------------------------72
c
      call MultAxB2(Q, A, tmp, n)       ! tmp = Q*A
      call MultAxBT2(tmp, Q, B, n)      ! B = tmp*QT

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MultQTAQ2(
     &   Q, A, B, n
     &   )

      implicit none

      integer n
      real*8  Q(n,n), A(n,n), B(n,n)

      real*8  tmp(n,n)
c
c---------------------------------------------------------------------72
c
      call MultATxB2(Q, A, tmp, n)      ! tmp = QT*A
      call MultAxB2(tmp, Q, B, n)       ! B = tmp*Q

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE EqualTensors2(
     &   tensA, tensB, n
     &   )

      implicit none

      integer n
      real*8  tensA(n), tensB(n)

      integer i
c
c---------------------------------------------------------------------72
c
      do i = 1, n
         tensB(i) = tensA(i)
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE AddTensors2(
     &   coef_A, tensA, coef_B, tensB, tensC, n
     &   )

      implicit none

      integer n
      real*8  coef_A, coef_B
      real*8  tensA(n), tensB(n), tensC(n)

      integer i
c
c---------------------------------------------------------------------72
c
       do i = 1, n
          tensC(i) = coef_A * tensA(i) + coef_B * tensB(i) 
       enddo

       return
       END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE SetToScaledTensor2(
     &   fac_A, tensA, tensB, n
     &   )

      implicit none

      integer n
      real*8  fac_A
      real*8  tensA(n), tensB(n)

      integer i
c
c---------------------------------------------------------------------72
c
      do i = 1, n
          tensB(i) = fac_A * tensA(i)
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE AddScaledTensor2(
     &   fac_A, tensA, tensB, n
     &   )

      implicit none

      integer n
      real*8  fac_A
      real*8  tensA(n), tensB(n)

      integer i
c
c---------------------------------------------------------------------72
c
      do i = 1, n
          tensB(i) = tensB(i) + fac_A * tensA(i)
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      real*8 FUNCTION DetOfMat3x32(
     &   tensor
     &   )

      implicit none
      
      real*8  tensor(3, 3)

      real*8  det11, det12, det13, det21, det22, det23
c
c---------------------------------------------------------------------72
c
c  Determinant of a second order tensor (3x3 matrix)
 
      det11 = tensor(1, 1) * tensor(2, 2) * tensor(3, 3)
      det12 = tensor(1, 2) * tensor(2, 3) * tensor(3, 1)
      det13 = tensor(1, 3) * tensor(2, 1) * tensor(3, 2)
      det21 = tensor(1, 1) * tensor(2, 3) * tensor(3, 2)
      det22 = tensor(1, 2) * tensor(2, 1) * tensor(3, 3)
      det23 = tensor(1, 3) * tensor(2, 2) * tensor(3, 1)
      
      DetOfMat3x32 = det11 + det12 + det13 - det21 - det22 - det23

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE InverseOfMat3x32(
     &   tensor, tensor_inv, tensor_det
     &   )

      implicit 	none

      real*8  tensor_det
      real*8  tensor(3, 3), tensor_inv(3, 3)

      real*8  tinv11, tinv12, tinv13
      real*8  tinv21, tinv22, tinv23
      real*8  tinv31, tinv32, tinv33
     
      real*8  DetOfMat3x32
c
c---------------------------------------------------------------------72
c
      tensor_det = DetOfMat3x32(tensor)

      tinv11 = tensor(2, 2) * tensor(3, 3) - tensor(2, 3) * tensor(3, 2)
      tinv12 = tensor(1, 3) * tensor(3, 2) - tensor(1, 2) * tensor(3, 3)
      tinv13 = tensor(1, 2) * tensor(2, 3) - tensor(1, 3) * tensor(2, 2)
      tinv21 = tensor(2, 3) * tensor(3, 1) - tensor(2, 1) * tensor(3, 3)
      tinv22 = tensor(1, 1) * tensor(3, 3) - tensor(1, 3) * tensor(3, 1)
      tinv23 = tensor(1, 3) * tensor(2, 1) - tensor(1, 1) * tensor(2, 3)
      tinv31 = tensor(2, 1) * tensor(3, 2) - tensor(2, 2) * tensor(3, 1)
      tinv32 = tensor(3, 1) * tensor(1, 2) - tensor(1, 1) * tensor(3, 2)
      tinv33 = tensor(1, 1) * tensor(2, 2) - tensor(1, 2) * tensor(2, 1)

      tensor_inv(1, 1) = tinv11 / tensor_det
      tensor_inv(1, 2) = tinv12 / tensor_det
      tensor_inv(1, 3) = tinv13 / tensor_det
      tensor_inv(2, 1) = tinv21 / tensor_det
      tensor_inv(2, 2) = tinv22 / tensor_det
      tensor_inv(2, 3) = tinv23 / tensor_det
      tensor_inv(3, 1) = tinv31 / tensor_det
      tensor_inv(3, 2) = tinv32 / tensor_det
      tensor_inv(3, 3) = tinv33 / tensor_det

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      real*8 FUNCTION SignOf2(
     &   value
     &   )

      implicit none

      real*8  value
c
c---------------------------------------------------------------------72
c
c------- compute the sign of a number
c
      SignOf2 = 1.0
      if (value .lt. 0.0) SignOf2 = -1.0

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      real*8 FUNCTION Power2(
     &   x, y
     &   )

      implicit none

      real*8  x, y
c
c---------------------------------------------------------------------72
c
c---- evaluates  x^y
c
      if (x .eq. 0.0) then
         if (y .gt. 0.0) then
            Power2 = 0.d0
         elseif (y .lt. 0.0) then
            Power2 = 1.d+300
         else
            Power2 = 1.d0
         endif
      else
         Power2 = y * log10(dabs(x))
         if (Power2 .gt. 300.0) then
            Power2 = 1.d+300
         else
            Power2 = 10.d0 ** Power2
         endif
         if (x .lt. 0.0) Power2 = -Power2
      endif

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE UnitVector2(
     &   vector, unitV, n
     &   )
      implicit none
      include 'params_xtal.inc'

      integer n
      real*8  vector(n), unitV(n)

      real*8  magnitude, InnerProductVec2
c
c---------------------------------------------------------------------72
c
      magnitude = dsqrt(InnerProductVec2(vector, vector, n))
      if (magnitude .le. 0.d0) 
     &        call RunTimeError2(XTAL_O,'UnitVector: magnitude <= 0.0')
      call SetToScaledTensor2(1./magnitude, vector, unitV, n)

      return
      END
c
c=====================================================================72
c
c=====================================================================72
c
      integer FUNCTION IndexMaxAbsValueOfVec2(
     &   vector, n
     &   )

      implicit none

      integer n
      real*8  vector(n)

      integer indexVal, i
      real*8  maxValue
c
c---------------------------------------------------------------------72
c
      indexVal = 1
      maxValue = dabs(vector(1))
      do i = 2, n
         if (dabs(vector(i)) .gt. maxValue) then
            indexVal = i
            maxValue = dabs(vector(i))
         endif
      enddo

      IndexMaxAbsValueOfVec2 = indexVal

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      real*8 FUNCTION MaxAbsValueOfVec2(
     &   vector, n
     &   )

      implicit none

      integer n
      real*8  vector(n)

      integer i
      real*8  maxValue
c
c---------------------------------------------------------------------72
c
      maxValue = dabs(vector(1))
      do i = 2, n
         if (dabs(vector(i)) .gt. maxValue) maxValue = dabs(vector(i))
      enddo

      MaxAbsValueOfVec2 = maxValue

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MatrixInverse2(
     &   a, n, np, det
     &   )
c
c---- routine borrowed from Hammid Youssef
c---- input:  a   : nonsym matrix to invert
c----         n   : order of matrix to invert
c----         np  : dimension of matrix in calling routine
c---- output: a   : matrix inverse
c----         det : determinant of matrix
c---- local:  k   : VECTEUR DE TRAVAIL ENTIER 
c
      IMPLICIT DOUBLE  PRECISION (A-H,O-Z)
      DIMENSION A(NP,NP),K(6)
      DATA ZERO,UN,EPS/0.D0,1.D0,1.D-13/
c
c---------------------------------------------------------------------72
c
c-------  Initialization
c
      DET=UN
      DO I=1,N
         K(I)=I
      END DO
c
c-------  inversion de la matrix A
c
      DO II=1,N
c
c-------  search for non-zero pivot in column II
c
         DO I=II,N
            XPIV=A(I,II)
            IF(DABS(XPIV).GT.EPS) GO TO 10
         END DO
         DET=ZERO
         GOTO 1000
c
c-------  exchange rows II and I
c
 10      DET=DET*XPIV
         IF(I.EQ.II) GO TO 20
         I1=K(II)
         K(II)=K(I)
         K(I)=I1
         DO J=1,N
            C=A(I,J)
            A(I,J)=A(II,J)
            A(II,J)=C
         END DO
         DET=-DET
c
c-------  normalize the row of the pivot
c
 20      C=UN/XPIV
         A(II,II)=UN
         DO J=1,N
            A(II,J)=A(II,J)*C
         END DO
c
c-------  Elimination
c
         DO I=1,N
            IF(I.EQ.II) GO TO 30
            C=A(I,II)
            A(I,II)=ZERO
            DO J=1,N
               A(I,J)=A(I,J)-C*A(II,J)
            END DO
 30      END DO
      END DO
c
c-------  re-order columns of inverse
c
      DO J=1,N
c-------  search J1 until K(J1)=J
         DO J1=J,N
            JJ=K(J1)
            IF(JJ.EQ.J) GO TO 100
         END DO
100      IF(J.EQ.J1) GO TO 110
c-------  exchange columns J and J1
         K(J1)=K(J)
         DO I=1,N
            C=A(I,J)
            A(I,J)=A(I,J1)
            A(I,J1)=C
         END DO
110   END DO

1000  CONTINUE

      return
      END
c
c=====================================================================72
c
c     
c=====================================================================72
c     
      SUBROUTINE lsolve2(a, x, nrhs, mystatus, n)
c     
c---- routine borrowed from Paul Dawson (Cornell University)
c---- Solve symmetric positive definite 5x5 linear systems. By calling 
c----  with the identity as right hand sides, you can use this routine
c----  to invert the matrix.
c
c     a    -- "n x n" matrix                                   (in/out)     
c     x    -- "n x nrhs" matrix, the array of right hand sides (in/out)
c     nrhs -- number of right hand sides                       (in)
c     mystatus -- return status                                (out)
c     n    -- size of system (must be equal to 5)
c     
      implicit none
      include 'params_xtal.inc'

      integer nrhs, mystatus, n
      real*8  a(n,n), x(n,*)
     
      integer i
      real*8  v1, v2, v3, v4, sum, not_zero, amax
c     
c---------------------------------------------------------------------72
c     
c------- check size of matrix
c
      if (n .ne. 5)
     &   call RunTimeError2(XTAL_O, 'lsolve : n .ne. 5')
c     
c------- set the scale according to the largest entry of a.
c     
      amax = 0.0d0
      amax = max(a(1,1), a(2,1), a(3,1), a(4,1), a(5,1))
      amax = max(a(2,2), a(3,2), a(4,2), a(5,2), amax)
      amax = max(a(3,3), a(4,3), a(5,3), amax)
      amax = max(a(4,4), a(5,4), amax)
      amax = max(a(5,5), amax)
     
      not_zero = amax * TINY
c     
c------- A = LDL'.
c
c------- At each step, check the size of each pivot
c---------- j = 1.
      if (a(1,1) .LT. not_zero) then
        mystatus = XTAL_SING_JACOBIAN
        return
      endif
     
      sum = 1.0/a(1,1)

      a(2,1) = a(2,1) * sum
      a(3,1) = a(3,1) * sum
      a(4,1) = a(4,1) * sum
      a(5,1) = a(5,1) * sum
c
c---------- j = 2.
      v1 = a(2,1)*a(1,1)
      a(2,2) = a(2,2) - a(2,1)*v1
      
      a(3,2) = a(3,2)-a(3,1)*v1
      a(4,2) = a(4,2)-a(4,1)*v1
      a(5,2) = a(5,2)-a(5,1)*v1

      if (a(2,2) .LT. not_zero) then
        mystatus = XTAL_SING_JACOBIAN
        return
      endif

      sum = 1.0/a(2,2)

      a(3,2) = a(3,2) * sum
      a(4,2) = a(4,2) * sum
      a(5,2) = a(5,2) * sum
c
c---------- j = 3.
      v1 = a(3,1)*a(1,1)
      v2 = a(3,2)*a(2,2)
      sum = a(3,1)*v1 + a(3,2)*v2

      a(3,3) = a(3,3) - sum
      
      if (a(3,3) .LT. not_zero) then
        mystatus = XTAL_SING_JACOBIAN
        return
      endif

      sum = 1.0/a(3,3)

      a(4,3) = a(4,3)-a(4,1)*v1-a(4,2)*v2
      a(5,3) = a(5,3)-a(5,1)*v1-a(5,2)*v2
      a(4,3) = a(4,3) * sum
      a(5,3) = a(5,3) * sum
c
c---------- j = 4.
      v1 = a(4,1)*a(1,1)
      v2 = a(4,2)*a(2,2)
      v3 = a(4,3)*a(3,3)
      sum = a(4,1)*v1+a(4,2)*v2+a(4,3)*v3

      a(4,4) = a(4,4) - sum

      if (a(4,4) .LT. not_zero) then
        mystatus = XTAL_SING_JACOBIAN
        return
      endif

      a(5,4) = a(5,4)-a(5,1)*v1-a(5,2)*v2-a(5,3)*v3
      a(5,4) = a(5,4)/a(4,4)
c
c---------- j = 5.
      v1 = a(5,1)*a(1,1)
      v2 = a(5,2)*a(2,2)
      v3 = a(5,3)*a(3,3)
      v4 = a(5,4)*a(4,4)
      sum = a(5,1)*v1+a(5,2)*v2+a(5,3)*v3+a(5,4)*v4

      a(5,5) = a(5,5) - sum

      if (a(5,5) .LT. not_zero) then
        mystatus = XTAL_SING_JACOBIAN
        return
      endif
c     
c------- solve for RHS
c
      do i=1, nrhs
c
c---------- Ly=b. 
        x(2, i) = x(2, i)
     &     - a(2,1)*x(1, i)
        x(3, i) = x(3, i)
     &     - a(3,1)*x(1, i) - a(3,2)*x(2, i)
        x(4, i) = x(4, i)
     &     - a(4,1)*x(1, i) - a(4,2)*x(2, i) - a(4,3)*x(3, i)
        x(5, i) = x(5, i)
     &     - a(5,1)*x(1, i) - a(5,2)*x(2, i) - a(5,3)*x(3, i)
     &     - a(5,4)*x(4, i)
c
c---------- Dz=y.
        x(1, i) = x(1, i)/a(1,1)
        x(2, i) = x(2, i)/a(2,2)
        x(3, i) = x(3, i)/a(3,3)
        x(4, i) = x(4, i)/a(4,4)
        x(5, i) = x(5, i)/a(5,5)
c
c---------- L'x=z.
        x(4, i) = x(4, i)
     &     - a(5,4)*x(5, i)
        x(3, i) = x(3, i)
     &     - a(4,3)*x(4, i) - a(5,3)*x(5, i)
        x(2, i) = x(2, i)
     &     - a(3,2)*x(3, i) - a(4,2)*x(4, i) - a(5,2)*x(5, i)
        x(1, i) = x(1, i)
     &     - a(2,1)*x(2, i) - a(3,1)*x(3, i) - a(4,1)*x(4, i)
     &     - a(5,1)*x(5, i)

      enddo

      return
      END
c     
c=====================================================================72
c     
c
c=====================================================================72
c
      SUBROUTINE Vec5x1ToVec6x12(
     &   vec5x1, vec6x1
     &   )

      implicit none
      include 'numbers.inc'

      real*8  vec5x1(5), vec6x1(6)
c
c---------------------------------------------------------------------72
c
      vec6x1(1) = 0.5 * (sqr2*vec5x1(1) - sqr23*vec5x1(2)) 
      vec6x1(2) =-0.5 * (sqr2*vec5x1(1) + sqr23*vec5x1(2)) 
      vec6x1(3) = vec5x1(2) * sqr23
      vec6x1(4) = vec5x1(3) / sqr2
      vec6x1(5) = vec5x1(4) / sqr2
      vec6x1(6) = vec5x1(5) / sqr2

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MultAxB_G2(
     &   A, B, C, m1, m2, m3
     &   )

      implicit none
      include 'numbers.inc'
     
      integer m1, m2, m3
      real*8  A(m1, m2), B(m2, m3), C(m1, m3)

      integer i, j, k
c
c---------------------------------------------------------------------72
c
c---- Dimensions: A is m1 x m2; B is m2 x m3; C is m1 x m3 
c
      call SetTensor2(C, pzero, m1*m3)

      do i = 1, m1 
         do j = 1, m3
            do k = 1, m2
               C(i,j) = C(i,j) + A(i,k) * B(k,j)
            enddo
         enddo
      enddo
   
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MultAxu_G2(
     &   A, u, v, m1, m2
     &   )

      implicit none
      include 'numbers.inc'

      integer m1, m2
      real*8  A(m1,m2), u(m2), v(m1)

      integer i, j
c
c---------------------------------------------------------------------72
c
c---- Dimensions: A is m1 x m2; u is m2; v is m1
c
      call SetTensor2(v, pzero, m1)

      do i = 1, m1
         do j = 1, m2
            v(i) = v(i) + A(i,j) * u(j)
         end do
      end do

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE MultATxu_G2(
     &   A, u, v, m1, m2
     &   )

      implicit none
      include 'numbers.inc'

      integer m1, m2
      real*8  A(m1,m2), u(m1), v(m2)

      integer i, j
c
c---------------------------------------------------------------------72
c
c---- Dimensions: A is m2 x m1; u is m1; v is m2
c
      call SetTensor2(v, pzero, m2)

      do i = 1, m2
         do j = 1, m1
            v(i) = v(i) + A(j,i) * u(j)
         end do
      end do

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE OuterProductVec_G2(
     &   vecU, vecV, outer, m1, m2
     &   )

      implicit none
      include 'numbers.inc'

      integer m1, m2
      real*8  vecU(m1), vecV(m2)
      real*8  outer(m1, m2)
     
      integer i, j
c
c---------------------------------------------------------------------72
c
      do i = 1, m1
         do j = 1, m2
            outer(i,j) = vecU(i) * vecV(j)
         enddo
      enddo

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      real*8 FUNCTION NormOfVector2(
     &   vec,  n
     &   )

      implicit none
      include 'numbers.inc'

      integer n
      real*8  vec(n)

      real*8  max_vec
      real*8  MaxAbsValueOfVec2, InnerProductVec2
c
c---------------------------------------------------------------------72
c
      max_vec = MaxAbsValueOfVec2(vec, n)
      if (max_vec .gt. pzero)
     &    call SetToScaledTensor2(pone/max_vec, vec, vec, n)

      NormOfVector2 = dsqrt(InnerProductVec2(vec, vec, n))

      if (max_vec .gt. pzero) then
         NormOfVector2 = max_vec*NormOfVector2
         call SetToScaledTensor2(max_vec, vec, vec, n)
      endif

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
c---- random number generator from numerical recipes
      FUNCTION ran1_nr2(idum)
      IMPLICIT NONE
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1_nr2,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     &    NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy/0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do j=NTAB+8,1,-1
           k=idum/IQ
           idum=IA*(idum-k*IQ)-IR*k
           if (idum.lt.0) idum=idum+IM
           if (j.le.NTAB) iv(j)=idum
        enddo
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1_nr2=min(AM*iy,RNMX)
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE AnglesToRotMatrix2(
     &   angle, crot, n
     &   )

      implicit none

      integer n
      real*8  angle(n)
      real*8  crot(n,n)

      real*8  sps, cps, sth, cth, sph, cph
c
c---------------------------------------------------------------------72
c
c------- Construct [C] matrix from euler angles (in radians)
c-------     {a}_sm = [C] {a}_xtal
c
      sps = dsin(angle(1))
      cps = dcos(angle(1))
      sth = dsin(angle(2))
      cth = dcos(angle(2))
      sph = dsin(angle(3))
      cph = dcos(angle(3))

      crot(1,1) = -sps * sph - cps * cph * cth
      crot(2,1) =  cps * sph - sps * cph * cth
      crot(3,1) =  cph * sth
      crot(1,2) =  cph * sps - sph * cps * cth
      crot(2,2) = -cps * cph - sps * sph * cth
      crot(3,2) =  sph * sth
      crot(1,3) =  cps * sth
      crot(2,3) =  sps * sth
      crot(3,3) =  cth

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE RotMatrixToAngles2(
     &   crot, angle, n
     &   )

      implicit none

      integer n
      real*8  angle(n)
      real*8  crot(n,n)
      
      real*8  sth
c
c---------------------------------------------------------------------72
c
c------- compute euler angles from [C] (in radians)
c
      angle(2) = acos(crot(3,3))
      if (dabs(crot(3,3)) .ne. 1.0) then
         sth = dsin(angle(2))
         angle(1) = datan2(crot(2,3)/sth, crot(1,3)/sth)
         angle(3) = datan2(crot(3,2)/sth, crot(3,1)/sth)
      else
         angle(1) = 0.
         angle(3) = datan2(-crot(1,2), -crot(1,1))
      endif

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      real*8 FUNCTION SSKineticEqn2(
     &   rss, crss, matProp, kflag
     &   )

      implicit none
      include 'params_xtal.inc'

      integer kflag
      real*8  crss
      real*8  matProp(NPROPS)

      real*8  xm, gam0, arg, pow, rss
      real*8  CheckArgPowLaw2, SignOf2, Power2
c
c---------------------------------------------------------------------72
c
c------- recover power law material constants
c
      xm   = matProp(3) 
      gam0 = matProp(4)
c
c------- check argument of power law
c
      arg = rss/crss
      if (dabs(arg) .ne. 0) arg = CheckArgPowLaw2(arg, SignOf2(arg))
c
c------- evaluate slip system kinetic equation 
c
      pow = Power2(dabs(arg), 1./xm-1.)
      if (kflag .eq. kGAMDOT) then 
         SSKineticEqn2 = gam0 * arg * pow
      else if (kflag .eq. kdGAMdTAU) then
         SSKineticEqn2 = gam0 / (xm*crss) * pow
      else if (kflag .eq. kdGAMdKAPP) then
         SSKineticEqn2 = -gam0 / (xm*crss) * arg * pow
      else
         call RunTimeError2(XTAL_O, 'SSKineticEqn: unkown kflag')
      endif

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      real*8 FUNCTION SSKineticEqnDrag2(
     &   rss, crss, matProp, kflag
     &   )

      implicit none
      include 'params_xtal.inc'

      integer kflag
      real*8  rss, crss
      real*8  matProp(NPROPS)

      real*8  xm, gam0, arg, pow, bdrag, ratio_drag, tmp1, tmp2
      real*8  gamdot_1, gamdot_2, gamdot_e, dgamdot_1, dgamdot_2
      real*8  CheckArgPowLaw2, SignOf2, Power2
c
c---------------------------------------------------------------------72
c
c------- recover power law material constants
c
      xm    = matProp(3) 
      gam0  = matProp(4)
      bdrag = matProp(11)
      ratio_drag = matProp(12)

      arg = rss/crss
      if (dabs(arg) .ge. ratio_drag) then
c
c------- flow rule based on dislocation drag kinetics
c
c---------- evaluate slip system kinetic equation 
         if (kflag .eq. kGAMDOT) then 
            SSKineticEqnDrag2 = arg / bdrag
         else if (kflag .eq. kdGAMdTAU) then
            SSKineticEqnDrag2 = 1.0 / (bdrag*crss)
         else if (kflag .eq. kdGAMdKAPP) then
            SSKineticEqnDrag2 = -1.0 / (bdrag*crss) * arg
         else
            call RunTimeError2(XTAL_O, 'SSKineticEqnDrag: unkown kflag')
         endif

      else
c
c------- flow rule based on power law & dislocation drag kinetics
c
c---------- check argument of power law
         if (dabs(arg) .ne. 0) arg = CheckArgPowLaw2(arg, SignOf2(arg))
c
c---------- evaluate slip system kinetic equation 
         pow = Power2(dabs(arg), 1./xm-1.)
         gamdot_1 = gam0 * arg * pow
         gamdot_2 = arg / bdrag
         gamdot_e =  gamdot_1*gamdot_2 / (gamdot_1 + gamdot_2)
         if (kflag .eq. kGAMDOT) then 
            SSKineticEqnDrag2 = gamdot_e
         else if (kflag .eq. kdGAMdTAU) then
            dgamdot_1 = gam0 / (xm*crss) * pow
            dgamdot_2 = 1.0 / (bdrag*crss)
            tmp1 = dgamdot_1/gamdot_1*gamdot_e/gamdot_1
            tmp2 = dgamdot_2/gamdot_2*gamdot_e/gamdot_2
            SSKineticEqnDrag2 = gamdot_e * (tmp1 + tmp2)
         else if (kflag .eq. kdGAMdKAPP) then
            dgamdot_1 = -gam0 / (xm*crss) * arg * pow
            dgamdot_2 = -1.0 / (bdrag*crss) * arg
            tmp1 = dgamdot_1/gamdot_1*gamdot_e/gamdot_1
            tmp2 = dgamdot_2/gamdot_2*gamdot_e/gamdot_2
            SSKineticEqnDrag2 = gamdot_e * (tmp1 + tmp2)
         else
            call RunTimeError2(XTAL_O, 'SSKineticEqnDrag: unkown kflag')
         endif
      endif

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE BoundForArgPowLaw2(
     &   xm
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      real*8  xm

      real*8  xmm

      real*8  argMin, argMax
      common /POWERLAWARG/ argMin, argMax

      real*8  EXPON
      data EXPON /280.d0/
c
c---------------------------------------------------------------------72
c
c---- limits (bounds) on the value of the argument for power law
c---- note: In general: xm <= 1.0 (1/xm >= 1.0)
c
      xmm = pone/xm - pone
      if (dabs(xm-pone) .lt. TINY) xmm = pone

      argMin = dexp(-EXPON/(xmm)*dlog(10.d0))
      argMax = dexp( EXPON/(xmm)*dlog(10.d0))

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      real*8 FUNCTION CheckArgPowLaw2(
     &   arg, sgn
     &   )
 
      implicit none
      
      real*8  arg, sgn

      real*8  argMin, argMax                 ! values assigned in
      common /POWERLAWARG/ argMin, argMax    ! BoundForArgPowerLaw
c
c---------------------------------------------------------------------72
c
c---- check range of argument for power law
c
      if (dabs(arg) .lt. argMin) then
         CheckArgPowLaw2 = argMin * sgn
      else if (dabs(arg) .ge. argMax) then
         CheckArgPowLaw2 = argMax * sgn
      else
         CheckArgPowLaw2 = arg
      endif

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE LimitRatioForDrag2(
     &   matProp
     &   )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      real*8  matProp(NPROPS)

      real*8  xm, gam0, bdrag, arg, argMax

      real*8  FACTOR
      data FACTOR /1.0d+10/

      real*8  EXPON
      data EXPON /280.d0/
c
c---------------------------------------------------------------------72
c
c------ (rss/crss)_limit to switch from (powerLaw & drag) to pure drag
c
      xm    = matProp(3)
      gam0  = matProp(4)
      bdrag = matProp(11)

      if (dabs(xm-pone) .lt. TINY) then
         matProp(12)=1.0d+300
         return
      endif

      arg = FACTOR/(bdrag*gam0)

      if (xm .gt. 0.700) then
         argMax = dexp(EXPON*(pone/xm-pone)*dlog(10.0d0))
         if (arg .gt. argMax) then
            matProp(12)=1.0d+300
            return
         endif
      endif

      matProp(12) = dexp(xm/(pone-xm)*dlog(arg))

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE RunTimeError2(
     &   io, message
     &   )

      implicit none
      
      character*(*) message
      integer io
c
c---------------------------------------------------------------------72
c
      write(io, 1000) message
      write(*, 1000) message
c      call closfl( )
      call CrystalCloseIOFiles2( )
      call CrystalCloseIOFiles_22( )
      stop

1000  format(/,'***ERROR Message: '/, 3x, a)

      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE WriteWarning2(
     &   io, message
     &   )

      implicit none
      
      character*(*) message
      integer io
c
c---------------------------------------------------------------------72
c
      write(io, 1000) message

1000  format(/,'***WARNING Message: '/, 3x, a)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE WriteMessage2(
     &   io, message
     &   )

      implicit none
      
      character*(*) message
      integer io
c
c---------------------------------------------------------------------72
c
      write(io, 1000) message

1000  format('***Message: ', a)

      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE CalculateSSD(gamdot, dtime, numslip, forden, forden_n,
     &                        incr, sub, sub_n)
     
      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'
      
      integer numslip, incr
      real*8  dtime
      real*8  gamdot(MAX_SLIP)
      real*8  forden(MAX_SLIP), forden_n(MAX_SLIP)
      real*8  sub(MAX_SLIP), sub_n(MAX_SLIP)
      
      real*8  k1, g, D, b, eps0
      real*8  T, k
      real*8  q, q0, q1
      
      real*8  k2, x
      real*8  gam_tot(MAX_SLIP), temp(MAX_SLIP)
      real*8  interim(MAX_SLIP)
      integer is
      real*8  fCeDev(NVECS, NVECS), fCeiDev(NVECS, NVECS)
      real*8  matProp(NPROPS, MAX_SLIP)
      real*8  fCeDevVol(NVECS), fCeVol, hardmtx(MAX_SLIP, MAX_SLIP)
      common /XtalProps / fCeDev, fCeiDev, matProp,
     &                    fCeDevVol, fCeVol, hardmtx 
      real*8  d_aba(3,3)
      common  /MacroStrainRate/ d_aba

      PARAMETER ( x=0.9 )   ! dislocation interaction factor, changed from 0.9
      PARAMETER ( k1=3.3d8 ) ! changed from  3.3d7
      PARAMETER ( g=1.25d-2 )  !an effective activation enthalpy, changed from  1.25d-2
      PARAMETER ( D=690.d6 )   ! a drag stress,changed from 900.d6
      PARAMETER ( b=2.56d-10 )
      PARAMETER ( eps0=1.d7 )    ! reference strain rate
      PARAMETER ( k=1.3806488d-23 )
      PARAMETER ( T=298.0 )
      PARAMETER ( q0=47.22 )
      PARAMETER ( q1=800 )
c
c---------------------------------------------------------------------72
c
      call SetTensor2(gam_tot, pzero, MAX_SLIP)
      call SetTensor2(temp, pzero, MAX_SLIP)

      q = q0*log(1+T/q1)
      if (incr .eq. 1) then
         k2= 1.5 !(k1*x*b/g)*(1.0-k*T/(D*(b**3))*log(1.d-4/eps0))
		 
      else
         k2= 1.5 !(k1*x*b/g)*(1.0-k*T/(D*(b**3))*log(matProp(4,1)/eps0))
		 
      endif
      
      do is=1, numslip
         gam_tot(is) = abs(gamdot(is))*dtime
         temp(is) = k1*sqrt(forden_n(is)) - k2*forden_n(is)
         forden(is) = forden_n(is) + temp(is) * gam_tot(is)
         interim(is) = q*b*sqrt(sub_n(is))*k2*forden_n(is)
         sub(is) = sub_n(is) + interim(is)*gam_tot(is)
      enddo

      
      return
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE UpdateCoordinates(coords, npt, noel)
      
      implicit none
      include 'params_xtal.inc'
      
      real*8  coords(3)
      integer npt, noel
      
      real*8  globalcoords1(3, 1, MAX_GRN)
      common  /IntPointCoords1/ globalcoords1
      save    /IntPointCoords1/
c
c---------------------------------------------------------------------72
c
      call EqualTensors2(coords, globalcoords1(1, 1, noel), 3)
      
      RETURN
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE CalculateGND(gamdot, numslip, npt, noel, gndden, 
     &                        gndden_n)
     
      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'
      
      real*8  gamdot(MAX_SLIP), gndden(MAX_SLIP), gndden_n(MAX_SLIP)
      integer numslip, npt, noel
      
      real*8  plasdef0(3, 3, 1, MAX_GRN)
      common  /PlasticDefGrad0/ plasdef0
      
      real*8  globalcoords0(3, 1, MAX_GRN)
      common  /IntPointCoords0/ globalcoords0
      
      integer neighborinfo(3, 1, MAX_GRN)
      common  /NeighInfo/ neighborinfo
      
      real*8  vecM(DIMS, MAX_SLIP), vecS(DIMS, MAX_SLIP)
      common  /SlipSystem/ vecM, vecS
      
      real*8  vecT(DIMS, MAX_SLIP)
      
      real*8  Fp_0(3,3), Fp_x(3,3), Fp_y(3,3), Fp_z(3,3)
      real*8  vec_0(DIMS), vec_x(DIMS), vec_y(DIMS), vec_z(DIMS)
      real*8  d_x, d_y, d_z
      real*8  tempv(DIMS), curl(DIMS)
      integer is, x, y, z
      real*8  rho_GNDs(MAX_SLIP), rho_GNDet(MAX_SLIP)
      real*8  rho_GNDen(MAX_SLIP), rho_GND(MAX_SLIP)
      real*8  b
      
      real*8 InnerProductVec2
      
      PARAMETER ( b=2.56d-10 )
c
c---------------------------------------------------------------------72
c 
c---------- find neighbor info in x, y and z directions
      x = neighborinfo(1, npt, noel)
      y = neighborinfo(2, npt, noel)
      z = neighborinfo(3, npt, noel)
      
      call EqualTensors2(plasdef0(1, 1, npt, noel), Fp_0(1,1), 9)
      call EqualTensors2(plasdef0(1, 1, x, noel), Fp_x(1,1), 9)
      call EqualTensors2(plasdef0(1, 1, y, noel), Fp_y(1,1), 9)
      call EqualTensors2(plasdef0(1, 1, z, noel), Fp_z(1,1), 9)
      
      d_x = ( globalcoords0(1, npt, noel) - globalcoords0(1, x, noel) )
     &      *1.d-3
      d_y = ( globalcoords0(2, npt, noel) - globalcoords0(2, y, noel) )
     &      *1.d-3
      d_z = ( globalcoords0(3, npt, noel) - globalcoords0(3, z, noel) )
     &      *1.d-3
      
      do is = 1, numslip
          call SetTensor2(vecT(1,is), pzero, DIMS)
          
          vecT(1,is) = vecM(2,is)*vecS(3,is) - vecM(3,is)*vecS(2,is)
          vecT(2,is) = vecM(3,is)*vecS(1,is) - vecM(1,is)*vecS(3,is)
          vecT(3,is) = vecM(1,is)*vecS(2,is) - vecM(2,is)*vecS(1,is)
      enddo
      
      do is = 1, numslip
          call SetTensor2(curl, pzero, DIMS)
          
          call MultATxu2(Fp_0,vecM(1,is),tempv,DIMS)
          call SetToScaledTensor2(gamdot(is),tempv,vec_0,DIMS)
          
          call MultATxu2(Fp_x,vecM(1,is),tempv,DIMS)
          call SetToScaledTensor2(gamdot(is),tempv,vec_x,DIMS)
          
          call MultATxu2(Fp_y,vecM(1,is),tempv,DIMS)
          call SetToScaledTensor2(gamdot(is),tempv,vec_y,DIMS)
          
          call MultATxu2(Fp_z,vecM(1,is),tempv,DIMS)
          call SetToScaledTensor2(gamdot(is),tempv,vec_z,DIMS)
          
          call CalculateCurl(vec_0, vec_x, vec_y, vec_z, d_x, d_y, d_z, 
     &                       curl)
          
          rho_GNDs(is) = (1.0/b)*InnerProductVec2(curl, vecS(1,is), 
     &                   DIMS)
          rho_GNDet(is) = (1.0/b)*InnerProductVec2(curl, vecT(1,is), 
     &                    DIMS)
          rho_GNDen(is) = (1.0/b)*InnerProductVec2(curl, vecM(1,is), 
     &                    DIMS)
          rho_GND(is) = sqrt(rho_GNDs(is)**2 + rho_GNDet(is)**2 + 
     &                       rho_GNDen(is)**2)
          
          gndden(is) = gndden_n(is) + rho_GND(is)
      enddo
      
      RETURN
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE UpdateFPandCoords()
      
      implicit none
      include 'params_xtal.inc'
      
      real*8  plasdef1(3, 3, 1, MAX_GRN)
      common  /PlasticDefGrad1/ plasdef1
      
      real*8  globalcoords1(3, 1, MAX_GRN)
      common  /IntPointCoords1/ globalcoords1
      
      real*8  plasdef0(3, 3, 1, MAX_GRN)
      common  /PlasticDefGrad0/ plasdef0
      save    /PlasticDefGrad0/
      
      real*8  globalcoords0(3, 1, MAX_GRN)
      common  /IntPointCoords0/ globalcoords0
      save    /IntPointCoords0/
      
      integer ig, ip
c
c---------------------------------------------------------------------72
c
      do ig = 1, MAX_GRN
c          do ip = 1, 8
              call EqualTensors2(plasdef1(1, 1, 1, ig), 
     &                           plasdef0(1, 1, 1, ig), 9)
              call EqualTensors2(globalcoords1(1, 1, ig),
     &                           globalcoords0(1, 1, ig), 3)
c          enddo
      enddo
      
      RETURN
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE UpdateFPandCoords2()
      
      implicit none
      include 'params_xtal.inc'
      
      real*8  plasdef1(3, 3, 1, MAX_GRN)
      common  /PlasticDefGrad1/ plasdef1
      
      real*8  globalcoords1(3, 1, MAX_GRN)
      common  /IntPointCoords1/ globalcoords1
      
      real*8  plasdef0(3, 3, 1, MAX_GRN)
      common  /PlasticDefGrad0/ plasdef0
      save    /PlasticDefGrad0/
      
      real*8  globalcoords0(3, 1, MAX_GRN)
      common  /IntPointCoords0/ globalcoords0
      save    /IntPointCoords0/
      
      integer ig, ip
c
c---------------------------------------------------------------------72
c
      do ig = 1, MAX_GRN
c          do ip = 1, 8
              call EqualTensors(plasdef1(1, 1, 1, ig), 
     &                           plasdef0(1, 1, 1, ig), 9)
              call EqualTensors(globalcoords1(1, 1, ig),
     &                           globalcoords0(1, 1, ig), 3)
c          enddo
      enddo
      
      RETURN
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE FindNeighbors()
      
      implicit none
      include 'params_xtal.inc'
      
      real*8  globalcoords1(3, 1, MAX_GRN)
      common  /IntPointCoords1/ globalcoords1
      
      integer neighborinfo(3, 8, MAX_GRN)
      common  /NeighInfo/ neighborinfo
      save    /NeighInfo/
      
      integer numel_aba, numqpt_aba
      common /data_aba/ numel_aba, numqpt_aba
      
      integer ig, ip, it, flag
      real*8  dx, dy, dz, toler
      PARAMETER ( toler=1.d-4 )
c
c---------------------------------------------------------------------72
c
      do ig = 801, numel_aba
          do ip = 1, 8
              flag = 0
              do it = 1, 8
                  if (ip .eq. it) CYCLE
                  dx = abs( (globalcoords1(1, it, ig) - 
     &                       globalcoords1(1, ip, ig)) / 
     &                       globalcoords1(1, ip, ig) )
                  dy = abs( (globalcoords1(2, it, ig) - 
     &                       globalcoords1(2, ip, ig)) / 
     &                       globalcoords1(2, ip, ig) )
                  dz = abs( (globalcoords1(3, it, ig) - 
     &                       globalcoords1(3, ip, ig)) / 
     &                       globalcoords1(3, ip, ig) )
                  if ( (dx .gt. toler) .and. (dy .lt. toler) .and. 
     &                 (dz .lt. toler) ) then
                      neighborinfo(1, ip, ig) = it
                      flag = flag + 1
                  endif
                  if ( (dx .lt. toler) .and. (dy .gt. toler) .and. 
     &                 (dz .lt. toler) ) then
                      neighborinfo(2, ip, ig) = it
                      flag = flag + 1
                  endif
                  if ( (dx .lt. toler) .and. (dy .lt. toler) .and. 
     &                 (dz .gt. toler) ) then
                      neighborinfo(3, ip, ig) = it
                      flag = flag + 1
                  endif
              enddo
              if (flag .ne. 3) then
                  write(XTAL_O, 1000) ip, ig
                  call RunTimeError2(XTAL_O, '')
              endif
          enddo
      enddo
      
1000  format('Did not find neighboring IP info at IP = ',i1,2x,
     &       'Ielem = ',i5)
      
      RETURN
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE CalculateCurl(vec_0, vec_x1, vec_y1, vec_z1, vec_x2,
     &                      vec_y2, vec_z2, d_x, d_y, d_z, curl)
      
      implicit none
      include 'params_xtal.inc'
      
      real*8  vec_0(DIMS), vec_x1(DIMS), vec_y1(DIMS), vec_z1(DIMS)
      real*8  vec_x2(DIMS), vec_y2(DIMS), vec_z2(DIMS)
      real*8  d_x, d_y, d_z
      real*8  curl(DIMS)
c
c---------------------------------------------------------------------72
c
      curl(1) = ((vec_y1(3)-vec_y2(3))/d_y)-((vec_z1(2)-vec_z2(2))/d_z)
      curl(2) = ((vec_z1(1)-vec_z2(1))/d_z)-((vec_x1(3)-vec_x2(3))/d_x)
      curl(3) = ((vec_x1(2)-vec_x2(2))/d_x)-((vec_y1(1)-vec_y2(1))/d_y)
      
      
      RETURN
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE UpdateKappa(kappa, forden, gndden, numslip, sub, incr, numgrn)
      
      implicit none
      include 'params_xtal.inc'
      
      real*8  kappa(NKAPP), forden(MAX_SLIP), gndden(MAX_SLIP)
      real*8  sub(MAX_SLIP), tau_for(MAX_SLIP)
      integer numslip, incr, numgrn, ig
      
      integer is, GIndex
      common  /newAddGIndex/ GIndex
      real*8  globalgrainsize(MAX_GRN)
      common  /newAddParameters/ globalgrainsize 
	  save    /newAddParameters/
      real*8  alpha, G, b, D, k, kHP, alpha_ref, den_ref, p,q   ! p, q are for alpha for the eq. 5 in the paper Devincre-2019
    
      
      PARAMETER ( alpha=0.5)    ! Changed from 0.3
      PARAMETER ( G=7.6d4 )
      PARAMETER ( b=2.56d-10 )
      PARAMETER ( k=0.086 )
      PARAMETER ( kHP=0.03 )
      PARAMETER ( D=50.d-9 )
	  PARAMETER ( alpha_ref = 0.405 )
	  PARAMETER ( den_ref = 10.d12)
	  
	 ! Initial dislocation density consists of forden(is) + gndden(is)
c
c---------------------------------------------------------------------72
c     
	  do is = 1, numslip
		!p = (log10(1/(alpha_ref*b*sqrt(forden(is) + gndden(is)))))
	    !q = (log10(1/(alpha_ref*b*sqrt(den_ref))))
	    !alpha = (p/q)*alpha_ref
!		alpha = (log10(1/(alpha_ref*b*sqrt(forden(is) + gndden(is) 
 !    &         )))/(log10(1/(alpha_ref*b*sqrt(den_ref))))*alpha_ref)
		tau_for(is) = alpha*G*b*sqrt(forden(is) + gndden(is))
     &                + kHP/sqrt(globalgrainsize(GIndex)) 
           
      enddo
      
cc---------update kappa   
      do is=1, numslip
         kappa(is) = tau_for(is)
      enddo
      
      RETURN
      END
c
c=====================================================================72
c=====================================================================72
c
      SUBROUTINE CalculateGND2(gamdot, numslip, npt, noel, gndden, 
     &                        gndden_n, incr, dtime)
     
      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'
      
      real*8  gamdot(MAX_SLIP), gndden(MAX_SLIP), gndden_n(MAX_SLIP)
      integer numslip, npt, noel, incr
      
      real*8  dtime
      
      real*8  plasdef0(3, 3, 1, MAX_GRN)
      common  /PlasticDefGrad0/ plasdef0
      
      real*8  globalcoords0(3, 1, MAX_GRN)
      common  /IntPointCoords0/ globalcoords0
      
      real*8  globalcoords1(3, 1, MAX_GRN)
      common  /IntPointCoords1/ globalcoords1
      
      integer neighborinfo(3, 1, MAX_GRN)
      common  /NeighInfo/ neighborinfo
      
      real*8  vecM(DIMS, MAX_SLIP), vecS(DIMS, MAX_SLIP)
      common  /SlipSystem/ vecM, vecS
      
      real*8  vecT(DIMS, MAX_SLIP)
      
      real*8  Fp_0(3,3), Fp_x1(3,3), Fp_y1(3,3), Fp_z1(3,3)
      real*8  Fp_x2(3,3), Fp_y2(3,3), Fp_z2(3,3)
      real*8  vec_0(DIMS), vec_x1(DIMS), vec_y1(DIMS), vec_z1(DIMS)
      real*8  vec_x2(DIMS), vec_y2(DIMS), vec_z2(DIMS)
      real*8  d_x, d_y, d_z
      real*8  tempv(DIMS), curl(DIMS)
      integer is, x, y, z
      real*8  rho_GNDs(MAX_SLIP), rho_GNDet(MAX_SLIP)
      real*8  rho_GNDen(MAX_SLIP), rho_GND(MAX_SLIP)
      real*8  b
      
      real*8 InnerProductVec2
      
      PARAMETER ( b=2.56d-10 )
      integer ms
      PARAMETER ( ms=6 )                                 ! The number of elements per edge
      real*8 uu,LL
      PARAMETER ( uu=55.4)
      PARAMETER ( LL=1.0)                                 !The length of the edge for the box
      
      call EqualTensors2(plasdef0(1, 1, 1, noel), Fp_0(1,1), 9)
	  
	  
	  
c---------------------------------------------------------------------72
c---------- find neighbor info in X direction
c---------------------------------------------------------------------72
      if ( (noel .GE. 1) .and. (noel .LE. ms*ms) ) then !X positive surface
        d_x = 2.0 * (globalcoords0(1, 1, noel) -
     & globalcoords0(1, 1, (noel+ms**2)) )*uu
      call EqualTensors2(plasdef0(1, 1, 1, (noel+ms**2)), Fp_x1(1,1), 9)
      call EqualTensors2(plasdef0(1, 1, 1, noel+ms**3-ms**2), Fp_x2(1,1)
     &, 9)  
     
      elseif ((noel .GT. (ms**3-ms**2)) .and. (noel .LE. ms**3) ) then !X negative surface
        d_x = 2.0 * ( globalcoords0(1, 1, (noel-ms**2)) -  
     &  globalcoords0(1, 1, noel) )*uu
      call EqualTensors2(plasdef0(1, 1, 1, noel-ms**3+ms**2), Fp_x1(1,1)
     & , 9)
      call EqualTensors2(plasdef0(1, 1, 1, (noel-ms**2)), Fp_x2(1,1), 9)
     
      else
        d_x = (globalcoords0(1, 1, (noel-ms**2)) -
     &   globalcoords0(1, 1, (noel+ms**2)) )*uu
      call EqualTensors2(plasdef0(1, 1, 1, (noel+ms**2)), Fp_x1(1,1), 9)
      call EqualTensors2(plasdef0(1, 1, 1, (noel-ms**2)), Fp_x2(1,1), 9)
     
      endif
     
c---------------------------------------------------------------------72
c---------- find neighbor info in Y direction
c---------------------------------------------------------------------72    
     
      if ( mod(noel,ms) .EQ. 1 ) then     !Y positive surface
        d_y = 2.0 * (globalcoords0(2, 1, noel) -  
     &  globalcoords0(2, 1, (noel+1)))*uu
      call EqualTensors2(plasdef0(1, 1, 1, (noel+1)), Fp_y1(1,1), 9)
      call EqualTensors2(plasdef0(1, 1, 1, noel+ms-1), Fp_y2(1,1), 9)
     
      elseif ( mod(noel,ms) .EQ. 0) then      !Y negative surface
        d_y = 2.0 * ( globalcoords0(2, 1, noel-1) -  
     &  globalcoords0(2, 1, noel))*uu
      call EqualTensors2(plasdef0(1, 1, 1, noel-ms+1),Fp_y1(1,1), 9)
      call EqualTensors2(plasdef0(1, 1, 1, (noel-1)), Fp_y2(1,1), 9)
     
      else 
        d_y = (globalcoords0(2, 1, (noel-1)) -
     &   globalcoords0(2, 1, (noel+1)))*uu    
      call EqualTensors2(plasdef0(1, 1, 1, (noel+1)), Fp_y1(1,1), 9)
      call EqualTensors2(plasdef0(1, 1, 1, (noel-1)), Fp_y2(1,1), 9)
     
      endif
     
c---------------------------------------------------------------------72
c---------- find neighbor info in Z direction
c---------------------------------------------------------------------72    
     
      if ( (noel .GE. 1) .and. (noel .LE. ms**3-ms**2+2*ms) .and. (mod(noel,ms**2) .GT. 0).and. (mod(noel,ms**2) .LE. ms)) then     !Z positive surface
        d_z = 2.0 * (globalcoords0(3, 1, noel) -  
     &   globalcoords0(3, 1, (noel+ms)) )*uu
      call EqualTensors2(plasdef0(1, 1, 1, (noel+ms)), Fp_z1(1,1), 9)
      call EqualTensors2(plasdef0(1, 1, 1, noel+ms**2-ms), Fp_z2(1,1),
     & 9)
     
      elseif ( (mod(noel,ms*ms) .GT. (ms*(ms-1)))
     & .OR. (mod(noel,ms*ms) .EQ. 0) ) then  !Z negative surface
        d_z = 2.0 * (  globalcoords0(3, 1, (noel-ms)) -  
     &  globalcoords0(3, 1, noel) )*uu
      call EqualTensors2(plasdef0(1, 1, 1, noel-ms**2+ms), Fp_z1(1,1),
     & 9)
      call EqualTensors2(plasdef0(1, 1, 1, (noel-ms)), Fp_z2(1,1), 9)
     
      else 
        d_z = (globalcoords0(3, 1, (noel-ms)) -
     &   globalcoords0(3, 1, (noel+ms)))*uu
      call EqualTensors2(plasdef0(1, 1, 1, (noel+ms)), Fp_z1(1,1), 9)
      call EqualTensors2(plasdef0(1, 1, 1, (noel-ms)), Fp_z2(1,1), 9)
     
      endif
      d_x=2*uu
	  d_y=2*uu
	  d_z=2*uu
c---------------------------------------------------------------------72
c---------------------------------------------------------------------72
      do is = 1, numslip
          call SetTensor2(vecT(1,is), pzero, DIMS)
          
          vecT(1,is) = vecM(2,is)*vecS(3,is) - vecM(3,is)*vecS(2,is)
          vecT(2,is) = vecM(3,is)*vecS(1,is) - vecM(1,is)*vecS(3,is)
          vecT(3,is) = vecM(1,is)*vecS(2,is) - vecM(2,is)*vecS(1,is)
      enddo
      
      do is = 1, numslip
          call SetTensor2(curl, pzero, DIMS)
          
          call Transpose3x3(Fp_0)
          call Transpose3x3(Fp_x1)
          call Transpose3x3(Fp_x2)
          call Transpose3x3(Fp_y1)
          call Transpose3x3(Fp_y2)
          call Transpose3x3(Fp_z1)
          call Transpose3x3(Fp_z2)
          
          call MultATxu2(Fp_0,vecM(1,is),tempv,DIMS)
          call SetToScaledTensor2(gamdot(is),tempv,vec_0,DIMS)
          
          call MultATxu2(Fp_x1,vecM(1,is),tempv,DIMS)
          call SetToScaledTensor2(gamdot(is),tempv,vec_x1,DIMS)
          
          call MultATxu2(Fp_x2,vecM(1,is),tempv,DIMS)
          call SetToScaledTensor2(gamdot(is),tempv,vec_x2,DIMS)
          
          call MultATxu2(Fp_y1,vecM(1,is),tempv,DIMS)
          call SetToScaledTensor2(gamdot(is),tempv,vec_y1,DIMS)
          
          call MultATxu2(Fp_y2,vecM(1,is),tempv,DIMS)
          call SetToScaledTensor2(gamdot(is),tempv,vec_y2,DIMS)
          
          call MultATxu2(Fp_z1,vecM(1,is),tempv,DIMS)
          call SetToScaledTensor2(gamdot(is),tempv,vec_z1,DIMS)
          
          call MultATxu2(Fp_z2,vecM(1,is),tempv,DIMS)
          call SetToScaledTensor2(gamdot(is),tempv,vec_z2,DIMS)
          
          call CalculateCurl(vec_0, vec_x1, vec_y1, vec_z1, vec_x2,
     &                      vec_y2, vec_z2, d_x, d_y, d_z, curl, incr)
          
          rho_GNDs(is) = (1.0/b)*InnerProductVec2(curl, vecS(1,is), 
     &                   DIMS)
          rho_GNDet(is) = (1.0/b)*InnerProductVec2(curl, vecT(1,is), 
     &                    DIMS)
          rho_GNDen(is) = (1.0/b)*InnerProductVec2(curl, vecM(1,is), 
     &                    DIMS)
          rho_GND(is) = sqrt(rho_GNDs(is)**2 + rho_GNDet(is)**2 + 
     &                       rho_GNDen(is)**2)
          
		  rho_GND(is) = 1d6*rho_GND(is)
          gndden(is) = gndden_n(is) + rho_GND(is)*dtime

      enddo
      
      
      RETURN
      END
c
c=====================================================================72
c
c
c=====================================================================72
c
      SUBROUTINE Initialelementcoordinate2(coords, noel)
      
      implicit none
      include 'params_xtal.inc'
      

      real*8  coords(3)
      integer noel
c
c---------------------------------------------------------------------72
c

      write(XTAL_G1, *)coords


      RETURN
      END
c
c=====================================================================72
c
c
c
c=====================================================================72
      SUBROUTINE CalculateForestDensity(density, numslip, forden, incr,
     &                                  noel, npt)
      
      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'
      
      integer is, numslip, incr, noel, npt, totEle, NincrF, ms
      real*8  forden(MAX_SLIP)
      real*8  density, totdensity, dentot, avgtot
      
      save    dentot, avgtot, NincrF
	  PARAMETER ( ms=6 )
c
c---------------------------------------------------------------------72
c
      do is = 1, numslip
        totdensity = totdensity + forden(is)
      enddo
      
      density = abs(totdensity / 12)     ! It was 11 before, I do not know why
c----------------------------------------------------------------------    
      if(noel .eq. 1) then
        dentot = 0
      endif
      
      dentot = dentot + density
      
      totEle = ms**3                        !!!! Make change for different mesh
	  
      if(noel .eq. totEle) then  
		avgtot = dentot / real(totEle)
		!write(XTAL_G1,*) 'Increment no.    =',incr
		!write(XTAL_G1,*) 'NincrF before if statement    =',NincrF
		if(incr .gt. NincrF) then
			NincrF = incr
			!write(XTAL_G1,*) 'NincrF     =',NincrF
			write(XTAL_G1,*) avgtot
		endif
      endif
c----------------------------------------------------------------------     

      return
      END
c
c=====================================================================72
c
c
c
c=====================================================================72
      SUBROUTINE CalculateGNDDensity(fdensity, numslip, gndden, incr,
     &                            noel, npt)
      
      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'
      
      integer is, numslip, incr, noel, npt, totEle, NincrG, ms 
      real*8  gndden(MAX_SLIP)
      real*8  fdensity, totfdensity, fdentot, favgtot
      
      save    fdentot, favgtot, NincrG
	  PARAMETER ( ms=6 )
c
c---------------------------------------------------------------------72
c
      do is = 1, numslip
        totfdensity = totfdensity + gndden(is)
      enddo
      
      fdensity = abs(totfdensity / 12)
c----------------------------------------------------------------------    
      if(noel .eq. 1) then
        fdentot = 0
      endif
      
      fdentot = fdentot + fdensity      
	  
	  totEle = ms**3                                    !!!! Make change for different mesh
      
      if(noel .eq. totEle) then  
        favgtot = fdentot / real(totEle)
		!write(XTAL_G1,*) ' NincrG    =',NincrG
		if(incr .gt. NincrG) then			
			NincrG = incr
			!write(XTAL_G1,*) ' NincrG    =',NincrG
			write(XTAL_G2,*) favgtot
        endif
      endif
c----------------------------------------------------------------------     

      return
      END
c
c=====================================================================72
c
      SUBROUTINE Initializegrain(noel)
     
      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'
      
      integer noel, ms, ig   
      real*8  globalgrainsize(MAX_GRN)
	  PARAMETER ( ms=6 )
	  
	 
	  
      common  /newAddParameters/ globalgrainsize  
      save    /newAddParameters/
	  
	  !f ((noel .GE. 1).and. (noel .LE. 500 then ! .and. (mod(noel, 100) .NE. 0) ) then
	  if ((noel .GE. 1) .and. (noel .LE.ms**3)) then
	  !if ((mod(noel, 100) .LT. 51).and. (mod(noel, 100) .GT. 0) ) then ! .and. (noel .LE.500)) then  ! works
		globalgrainsize(noel) =55.4*1.d-6
        !write(XTAL_O, *) 'noel  =  ', noel	  
		!write(XTAL_O, *) 'globalgrainsize(noel)  =  ', globalgrainsize(noel)
	  !elseif ((noel .GE. 15) .and. (noel .LE.27)) then
	  !elseif  ((mod(noel, 100) .GE. 51) .OR.(mod(noel, 100) .EQ. 0) )  then   ! works
	  ! else
	 !globalgrainsize(noel) = 1000.d-9
		!write(XTAL_O, *) 'noel  =  ', noel	  
		!write(XTAL_O, *) 'globalgrainsize(noel)  =  ', globalgrainsize(noel)
      endif

	  
	
        
      return
      END
      
      
      SUBROUTINE Transpose3x3(M3x3)
     
      implicit none
      include 'numbers.inc'

      real*8  M3x3(3,3), T

      integer i, j
c
c---------------------------------------------------------------------72

      do i = 1, 3
         do j = 1, 3
            T = M3x3(i,j)
            M3x3(i,j) = M3x3(j,i) 
            M3x3(j,i) = T
         end do
      end do

      return
      END
      
