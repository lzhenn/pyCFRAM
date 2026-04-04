  module rrtmg_lw_rad
      ! Modified from rrtmg_lw.1col.f90

      use parkind, only : im => kind_im, rb => kind_rb
      use parrrtm, only : mxlay, nbndlw, ngptlw, maxxsec, mxmol
      use rrlw_con, only: fluxfac, heatfac, oneminus, pi
      use rrlw_wvn, only: ng, ngb, nspa, nspb, wavenum1, wavenum2, delwave
      use rrlw_vsn
      use mcica_subcol_gen_lw, only: get_alpha, mcica_subcol_lw
      use rrtmg_lw_cldprop, only: cldprop
      use rrtmg_lw_cldprmc, only: cldprmc
      use rrtmg_lw_init, only: rrtmg_lw_ini
      use rrtmg_lw_rtrn, only: rtrn
      use rrtmg_lw_rtrnmr, only: rtrnmr
      use rrtmg_lw_rtrnmc, only: rtrnmc
      use rrtmg_lw_setcoef, only: setcoef
      use rrtmg_lw_taumol, only: taumol

      implicit none

! public interfaces
        public :: rrtmg_lw
        contains

      subroutine rrtmg_lw &
              (iaer, icld, nlayers, pavel, tavel, pz, tz, tbound, coldry, &
              wbrodl, wkl, pwvcm, semiss, cldfrac, ciwp, clwp, rel, rei, tauaer, &
              uflxsum, dflxsum, fnetsum, htrsum)
          
! ------- Declarations -------

! ----- Local -----

! Control
      integer(kind=im),intent(in) :: nlayers             ! total number of layers
      
      integer(kind=im) :: istart              ! beginning band of calculation
      integer(kind=im) :: iend                ! ending band of calculation
      integer(kind=im),intent(in) :: icld                ! clear/cloud flag
      integer(kind=im) :: iout                ! output option flag
      integer(kind=im),intent(in) :: iaer                ! aerosol option flag
      integer(kind=im) :: i                   ! output index
      integer(kind=im) :: ig                  ! g-point index
      integer(kind=im) :: iplon               ! column loop index
      integer(kind=im) :: ims                 ! mcica statistical loop index
      integer(kind=im) :: imca                ! flag for mcica [0=off, 1=on]
      integer(kind=im) :: nmca                ! number of mcica samples (mcica mode)
      integer(kind=im) :: irng                ! flag for random number generator
                                              ! [0=kissvec, 1=mersenne twister (default)]
      integer(kind=im) :: idrv                ! flag for calculation of dFdT, the change
                                              ! in upward flux as a function of surface 
                                              ! temperature [0=off, 1=on]
      integer(kind=im) :: lev, l              ! level indices
      integer(kind=im), parameter :: ncol = 1 ! total number of columns
      character page 

! Atmosphere
      real(kind=rb),intent(in) :: pavel(mxlay)           ! layer pressures (mb) 
      real(kind=rb),intent(in) :: tavel(mxlay)           ! layer temperatures (K)
      real(kind=rb),intent(in) :: pz(0:mxlay)            ! level (interface) pressures (hPa, mb)
      real(kind=rb),intent(in) :: tz(0:mxlay)            ! level (interface) temperatures (K)
      real(kind=rb),intent(in) :: tbound                 ! surface temperature (K)
      real(kind=rb),intent(in) :: coldry(mxlay)          ! dry air column density (mol/cm2)
      real(kind=rb),intent(in) :: wbrodl(mxlay)          ! broadening gas column density (mol/cm2)
      real(kind=rb),intent(in) :: wkl(mxmol,mxlay)       ! molecular amounts (mol/cm-2)
      real(kind=rb) :: wx(maxxsec,mxlay)      ! cross-section amounts (mol/cm-2)
      real(kind=rb),intent(in) :: pwvcm                  ! precipitable water vapor (cm)
      real(kind=rb),intent(in) :: semiss(nbndlw)         ! lw surface emissivity
      real(kind=rb) :: fracs(mxlay,ngptlw)    ! 
      real(kind=rb) :: taug(mxlay,ngptlw)     ! gaseous optical depths
      real(kind=rb) :: taut(mxlay,ngptlw)     ! gaseous + aerosol optical depths

      real(kind=rb),intent(in) :: tauaer(mxlay,nbndlw)   ! aerosol optical depth
!      real(kind=rb) :: ssaaer(mxlay,nbndlw)  ! aerosol single scattering albedo
                                              ! for future expansion 
                                              !   (lw aerosol scattering not yet available)
!      real(kind=rb) :: asmaer(mxlay,nbndlw)  ! aerosol asymmetry parameter
                                              ! for future expansion 
                                              !   (lw aerosol scattering not yet available)

! Atmosphere - setcoef
      integer(kind=im) :: laytrop             ! tropopause layer index
      integer(kind=im) :: jp(mxlay)           ! 
      integer(kind=im) :: jt(mxlay)           !
      integer(kind=im) :: jt1(mxlay)          !
      real(kind=rb) :: planklay(mxlay,nbndlw)   ! 
      real(kind=rb) :: planklev(0:mxlay,nbndlw) ! 
      real(kind=rb) :: plankbnd(nbndlw)       ! 
      real(kind=rb) :: dplankbnd_dt(nbndlw)   ! 

      real(kind=rb) :: colh2o(mxlay)          ! column amount (h2o)
      real(kind=rb) :: colco2(mxlay)          ! column amount (co2)
      real(kind=rb) :: colo3(mxlay)           ! column amount (o3)
      real(kind=rb) :: coln2o(mxlay)          ! column amount (n2o)
      real(kind=rb) :: colco(mxlay)           ! column amount (co)
      real(kind=rb) :: colch4(mxlay)          ! column amount (ch4)
      real(kind=rb) :: colo2(mxlay)           ! column amount (o2)
      real(kind=rb) :: colbrd(mxlay)          ! column amount (broadening gases)

      integer(kind=im) :: indself(mxlay)
      integer(kind=im) :: indfor(mxlay)
      real(kind=rb) :: selffac(mxlay)
      real(kind=rb) :: selffrac(mxlay)
      real(kind=rb) :: forfac(mxlay)
      real(kind=rb) :: forfrac(mxlay)

      integer(kind=im) :: indminor(mxlay)
      real(kind=rb) :: minorfrac(mxlay)
      real(kind=rb) :: scaleminor(mxlay)
      real(kind=rb) :: scaleminorn2(mxlay)

      real(kind=rb) :: &                      !
                         fac00(mxlay), fac01(mxlay), &
                         fac10(mxlay), fac11(mxlay) 
      real(kind=rb) :: &                      !
                         rat_h2oco2(mxlay),rat_h2oco2_1(mxlay), &
                         rat_h2oo3(mxlay),rat_h2oo3_1(mxlay), &
                         rat_h2on2o(mxlay),rat_h2on2o_1(mxlay), &
                         rat_h2och4(mxlay),rat_h2och4_1(mxlay), &
                         rat_n2oco2(mxlay),rat_n2oco2_1(mxlay), &
                         rat_o3co2(mxlay),rat_o3co2_1(mxlay)

! Atmosphere/clouds - cldprop
      integer(kind=im) :: ncbands             ! number of cloud spectral bands
      integer(kind=im) :: inflag              ! flag for cloud property method
      integer(kind=im) :: iceflag             ! flag for ice cloud properties
      integer(kind=im) :: liqflag             ! flag for liquid cloud properties

      real(kind=rb),intent(in) :: cldfrac(mxlay)         ! layer cloud fraction
      real(kind=rb) :: tauc(nbndlw,mxlay)     ! in-cloud optical depth (non-delta scaled)
!      real(kind=rb) :: ssac(nbndlw,mxlay)    ! in-cloud single scattering albedo (non-delta scaled)
                                              ! for future expansion 
                                              !   (lw scattering not yet available)
!      real(kind=rb) :: asmc(nbndlw,mxlay)    ! in-cloud asymmetry parameter (non-delta scaled)
                                              ! for future expansion 
                                              !   (lw scattering not yet available)
      real(kind=rb),intent(in) :: ciwp(mxlay)            ! in-cloud ice water path
      real(kind=rb),intent(in) :: clwp(mxlay)            ! in-cloud liquid water path
      real(kind=rb),intent(in) :: rei(mxlay)             ! cloud ice particle effective size (microns)
                                              ! specific definition of rei depends on setting of iceflag:
                                              ! iceflag = 0: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                              !              r_ec must be >= 10.0 microns
                                              ! iceflag = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                              !              r_ec range is limited to 13.0 to 130.0 microns
                                              ! iceflag = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996)
                                              !              r_k range is limited to 5.0 to 131.0 microns
                                              ! iceflag = 3: generalized effective size, dge, (Fu, 1996),
                                              !              dge range is limited to 5.0 to 140.0 microns
                                              !              [dge = 1.0315 * r_ec]
      real(kind=rb),intent(in) :: rel(mxlay)             ! cloud liquid particle effective radius (microns)

      real(kind=rb) :: taucloud(mxlay,nbndlw) ! in-cloud optical depth; delta scaled
!      real(kind=rb) :: ssacloud(mxlay,nbndlw)! in-cloud single scattering albedo; delta scaled
                                              ! for future expansion 
                                              !   (lw scattering not yet available)
!      real(kind=rb) :: asmcloud(mxlay,nbndlw)! in-cloud asymmetry parameter; delta scaled
                                              ! for future expansion 
                                              !   (lw scattering not yet available)

! Atmosphere/clouds - cldprmc [mcica]
      real(kind=rb) :: cldfmc(ngptlw,mxlay)   ! cloud fraction [mcica]
      real(kind=rb) :: ciwpmc(ngptlw,mxlay)   ! in-cloud ice water path [mcica]
      real(kind=rb) :: clwpmc(ngptlw,mxlay)   ! in-cloud liquid water path [mcica]
      real(kind=rb) :: relqmc(mxlay)          ! liquid particle effective radius (microns)
      real(kind=rb) :: reicmc(mxlay)          ! ice particle effective radius (microns)
      real(kind=rb) :: taucmc(ngptlw,mxlay)   ! in-cloud optical depth [mcica]
!      real(kind=rb) :: ssacmc(ngptlw,mxlay)  ! in-cloud single scattering albedo [mcica]
                                              ! for future expansion 
                                              !   (lw scattering not yet available)
!      real(kind=rb) :: asmcmc(ngptlw,mxlay)  ! in-cloud asymmetry parameter [mcica]
                                              ! for future expansion 

                                              !   (lw scattering not yet available)

      real(kind=rb) :: dtotuflux_dt(0:mxlay)  ! change in upward longwave flux (w/m2/k)
                                              ! with respect to surface temperature
      real(kind=rb) :: dtotuclfl_dt(0:mxlay)  ! change in clear sky upward longwave flux (w/m2/k)
                                              ! with respect to surface temperature
      real(kind=rb) :: alpha(mxlay)           ! vertical cloud fraction correlation parameter

! Parameters
      real(kind=rb), parameter :: cpdair = 1.004e3_rb  ! Specific heat capacity of dry air
                                                       ! at constant pressure at 273 K
                                                       ! (J kg-1 K-1)
! Output
      real(kind=rb) :: totuflux(0:mxlay)      ! upward longwave flux (w/m2)
      real(kind=rb) :: totdflux(0:mxlay)      ! downward longwave flux (w/m2)
      real(kind=rb) :: fnet(0:mxlay)          ! net longwave flux (w/m2)
      real(kind=rb) :: htr(0:mxlay)           ! longwave heating rate (k/day)
      real(kind=rb) :: totuclfl(0:mxlay)      ! clear sky upward longwave flux (w/m2)
      real(kind=rb) :: totdclfl(0:mxlay)      ! clear sky downward longwave flux (w/m2)
      real(kind=rb) :: fnetc(0:mxlay)         ! clear sky net longwave flux (w/m2)
      real(kind=rb) :: htrc(0:mxlay)          ! clear sky longwave heating rate (k/day)
! Output (mean output for McICA calculation)
      real(kind=rb),intent(out) :: uflxsum(0:mxlay)       ! upward longwave flux (w/m2)
      real(kind=rb),intent(out) :: dflxsum(0:mxlay)       ! downward longwave flux (w/m2)
      real(kind=rb),intent(out) :: fnetsum(0:mxlay)       ! net longwave flux (w/m2)
      real(kind=rb),intent(out) :: htrsum(0:mxlay)        ! longwave heating rate (k/day)

!
! Initializations

      oneminus = 1._rb - 1.e-6_rb
      pi = 2._rb * asin(1._rb)
      fluxfac = pi * 2.e4_rb                ! orig:   fluxfac = pi * 2.d4  
      page = char(12)

      uflxsum(0:) = 0._rb
      dflxsum(0:) = 0._rb
      fnetsum(0:) = 0._rb
      htrsum(0:) = 0._rb

! Set imca to select calculation type
!  (read by subroutine readprof from input file INPUT_RRTM):
! imca = 0, use standard forward model calculation
! imca = 1, use McICA for Monte Carlo treatment of sub-grid cloud variability

! Set irng to select random number generator for McICA (used when imca = 1)
! irng = 0, KISSVEC
! irng = 1, Mersenne Twister
!      irng = 0
      irng = 1

      call rrtmg_lw_ini(cpdair)

! This is the main longitude/column loop within rrtmg.

      do iplon = 1, ncol

! Input atmospheric profile from INPUT_RRTM.
         
         imca = 1
         iout = 0
         idrv = 0
         do l = 1,mxlay
            do i = 1, maxxsec
               wx(i,l) = 0.0
            end do
         end do
        inflag  = 2
        iceflag = 2
        liqflag = 1
        tauc    = 0

         istart = 1
         iend = 16

! Set nmca to sample size for Monte Carlo calculation
         if (imca.eq.1) nmca = 200

! This is the statistical sampling loop for McICA

         do ims = 1, nmca

! Call sub-colum cloud generator for McICA calculations
! Output will be written for all nmca samples.  This will be excessive if
! band output (iout=99) option is selected. 

            alpha(:) = 0.0_rb
            call mcica_subcol_lw(iplon, nlayers, icld, ims, irng, pavel, &
                    cldfrac, ciwp, clwp, rei, rel, tauc, alpha, &
                    cldfmc, ciwpmc, clwpmc, reicmc, relqmc, taucmc)

!  For cloudy atmosphere, use cldprop to set cloud optical properties based on
!  input cloud physical properties.  Select method based on choices described
!  in cldprop.  Cloud fraction, water path, liquid droplet and ice particle
!  effective radius must be passed into cldprop.  Cloud fraction and cloud
!  optical depth are transferred to rrtmg_lw arrays in cldprop.  

!  If McICA is requested use cloud fraction and cloud physical properties 
!  generated by sub-column cloud generator above. 

            call cldprmc(nlayers, inflag, iceflag, liqflag, cldfmc, &
                    ciwpmc, clwpmc, reicmc, relqmc, ncbands, taucmc)

! Calculate information needed by the radiative transfer routine
! that is specific to this atmosphere, especially some of the 
! coefficients and indices needed to compute the optical depths
! by interpolating data from stored reference atmospheres. 

            call setcoef(nlayers, istart, pavel, tavel, tz, tbound, semiss, &
                         coldry, wkl, wbrodl, &
                         laytrop, jp, jt, jt1, planklay, planklev, plankbnd, &
                         idrv, dplankbnd_dt, &
                         colh2o, colco2, colo3, coln2o, colco, colch4, colo2, &
                         colbrd, fac00, fac01, fac10, fac11, &
                         rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, &
                         rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, &
                         rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, &
                         selffac, selffrac, indself, forfac, forfrac, indfor, &
                         minorfrac, scaleminor, scaleminorn2, indminor)

!  Calculate the gaseous optical depths and Planck fractions for 
!  each longwave spectral band.

            call taumol(nlayers, pavel, wx, coldry, &
                        laytrop, jp, jt, jt1, planklay, planklev, plankbnd, &
                        colh2o, colco2, colo3, coln2o, colco, colch4, colo2, &
                        colbrd, fac00, fac01, fac10, fac11, &
                        rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, &
                        rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, &
                        rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, &
                        selffac, selffrac, indself, forfac, forfrac, indfor, &
                        minorfrac, scaleminor, scaleminorn2, indminor, &
                        fracs, taug)

! Combine gaseous and aerosol optical depths, if aerosol active
            if (iaer .eq. 0) then
               do i = 1, nlayers
                  do ig = 1, ngptlw
                     taut(i,ig) = taug(i,ig)
                  enddo
               enddo
            elseif (iaer .eq. 10) then
               do i = 1, nlayers
                  do ig = 1, ngptlw
                     taut(i,ig) = taug(i,ig) + tauaer(i,ngb(ig))
                  enddo
               enddo
            endif

! Call the radiative transfer routine.
! Either routine can be called to do clear sky calculation.  If clouds
! are present, then select routine based on cloud overlap assumption
! to be used.  Clear sky calculation is done simultaneously.
! For McICA, only RTRN is called for clear and cloudy calculations.

             call rtrnmc(nlayers, istart, iend, iout, pz, semiss, ncbands, &
                           cldfmc, taucmc, planklay, planklev, plankbnd, &
                           pwvcm, fracs, taut, &
                           totuflux, totdflux, fnet, htr, &
                           totuclfl, totdclfl, fnetc, htrc, &
                           idrv, dplankbnd_dt, dtotuflux_dt, dtotuclfl_dt )

               do i = nlayers, 0, -1
                  uflxsum(i) = uflxsum(i) + totuflux(i)
                  dflxsum(i) = dflxsum(i) + totdflux(i)
                  fnetsum(i) = fnetsum(i) + fnet(i)
                  htrsum(i) = htrsum(i) + htr(i)
               enddo

! Output average over samples when last sample reached.  Comment this if-check
! and swap flux output write statements below to output all McICA samples.
               if (ims .eq. nmca) then
                  do i = nlayers, 0, -1
                     uflxsum(i) = uflxsum(i)/nmca
                     dflxsum(i) = dflxsum(i)/nmca
                     fnetsum(i) = fnetsum(i)/nmca
                     htrsum(i) = htrsum(i)/nmca
                  enddo
               endif
! End statistical loop for McICA
         enddo
! End longitude/column loop
      enddo

     end subroutine rrtmg_lw
    
 end module rrtmg_lw_rad

