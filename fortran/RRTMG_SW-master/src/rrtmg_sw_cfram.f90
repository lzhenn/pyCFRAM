  module rrtmg_sw_rad
 ! Modified from rrtmg_sw.1col.f90
  
      use parkind, only : im => kind_im, rb => kind_rb
      use parrrsw, only : mxlay, nbndsw, ngptsw, naerec, nstr, nmol, mxmol, &
                          jpband, jpb1, jpb2
      use rrsw_aer, only : rsrtaua, rsrpiza, rsrasya
      use rrsw_con, only : heatfac, oneminus, pi
      use rrsw_wvn, only : wavenum1, wavenum2
      use rrsw_vsn
      use mcica_subcol_gen_sw, only: get_alpha, mcica_subcol_sw
      use rrtmg_sw_cldprop, only: cldprop_sw
      use rrtmg_sw_cldprmc, only: cldprmc_sw
      use rrtmg_sw_init, only: rrtmg_sw_ini
      use rrtmg_sw_setcoef, only: setcoef_sw
      use rrtmg_sw_spcvrt, only: spcvrt_sw
      use rrtmg_sw_spcvmc, only: spcvmc_sw

      implicit none
      public :: rrtmg_sw

      contains

      subroutine rrtmg_sw(iaer, icld, nlayers, pdp, pavel, tavel, pz, tz, tbound, semiss,&
              zenith, coldry, wbrodl, wkl, cldfrac, ciwp, clwp, rei, rel, tauaer, ssaaer, asmaer,&
              uflxsum, dflxsum, fnetsum, htrsum)

      integer(kind=im),intent(in) :: nlayers             ! total number of layers
      integer(kind=im) :: istart              ! beginning band of calculation
      integer(kind=im) :: iend                ! ending band of calculation
      integer(kind=im),intent(in) :: icld                ! clear/cloud and cloud overlap flag
      integer(kind=im) :: icpr                ! cldprop/cldprmc use flag
      integer(kind=im) :: iout                ! output option flag
      integer(kind=im),intent(in) :: iaer                ! aerosol option flag
      integer(kind=im) :: idelm               ! delta-m scaling flag
                                              ! [0 = direct and diffuse fluxes are unscaled]
                                              ! [1 = direct and diffuse fluxes are scaled]
      integer(kind=im) :: isccos              ! instrumental cosine response flag
      integer(kind=im) :: i                   ! layer loop index                      ! jk
      integer(kind=im) :: ib                  ! band loop index                       ! jsw
      integer(kind=im) :: ia, ig              ! indices
      integer(kind=im) :: iplon               ! column loop index                     ! jl
      integer(kind=im) :: ims                 ! mcica statistical loop index
      integer(kind=im) :: imca                ! flag for mcica [0=off, 1=on]
      integer(kind=im) :: nmca                ! number of mcica samples (mcica mode)
      integer(kind=im) :: irng                ! flag for random number generator
                                              ! [0=kissvec, 1=mersenne twister (default)]
      integer(kind=im), parameter :: ncol = 1 ! total number of columns

      real(kind=rb) :: zepsec, zepzen         ! epsilon
      real(kind=rb) :: zdpgcp                 ! flux to heating conversion ratio


! Atmosphere
      real(kind=rb),intent(in) :: pavel(mxlay)           ! layer pressures (mb) 
      real(kind=rb),intent(in) :: tavel(mxlay)           ! layer temperatures (K)
      real(kind=rb),intent(in) :: pz(0:mxlay)            ! level (interface) pressures (hPa, mb)
      real(kind=rb),intent(in) :: tz(0:mxlay)            ! level (interface) temperatures (K)
      real(kind=rb),intent(in) :: tbound                 ! surface temperature (K)
      real(kind=rb),intent(in) :: pdp(mxlay)             ! layer pressure thickness (hPa, mb)
      real(kind=rb),intent(in) :: coldry(mxlay)          ! 
      real(kind=rb),intent(in) :: wbrodl(mxlay)          !
      real(kind=rb),intent(in) :: wkl(mxmol,mxlay)       ! molecular amounts (mol/cm-2)

      ! real(kind=rb) :: cossza, zenith         ! cosine of solar zenith angle 
      real(kind=rb) :: cossza
      real(kind=rb),intent(in) :: zenith
      !      real(kind=rb) :: earth_sun              ! function for Earth/Sun distance factor
      real(kind=rb) :: adjflux(jpband)        ! adjustment for current Earth/Sun distance
      real(kind=rb) :: solvar(jpband)         ! solar constant scaling factor from rrtmg_sw
                                              !  default value of 1368.22 Wm-2 at 1 AU
      real(kind=rb),intent(in) :: semiss(jpband)         ! surface emissivity
      real(kind=rb) :: albdir(nbndsw)         ! surface albedo, direct          ! zalbp
      real(kind=rb) :: albdif(nbndsw)         ! surface albedo, diffuse         ! zalbd

      real(kind=rb),intent(in) :: tauaer(mxlay,jpband)   ! aerosol optical depth (iaer=10 only)
                                              ! (non-delta scaled)      
      real(kind=rb),intent(in) :: ssaaer(mxlay,jpband)   ! aerosol single scattering albedo (iaer=10 only)
                                              ! (non-delta scaled)      
      real(kind=rb),intent(in) :: asmaer(mxlay,jpband)   ! aerosol asymmetry parameter (iaer=10 only)
                                              ! (non-delta scaled)      
                                              !   first moment of input phase function
      real(kind=rb) :: ecaer(mxlay,naerec)    ! aerosol optical thickness at 0.55 micron (iaer=6 only)
                                              ! (non-delta scaled)      

! Atmosphere - setcoef
      integer(kind=im) :: laytrop             ! tropopause layer index
      integer(kind=im) :: layswtch            ! tropopause layer index
      integer(kind=im) :: laylow              ! tropopause layer index
      integer(kind=im) :: jp(mxlay)           ! 
      integer(kind=im) :: jt(mxlay)           !
      integer(kind=im) :: jt1(mxlay)          !

      real(kind=rb) :: colh2o(mxlay)          ! column amount (h2o)
      real(kind=rb) :: colco2(mxlay)          ! column amount (co2)
      real(kind=rb) :: colo3(mxlay)           ! column amount (o3)
      real(kind=rb) :: coln2o(mxlay)          ! column amount (n2o)
      real(kind=rb) :: colch4(mxlay)          ! column amount (ch4)
      real(kind=rb) :: colo2(mxlay)           ! column amount (o2)
      real(kind=rb) :: colmol(mxlay)          ! column amount
      real(kind=rb) :: co2mult(mxlay)         ! column amount 

      integer(kind=im) :: indself(mxlay)
      integer(kind=im) :: indfor(mxlay)
      real(kind=rb) :: selffac(mxlay)
      real(kind=rb) :: selffrac(mxlay)
      real(kind=rb) :: forfac(mxlay)
      real(kind=rb) :: forfrac(mxlay)

      real(kind=rb) :: &                      !
                         fac00(mxlay), fac01(mxlay), &
                         fac10(mxlay), fac11(mxlay) 

! Atmosphere/clouds - cldprop
      integer(kind=im) :: ncbands             ! number of cloud spectral bands
      integer(kind=im) :: inflag              ! flag for cloud property method
      integer(kind=im) :: iceflag             ! flag for ice cloud properties
      integer(kind=im) :: liqflag             ! flag for liquid cloud properties

      real(kind=rb),intent(in) :: cldfrac(mxlay)         ! layer cloud fraction
      real(kind=rb) :: tauc(nbndsw,mxlay)     ! in-cloud optical depth (non-delta scaled)
      real(kind=rb) :: ssac(nbndsw,mxlay)     ! in-cloud single scattering albedo (non-delta scaled)
      real(kind=rb) :: asmc(nbndsw,mxlay)     ! in-cloud asymmetry parameter (non-delta scaled)
      real(kind=rb) :: fsfc(nbndsw,mxlay)     ! in-cloud forward scattering fraction (non-delta scaled)
      real(kind=rb),intent(in) :: ciwp(mxlay)            ! in-cloud ice water path
      real(kind=rb),intent(in) :: clwp(mxlay)            ! in-cloud liquid water path
      real(kind=rb),intent(in) :: rei(mxlay)             ! cloud ice particle effective size (microns)
      real(kind=rb),intent(in) :: rel(mxlay)             ! cloud liquid particle effective radius (microns)

      real(kind=rb) :: taucloud(mxlay,jpband) ! in-cloud optical depth
      real(kind=rb) :: taucldorig(mxlay,jpband)! in-cloud optical depth (non-delta scaled)
      real(kind=rb) :: ssacloud(mxlay,jpband) ! in-cloud single scattering albedo
      real(kind=rb) :: asmcloud(mxlay,jpband) ! in-cloud asymmetry parameter

! Atmosphere/clouds - cldprmc [mcica]
      real(kind=rb) :: cldfmc(ngptsw,mxlay)   ! cloud fraction [mcica]
      real(kind=rb) :: ciwpmc(ngptsw,mxlay)   ! in-cloud ice water path [mcica]
      real(kind=rb) :: clwpmc(ngptsw,mxlay)   ! in-cloud liquid water path [mcica]
      real(kind=rb) :: relqmc(mxlay)          ! liquid particle effective radius (microns)
      real(kind=rb) :: reicmc(mxlay)          ! ice particle effective radius (microns)
      real(kind=rb) :: taucmc(ngptsw,mxlay)   ! in-cloud optical depth [mcica]
      real(kind=rb) :: taormc(ngptsw,mxlay)   ! unscaled in-cloud optical depth [mcica]
      real(kind=rb) :: ssacmc(ngptsw,mxlay)   ! in-cloud single scattering albedo [mcica]
      real(kind=rb) :: asmcmc(ngptsw,mxlay)   ! in-cloud asymmetry parameter [mcica]
      real(kind=rb) :: fsfcmc(ngptsw,mxlay)   ! in-cloud forward scattering fraction [mcica]

! Atmosphere/clouds/aerosol - spcvrt,spcvmc
      real(kind=rb) :: ztauc(mxlay,nbndsw)    ! cloud optical depth
      real(kind=rb) :: ztaucorig(mxlay,nbndsw)! unscaled cloud optical depth
      real(kind=rb) :: zasyc(mxlay,nbndsw)    ! cloud asymmetry parameter 
                                              !  (first moment of phase function)
      real(kind=rb) :: zomgc(mxlay,nbndsw)    ! cloud single scattering albedo
      real(kind=rb) :: ztaua(mxlay,nbndsw)    ! total aerosol optical depth
      real(kind=rb) :: zasya(mxlay,nbndsw)    ! total aerosol asymmetry parameter 
      real(kind=rb) :: zomga(mxlay,nbndsw)    ! total aerosol single scattering albedo

      real(kind=rb) :: zcldfmc(mxlay,ngptsw)  ! cloud fraction [mcica]
      real(kind=rb) :: ztaucmc(mxlay,ngptsw)  ! cloud optical depth [mcica]
      real(kind=rb) :: ztaormc(mxlay,ngptsw)  ! unscaled cloud optical depth [mcica]
      real(kind=rb) :: zasycmc(mxlay,ngptsw)  ! cloud asymmetry parameter [mcica] 
      real(kind=rb) :: zomgcmc(mxlay,ngptsw)  ! cloud single scattering albedo [mcica]

      real(kind=rb) :: zbbfu(mxlay+1)         ! temporary upward shortwave flux (w/m2)
      real(kind=rb) :: zbbfd(mxlay+1)         ! temporary downward shortwave flux (w/m2)
      real(kind=rb) :: zbbcu(mxlay+1)         ! temporary clear sky upward shortwave flux (w/m2)
      real(kind=rb) :: zbbcd(mxlay+1)         ! temporary clear sky downward shortwave flux (w/m2)
      real(kind=rb) :: zbbfddir(mxlay+1)      ! temporary downward direct shortwave flux (w/m2)
      real(kind=rb) :: zbbcddir(mxlay+1)      ! temporary clear sky downward direct shortwave flux (w/m2)
      real(kind=rb) :: zuvfd(mxlay+1)         ! temporary UV downward shortwave flux (w/m2)
      real(kind=rb) :: zuvcd(mxlay+1)         ! temporary clear sky UV downward shortwave flux (w/m2)
      real(kind=rb) :: zuvfddir(mxlay+1)      ! temporary UV downward direct shortwave flux (w/m2)
      real(kind=rb) :: zuvcddir(mxlay+1)      ! temporary clear sky UV downward direct shortwave flux (w/m2)
      real(kind=rb) :: znifd(mxlay+1)         ! temporary near-IR downward shortwave flux (w/m2)
      real(kind=rb) :: znicd(mxlay+1)         ! temporary clear sky near-IR downward shortwave flux (w/m2)
      real(kind=rb) :: znifddir(mxlay+1)      ! temporary near-IR downward direct shortwave flux (w/m2)
      real(kind=rb) :: znicddir(mxlay+1)      ! temporary clear sky near-IR downward direct shortwave flux (w/m2)

      real(kind=rb) :: alpha(mxlay)           ! vertical cloud fraction correlation parameter

! Parameters
      real(kind=rb), parameter :: cpdair = 1.004e3_rb  ! Specific heat capacity of dry air
                                                       ! at constant pressure at 273 K (J kg-1 K-1)
! Output
      real(kind=rb) :: totuflux(0:mxlay)      ! upward shortwave flux (w/m2)                  ! pfup
      real(kind=rb) :: totdflux(0:mxlay)      ! downward shortwave flux (w/m2)                ! pfdown
      real(kind=rb) :: fnet(0:mxlay)          ! net shortwave flux (w/m2)                     ! pfls
      real(kind=rb) :: htr(0:mxlay)           ! shortwave heating rate (k/day)                ! pheat
      real(kind=rb) :: totuclfl(0:mxlay)      ! clear sky upward shortwave flux (w/m2)        ! pcup 
      real(kind=rb) :: totdclfl(0:mxlay)      ! clear sky downward shortwave flux (w/m2)      ! pcdown 
      real(kind=rb) :: fnetc(0:mxlay)         ! clear sky net shortwave flux (w/m2)           ! pfcs
      real(kind=rb) :: htrc(0:mxlay)          ! clear sky shortwave heating rate (k/day)      ! pheac

      real(kind=rb) :: dirdflux(0:mxlay)      ! direct downward shortwave flux (w/m2)         ! dirdownflux
      real(kind=rb) :: difdflux(0:mxlay)      ! diffuse downward shortwave flux (w/m2)        ! difdownflux
      real(kind=rb) :: dflxuv(0:mxlay)        ! Total sky downward shortwave flux, UV/vis     ! pfdnuv
      real(kind=rb) :: dflxir(0:mxlay)        ! Total sky downward shortwave flux, near-IR    ! pfdnir 
      real(kind=rb) :: dirdnuv(0:mxlay)       ! Direct downward shortwave surface flux, UV/vis
      real(kind=rb) :: difdnuv(0:mxlay)       ! Diffuse downward shortwave surface flux, UV/vis
      real(kind=rb) :: dirdnir(0:mxlay)       ! Direct downward shortwave surface flux, near-IR
      real(kind=rb) :: difdnir(0:mxlay)       ! Diffuse downward shortwave surface flux, near-IR

! Output (mean output for McICA calculation)
      real(kind=rb) :: dirdsum(0:mxlay)       ! direct downward shortwave flux (w/m2)
      real(kind=rb) :: difdsum(0:mxlay)       ! diffuse downward shortwave flux (w/m2)
      real(kind=rb), intent(out) :: uflxsum(0:mxlay)       ! upward shortwave flux (w/m2)
      real(kind=rb), intent(out) :: dflxsum(0:mxlay)       ! downward shortwave flux (w/m2)
      real(kind=rb), intent(out) :: fnetsum(0:mxlay)       ! net shortwave flux (w/m2)
      real(kind=rb), intent(out) :: htrsum(0:mxlay)        ! shortwave heating rate (k/day)

! Solar variability
      integer(kind=im) :: isolvar             ! Flag for solar variability method
      real(kind=rb) :: svar_f                 ! Solar variability facular multiplier
      real(kind=rb) :: svar_s                 ! Solar variability sunspot multiplier
      real(kind=rb) :: svar_i                 ! Solar variability baseline irradiance multiplier
      real(kind=rb) :: svar_f_bnd(jpband)     ! Solar variability facular multiplier (by band)
      real(kind=rb) :: svar_s_bnd(jpband)     ! Solar variability sunspot multiplier (by band)
      real(kind=rb) :: svar_i_bnd(jpband)     ! Solar variability baseline irradiance multiplier (by band)

! Initializations

      zepsec = 1.e-06_rb
      zepzen = 1.e-10_rb
      oneminus = 1.0_rb - zepsec
      pi = 2._rb * asin(1._rb)

      icpr = 1
      
      uflxsum(0:) = 0._rb
      dflxsum(0:) = 0._rb
      dirdsum(0:) = 0._rb
      difdsum(0:) = 0._rb
      fnetsum(0:) = 0._rb
      htrsum(0:) = 0._rb

! Set irng to select random number generator for McICA (use when imca = 1)
! irng = 0, KISSVEC
! irng = 1, Mersenne Twister
!      irng = 0
      irng = 1

! Call model and data initialization, compute lookup tables, perform
! reduction of g-points from 224 to 112 for input absorption
! coefficient data and other arrays.

      call rrtmg_sw_ini(cpdair)

! This is the main longitude/column loop within rrtmg.
      do iplon = 1, ncol

!--------------------------
! Pass in from driver
! Namelist
          iout = 0 
          imca = 1
          isccos = 0
          idelm = 1 
          adjflux = 1.0
          inflag = 2
          iceflag = 2
          liqflag = 1
          tauc = 0.0
          ssac = 0.0
          asmc = 0.0
          fsfc = 0.0
          isolvar = 0
          svar_f = 1
          svar_s = 1
          svar_i = 1
!-----------------------------

         istart = jpb1                 ! jpb1 = 16
         iend = jpb2                   ! jpb2 = 29
! Set nmca to sample size for Monte Carlo calculation
         if (imca.eq.1) nmca = 200
! This is the statistical sampling loop for McICA

         do ims = 1, nmca

! Call sub-colum cloud generator for McICA calculations.
! Output will be averaged over all nmca samples.  The code can be modified to
! write output for each individual sample (this will be excessive if output
! is also requested for each spectral band).  

               alpha(:) = 0.0_rb
               call mcica_subcol_sw(iplon, nlayers, icld, ims, irng, pavel, &
                          cldfrac, ciwp, clwp, rei, rel, tauc, ssac, asmc, fsfc, &
                          alpha, cldfmc, ciwpmc, clwpmc, reicmc, relqmc, taucmc, &
                          ssacmc, asmcmc, fsfcmc)

!  For cloudy atmosphere, use cldprop to set cloud optical properties based on
!  input cloud physical properties.  Select method based on choices described
!  in cldprop.  Cloud fraction, water path, liquid droplet and ice particle
!  effective radius must be passed in cldprop.  Cloud fraction and cloud
!  optical properties are transferred to rrtmg_sw arrays in cldprop.  
!  Note: Model will be stopped if partial cloud present without McICA.

!  If McICA is requested use cloud fraction and cloud physical properties 
!  generated by sub-column cloud generator above. 

               call cldprmc_sw(nlayers, inflag, iceflag, liqflag, cldfmc, &
                               ciwpmc, clwpmc, reicmc, relqmc, &
                               taormc, taucmc, ssacmc, asmcmc, fsfcmc)

! Calculate coefficients for the temperature and pressure dependence of the 
! molecular absorption coefficients by interpolating data from stored
! reference atmospheres.

            call setcoef_sw(nlayers, pavel, tavel, pz, tz, tbound, coldry, wkl, &
                            laytrop, layswtch, laylow, jp, jt, jt1, &
                            co2mult, colch4, colco2, colh2o, colmol, coln2o, &
                            colo2, colo3, fac00, fac01, fac10, fac11, &
                            selffac, selffrac, indself, forfac, forfrac, indfor)

! Cosine of the solar zenith angle 
!  Prevent using value of zero;

            cossza = zenith
            if (cossza.eq.0._rb) cossza = zepzen

! Transfer albedo, cloud and aerosol properties into arrays for 2-stream radiative transfer 
  
! Surface albedo
            do ib=1,nbndsw
               albdif(ib) = 1._rb - semiss(jpb1-1+ib)
               albdir(ib) = 1._rb - semiss(jpb1-1+ib)
            enddo

! Clouds
            if (icld.eq.0) then

               ztauc(:,:) = 0._rb
               ztaucorig(:,:) = 0._rb
               zasyc(:,:) = 0._rb
               zomgc(:,:) = 1._rb
               zcldfmc(:,:) = 0._rb
               ztaucmc(:,:) = 0._rb
               ztaormc(:,:) = 0._rb
               zasycmc(:,:) = 0._rb
               zomgcmc(:,:) = 1._rb

            elseif (icld.ge.1) then
                do i=1,nlayers
                    do ig=1,ngptsw
                        zcldfmc(i,ig) = cldfmc(ig,i)
                        ztaucmc(i,ig) = taucmc(ig,i)
                        ztaormc(i,ig) = taormc(ig,i)
                        zasycmc(i,ig) = asmcmc(ig,i)
                        zomgcmc(i,ig) = ssacmc(ig,i)
                     enddo
                enddo
            endif   

! Aerosol
! IAER = 0: no aerosols
            if (iaer.eq.0) then

               ztaua(:,:) = 0._rb
               zasya(:,:) = 0._rb
               zomga(:,:) = 1._rb

! IAER=10: Direct specification of aerosol properties from IN_AER_RRTM.
            elseif (iaer.eq.10) then

               do i = 1 ,nlayers
                  do ib = 1 ,nbndsw
                     ztaua(i,ib) = tauaer(i,jpb1-1+ib)
                     zasya(i,ib) = asmaer(i,jpb1-1+ib)
                     zomga(i,ib) = ssaaer(i,jpb1-1+ib)
                  enddo
               enddo

            endif

! Call the 2-stream radiation transfer model

            do i=1,nlayers+1
               zbbcu(i) = 0._rb
               zbbcd(i) = 0._rb
               zbbfu(i) = 0._rb
               zbbfd(i) = 0._rb
               zbbcddir(i) = 0._rb
               zbbfddir(i) = 0._rb
               zuvcd(i) = 0._rb
               zuvfd(i) = 0._rb
               zuvcddir(i) = 0._rb
               zuvfddir(i) = 0._rb
               znicd(i) = 0._rb
               znifd(i) = 0._rb
               znicddir(i) = 0._rb
               znifddir(i) = 0._rb
            enddo

            call spcvmc_sw &
                (nlayers, istart, iend, icpr, idelm, iout, &
                pavel, tavel, pz, tz, tbound, albdif, albdir, &
                zcldfmc, ztaucmc, zasycmc, zomgcmc, ztaormc, &
                ztaua, zasya, zomga, cossza, coldry, wkl, adjflux, &	 
                isolvar, svar_f, svar_s, svar_i, &
                svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                laytrop, layswtch, laylow, jp, jt, jt1, &
                co2mult, colch4, colco2, colh2o, colmol, coln2o, colo2, colo3, &
                fac00, fac01, fac10, fac11, &
                selffac, selffrac, indself, forfac, forfrac, indfor, &
                zbbfd, zbbfu, zbbcd, zbbcu, zuvfd, zuvcd, znifd, znicd, &
                zbbfddir, zbbcddir, zuvfddir, zuvcddir, znifddir, znicddir)

! Prepare output up and down, clear and total flux output
            do i = 1, nlayers+1
               totuclfl(i-1) = zbbcu(i)
               totdclfl(i-1) = zbbcd(i)
               totuflux(i-1) = zbbfu(i)
               totdflux(i-1) = zbbfd(i)
! Prepare direct/diffuse flux output
               dirdflux(i-1) = zbbfddir(i)
               difdflux(i-1) = totdflux(i-1) - dirdflux(i-1)
            enddo

! Prepare net clear and total flux output
            do i = 1, nlayers+1
               fnetc(i-1) = totdclfl(i-1) - totuclfl(i-1)
               fnet(i-1) = totdflux(i-1) - totuflux(i-1)
            enddo

! Output clear and total heating rates
            do i = 1, nlayers
               zdpgcp = heatfac / pdp(i)
               htrc(i-1) = (fnetc(i) - fnetc(i-1)) * zdpgcp
               htr(i-1) = (fnet(i) - fnet(i-1)) * zdpgcp
            enddo
            htr(nlayers) = 0._rb
            htrc(nlayers) = 0._rb

! Process output.

                do i = nlayers, 0, -1
                  uflxsum(i) = uflxsum(i) + totuflux(i)
                  dflxsum(i) = dflxsum(i) + totdflux(i)
                  dirdsum(i) = dirdsum(i) + dirdflux(i)
                  difdsum(i) = difdsum(i) + difdflux(i)
                  fnetsum(i) = fnetsum(i) + fnet(i)
                  htrsum(i) = htrsum(i) + htr(i)
               enddo

! Output average over samples when last sample reached. Comment this if-check 
! and swap flux output write statement below to output all McICA samples. 
               if (ims .eq. nmca) then
                  do i = nlayers, 0, -1
                     uflxsum(i) = uflxsum(i)/nmca
                     dflxsum(i) = dflxsum(i)/nmca
                     dirdsum(i) = dirdsum(i)/nmca
                     difdsum(i) = difdsum(i)/nmca
                     fnetsum(i) = fnetsum(i)/nmca
                     htrsum(i) = htrsum(i)/nmca
                  enddo
               endif
! End statistical loop for McICA
         enddo
! End longitude/column loop
      enddo

     end subroutine rrtmg_sw
  end module rrtmg_sw_rad


