subroutine rad_driver_lw(nlayer,iaer, icld,  co2, ch4, n2o, ps, ts, &
    albedo_lw, plev, t, q, o3, cldfrac_cfram, cldlwc, cldiwc, tauaer_lw_cfram,&
    fdl, ful, htr_lw, lw_base)

    use parkind,      only: im=>kind_im, rb=>kind_rb
    use parrrtm,      only: mxlay, nbndlw, mxmol
    use parrrsw,      only: jpband,jpb1,jpb2,nmol
    use rrtmg_lw_rad, only: rrtmg_lw

    implicit none

    integer(kind=im), intent(in) :: nlayer   ! number of layers
    integer(kind=im), intent(in) :: iaer     ! 0-no aerosol, 10-aerosol
    integer(kind=im), intent(in) :: icld     ! 0-no cloud, 2-cloud
    real(kind=rb), intent(in)    :: co2          ! co2 ppmv
    real(kind=rb), intent(in)    :: ch4          ! ch4 ppmv
    real(kind=rb), intent(in)    :: n2o          ! n2o ppmv
    real(kind=rb), intent(in)    :: ps           ! surface pressure (hPa)
    real(kind=rb), intent(in)    :: ts           ! surface temperature
    real(kind=rb), intent(in)    :: albedo_lw    ! surface albedo for lw
    real(kind=rb), intent(in)    :: plev(nlayer) ! pressures of ERA
    real(kind=rb), intent(in)    :: t(nlayer)    ! atmospheric temperature
    real(kind=rb), intent(in)    :: q(nlayer)    ! specific humidity (kg/kg)
    real(kind=rb), intent(in)    :: o3(nlayer)   ! ozone ppmv
    real(kind=rb), intent(in)    :: cldfrac_cfram(nlayer)! cloud fraction (0-1)
    real(kind=rb), intent(in)    :: cldlwc(nlayer)       ! cloud liquid water path (kg/kg)
    real(kind=rb), intent(in)    :: cldiwc(nlayer)       ! cloud ice water path (kg/kg)
    real(kind=rb), intent(in)    :: tauaer_lw_cfram(nlayer,nbndlw)    ! optical depth for longwave each band
    real(kind=rb), intent(out)   :: fdl(nlayer+1)
    real(kind=rb), intent(out)   :: ful(nlayer+1)
    real(kind=rb), intent(out)   :: htr_lw(nlayer)
    real(kind=rb), intent(out)   :: lw_base(nlayer)

    integer(kind=im) :: ilayer, imol, ib
    real(kind=rb) :: amm                 ! moist air molecular weight
    real(kind=rb) :: summol              ! summation of absoption gases
    real(kind=rb) :: coef, tbound

    real(kind=rb) :: coldry_cfram(mxlay) ! column density of dry air
    real(kind=rb) :: pint_cfram(0:mxlay)
    real(kind=rb) :: tint_cfram(0:mxlay)
   ! column density for broadening gases(molecules/cm**2)
    real(kind=rb) :: wkl_cfram(mxmol,mxlay)! column volumn mixing ratio for 7 gases
    real(kind=rb) :: wbrodl_cfram(mxlay)
    real(kind=rb) :: ciwp_cfram(mxlay)
    real(kind=rb) :: clwp_cfram(mxlay)
    real(kind=rb) :: dp_cfram(mxlay)

    real(kind=rb) :: pave(mxlay)  ! mid-layer pressure
    real(kind=rb) :: tave(mxlay)  ! mid-layer temperature
    real(kind=rb) :: pz(0:mxlay)  ! interface pressure
    real(kind=rb) :: tz(0:mxlay)  ! interface temperature
    real(kind=rb) :: pdp(mxlay)   ! layer thickness
    real(kind=rb) :: semiss_lw(nbndlw) ! lw surface emissivity for each band
    real(kind=rb) :: wkl(mxmol,mxlay)! column volumn mixing ratio for 7 gases
    real(kind=rb) :: coldry(mxlay)  ! column density of dry air
    real(kind=rb) :: wbrodl(mxlay)
   ! column density for broadening gases(molecules/cm**2)
    real(kind=rb) :: cldfrac(mxlay)! cloud fraction (0.XX)
    real(kind=rb) :: ciwp(mxlay) ! cloud ice water path for the layer
    ! fraction of layer's cloud water path in the form of ice particles
    real(kind=rb) :: clwp(mxlay) ! cloud liquid water path for the layer (g/m2)
    real(kind=rb) :: tauaer_lw(mxlay,nbndlw)
    ! aerosol optical depth for each layer and each waveband
    ! refer rrtmg documentation for the waveband details
    real(kind=rb) :: rei(mxlay) ! cloud ice particle effective size (microns)
    real(kind=rb) :: rel(mxlay) ! cloud liquid partical effective size (microns)

    real(kind=rb) :: uflxsum_lw(0:mxlay) ! lw upward flux (W/m**2)
    real(kind=rb) :: dflxsum_lw(0:mxlay) ! lw downward flux
    real(kind=rb) :: fnetsum_lw(0:mxlay) ! lw net flux
    real(kind=rb) :: htrsum_lw(0:mxlay)  ! lw heating rate (K/day)

    real(kind=rb) :: amw                 ! Molecular weight of water vapor (g/mol)
    real(kind=rb) :: amd                 ! Effective molecular weight of dry air (g/mol)
    real(kind=rb) :: grav                ! acceleration of gravity
    real(kind=rb) :: pwvcm               ! precipitable water vapor (cm)
    real(kind=rb) :: wvsh                ! water vapor vertical total specific humitidy
    real(kind=rb) :: amttl               ! moist air vertical sum (molecular amount)
    real(kind=rb) :: wvttl               ! water vapor vertical sum (molecular amount)


!------------------------------------------------------------------------------------
    ! calculate pint and tint
    ! interface top level is 0.5 hPa
    pint_cfram(0) = 0.5
    tint_cfram(0) = t(1)
    ! interface bottom level is the surface
    pint_cfram(nlayer) = ps
    tint_cfram(nlayer) = ts

    ! interface is the average of mid-layer
    do ilayer = 1, nlayer-1
        pint_cfram(ilayer) = (plev(ilayer) + plev(ilayer+1))/2.0
        tint_cfram(ilayer) = (t(ilayer) + t(ilayer+1))/2.0
    end do
    
    do ilayer =1 ,nlayer
        wkl_cfram(1,ilayer) = q(ilayer)* 28.966/18.016
        wkl_cfram(2,ilayer) = co2*1.e-6
        wkl_cfram(3,ilayer) = o3(ilayer) * 28.966/48.0
        wkl_cfram(4,ilayer) = n2o*1.e-6
        wkl_cfram(6,ilayer) = ch4*1.e-6
        wkl_cfram(7,ilayer) = 0.209
        ! the molecular weight of moist air
        amm = (1. - wkl_cfram(1,ilayer)) * 28.966 + wkl_cfram(1,ilayer) * 18.016
        ! unit: molecules/cm**2
        coldry_cfram(ilayer) = (pint_cfram(ilayer) - pint_cfram(ilayer-1)) * 1.e3 * 6.02e23 / &
            (1.e2 * 9.8 * amm * (1.+ wkl_cfram(1,ilayer)))
        summol = wkl_cfram(2,ilayer)+wkl_cfram(3,ilayer)+wkl_cfram(4,ilayer)+wkl_cfram(6,ilayer)+wkl_cfram(7,ilayer)
        ! column density for broadening gases
        wbrodl_cfram(ilayer) = coldry_cfram(ilayer) * (1-summol)
        dp_cfram(ilayer) = pint_cfram(ilayer) - pint_cfram(ilayer-1)
   end do

    ! upward: from top to bottom
    do ilayer = 1, nlayer
        if (cldfrac_cfram(ilayer) .gt. 1.e-5) then
            coef = (1.0/9.8) * 1.e2 * 1.e3
            !fracice_cfram(ilayer) = cldiwc(ilayer)/(cldiwc(ilayer) + cldlwc(ilayer))
            !cwp_cfram(ilayer) = coef*(cldiwc(ilayer)+cldlwc(ilayer))*dp_cfram(ilayer)
            ciwp_cfram(ilayer) = coef * cldiwc(ilayer) * dp_cfram(ilayer)
            clwp_cfram(ilayer) = coef * cldlwc(ilayer) * dp_cfram(ilayer)
        end if
    end do

    ! convert top to bottom To bottom to top 
    do ilayer = 0,nlayer
        pz(ilayer) = pint_cfram(nlayer - ilayer)
        tz(ilayer) = tint_cfram(nlayer - ilayer)
    end do

    do ilayer = 1, nlayer
        pave(nlayer - ilayer +1)      = plev(ilayer)
        pdp(nlayer - ilayer +1)       = dp_cfram(ilayer)
        tave(nlayer - ilayer +1)      = t(ilayer)
        wkl(:,nlayer -ilayer+1)       = wkl_cfram(:,ilayer)
        coldry(nlayer - ilayer+1)     = coldry_cfram(ilayer)
        wbrodl(nlayer - ilayer+1)     = wbrodl_cfram(ilayer)
    end do
   
    ! Convert wkl from mixing ratios to column densities
    do ilayer = 1,nlayer
        do imol = 1, nmol
            wkl(imol,ilayer) = coldry(ilayer)*wkl(imol,ilayer)
        end do
    end do

    amw = 18.016
    amd = 28.966
    grav= 9.8066
    amttl = 0.0
    wvttl = 0.0
    do ilayer = 1,nlayer
       amttl = amttl + coldry(ilayer)+wkl(1,ilayer)
       wvttl = wvttl + wkl(1,ilayer)
   end do
    !  Calculate total precipitable water 
    wvsh = (amw * wvttl) / (amd * amttl)
    pwvcm = wvsh * (1.e3_rb * pz(0)) / (1.e2_rb * grav)

    ! cloud partical size: ice 5 microns, liquid 20 microns
    tbound       = ts   ! tbound < 0: tbound = tz(0)
    semiss_lw    = 1 - albedo_lw
   
    do ib = 1,nbndlw
        semiss_lw(ib) = 1 - albedo_lw
    end do

   if (icld .ne. 0) then
       do ilayer = 1, nlayer
           cldfrac(nlayer - ilayer +1) = cldfrac_cfram(ilayer)
           ciwp(nlayer - ilayer +1)    = ciwp_cfram(ilayer)
           clwp(nlayer - ilayer +1)    = clwp_cfram(ilayer)
           rei(ilayer) = 5
           rel(ilayer) = 20
       end do
   end if

   if (iaer .ne. 0) then
       do ilayer = 1,nlayer
           tauaer_lw(nlayer-ilayer+1,:) = tauaer_lw_cfram(ilayer,:)
       end do
   end if

    call rrtmg_lw(iaer, icld, nlayer, pave, tave, pz, tz, tbound, coldry,&
        wbrodl, wkl, pwvcm, semiss_lw, cldfrac, ciwp, clwp, rel, rei, tauaer_lw,&
        uflxsum_lw, dflxsum_lw, fnetsum_lw, htrsum_lw)

    do ilayer = 1,nlayer+1
        fdl(ilayer) = dflxsum_lw(nlayer+1-ilayer)
        ful(ilayer) = uflxsum_lw(nlayer+1-ilayer)
    end do

    do ilayer = 1,nlayer
        htr_lw(ilayer) = htrsum_lw(nlayer-ilayer)
    end do

    do ilayer = 1, nlayer
       lw_base(ilayer) = fdl(ilayer) + ful(ilayer+1) - fdl(ilayer+1) - ful(ilayer)
    end do

    return
    end subroutine rad_driver_lw
    
