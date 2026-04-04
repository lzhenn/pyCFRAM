subroutine calc_drdt(nlayer, iaer, icld, co2, ch4, n2o, ps, ts, &
        albedo_lw, plev, t, q, o3, cldfrac, cldlwc,cldiwc, tauaer_lw, lw_base, &
        fdl_base, ful_base, drdt)

    ! Planck Jacobian
    ! Dimension: (nlayer+1)x(nlayer+1)
    ! Row: R change; Column: T change
    ! Atm: convergence of LW fluxes; Surf: LW towards the surface
    
    use parkind, only : im=>kind_im, rb=>kind_rb
    use parrrtm, only:nbndlw
    use parrrsw, only:jpband, jpb1, jpb2
    
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
    real(kind=rb), intent(in)    :: cldfrac(nlayer)! cloud fraction (0-1)
    real(kind=rb), intent(in)    :: cldlwc(nlayer) ! cloud liquid water path (kg/kg)
    real(kind=rb), intent(in)    :: cldiwc(nlayer) ! cloud ice water path (kg/kg)
    real(kind=rb), intent(in)    :: tauaer_lw(nlayer,nbndlw)    ! optical depth for longwave each band
    real(kind=rb), intent(in)    :: lw_base(nlayer)
    real(kind=rb), intent(in)    :: fdl_base(nlayer+1)
    real(kind=rb), intent(in)    :: ful_base(nlayer+1)
    real(kind=rb), intent(out)   :: drdt(nlayer+1,nlayer+1)

    integer(kind=im):: ilayer
    real(kind=rb)   :: fdl_1k(nlayer+1)
    real(kind=rb)   :: ful_1k(nlayer+1)
    real(kind=rb)   :: htr_lw_1k(nlayer)
    real(kind=rb)   :: lw_1k(nlayer)
    real(kind=rb)   :: t_1k(nlayer), ts_1k

    do ilayer = 1,nlayer
        t_1k(1:nlayer) = t(1:nlayer)
        t_1k(ilayer) = t(ilayer) + 1.0
        
        call rad_driver_lw(nlayer, iaer, icld,  co2, ch4, n2o, ps, ts, albedo_lw, &
        plev, t_1k, q, o3, cldfrac, cldlwc,cldiwc, tauaer_lw,  &
        fdl_1k, ful_1k, htr_lw_1k, lw_1k)
        
        drdt(ilayer,1:nlayer) = lw_1k(1:nlayer) - lw_base(1:nlayer)
        drdt(ilayer,nlayer+1) = fdl_1k(nlayer+1) - ful_1k(nlayer+1) - (fdl_base(nlayer+1) - ful_base(nlayer+1))

    end do

    ts_1k = ts + 1.0

    call rad_driver_lw(nlayer, iaer, icld,  co2, ch4, n2o, ps, ts_1k,  albedo_lw, &
        plev, t, q, o3, cldfrac, cldlwc,cldiwc, tauaer_lw,  &
        fdl_1k, ful_1k, htr_lw_1k, lw_1k)

        drdt(nlayer+1,1:nlayer) = lw_1k(1:nlayer) - lw_base(1:nlayer)
        drdt(nlayer+1, nlayer+1) = fdl_1k(nlayer+1) - ful_1k(nlayer+1) - (fdl_base(nlayer+1) - ful_base(nlayer+1))

end subroutine
     


