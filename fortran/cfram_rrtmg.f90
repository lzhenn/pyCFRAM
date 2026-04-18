 Program CFRAM

    use parkind,   only:im=>kind_im, rb=>kind_rb
    use parrrtm,   only:nbndlw
    use parrrsw,   only:jpband, jpb1, jpb2
    use output,    only:write_out_3d
    use math,      only:inv

    implicit none

    integer(kind=im),parameter :: nlat=1,nlon=1,nlev=37,ntime=1
    integer(kind=im),parameter :: icld=2   ! 0-no cloud, 2-cloud
    integer(kind=im),parameter :: iaer=10   ! 0-no aerosol, 10-aerosol
    integer(kind=im),parameter :: nspecies=6 ! per-species: bc,ocphi,ocpho,sulf,ss,dust
    real(kind=rb),   parameter :: scon=1360.98 ! solar constant
    integer(kind=im)           :: ilat, ilon, ilev, itime, ib, i, j
    integer(kind=im)           :: ilayer,nlayer,isp

    real(kind=rb) :: co2ppmv_base,co2ppmv_warm,ch4ppmv,n2oppmv
    real(kind=rb) :: plev(nlev)
    real(kind=rb) :: ts_base(nlat,nlon),ts_warm(nlat,nlon)
    real(kind=rb) :: ps_base(nlat,nlon),ps_warm(nlat,nlon)
    real(kind=rb) :: albedo_lw_base(nlat,nlon),albedo_lw_warm(nlat,nlon)
    real(kind=rb) :: albedo_sw_base(nlat,nlon),albedo_sw_warm(nlat,nlon)
    real(kind=rb) :: solin_base(nlat,nlon), solin_warm(nlat,nlon)
    real(kind=rb) :: zenith_base(nlat,nlon), zenith_warm(nlat,nlon)
    real(kind=rb) :: swdn_surf_base(nlat,nlon),swdn_surf_warm(nlat,nlon)
    real(kind=rb) :: swup_surf_base(nlat,nlon),swup_surf_warm(nlat,nlon)
    real(kind=rb) :: t_base(nlev,nlat,nlon), t_warm(nlev,nlat,nlon)
    real(kind=rb) :: q_base(nlev,nlat,nlon), q_warm(nlev,nlat,nlon)
    real(kind=rb) :: o3_base(nlev,nlat,nlon),o3_warm(nlev,nlat,nlon)
    real(kind=rb) :: cldfrac_base(nlev,nlat,nlon), cldfrac_warm(nlev,nlat,nlon)
    real(kind=rb) :: cldlwc_base(nlev,nlat,nlon), cldlwc_warm(nlev,nlat,nlon)
    real(kind=rb) :: cldiwc_base(nlev,nlat,nlon), cldiwc_warm(nlev,nlat,nlon)
    real(kind=rb) :: tauaer_sw_base(nlev,nlat,nlon,jpband),tauaer_sw_warm(nlev,nlat,nlon,jpband)
    real(kind=rb) :: ssaaer_sw_base(nlev,nlat,nlon,jpband),ssaaer_sw_warm(nlev,nlat,nlon,jpband)
    real(kind=rb) :: asmaer_sw_base(nlev,nlat,nlon,jpband),asmaer_sw_warm(nlev,nlat,nlon,jpband)
    real(kind=rb) :: tauaer_lw_base(nlev,nlat,nlon,nbndlw),tauaer_lw_warm(nlev,nlat,nlon,jpband)

    ! Per-species aerosol optical properties. Species order must match Python writer
    ! (SPECIES_ORDER in run_parallel_python.py / extract_full_field.py):
    ! 1=bc, 2=ocphi, 3=ocpho, 4=sulf, 5=ss, 6=dust.
    ! File layout on disk (from numpy tofile): C-order bytes with species fastest.
    ! Fortran array declared with species as LAST dimension so READ loop reads
    ! species innermost (matches disk layout).
    real(kind=rb) :: tauaer_sw_base_spc(nlev,nlat,nlon,jpband,nspecies)
    real(kind=rb) :: tauaer_sw_warm_spc(nlev,nlat,nlon,jpband,nspecies)
    real(kind=rb) :: ssaaer_sw_base_spc(nlev,nlat,nlon,jpband,nspecies)
    real(kind=rb) :: ssaaer_sw_warm_spc(nlev,nlat,nlon,jpband,nspecies)
    real(kind=rb) :: asmaer_sw_base_spc(nlev,nlat,nlon,jpband,nspecies)
    real(kind=rb) :: asmaer_sw_warm_spc(nlev,nlat,nlon,jpband,nspecies)
    real(kind=rb) :: tauaer_lw_base_spc(nlev,nlat,nlon,nbndlw,nspecies)
    real(kind=rb) :: tauaer_lw_warm_spc(nlev,nlat,nlon,nbndlw,nspecies)

    !Output: forcing only (dT solved in Python)
    real(kind=rb) :: frc_warm_output(ntime,nlev+1,nlat,nlon), frc_co2_output(ntime,nlev+1,nlat,nlon)
    real(kind=rb) :: frc_q_output(ntime,nlev+1,nlat,nlon), frc_ts_output(ntime,nlev+1,nlat,nlon)
    real(kind=rb) :: frc_o3_output(ntime,nlev+1,nlat,nlon), frc_solar_output(ntime,nlev+1,nlat,nlon)
    real(kind=rb) :: frc_albedo_output(ntime,nlev+1,nlat,nlon), frc_cloud_output(ntime,nlev+1,nlat,nlon)
    real(kind=rb) :: frc_aerosol_output(ntime,nlev+1,nlat,nlon)
    ! Per-species aerosol forcing outputs (Phase 3).
    real(kind=rb) :: frc_bc_output(ntime,nlev+1,nlat,nlon)
    real(kind=rb) :: frc_ocphi_output(ntime,nlev+1,nlat,nlon)
    real(kind=rb) :: frc_ocpho_output(ntime,nlev+1,nlat,nlon)
    real(kind=rb) :: frc_sulf_output(ntime,nlev+1,nlat,nlon)
    real(kind=rb) :: frc_ss_output(ntime,nlev+1,nlat,nlon)
    real(kind=rb) :: frc_dust_output(ntime,nlev+1,nlat,nlon)
    ! Cloud LW/SW split forcing outputs
    real(kind=rb) :: frc_cloud_lw_output(ntime,nlev+1,nlat,nlon)
    real(kind=rb) :: frc_cloud_sw_output(ntime,nlev+1,nlat,nlon)

    ! Thread-private fixed-size work arrays (nlev+1 = 38 max)
    real(kind=rb) :: fds_1d(nlev+1), fus_1d(nlev+1), htr_sw_1d(nlev)
    real(kind=rb) :: fdl_1d(nlev+1), ful_1d(nlev+1), htr_lw_1d(nlev)
    real(kind=rb) :: htr_1d(nlev), lw_1d(nlev), sw_1d(nlev)
    real(kind=rb) :: rad_1d_base(nlev), fd_1d_base(nlev+1), fu_1d_base(nlev+1)
    real(kind=rb) :: rad_1d_warm(nlev), fd_1d_warm(nlev+1), fu_1d_warm(nlev+1)
    real(kind=rb) :: rad_1d_t(nlev),    fd_1d_t(nlev+1),    fu_1d_t(nlev+1)
    real(kind=rb) :: rad_1d_q(nlev),    fd_1d_q(nlev+1),    fu_1d_q(nlev+1)
    real(kind=rb) :: rad_1d_ts(nlev),   fd_1d_ts(nlev+1),   fu_1d_ts(nlev+1)
    real(kind=rb) :: rad_1d_o3(nlev),   fd_1d_o3(nlev+1),   fu_1d_o3(nlev+1)
    real(kind=rb) :: rad_1d_solar(nlev),fd_1d_solar(nlev+1),fu_1d_solar(nlev+1)
    real(kind=rb) :: rad_1d_albedo(nlev),fd_1d_albedo(nlev+1),fu_1d_albedo(nlev+1)
    real(kind=rb) :: rad_1d_cloud(nlev), fd_1d_cloud(nlev+1), fu_1d_cloud(nlev+1)
    real(kind=rb) :: rad_1d_aerosol(nlev),fd_1d_aerosol(nlev+1),fu_1d_aerosol(nlev+1)
    real(kind=rb) :: rad_1d_co2(nlev),   fd_1d_co2(nlev+1),   fu_1d_co2(nlev+1)
    real(kind=rb) :: lw_1d_base(nlev),   fdl_1d_base(nlev+1), ful_1d_base(nlev+1)

    real(kind=rb) :: frc_co2(nlev+1), frc_t(nlev+1), frc_q(nlev+1), frc_albedo(nlev+1)
    real(kind=rb) :: frc_ts(nlev+1), frc_o3(nlev+1), frc_solar(nlev+1), frc_cloud(nlev+1)
    real(kind=rb) :: frc_warm(nlev+1), frc_aerosol(nlev+1)
    real(kind=rb) :: frc_cloud_lw(nlev+1), frc_cloud_sw(nlev+1)
    ! Cloud-state LW/SW component fluxes, saved right after the cloud rad_driver call
    real(kind=rb) :: lw_1d_cloud(nlev), sw_1d_cloud(nlev)
    real(kind=rb) :: fdl_1d_cloud(nlev+1), ful_1d_cloud(nlev+1)
    real(kind=rb) :: fds_1d_cloud(nlev+1), fus_1d_cloud(nlev+1)

    ! Per-species perturbation work arrays (thread-private)
    real(kind=rb) :: tauaer_sw_mix(nlev,jpband), ssaaer_sw_mix(nlev,jpband), asmaer_sw_mix(nlev,jpband)
    real(kind=rb) :: tauaer_lw_mix(nlev,nbndlw)
    real(kind=rb) :: tau_ssa_base_sw(nlev,jpband), tau_ssa_g_base_sw(nlev,jpband)
    real(kind=rb) :: tau_ssa_mix(nlev,jpband), tau_ssa_g_mix(nlev,jpband)
    real(kind=rb) :: rad_1d_spc(nlev), fd_1d_spc(nlev+1), fu_1d_spc(nlev+1)
    real(kind=rb) :: frc_spc(nlev+1,nspecies)

    real(kind=rb) :: drdt(nlev+1,nlev+1), drdt_inv(nlev+1,nlev+1)

    plev = (/1.,2.,3.,5.,7.,10.,20.,30.,50.,70.,100.,125.,150.,175.,200.,&
        225.,250.,300.,350.,400.,450.,500.,550.,600.,650.,700.,750.,775.,&
        800.,825.,850.,875.,900.,925.,950.,975.,1000./)

!-------------------------------- Open Input Files-----------------------------------------------
    open(unit=101,file='data_prep/hus_base.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8)
    open(unit=102,file='data_prep/t_base.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8)
    open(unit=103,file='data_prep/solarin_base.dat',form='unformatted',access='direct',recl=nlat*nlon*8)
    open(unit=104,file='data_prep/ssrd_base.dat',form='unformatted', access='direct',recl=nlat*nlon*8)
    open(unit=105,file='data_prep/ssru_base.dat',form='unformatted', access='direct',recl=nlat*nlon*8)
    open(unit=106,file='data_prep/skt_base.dat',form='unformatted', access='direct',recl=nlat*nlon*8)
    open(unit=107,file='data_prep/sp_base.dat',form='unformatted', access='direct',recl=nlat*nlon*8)
    open(unit=108,file='data_prep/O3_base.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8)
    open(unit=201,file='data_prep/hus_warm.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8)
    open(unit=202,file='data_prep/t_warm.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8)
    open(unit=203,file='data_prep/solarin_warm.dat',form='unformatted',access='direct',recl=nlat*nlon*8)
    open(unit=204,file='data_prep/ssrd_warm.dat',form='unformatted', access='direct',recl=nlat*nlon*8)
    open(unit=205,file='data_prep/ssru_warm.dat',form='unformatted', access='direct',recl=nlat*nlon*8)
    open(unit=206,file='data_prep/skt_warm.dat',form='unformatted', access='direct',recl=nlat*nlon*8)
    open(unit=207,file='data_prep/sp_warm.dat',form='unformatted', access='direct',recl=nlat*nlon*8)
    open(unit=208,file='data_prep/O3_warm.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8)

    if (icld == 2) then
        print *, "Reading Cloud Profiles"
        open(unit=111,file='data_prep/cc_base.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8)
        open(unit=112,file='data_prep/clwc_base.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8)
        open(unit=113,file='data_prep/ciwc_base.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8)
        open(unit=211,file='data_prep/cc_warm.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8)
        open(unit=212,file='data_prep/clwc_warm.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8)
        open(unit=213,file='data_prep/ciwc_warm.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8)
    else
        print *, "Clear sky, no cloud"
    end if

    if (iaer == 10) then
        print *, "Reading Aerosol Profiles"
        open(unit=121,file='data_prep/aerosol_aod_lw_base.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8*nbndlw)
        open(unit=122,file='data_prep/aerosol_aod_sw_base.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8*jpband)
        open(unit=123,file='data_prep/aerosol_ssa_sw_base.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8*jpband)
        open(unit=124,file='data_prep/aerosol_g_sw_base.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8*jpband)
        open(unit=221,file='data_prep/aerosol_aod_lw_warm.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8*nbndlw)
        open(unit=222,file='data_prep/aerosol_aod_sw_warm.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8*jpband)
        open(unit=223,file='data_prep/aerosol_ssa_sw_warm.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8*jpband)
        open(unit=224,file='data_prep/aerosol_g_sw_warm.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8*jpband)
        ! Per-species optical properties (Phase 2: load only, not used yet)
        open(unit=125,file='data_prep/aerosol_aod_lw_base_spc.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8*nbndlw*nspecies)
        open(unit=126,file='data_prep/aerosol_aod_sw_base_spc.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8*jpband*nspecies)
        open(unit=127,file='data_prep/aerosol_ssa_sw_base_spc.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8*jpband*nspecies)
        open(unit=128,file='data_prep/aerosol_g_sw_base_spc.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8*jpband*nspecies)
        open(unit=225,file='data_prep/aerosol_aod_lw_warm_spc.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8*nbndlw*nspecies)
        open(unit=226,file='data_prep/aerosol_aod_sw_warm_spc.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8*jpband*nspecies)
        open(unit=227,file='data_prep/aerosol_ssa_sw_warm_spc.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8*jpband*nspecies)
        open(unit=228,file='data_prep/aerosol_g_sw_warm_spc.dat',form='unformatted', access='direct',recl=nlev*nlat*nlon*8*jpband*nspecies)
    else
        print *, "No aerosol"
    end if

    open(unit=131,file='data_prep/co2_b.dat',form='unformatted', access='direct',recl = 8 )
    open(unit=231,file='data_prep/co2_w.dat',form='unformatted', access='direct',recl = 8 )

!-------------------------------Read Input Data--------------------------------------------------
    do itime = 1,ntime
        read(131,rec=itime) co2ppmv_base
        read(101,rec=itime) (((q_base(ilev,ilat,ilon),ilon=1,nlon),ilat=1,nlat),ilev=1,nlev)
        read(102,rec=itime) (((t_base(ilev,ilat,ilon),ilon=1,nlon),ilat=1,nlat),ilev=1,nlev)
        read(103,rec=itime) ((solin_base(ilat,ilon),ilon=1,nlon),ilat=1,nlat)
        read(104,rec=itime) ((swdn_surf_base(ilat,ilon),ilon=1,nlon),ilat=1,nlat)
        read(105,rec=itime) ((swup_surf_base(ilat,ilon),ilon=1,nlon),ilat=1,nlat)
        read(106,rec=itime) ((ts_base(ilat,ilon),ilon=1,nlon),ilat=1,nlat)
        read(107,rec=itime) ((ps_base(ilat,ilon),ilon=1,nlon),ilat=1,nlat)
        read(108,rec=itime) (((o3_base(ilev,ilat,ilon),ilon=1,nlon),ilat=1,nlat),ilev=1,nlev)

        read(231,rec=itime) co2ppmv_warm
        read(201,rec=itime) (((q_warm(ilev,ilat,ilon),ilon=1,nlon),ilat=1,nlat),ilev=1,nlev)
        read(202,rec=itime) (((t_warm(ilev,ilat,ilon),ilon=1,nlon),ilat=1,nlat),ilev=1,nlev)
        read(203,rec=itime) ((solin_warm(ilat,ilon),ilon=1,nlon),ilat=1,nlat)
        read(204,rec=itime) ((swdn_surf_warm(ilat,ilon),ilon=1,nlon),ilat=1,nlat)
        read(205,rec=itime) ((swup_surf_warm(ilat,ilon),ilon=1,nlon),ilat=1,nlat)
        read(206,rec=itime) ((ts_warm(ilat,ilon),ilon=1,nlon),ilat=1,nlat)
        read(207,rec=itime) ((ps_warm(ilat,ilon),ilon=1,nlon),ilat=1,nlat)
        read(208,rec=itime) (((o3_warm(ilev,ilat,ilon),ilon=1,nlon),ilat=1,nlat),ilev=1,nlev)

        if (icld == 2) then
            read(111,rec=itime) (((cldfrac_base(ilev,ilat,ilon),ilon=1,nlon),ilat=1,nlat),ilev=1,nlev)
            read(112,rec=itime) (((cldlwc_base(ilev,ilat,ilon),ilon=1,nlon),ilat=1,nlat),ilev=1,nlev)
            read(113,rec=itime)(((cldiwc_base(ilev,ilat,ilon),ilon=1,nlon),ilat=1,nlat),ilev=1,nlev)
            read(211,rec=itime) (((cldfrac_warm(ilev,ilat,ilon),ilon=1,nlon),ilat=1,nlat),ilev=1,nlev)
            read(212,rec=itime) (((cldlwc_warm(ilev,ilat,ilon),ilon=1,nlon),ilat=1,nlat),ilev=1,nlev)
            read(213,rec=itime)(((cldiwc_warm(ilev,ilat,ilon),ilon=1,nlon),ilat=1,nlat),ilev=1,nlev)
        else
            cldfrac_base = 0.0; cldlwc_base = 0.0; cldiwc_base = 0.0
            cldfrac_warm = 0.0; cldlwc_warm = 0.0; cldiwc_warm = 0.0
        end if

        if (iaer == 10) then
            read(121,rec=itime) ((((tauaer_lw_base(ilev,ilat,ilon,ib),&
                ib=1,nbndlw), ilon=1,nlon), ilat=1,nlat),ilev=1,nlev)
            read(122,rec=itime) ((((tauaer_sw_base(ilev,ilat,ilon,ib),&
                ib=jpb1,jpb2),ilon=1,nlon), ilat=1,nlat),ilev=1,nlev)
            read(123,rec=itime) ((((ssaaer_sw_base(ilev,ilat,ilon,ib),&
                ib=jpb1,jpb2),ilon=1,nlon), ilat=1,nlat),ilev=1,nlev)
            read(124,rec=itime) ((((asmaer_sw_base(ilev,ilat,ilon,ib),&
                ib=jpb1,jpb2),ilon=1,nlon), ilat=1,nlat),ilev=1,nlev)
            read(221,rec=itime) ((((tauaer_lw_warm(ilev,ilat,ilon,ib),&
                ib=1,nbndlw), ilon=1,nlon), ilat=1,nlat),ilev=1,nlev)
            read(222,rec=itime) ((((tauaer_sw_warm(ilev,ilat,ilon,ib),&
                ib=jpb1,jpb2),ilon=1,nlon), ilat=1,nlat),ilev=1,nlev)
            read(223,rec=itime) ((((ssaaer_sw_warm(ilev,ilat,ilon,ib),&
                ib=jpb1,jpb2),ilon=1,nlon), ilat=1,nlat),ilev=1,nlev)
            read(224,rec=itime) ((((asmaer_sw_warm(ilev,ilat,ilon,ib),&
                ib=jpb1,jpb2),ilon=1,nlon), ilat=1,nlat),ilev=1,nlev)
            ! Per-species: on disk bytes are C-order with species fastest, so READ
            ! loop has isp innermost, then ib, ilon, ilat, ilev outermost.
            read(125,rec=itime) (((((tauaer_lw_base_spc(ilev,ilat,ilon,ib,isp),&
                isp=1,nspecies),ib=1,nbndlw),ilon=1,nlon),ilat=1,nlat),ilev=1,nlev)
            read(126,rec=itime) (((((tauaer_sw_base_spc(ilev,ilat,ilon,ib,isp),&
                isp=1,nspecies),ib=jpb1,jpb2),ilon=1,nlon),ilat=1,nlat),ilev=1,nlev)
            read(127,rec=itime) (((((ssaaer_sw_base_spc(ilev,ilat,ilon,ib,isp),&
                isp=1,nspecies),ib=jpb1,jpb2),ilon=1,nlon),ilat=1,nlat),ilev=1,nlev)
            read(128,rec=itime) (((((asmaer_sw_base_spc(ilev,ilat,ilon,ib,isp),&
                isp=1,nspecies),ib=jpb1,jpb2),ilon=1,nlon),ilat=1,nlat),ilev=1,nlev)
            read(225,rec=itime) (((((tauaer_lw_warm_spc(ilev,ilat,ilon,ib,isp),&
                isp=1,nspecies),ib=1,nbndlw),ilon=1,nlon),ilat=1,nlat),ilev=1,nlev)
            read(226,rec=itime) (((((tauaer_sw_warm_spc(ilev,ilat,ilon,ib,isp),&
                isp=1,nspecies),ib=jpb1,jpb2),ilon=1,nlon),ilat=1,nlat),ilev=1,nlev)
            read(227,rec=itime) (((((ssaaer_sw_warm_spc(ilev,ilat,ilon,ib,isp),&
                isp=1,nspecies),ib=jpb1,jpb2),ilon=1,nlon),ilat=1,nlat),ilev=1,nlev)
            read(228,rec=itime) (((((asmaer_sw_warm_spc(ilev,ilat,ilon,ib,isp),&
                isp=1,nspecies),ib=jpb1,jpb2),ilon=1,nlon),ilat=1,nlat),ilev=1,nlev)
        else
            tauaer_lw_base = 0.0; tauaer_sw_base = 0.0; ssaaer_sw_base = 0.0; asmaer_sw_base = 0.0
            tauaer_lw_warm = 0.0; tauaer_sw_warm = 0.0; ssaaer_sw_warm = 0.0; asmaer_sw_warm = 0.0
            tauaer_lw_base_spc = 0.0; tauaer_sw_base_spc = 0.0
            ssaaer_sw_base_spc = 0.0; asmaer_sw_base_spc = 0.0
            tauaer_lw_warm_spc = 0.0; tauaer_sw_warm_spc = 0.0
            ssaaer_sw_warm_spc = 0.0; asmaer_sw_warm_spc = 0.0
        end if

        ch4ppmv = 1.6
        n2oppmv = 0.28

        do ilat=1,nlat
           do ilon = 1,nlon
             if (swdn_surf_base(ilat,ilon) .eq. 0.) then
                albedo_sw_base(ilat,ilon) = 0.0
             else
                albedo_sw_base(ilat,ilon) = swup_surf_base(ilat,ilon)/swdn_surf_base(ilat,ilon)
             end if
             if (swdn_surf_warm(ilat,ilon) .eq. 0.) then
                albedo_sw_warm(ilat,ilon) = 0.0
             else
                albedo_sw_warm(ilat,ilon) = swup_surf_warm(ilat,ilon)/swdn_surf_warm(ilat,ilon)
             end if
         end do
        end do

        albedo_lw_base = 0.0
        albedo_lw_warm = 0.0
        zenith_base = solin_base/scon
        zenith_warm = solin_warm/scon

        print *, "Data loaded, starting computation..."

!--------------------RT Calculation-----------------------------------------
! Fortran computes: 9 radiative forcings + Planck inverse matrix
! Python handles: dT = -(dR/dT)^-1 * forcing for all terms
!$OMP PARALLEL DO PRIVATE(ilon, ilayer, nlayer, &
!$OMP   fds_1d, fus_1d, htr_sw_1d, fdl_1d, ful_1d, htr_lw_1d, &
!$OMP   htr_1d, lw_1d, sw_1d, &
!$OMP   fd_1d_base, fu_1d_base, rad_1d_base, &
!$OMP   fd_1d_warm, fu_1d_warm, rad_1d_warm, &
!$OMP   fd_1d_t, fu_1d_t, rad_1d_t, &
!$OMP   fd_1d_q, fu_1d_q, rad_1d_q, &
!$OMP   fd_1d_ts, fu_1d_ts, rad_1d_ts, &
!$OMP   fd_1d_o3, fu_1d_o3, rad_1d_o3, &
!$OMP   fd_1d_solar, fu_1d_solar, rad_1d_solar, &
!$OMP   fd_1d_albedo, fu_1d_albedo, rad_1d_albedo, &
!$OMP   fd_1d_cloud, fu_1d_cloud, rad_1d_cloud, &
!$OMP   fd_1d_aerosol, fu_1d_aerosol, rad_1d_aerosol, &
!$OMP   fd_1d_co2, fu_1d_co2, rad_1d_co2, &
!$OMP   lw_1d_base, fdl_1d_base, ful_1d_base, &
!$OMP   frc_co2, frc_t, frc_q, frc_ts, frc_o3, frc_solar, &
!$OMP   frc_cloud, frc_albedo, frc_warm, frc_aerosol, &
!$OMP   drdt, drdt_inv) &
!$OMP SCHEDULE(dynamic)
       do ilat = 1,nlat
           do ilon = 1,nlon

               ! find the index of the lowest layer
               do ilayer = 1, nlev
                   if(plev(ilayer) .lt. ps_base(ilat,ilon)/100.0) then
                         nlayer = ilayer
                   else
                        exit
                   end if
               end do

               ! Zero work arrays
               fds_1d = 0; fus_1d = 0; htr_sw_1d = 0
               fdl_1d = 0; ful_1d = 0; htr_lw_1d = 0
               htr_1d = 0; lw_1d = 0; sw_1d = 0
               drdt = 0; drdt_inv = 0

               ! Baseline
               call rad_driver(nlayer, iaer, icld, co2ppmv_base, ch4ppmv, n2oppmv, ps_base(ilat,ilon)/100.0,&
                   ts_base(ilat,ilon), zenith_base(ilat,ilon), albedo_sw_base(ilat,ilon), albedo_lw_base(ilat,ilon),&
                   plev, t_base(:,ilat,ilon), q_base(:,ilat,ilon), o3_base(:,ilat,ilon), cldfrac_base(:,ilat,ilon), cldlwc_base(:,ilat,ilon),&
                   cldiwc_base(:,ilat,ilon), tauaer_sw_base(:,ilat,ilon,:), ssaaer_sw_base(:,ilat,ilon,:), asmaer_sw_base(:,ilat,ilon,:),&
                   tauaer_lw_base(:,ilat,ilon,:), fds_1d, fus_1d, htr_sw_1d,&
                   fdl_1d_base, ful_1d_base, htr_lw_1d, fd_1d_base, fu_1d_base, htr_1d, rad_1d_base, lw_1d_base, sw_1d)

               ! Warm
               call rad_driver(nlayer, iaer, icld, co2ppmv_warm, ch4ppmv, n2oppmv, ps_base(ilat,ilon)/100.0,&
                   ts_warm(ilat,ilon), zenith_warm(ilat,ilon), albedo_sw_warm(ilat,ilon), albedo_lw_warm(ilat,ilon),&
                   plev, t_base(:,ilat,ilon), q_warm(:,ilat,ilon), o3_warm(:,ilat,ilon), cldfrac_warm(:,ilat,ilon), cldlwc_warm(:,ilat,ilon),&
                   cldiwc_warm(:,ilat,ilon), tauaer_sw_warm(:,ilat,ilon,:), ssaaer_sw_warm(:,ilat,ilon,:), asmaer_sw_warm(:,ilat,ilon,:),&
                   tauaer_lw_warm(:,ilat,ilon,:), fds_1d, fus_1d, htr_sw_1d,&
                   fdl_1d, ful_1d, htr_lw_1d, fd_1d_warm, fu_1d_warm, htr_1d, rad_1d_warm, lw_1d, sw_1d)

               ! Planck Matrix
               call calc_drdt(nlayer, iaer, icld, co2ppmv_base, ch4ppmv, n2oppmv, ps_base(ilat,ilon)/100.0,&
                   ts_base(ilat,ilon),  albedo_lw_base(ilat,ilon),&
                   plev, t_base(:,ilat,ilon), q_base(:,ilat,ilon), o3_base(:,ilat,ilon), cldfrac_base(:,ilat,ilon), cldlwc_base(:,ilat,ilon),&
                   cldiwc_base(:,ilat,ilon), tauaer_lw_base(:,ilat,ilon,:), lw_1d_base, fdl_1d_base, ful_1d_base, &
                   drdt(1:nlayer+1,1:nlayer+1))

                  drdt_inv(1:nlayer+1,1:nlayer+1) = inv(drdt(1:nlayer+1,1:nlayer+1))

               ! CO2
               call rad_driver(nlayer, iaer, icld, co2ppmv_warm, ch4ppmv, n2oppmv, ps_base(ilat,ilon)/100.0,&
                   ts_base(ilat,ilon), zenith_base(ilat,ilon), albedo_sw_base(ilat,ilon), albedo_lw_base(ilat,ilon),&
                   plev, t_base(:,ilat,ilon), q_base(:,ilat,ilon), o3_base(:,ilat,ilon), cldfrac_base(:,ilat,ilon), cldlwc_base(:,ilat,ilon),&
                   cldiwc_base(:,ilat,ilon), tauaer_sw_base(:,ilat,ilon,:), ssaaer_sw_base(:,ilat,ilon,:), asmaer_sw_base(:,ilat,ilon,:),&
                   tauaer_lw_base(:,ilat,ilon,:), fds_1d, fus_1d, htr_sw_1d,&
                   fdl_1d, ful_1d, htr_lw_1d, fd_1d_co2, fu_1d_co2, htr_1d, rad_1d_co2, lw_1d, sw_1d)

               ! Q
             call rad_driver(nlayer, iaer, icld, co2ppmv_base, ch4ppmv, n2oppmv, ps_base(ilat,ilon)/100.0,&
                   ts_base(ilat,ilon), zenith_base(ilat,ilon), albedo_sw_base(ilat,ilon), albedo_lw_base(ilat,ilon),&
                   plev, t_base(:,ilat,ilon), q_warm(:,ilat,ilon), o3_base(:,ilat,ilon), cldfrac_base(:,ilat,ilon), cldlwc_base(:,ilat,ilon),&
                   cldiwc_base(:,ilat,ilon), tauaer_sw_base(:,ilat,ilon,:), ssaaer_sw_base(:,ilat,ilon,:), asmaer_sw_base(:,ilat,ilon,:),&
                   tauaer_lw_base(:,ilat,ilon,:), fds_1d, fus_1d, htr_sw_1d,&
                   fdl_1d, ful_1d, htr_lw_1d, fd_1d_q, fu_1d_q, htr_1d, rad_1d_q, lw_1d, sw_1d)

               ! ts
             call rad_driver(nlayer, iaer, icld, co2ppmv_base, ch4ppmv, n2oppmv, ps_base(ilat,ilon)/100.0,&
                   ts_warm(ilat,ilon), zenith_base(ilat,ilon), albedo_sw_base(ilat,ilon), albedo_lw_base(ilat,ilon),&
                   plev, t_base(:,ilat,ilon), q_base(:,ilat,ilon), o3_base(:,ilat,ilon), cldfrac_base(:,ilat,ilon), cldlwc_base(:,ilat,ilon),&
                   cldiwc_base(:,ilat,ilon), tauaer_sw_base(:,ilat,ilon,:), ssaaer_sw_base(:,ilat,ilon,:), asmaer_sw_base(:,ilat,ilon,:),&
                   tauaer_lw_base(:,ilat,ilon,:), fds_1d, fus_1d, htr_sw_1d,&
                   fdl_1d, ful_1d, htr_lw_1d, fd_1d_ts, fu_1d_ts, htr_1d, rad_1d_ts, lw_1d, sw_1d)

               ! o3
             call rad_driver(nlayer, iaer, icld, co2ppmv_base, ch4ppmv, n2oppmv, ps_base(ilat,ilon)/100.0,&
                   ts_base(ilat,ilon), zenith_base(ilat,ilon), albedo_sw_base(ilat,ilon), albedo_lw_base(ilat,ilon),&
                   plev, t_base(:,ilat,ilon), q_base(:,ilat,ilon), o3_warm(:,ilat,ilon), cldfrac_base(:,ilat,ilon), cldlwc_base(:,ilat,ilon),&
                   cldiwc_base(:,ilat,ilon), tauaer_sw_base(:,ilat,ilon,:), ssaaer_sw_base(:,ilat,ilon,:), asmaer_sw_base(:,ilat,ilon,:),&
                   tauaer_lw_base(:,ilat,ilon,:), fds_1d, fus_1d, htr_sw_1d,&
                   fdl_1d, ful_1d, htr_lw_1d, fd_1d_o3, fu_1d_o3, htr_1d, rad_1d_o3, lw_1d, sw_1d)

               ! solar
             call rad_driver(nlayer, iaer, icld, co2ppmv_base, ch4ppmv, n2oppmv, ps_base(ilat,ilon)/100.0,&
                   ts_base(ilat,ilon), zenith_warm(ilat,ilon), albedo_sw_base(ilat,ilon), albedo_lw_base(ilat,ilon),&
                   plev, t_base(:,ilat,ilon), q_base(:,ilat,ilon), o3_base(:,ilat,ilon), cldfrac_base(:,ilat,ilon), cldlwc_base(:,ilat,ilon),&
                   cldiwc_base(:,ilat,ilon), tauaer_sw_base(:,ilat,ilon,:), ssaaer_sw_base(:,ilat,ilon,:), asmaer_sw_base(:,ilat,ilon,:),&
                   tauaer_lw_base(:,ilat,ilon,:), fds_1d, fus_1d, htr_sw_1d,&
                   fdl_1d, ful_1d, htr_lw_1d, fd_1d_solar, fu_1d_solar, htr_1d, rad_1d_solar, lw_1d, sw_1d)

               ! albedo
             call rad_driver(nlayer, iaer, icld, co2ppmv_base, ch4ppmv, n2oppmv, ps_base(ilat,ilon)/100.0,&
                   ts_base(ilat,ilon), zenith_base(ilat,ilon), albedo_sw_warm(ilat,ilon), albedo_lw_base(ilat,ilon),&
                   plev, t_base(:,ilat,ilon), q_base(:,ilat,ilon), o3_base(:,ilat,ilon), cldfrac_base(:,ilat,ilon), cldlwc_base(:,ilat,ilon),&
                   cldiwc_base(:,ilat,ilon), tauaer_sw_base(:,ilat,ilon,:), ssaaer_sw_base(:,ilat,ilon,:), asmaer_sw_base(:,ilat,ilon,:),&
                   tauaer_lw_base(:,ilat,ilon,:), fds_1d, fus_1d, htr_sw_1d,&
                   fdl_1d, ful_1d, htr_lw_1d, fd_1d_albedo, fu_1d_albedo, htr_1d, rad_1d_albedo, lw_1d, sw_1d)

               ! cloud
             call rad_driver(nlayer, iaer, icld, co2ppmv_base, ch4ppmv, n2oppmv, ps_base(ilat,ilon)/100.0,&
                   ts_base(ilat,ilon), zenith_base(ilat,ilon), albedo_sw_base(ilat,ilon), albedo_lw_base(ilat,ilon),&
                   plev, t_base(:,ilat,ilon), q_base(:,ilat,ilon), o3_base(:,ilat,ilon), cldfrac_warm(:,ilat,ilon), cldlwc_warm(:,ilat,ilon),&
                   cldiwc_warm(:,ilat,ilon), tauaer_sw_base(:,ilat,ilon,:), ssaaer_sw_base(:,ilat,ilon,:), asmaer_sw_base(:,ilat,ilon,:),&
                   tauaer_lw_base(:,ilat,ilon,:), fds_1d, fus_1d, htr_sw_1d,&
                   fdl_1d, ful_1d, htr_lw_1d, fd_1d_cloud, fu_1d_cloud, htr_1d, rad_1d_cloud, lw_1d, sw_1d)

               ! Save cloud-state LW/SW components for LW/SW forcing split.
               ! (Shared workspace fds/fus/fdl/ful/lw_1d/sw_1d is about to be
               ! overwritten by the aerosol call, so snapshot now.)
               lw_1d_cloud  = lw_1d
               sw_1d_cloud  = sw_1d
               fdl_1d_cloud = fdl_1d
               ful_1d_cloud = ful_1d
               fds_1d_cloud = fds_1d
               fus_1d_cloud = fus_1d

              ! Aerosol
               call rad_driver(nlayer, iaer, icld, co2ppmv_base, ch4ppmv, n2oppmv, ps_base(ilat,ilon)/100.0,&
                   ts_base(ilat,ilon), zenith_base(ilat,ilon), albedo_sw_base(ilat,ilon), albedo_lw_base(ilat,ilon),&
                   plev, t_base(:,ilat,ilon), q_base(:,ilat,ilon), o3_base(:,ilat,ilon), cldfrac_base(:,ilat,ilon), cldlwc_base(:,ilat,ilon),&
                   cldiwc_base(:,ilat,ilon), tauaer_sw_warm(:,ilat,ilon,:), ssaaer_sw_warm(:,ilat,ilon,:), asmaer_sw_warm(:,ilat,ilon,:),&
                   tauaer_lw_warm(:,ilat,ilon,:), fds_1d, fus_1d, htr_sw_1d,&
                   fdl_1d, ful_1d, htr_lw_1d, fd_1d_aerosol, fu_1d_aerosol, htr_1d, rad_1d_aerosol, lw_1d, sw_1d)

               ! Per-species aerosol perturbation (Phase 3).
               ! For each species isp, build a mixed optical state where only isp is
               ! at warm, all others stay at base. Use AOD-weighted combination:
               !   tau_mix      = tau_bulk_base - tau_base_spc(isp) + tau_warm_spc(isp)
               !   tau_ssa_mix  = sum_j (tau_j * ssa_j) with isp warm, others base
               !   ssa_mix      = tau_ssa_mix / tau_mix
               !   g_mix        = tau_ssa_g_mix / tau_ssa_mix
               ! Precompute base bulk products (reused across species):
               tau_ssa_base_sw(:,:)   = tauaer_sw_base(:,ilat,ilon,:) * ssaaer_sw_base(:,ilat,ilon,:)
               tau_ssa_g_base_sw(:,:) = tau_ssa_base_sw(:,:) * asmaer_sw_base(:,ilat,ilon,:)
               do isp = 1, nspecies
                   ! SW: swap species isp from base -> warm
                   tauaer_sw_mix(:,:) = tauaer_sw_base(:,ilat,ilon,:) &
                       - tauaer_sw_base_spc(:,ilat,ilon,:,isp) + tauaer_sw_warm_spc(:,ilat,ilon,:,isp)
                   tau_ssa_mix(:,:) = tau_ssa_base_sw(:,:) &
                       - tauaer_sw_base_spc(:,ilat,ilon,:,isp) * ssaaer_sw_base_spc(:,ilat,ilon,:,isp) &
                       + tauaer_sw_warm_spc(:,ilat,ilon,:,isp) * ssaaer_sw_warm_spc(:,ilat,ilon,:,isp)
                   tau_ssa_g_mix(:,:) = tau_ssa_g_base_sw(:,:) &
                       - tauaer_sw_base_spc(:,ilat,ilon,:,isp) * ssaaer_sw_base_spc(:,ilat,ilon,:,isp) * asmaer_sw_base_spc(:,ilat,ilon,:,isp) &
                       + tauaer_sw_warm_spc(:,ilat,ilon,:,isp) * ssaaer_sw_warm_spc(:,ilat,ilon,:,isp) * asmaer_sw_warm_spc(:,ilat,ilon,:,isp)
                   where (tauaer_sw_mix > 0.0)
                       ssaaer_sw_mix = tau_ssa_mix / tauaer_sw_mix
                   elsewhere
                       ssaaer_sw_mix = 0.0
                   end where
                   where (tau_ssa_mix > 0.0)
                       asmaer_sw_mix = tau_ssa_g_mix / tau_ssa_mix
                   elsewhere
                       asmaer_sw_mix = 0.0
                   end where
                   ! LW: tau only
                   tauaer_lw_mix(:,:) = tauaer_lw_base(:,ilat,ilon,:) &
                       - tauaer_lw_base_spc(:,ilat,ilon,:,isp) + tauaer_lw_warm_spc(:,ilat,ilon,:,isp)

                   call rad_driver(nlayer, iaer, icld, co2ppmv_base, ch4ppmv, n2oppmv, ps_base(ilat,ilon)/100.0,&
                       ts_base(ilat,ilon), zenith_base(ilat,ilon), albedo_sw_base(ilat,ilon), albedo_lw_base(ilat,ilon),&
                       plev, t_base(:,ilat,ilon), q_base(:,ilat,ilon), o3_base(:,ilat,ilon), cldfrac_base(:,ilat,ilon), cldlwc_base(:,ilat,ilon),&
                       cldiwc_base(:,ilat,ilon), tauaer_sw_mix, ssaaer_sw_mix, asmaer_sw_mix,&
                       tauaer_lw_mix, fds_1d, fus_1d, htr_sw_1d,&
                       fdl_1d, ful_1d, htr_lw_1d, fd_1d_spc, fu_1d_spc, htr_1d, rad_1d_spc, lw_1d, sw_1d)

                   frc_spc(1:nlayer, isp) = rad_1d_spc(1:nlayer) - rad_1d_base(1:nlayer)
                   frc_spc(nlayer+1, isp) = fd_1d_spc(nlayer+1) - fu_1d_spc(nlayer+1) &
                       - (fd_1d_base(nlayer+1) - fu_1d_base(nlayer+1))
               end do

               ! Compute forcing (W/m2)
               frc_warm(1:nlayer)    = rad_1d_warm(1:nlayer) - rad_1d_base(1:nlayer)
               frc_warm(nlayer+1)    = fd_1d_warm(nlayer+1)-fu_1d_warm(nlayer+1)-(fd_1d_base(nlayer+1)-fu_1d_base(nlayer+1))
               frc_co2(1:nlayer)     = rad_1d_co2(1:nlayer) - rad_1d_base(1:nlayer)
               frc_co2(nlayer+1)     = fd_1d_co2(nlayer+1)-fu_1d_co2(nlayer+1)-(fd_1d_base(nlayer+1)-fu_1d_base(nlayer+1))
               frc_q(1:nlayer)       = rad_1d_q(1:nlayer) - rad_1d_base(1:nlayer)
               frc_q(nlayer+1)       = fd_1d_q(nlayer+1)-fu_1d_q(nlayer+1)-(fd_1d_base(nlayer+1)-fu_1d_base(nlayer+1))
               frc_ts(1:nlayer)      = rad_1d_ts(1:nlayer) - rad_1d_base(1:nlayer)
               frc_ts(nlayer+1)      = fd_1d_ts(nlayer+1)-fu_1d_ts(nlayer+1)-(fd_1d_base(nlayer+1)-fu_1d_base(nlayer+1))
               frc_o3(1:nlayer)      = rad_1d_o3(1:nlayer) - rad_1d_base(1:nlayer)
               frc_o3(nlayer+1)      = fd_1d_o3(nlayer+1)-fu_1d_o3(nlayer+1)-(fd_1d_base(nlayer+1)-fu_1d_base(nlayer+1))
               frc_solar(1:nlayer)   = rad_1d_solar(1:nlayer) - rad_1d_base(1:nlayer)
               frc_solar(nlayer+1)   = fd_1d_solar(nlayer+1)-fu_1d_solar(nlayer+1)-(fd_1d_base(nlayer+1)-fu_1d_base(nlayer+1))
               frc_albedo(1:nlayer)  = rad_1d_albedo(1:nlayer) - rad_1d_base(1:nlayer)
               frc_albedo(nlayer+1)  = fd_1d_albedo(nlayer+1)-fu_1d_albedo(nlayer+1)-(fd_1d_base(nlayer+1)-fu_1d_base(nlayer+1))
               frc_cloud(1:nlayer)   = rad_1d_cloud(1:nlayer) - rad_1d_base(1:nlayer)
               frc_cloud(nlayer+1)   = fd_1d_cloud(nlayer+1)-fu_1d_cloud(nlayer+1)-(fd_1d_base(nlayer+1)-fu_1d_base(nlayer+1))
               ! Cloud LW component: lw_1d_base is already saved; LW surface fluxes
               ! fdl/ful_1d_base are also saved. Direct subtraction.
               frc_cloud_lw(1:nlayer) = lw_1d_cloud(1:nlayer) - lw_1d_base(1:nlayer)
               frc_cloud_lw(nlayer+1) = (fdl_1d_cloud(nlayer+1) - ful_1d_cloud(nlayer+1)) &
                   - (fdl_1d_base(nlayer+1) - ful_1d_base(nlayer+1))
               ! Cloud SW component: SW base derived from (total - LW) base since
               ! sw_1d_base / fds_1d_base / fus_1d_base are not stored separately.
               frc_cloud_sw(1:nlayer) = sw_1d_cloud(1:nlayer) &
                   - (rad_1d_base(1:nlayer) - lw_1d_base(1:nlayer))
               frc_cloud_sw(nlayer+1) = (fds_1d_cloud(nlayer+1) - fus_1d_cloud(nlayer+1)) &
                   - ((fd_1d_base(nlayer+1) - fu_1d_base(nlayer+1)) &
                    - (fdl_1d_base(nlayer+1) - ful_1d_base(nlayer+1)))
               frc_aerosol(1:nlayer) = rad_1d_aerosol(1:nlayer) - rad_1d_base(1:nlayer)
               frc_aerosol(nlayer+1) = fd_1d_aerosol(nlayer+1)-fu_1d_aerosol(nlayer+1)-(fd_1d_base(nlayer+1)-fu_1d_base(nlayer+1))

               ! Store forcing output
               do ilayer = 1, nlayer
                  frc_warm_output(itime,ilayer,ilat,ilon) = frc_warm(ilayer)
                  frc_co2_output(itime,ilayer,ilat,ilon)  = frc_co2(ilayer)
                  frc_q_output(itime,ilayer,ilat,ilon)    = frc_q(ilayer)
                  frc_ts_output(itime,ilayer,ilat,ilon)   = frc_ts(ilayer)
                  frc_o3_output(itime,ilayer,ilat,ilon)   = frc_o3(ilayer)
                  frc_solar_output(itime,ilayer,ilat,ilon) = frc_solar(ilayer)
                  frc_albedo_output(itime,ilayer,ilat,ilon)= frc_albedo(ilayer)
                  frc_cloud_output(itime,ilayer,ilat,ilon) = frc_cloud(ilayer)
                  frc_cloud_lw_output(itime,ilayer,ilat,ilon) = frc_cloud_lw(ilayer)
                  frc_cloud_sw_output(itime,ilayer,ilat,ilon) = frc_cloud_sw(ilayer)
                  frc_aerosol_output(itime,ilayer,ilat,ilon) = frc_aerosol(ilayer)
                  frc_bc_output(itime,ilayer,ilat,ilon)    = frc_spc(ilayer,1)
                  frc_ocphi_output(itime,ilayer,ilat,ilon) = frc_spc(ilayer,2)
                  frc_ocpho_output(itime,ilayer,ilat,ilon) = frc_spc(ilayer,3)
                  frc_sulf_output(itime,ilayer,ilat,ilon)  = frc_spc(ilayer,4)
                  frc_ss_output(itime,ilayer,ilat,ilon)    = frc_spc(ilayer,5)
                  frc_dust_output(itime,ilayer,ilat,ilon)  = frc_spc(ilayer,6)
               end do

               if (nlayer+1 .le. nlev) then
                  frc_warm_output(itime,nlayer+1:nlev,ilat,ilon) = -999.0
                  frc_co2_output(itime,nlayer+1:nlev,ilat,ilon)  = -999.0
                  frc_q_output(itime,nlayer+1:nlev,ilat,ilon)    = -999.0
                  frc_ts_output(itime,nlayer+1:nlev,ilat,ilon)   = -999.0
                  frc_o3_output(itime,nlayer+1:nlev,ilat,ilon)   = -999.0
                  frc_solar_output(itime,nlayer+1:nlev,ilat,ilon)= -999.0
                  frc_albedo_output(itime,nlayer+1:nlev,ilat,ilon)= -999.0
                  frc_cloud_output(itime,nlayer+1:nlev,ilat,ilon)= -999.0
                  frc_cloud_lw_output(itime,nlayer+1:nlev,ilat,ilon) = -999.0
                  frc_cloud_sw_output(itime,nlayer+1:nlev,ilat,ilon) = -999.0
                  frc_aerosol_output(itime,nlayer+1:nlev,ilat,ilon) = -999.0
                  frc_bc_output(itime,nlayer+1:nlev,ilat,ilon)    = -999.0
                  frc_ocphi_output(itime,nlayer+1:nlev,ilat,ilon) = -999.0
                  frc_ocpho_output(itime,nlayer+1:nlev,ilat,ilon) = -999.0
                  frc_sulf_output(itime,nlayer+1:nlev,ilat,ilon)  = -999.0
                  frc_ss_output(itime,nlayer+1:nlev,ilat,ilon)    = -999.0
                  frc_dust_output(itime,nlayer+1:nlev,ilat,ilon)  = -999.0
               end if

               ! Surface forcing
               frc_warm_output(itime,nlev+1,ilat,ilon) = frc_warm(nlayer+1)
               frc_co2_output(itime,nlev+1,ilat,ilon)  = frc_co2(nlayer+1)
               frc_q_output(itime,nlev+1,ilat,ilon)    = frc_q(nlayer+1)
               frc_ts_output(itime,nlev+1,ilat,ilon)   = frc_ts(nlayer+1)
               frc_o3_output(itime,nlev+1,ilat,ilon)   = frc_o3(nlayer+1)
               frc_solar_output(itime,nlev+1,ilat,ilon) = frc_solar(nlayer+1)
               frc_albedo_output(itime,nlev+1,ilat,ilon)= frc_albedo(nlayer+1)
               frc_cloud_output(itime,nlev+1,ilat,ilon) = frc_cloud(nlayer+1)
               frc_cloud_lw_output(itime,nlev+1,ilat,ilon) = frc_cloud_lw(nlayer+1)
               frc_cloud_sw_output(itime,nlev+1,ilat,ilon) = frc_cloud_sw(nlayer+1)
               frc_aerosol_output(itime,nlev+1,ilat,ilon) = frc_aerosol(nlayer+1)
               frc_bc_output(itime,nlev+1,ilat,ilon)    = frc_spc(nlayer+1,1)
               frc_ocphi_output(itime,nlev+1,ilat,ilon) = frc_spc(nlayer+1,2)
               frc_ocpho_output(itime,nlev+1,ilat,ilon) = frc_spc(nlayer+1,3)
               frc_sulf_output(itime,nlev+1,ilat,ilon)  = frc_spc(nlayer+1,4)
               frc_ss_output(itime,nlev+1,ilat,ilon)    = frc_spc(nlayer+1,5)
               frc_dust_output(itime,nlev+1,ilat,ilon)  = frc_spc(nlayer+1,6)

               ! Write Planck inverse matrix for this grid point
               open(98,file='data_output/drdt_inv.dat',access='sequential',form='unformatted',status='replace')
               write(98) nlayer
               write(98) ((drdt_inv(i,j), j=1,nlayer+1), i=1,nlayer+1)
               close(98)

           end do
       end do
!$OMP END PARALLEL DO
       print*,'Finish time slice', itime
   end do

    ! Write forcing only (dT solved in Python)
    call write_out_3d(frc_co2_output,    "frc_co2.dat     ",      ntime,nlev,nlat,nlon)
    call write_out_3d(frc_o3_output,     "frc_o3.dat      ",      ntime,nlev,nlat,nlon)
    call write_out_3d(frc_q_output,      "frc_q.dat       ",      ntime,nlev,nlat,nlon)
    call write_out_3d(frc_albedo_output, "frc_albedo.dat  ",      ntime,nlev,nlat,nlon)
    call write_out_3d(frc_ts_output,     "frc_ts.dat      ",      ntime,nlev,nlat,nlon)
    call write_out_3d(frc_cloud_output,  "frc_cloud.dat   ",      ntime,nlev,nlat,nlon)
    call write_out_3d(frc_cloud_lw_output,"frc_cloud_lw.dat",     ntime,nlev,nlat,nlon)
    call write_out_3d(frc_cloud_sw_output,"frc_cloud_sw.dat",     ntime,nlev,nlat,nlon)
    call write_out_3d(frc_aerosol_output,"frc_aerosol.dat ",      ntime,nlev,nlat,nlon)
    call write_out_3d(frc_warm_output,   "frc_warm.dat    ",      ntime,nlev,nlat,nlon)
    call write_out_3d(frc_solar_output,  "frc_solar.dat   ",      ntime,nlev,nlat,nlon)
    ! Per-species aerosol forcing (Phase 3)
    call write_out_3d(frc_bc_output,     "frc_bc.dat      ",      ntime,nlev,nlat,nlon)
    call write_out_3d(frc_ocphi_output,  "frc_ocphi.dat   ",      ntime,nlev,nlat,nlon)
    call write_out_3d(frc_ocpho_output,  "frc_ocpho.dat   ",      ntime,nlev,nlat,nlon)
    call write_out_3d(frc_sulf_output,   "frc_sulf.dat    ",      ntime,nlev,nlat,nlon)
    call write_out_3d(frc_ss_output,     "frc_ss.dat      ",      ntime,nlev,nlat,nlon)
    call write_out_3d(frc_dust_output,   "frc_dust.dat    ",      ntime,nlev,nlat,nlon)

 End Program CFRAM
