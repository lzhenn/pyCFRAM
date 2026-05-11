 Program CFRAM

    use parkind,   only:im=>kind_im, rb=>kind_rb
    use parrrtm,   only:nbndlw
    use parrrsw,   only:jpband, jpb1, jpb2
    use output,    only:write_out_3d
    use math,      only:inv

    implicit none

    integer(kind=im),parameter :: nlat=1,nlon=1,ntime=1
    integer(kind=im),parameter :: icld=2   ! 0-no cloud, 2-cloud
    integer(kind=im),parameter :: iaer=10   ! 0-no aerosol, 10-aerosol
    integer(kind=im),parameter :: nspecies=6 ! per-species: bc,ocphi,ocpho,sulf,ss,dust
    real(kind=rb),   parameter :: scon=1360.98 ! solar constant
    integer(kind=im)           :: nlev      ! determined at runtime from plev.dat size
    integer(kind=im)           :: plev_fsize
    integer(kind=im)           :: ilat, ilon, ilev, itime, ib, i, j
    integer(kind=im)           :: ilayer,nlayer,isp

    real(kind=rb) :: co2ppmv_base,co2ppmv_warm,ch4ppmv,n2oppmv
    real(kind=rb), allocatable :: plev(:)
    real(kind=rb) :: ts_base(nlat,nlon),ts_warm(nlat,nlon)
    real(kind=rb) :: ps_base(nlat,nlon),ps_warm(nlat,nlon)
    real(kind=rb) :: albedo_lw_base(nlat,nlon),albedo_lw_warm(nlat,nlon)
    real(kind=rb) :: albedo_sw_base(nlat,nlon),albedo_sw_warm(nlat,nlon)
    real(kind=rb) :: solin_base(nlat,nlon), solin_warm(nlat,nlon)
    real(kind=rb) :: zenith_base(nlat,nlon), zenith_warm(nlat,nlon)
    real(kind=rb) :: swdn_surf_base(nlat,nlon),swdn_surf_warm(nlat,nlon)
    real(kind=rb) :: swup_surf_base(nlat,nlon),swup_surf_warm(nlat,nlon)
    real(kind=rb), allocatable :: t_base(:,:,:), t_warm(:,:,:)
    real(kind=rb), allocatable :: q_base(:,:,:), q_warm(:,:,:)
    real(kind=rb), allocatable :: o3_base(:,:,:),o3_warm(:,:,:)
    real(kind=rb), allocatable :: cldfrac_base(:,:,:), cldfrac_warm(:,:,:)
    real(kind=rb), allocatable :: cldlwc_base(:,:,:), cldlwc_warm(:,:,:)
    real(kind=rb), allocatable :: cldiwc_base(:,:,:), cldiwc_warm(:,:,:)
    real(kind=rb), allocatable :: tauaer_sw_base(:,:,:,:),tauaer_sw_warm(:,:,:,:)
    real(kind=rb), allocatable :: ssaaer_sw_base(:,:,:,:),ssaaer_sw_warm(:,:,:,:)
    real(kind=rb), allocatable :: asmaer_sw_base(:,:,:,:),asmaer_sw_warm(:,:,:,:)
    real(kind=rb), allocatable :: tauaer_lw_base(:,:,:,:),tauaer_lw_warm(:,:,:,:)

    ! Per-species aerosol optical properties. Species order must match Python writer
    ! (SPECIES_ORDER in run_parallel_python.py): 1=bc, 2=ocphi, 3=ocpho, 4=sulf, 5=ss, 6=dust.
    real(kind=rb), allocatable :: tauaer_sw_base_spc(:,:,:,:,:),tauaer_sw_warm_spc(:,:,:,:,:)
    real(kind=rb), allocatable :: ssaaer_sw_base_spc(:,:,:,:,:),ssaaer_sw_warm_spc(:,:,:,:,:)
    real(kind=rb), allocatable :: asmaer_sw_base_spc(:,:,:,:,:),asmaer_sw_warm_spc(:,:,:,:,:)
    real(kind=rb), allocatable :: tauaer_lw_base_spc(:,:,:,:,:),tauaer_lw_warm_spc(:,:,:,:,:)

    !Output: forcing only (dT solved in Python). All (ntime, nlev+1, nlat, nlon).
    real(kind=rb), allocatable :: frc_warm_output(:,:,:,:), frc_co2_output(:,:,:,:)
    real(kind=rb), allocatable :: frc_q_output(:,:,:,:), frc_ts_output(:,:,:,:)
    real(kind=rb), allocatable :: frc_o3_output(:,:,:,:), frc_solar_output(:,:,:,:)
    real(kind=rb), allocatable :: frc_albedo_output(:,:,:,:), frc_cloud_output(:,:,:,:)
    real(kind=rb), allocatable :: frc_aerosol_output(:,:,:,:)
    real(kind=rb), allocatable :: frc_bc_output(:,:,:,:), frc_ocphi_output(:,:,:,:)
    real(kind=rb), allocatable :: frc_ocpho_output(:,:,:,:), frc_sulf_output(:,:,:,:)
    real(kind=rb), allocatable :: frc_ss_output(:,:,:,:), frc_dust_output(:,:,:,:)
    real(kind=rb), allocatable :: frc_cloud_lw_output(:,:,:,:), frc_cloud_sw_output(:,:,:,:)

    ! Thread-private work arrays (size nlev or nlev+1)
    real(kind=rb), allocatable :: fds_1d(:), fus_1d(:), htr_sw_1d(:)
    real(kind=rb), allocatable :: fdl_1d(:), ful_1d(:), htr_lw_1d(:)
    real(kind=rb), allocatable :: htr_1d(:), lw_1d(:), sw_1d(:)
    real(kind=rb), allocatable :: rad_1d_base(:), fd_1d_base(:), fu_1d_base(:)
    real(kind=rb), allocatable :: rad_1d_warm(:), fd_1d_warm(:), fu_1d_warm(:)
    real(kind=rb), allocatable :: rad_1d_t(:),    fd_1d_t(:),    fu_1d_t(:)
    real(kind=rb), allocatable :: rad_1d_q(:),    fd_1d_q(:),    fu_1d_q(:)
    real(kind=rb), allocatable :: rad_1d_ts(:),   fd_1d_ts(:),   fu_1d_ts(:)
    real(kind=rb), allocatable :: rad_1d_o3(:),   fd_1d_o3(:),   fu_1d_o3(:)
    real(kind=rb), allocatable :: rad_1d_solar(:),fd_1d_solar(:),fu_1d_solar(:)
    real(kind=rb), allocatable :: rad_1d_albedo(:),fd_1d_albedo(:),fu_1d_albedo(:)
    real(kind=rb), allocatable :: rad_1d_cloud(:), fd_1d_cloud(:), fu_1d_cloud(:)
    real(kind=rb), allocatable :: rad_1d_aerosol(:),fd_1d_aerosol(:),fu_1d_aerosol(:)
    real(kind=rb), allocatable :: rad_1d_co2(:),   fd_1d_co2(:),   fu_1d_co2(:)
    real(kind=rb), allocatable :: lw_1d_base(:),   fdl_1d_base(:), ful_1d_base(:)

    real(kind=rb), allocatable :: frc_co2(:), frc_t(:), frc_q(:), frc_albedo(:)
    real(kind=rb), allocatable :: frc_ts(:), frc_o3(:), frc_solar(:), frc_cloud(:)
    real(kind=rb), allocatable :: frc_warm(:), frc_aerosol(:)
    real(kind=rb), allocatable :: frc_cloud_lw(:), frc_cloud_sw(:)
    ! Cloud-state LW/SW snapshot
    real(kind=rb), allocatable :: lw_1d_cloud(:), sw_1d_cloud(:)
    real(kind=rb), allocatable :: fdl_1d_cloud(:), ful_1d_cloud(:)
    real(kind=rb), allocatable :: fds_1d_cloud(:), fus_1d_cloud(:)

    ! Lu/Cai full-state perturbation arrays
    real(kind=rb), allocatable :: rad_1d_full(:), fd_1d_full(:), fu_1d_full(:)
    real(kind=rb), allocatable :: lw_1d_full(:), sw_1d_full(:)
    real(kind=rb), allocatable :: fdl_1d_full(:), ful_1d_full(:)
    real(kind=rb), allocatable :: fds_1d_full(:), fus_1d_full(:)
    real(kind=rb), allocatable :: frc_full(:)
    real(kind=rb), allocatable :: frc_full_output(:,:,:,:)

    ! Skip flags (read from data_prep/skip_flags.txt at startup; default = run all)
    logical :: skip_aerosol

    ! Per-species perturbation work arrays
    real(kind=rb), allocatable :: tauaer_sw_mix(:,:), ssaaer_sw_mix(:,:), asmaer_sw_mix(:,:)
    real(kind=rb), allocatable :: tauaer_lw_mix(:,:)
    real(kind=rb), allocatable :: tau_ssa_base_sw(:,:), tau_ssa_g_base_sw(:,:)
    real(kind=rb), allocatable :: tau_ssa_mix(:,:), tau_ssa_g_mix(:,:)
    real(kind=rb), allocatable :: rad_1d_spc(:), fd_1d_spc(:), fu_1d_spc(:)
    real(kind=rb), allocatable :: frc_spc(:,:)

    real(kind=rb), allocatable :: drdt(:,:), drdt_inv(:,:)

    ! Second-order CFRAM (midstate-Planck) option. When data_prep/drdt_midstate.flag
    ! exists, the Planck matrix is evaluated at the *midstate* (base+warm)/2
    ! rather than the base state — equivalent to a 2-term Taylor expansion
    ! centered at the midpoint, reducing the linearisation residual from
    ! O(ΔT²) (CFRAM-1 around base) to O((ΔT/2)²) — ~4× improvement.
    ! Cost: 1 extra rad_driver call per cell (negligible vs ~28 baseline).
    logical :: use_midstate_planck
    real(kind=rb), allocatable :: t_mid(:), q_mid(:), o3_mid(:)
    real(kind=rb), allocatable :: cldfrac_mid(:), cldlwc_mid(:), cldiwc_mid(:)
    real(kind=rb), allocatable :: tauaer_lw_mid(:,:)
    real(kind=rb), allocatable :: tauaer_sw_mid(:,:), ssaaer_sw_mid(:,:), asmaer_sw_mid(:,:)
    real(kind=rb), allocatable :: fdl_1d_mid(:), ful_1d_mid(:)
    real(kind=rb), allocatable :: fds_1d_mid(:), fus_1d_mid(:)
    real(kind=rb), allocatable :: fd_1d_mid(:), fu_1d_mid(:)
    real(kind=rb), allocatable :: rad_1d_mid(:), lw_1d_mid(:), sw_1d_mid(:)
    real(kind=rb) :: ts_mid, ps_mid_hPa, albedo_lw_mid, albedo_sw_mid
    real(kind=rb) :: zenith_mid, co2_mid

    ! Manabe q-feedback option (RRTMG only). When data_prep/q_feedback.flag
    ! exists, the q partial is recomputed in the warm-T atmosphere:
    !   frc_q = R(q_warm + T_warm + others_base) − R(q_base + T_warm + others_base)
    ! Physically meaningful q decomposition when q is a feedback tied to T
    ! (Manabe RH-fixed). Avoids the supersaturation artifact in the default
    ! formulation (q_warm in T_base atmosphere). Cost: +1 rad_driver call per
    ! cell (the warm-T q-base reference state).
    logical :: use_q_feedback
    real(kind=rb), allocatable :: rad_1d_q_ref(:), fd_1d_q_ref(:), fu_1d_q_ref(:)

    !-------------------------------- Determine nlev from plev.dat size --------
    ! plev.dat is nlev × float64 (8 bytes). Inferring nlev from file size lets a
    ! single binary handle any vertical grid (per-case override via plev.dat).
    inquire(file='data_prep/plev.dat', size=plev_fsize)
    if (plev_fsize <= 0) then
        print *, 'ERROR: data_prep/plev.dat missing or empty'
        stop 1
    end if
    nlev = plev_fsize / 8
    print *, 'Runtime nlev =', nlev, ' (from plev.dat size =', plev_fsize, 'bytes)'

    !-------------------------------- Allocate all nlev-dependent arrays -------
    allocate(plev(nlev))
    allocate(t_base(nlev,nlat,nlon), t_warm(nlev,nlat,nlon))
    allocate(q_base(nlev,nlat,nlon), q_warm(nlev,nlat,nlon))
    allocate(o3_base(nlev,nlat,nlon), o3_warm(nlev,nlat,nlon))
    allocate(cldfrac_base(nlev,nlat,nlon), cldfrac_warm(nlev,nlat,nlon))
    allocate(cldlwc_base(nlev,nlat,nlon), cldlwc_warm(nlev,nlat,nlon))
    allocate(cldiwc_base(nlev,nlat,nlon), cldiwc_warm(nlev,nlat,nlon))
    allocate(tauaer_sw_base(nlev,nlat,nlon,jpband), tauaer_sw_warm(nlev,nlat,nlon,jpband))
    allocate(ssaaer_sw_base(nlev,nlat,nlon,jpband), ssaaer_sw_warm(nlev,nlat,nlon,jpband))
    allocate(asmaer_sw_base(nlev,nlat,nlon,jpband), asmaer_sw_warm(nlev,nlat,nlon,jpband))
    allocate(tauaer_lw_base(nlev,nlat,nlon,nbndlw), tauaer_lw_warm(nlev,nlat,nlon,nbndlw))
    allocate(tauaer_sw_base_spc(nlev,nlat,nlon,jpband,nspecies))
    allocate(tauaer_sw_warm_spc(nlev,nlat,nlon,jpband,nspecies))
    allocate(ssaaer_sw_base_spc(nlev,nlat,nlon,jpband,nspecies))
    allocate(ssaaer_sw_warm_spc(nlev,nlat,nlon,jpband,nspecies))
    allocate(asmaer_sw_base_spc(nlev,nlat,nlon,jpband,nspecies))
    allocate(asmaer_sw_warm_spc(nlev,nlat,nlon,jpband,nspecies))
    allocate(tauaer_lw_base_spc(nlev,nlat,nlon,nbndlw,nspecies))
    allocate(tauaer_lw_warm_spc(nlev,nlat,nlon,nbndlw,nspecies))
    allocate(frc_warm_output(ntime,nlev+1,nlat,nlon), frc_co2_output(ntime,nlev+1,nlat,nlon))
    allocate(frc_q_output(ntime,nlev+1,nlat,nlon), frc_ts_output(ntime,nlev+1,nlat,nlon))
    allocate(frc_o3_output(ntime,nlev+1,nlat,nlon), frc_solar_output(ntime,nlev+1,nlat,nlon))
    allocate(frc_albedo_output(ntime,nlev+1,nlat,nlon), frc_cloud_output(ntime,nlev+1,nlat,nlon))
    allocate(frc_aerosol_output(ntime,nlev+1,nlat,nlon))
    allocate(frc_bc_output(ntime,nlev+1,nlat,nlon), frc_ocphi_output(ntime,nlev+1,nlat,nlon))
    allocate(frc_ocpho_output(ntime,nlev+1,nlat,nlon), frc_sulf_output(ntime,nlev+1,nlat,nlon))
    allocate(frc_ss_output(ntime,nlev+1,nlat,nlon), frc_dust_output(ntime,nlev+1,nlat,nlon))
    allocate(frc_cloud_lw_output(ntime,nlev+1,nlat,nlon), frc_cloud_sw_output(ntime,nlev+1,nlat,nlon))
    allocate(fds_1d(nlev+1), fus_1d(nlev+1), htr_sw_1d(nlev))
    allocate(fdl_1d(nlev+1), ful_1d(nlev+1), htr_lw_1d(nlev))
    allocate(htr_1d(nlev), lw_1d(nlev), sw_1d(nlev))
    allocate(rad_1d_base(nlev), fd_1d_base(nlev+1), fu_1d_base(nlev+1))
    allocate(rad_1d_warm(nlev), fd_1d_warm(nlev+1), fu_1d_warm(nlev+1))
    allocate(rad_1d_t(nlev),    fd_1d_t(nlev+1),    fu_1d_t(nlev+1))
    allocate(rad_1d_q(nlev),    fd_1d_q(nlev+1),    fu_1d_q(nlev+1))
    allocate(rad_1d_ts(nlev),   fd_1d_ts(nlev+1),   fu_1d_ts(nlev+1))
    allocate(rad_1d_o3(nlev),   fd_1d_o3(nlev+1),   fu_1d_o3(nlev+1))
    allocate(rad_1d_solar(nlev),fd_1d_solar(nlev+1),fu_1d_solar(nlev+1))
    allocate(rad_1d_albedo(nlev),fd_1d_albedo(nlev+1),fu_1d_albedo(nlev+1))
    allocate(rad_1d_cloud(nlev), fd_1d_cloud(nlev+1), fu_1d_cloud(nlev+1))
    allocate(rad_1d_aerosol(nlev),fd_1d_aerosol(nlev+1),fu_1d_aerosol(nlev+1))
    allocate(rad_1d_co2(nlev),   fd_1d_co2(nlev+1),   fu_1d_co2(nlev+1))
    allocate(lw_1d_base(nlev),   fdl_1d_base(nlev+1), ful_1d_base(nlev+1))
    allocate(frc_co2(nlev+1), frc_t(nlev+1), frc_q(nlev+1), frc_albedo(nlev+1))
    allocate(frc_ts(nlev+1), frc_o3(nlev+1), frc_solar(nlev+1), frc_cloud(nlev+1))
    allocate(frc_warm(nlev+1), frc_aerosol(nlev+1))
    allocate(frc_cloud_lw(nlev+1), frc_cloud_sw(nlev+1))
    allocate(lw_1d_cloud(nlev), sw_1d_cloud(nlev))
    allocate(fdl_1d_cloud(nlev+1), ful_1d_cloud(nlev+1))
    allocate(fds_1d_cloud(nlev+1), fus_1d_cloud(nlev+1))
    allocate(rad_1d_full(nlev), fd_1d_full(nlev+1), fu_1d_full(nlev+1))
    allocate(lw_1d_full(nlev), sw_1d_full(nlev))
    allocate(fdl_1d_full(nlev+1), ful_1d_full(nlev+1))
    allocate(fds_1d_full(nlev+1), fus_1d_full(nlev+1))
    allocate(frc_full(nlev+1))
    allocate(frc_full_output(ntime,nlev+1,nlat,nlon))
    allocate(tauaer_sw_mix(nlev,jpband), ssaaer_sw_mix(nlev,jpband), asmaer_sw_mix(nlev,jpband))
    allocate(tauaer_lw_mix(nlev,nbndlw))
    allocate(tau_ssa_base_sw(nlev,jpband), tau_ssa_g_base_sw(nlev,jpband))
    allocate(tau_ssa_mix(nlev,jpband), tau_ssa_g_mix(nlev,jpband))
    allocate(rad_1d_spc(nlev), fd_1d_spc(nlev+1), fu_1d_spc(nlev+1))
    allocate(frc_spc(nlev+1,nspecies))
    allocate(drdt(nlev+1,nlev+1), drdt_inv(nlev+1,nlev+1))

    ! Midstate-Planck buffers — only really used when use_midstate_planck=.true.,
    ! but allocate unconditionally so the code paths can share variables.
    allocate(t_mid(nlev), q_mid(nlev), o3_mid(nlev))
    allocate(cldfrac_mid(nlev), cldlwc_mid(nlev), cldiwc_mid(nlev))
    allocate(tauaer_lw_mid(nlev,nbndlw))
    allocate(tauaer_sw_mid(nlev,jpband), ssaaer_sw_mid(nlev,jpband), asmaer_sw_mid(nlev,jpband))
    allocate(fdl_1d_mid(nlev+1), ful_1d_mid(nlev+1))
    allocate(fds_1d_mid(nlev+1), fus_1d_mid(nlev+1))
    allocate(fd_1d_mid(nlev+1), fu_1d_mid(nlev+1))
    allocate(rad_1d_mid(nlev), lw_1d_mid(nlev), sw_1d_mid(nlev))
    allocate(rad_1d_q_ref(nlev), fd_1d_q_ref(nlev+1), fu_1d_q_ref(nlev+1))

    ! Detect midstate-Planck mode via flag file. Python writes
    ! data_prep/drdt_midstate.flag when case.yaml has radiation.drdt_eval=midstate.
    inquire(file='data_prep/drdt_midstate.flag', exist=use_midstate_planck)
    inquire(file='data_prep/q_feedback.flag',    exist=use_q_feedback)
    if (use_midstate_planck) then
        print *, '  Planck matrix: midstate (2nd-order CFRAM, drdt at (T_base+T_warm)/2)'
    else
        print *, '  Planck matrix: base state (1st-order CFRAM, drdt at T_base)'
    end if
    if (use_q_feedback) then
        print *, '  q partial: feedback (Manabe RH-fixed: frc_q computed in T_warm atmosphere)'
    end if

    !-------------------------------- Read plev values --------
    open(unit=99, file='data_prep/plev.dat', form='unformatted', access='direct', recl=nlev*8)
    read(99, rec=1) plev
    close(99)

!-------------------------------- Read Skip Flags (optional) ----------------------------------
    ! data_prep/skip_flags.txt may contain lines like "skip_aerosol=1" to suppress
    ! perturbation calls when the variable is identical between base and warm states.
    ! Default (file absent or flag not set): run all perturbations.
    skip_aerosol = .false.
    block
        logical :: flag_exists
        character(len=128) :: line
        integer :: ios
        inquire(file='data_prep/skip_flags.txt', exist=flag_exists)
        if (flag_exists) then
            open(unit=97, file='data_prep/skip_flags.txt', status='old', action='read')
            do
                read(97, '(A)', iostat=ios) line
                if (ios /= 0) exit
                if (index(line, 'skip_aerosol=1') > 0) skip_aerosol = .true.
            end do
            close(97)
            print *, "skip_flags loaded: skip_aerosol=", skip_aerosol
        else
            print *, "no skip_flags.txt; running all perturbations"
        end if
    end block

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
!$OMP   rad_1d_full, lw_1d_full, sw_1d_full, &
!$OMP   fdl_1d_full, ful_1d_full, fds_1d_full, fus_1d_full, &
!$OMP   fd_1d_full, fu_1d_full, frc_full, &
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

               ! Planck Matrix — base state (default, 1st-order CFRAM) OR
               ! midstate (= (base+warm)/2, 2nd-order Taylor centered at midpoint).
               if (use_midstate_planck) then
                   ! Build midstate atmospheric profile + surface scalars.
                   t_mid(:)        = 0.5 * (t_base(:,ilat,ilon)       + t_warm(:,ilat,ilon))
                   q_mid(:)        = 0.5 * (q_base(:,ilat,ilon)       + q_warm(:,ilat,ilon))
                   o3_mid(:)       = 0.5 * (o3_base(:,ilat,ilon)      + o3_warm(:,ilat,ilon))
                   cldfrac_mid(:)  = 0.5 * (cldfrac_base(:,ilat,ilon) + cldfrac_warm(:,ilat,ilon))
                   cldlwc_mid(:)   = 0.5 * (cldlwc_base(:,ilat,ilon)  + cldlwc_warm(:,ilat,ilon))
                   cldiwc_mid(:)   = 0.5 * (cldiwc_base(:,ilat,ilon)  + cldiwc_warm(:,ilat,ilon))
                   tauaer_lw_mid(:,:) = 0.5 * (tauaer_lw_base(:,ilat,ilon,:) + tauaer_lw_warm(:,ilat,ilon,:))
                   tauaer_sw_mid(:,:) = 0.5 * (tauaer_sw_base(:,ilat,ilon,:) + tauaer_sw_warm(:,ilat,ilon,:))
                   ssaaer_sw_mid(:,:) = 0.5 * (ssaaer_sw_base(:,ilat,ilon,:) + ssaaer_sw_warm(:,ilat,ilon,:))
                   asmaer_sw_mid(:,:) = 0.5 * (asmaer_sw_base(:,ilat,ilon,:) + asmaer_sw_warm(:,ilat,ilon,:))
                   ts_mid          = 0.5 * (ts_base(ilat,ilon)         + ts_warm(ilat,ilon))
                   ps_mid_hPa      = 0.5 * (ps_base(ilat,ilon)         + ps_warm(ilat,ilon)) / 100.0
                   albedo_lw_mid   = 0.5 * (albedo_lw_base(ilat,ilon)  + albedo_lw_warm(ilat,ilon))
                   albedo_sw_mid   = 0.5 * (albedo_sw_base(ilat,ilon)  + albedo_sw_warm(ilat,ilon))
                   zenith_mid      = 0.5 * (zenith_base(ilat,ilon)     + zenith_warm(ilat,ilon))
                   co2_mid         = 0.5 * (co2ppmv_base               + co2ppmv_warm)

                   ! Baseline radiation at the midstate (the +1K perturbation pivots around this).
                   call rad_driver(nlayer, iaer, icld, co2_mid, ch4ppmv, n2oppmv, ps_mid_hPa,&
                       ts_mid, zenith_mid, albedo_sw_mid, albedo_lw_mid,&
                       plev, t_mid, q_mid, o3_mid, cldfrac_mid, cldlwc_mid,&
                       cldiwc_mid, tauaer_sw_mid, ssaaer_sw_mid, asmaer_sw_mid,&
                       tauaer_lw_mid, fds_1d_mid, fus_1d_mid, htr_sw_1d,&
                       fdl_1d_mid, ful_1d_mid, htr_lw_1d, fd_1d_mid, fu_1d_mid, htr_1d,&
                       rad_1d_mid, lw_1d_mid, sw_1d_mid)

                   ! Planck Jacobian at midstate.
                   call calc_drdt(nlayer, iaer, icld, co2_mid, ch4ppmv, n2oppmv, ps_mid_hPa,&
                       ts_mid, albedo_lw_mid,&
                       plev, t_mid, q_mid, o3_mid, cldfrac_mid, cldlwc_mid,&
                       cldiwc_mid, tauaer_lw_mid, lw_1d_mid, fdl_1d_mid, ful_1d_mid, &
                       drdt(1:nlayer+1,1:nlayer+1))
               else
                   ! Default: 1st-order CFRAM, Planck at base state.
                   call calc_drdt(nlayer, iaer, icld, co2ppmv_base, ch4ppmv, n2oppmv, ps_base(ilat,ilon)/100.0,&
                       ts_base(ilat,ilon),  albedo_lw_base(ilat,ilon),&
                       plev, t_base(:,ilat,ilon), q_base(:,ilat,ilon), o3_base(:,ilat,ilon), cldfrac_base(:,ilat,ilon), cldlwc_base(:,ilat,ilon),&
                       cldiwc_base(:,ilat,ilon), tauaer_lw_base(:,ilat,ilon,:), lw_1d_base, fdl_1d_base, ful_1d_base, &
                       drdt(1:nlayer+1,1:nlayer+1))
               end if

                  drdt_inv(1:nlayer+1,1:nlayer+1) = inv(drdt(1:nlayer+1,1:nlayer+1))

               ! CO2
               call rad_driver(nlayer, iaer, icld, co2ppmv_warm, ch4ppmv, n2oppmv, ps_base(ilat,ilon)/100.0,&
                   ts_base(ilat,ilon), zenith_base(ilat,ilon), albedo_sw_base(ilat,ilon), albedo_lw_base(ilat,ilon),&
                   plev, t_base(:,ilat,ilon), q_base(:,ilat,ilon), o3_base(:,ilat,ilon), cldfrac_base(:,ilat,ilon), cldlwc_base(:,ilat,ilon),&
                   cldiwc_base(:,ilat,ilon), tauaer_sw_base(:,ilat,ilon,:), ssaaer_sw_base(:,ilat,ilon,:), asmaer_sw_base(:,ilat,ilon,:),&
                   tauaer_lw_base(:,ilat,ilon,:), fds_1d, fus_1d, htr_sw_1d,&
                   fdl_1d, ful_1d, htr_lw_1d, fd_1d_co2, fu_1d_co2, htr_1d, rad_1d_co2, lw_1d, sw_1d)

               ! Q — water vapour partial.
               ! Default (q_handling=independent): q_warm in T_base atmosphere
               !   (standard CFRAM). frc_q subtracts rad_1d_base.
               ! Manabe (q_handling=feedback): q_warm in T_warm atmosphere, with
               !   an extra reference call R(q_base + T_warm + others_base).
               !   frc_q subtracts rad_1d_q_ref. This isolates the q radiative
               !   impact within the warm-T atmosphere, avoiding the
               !   supersaturation artifact of q_warm in cold T_base.
               if (use_q_feedback) then
                   ! Perturbed: q_warm + T_warm + others_base
                   call rad_driver(nlayer, iaer, icld, co2ppmv_base, ch4ppmv, n2oppmv, ps_base(ilat,ilon)/100.0,&
                       ts_base(ilat,ilon), zenith_base(ilat,ilon), albedo_sw_base(ilat,ilon), albedo_lw_base(ilat,ilon),&
                       plev, t_warm(:,ilat,ilon), q_warm(:,ilat,ilon), o3_base(:,ilat,ilon), cldfrac_base(:,ilat,ilon), cldlwc_base(:,ilat,ilon),&
                       cldiwc_base(:,ilat,ilon), tauaer_sw_base(:,ilat,ilon,:), ssaaer_sw_base(:,ilat,ilon,:), asmaer_sw_base(:,ilat,ilon,:),&
                       tauaer_lw_base(:,ilat,ilon,:), fds_1d, fus_1d, htr_sw_1d,&
                       fdl_1d, ful_1d, htr_lw_1d, fd_1d_q, fu_1d_q, htr_1d, rad_1d_q, lw_1d, sw_1d)
                   ! Reference: q_base + T_warm + others_base (extra rad_driver call)
                   call rad_driver(nlayer, iaer, icld, co2ppmv_base, ch4ppmv, n2oppmv, ps_base(ilat,ilon)/100.0,&
                       ts_base(ilat,ilon), zenith_base(ilat,ilon), albedo_sw_base(ilat,ilon), albedo_lw_base(ilat,ilon),&
                       plev, t_warm(:,ilat,ilon), q_base(:,ilat,ilon), o3_base(:,ilat,ilon), cldfrac_base(:,ilat,ilon), cldlwc_base(:,ilat,ilon),&
                       cldiwc_base(:,ilat,ilon), tauaer_sw_base(:,ilat,ilon,:), ssaaer_sw_base(:,ilat,ilon,:), asmaer_sw_base(:,ilat,ilon,:),&
                       tauaer_lw_base(:,ilat,ilon,:), fds_1d, fus_1d, htr_sw_1d,&
                       fdl_1d, ful_1d, htr_lw_1d, fd_1d_q_ref, fu_1d_q_ref, htr_1d, rad_1d_q_ref, lw_1d, sw_1d)
               else
                   ! Standard CFRAM: q_warm in T_base atmosphere; reference = base
                   call rad_driver(nlayer, iaer, icld, co2ppmv_base, ch4ppmv, n2oppmv, ps_base(ilat,ilon)/100.0,&
                       ts_base(ilat,ilon), zenith_base(ilat,ilon), albedo_sw_base(ilat,ilon), albedo_lw_base(ilat,ilon),&
                       plev, t_base(:,ilat,ilon), q_warm(:,ilat,ilon), o3_base(:,ilat,ilon), cldfrac_base(:,ilat,ilon), cldlwc_base(:,ilat,ilon),&
                       cldiwc_base(:,ilat,ilon), tauaer_sw_base(:,ilat,ilon,:), ssaaer_sw_base(:,ilat,ilon,:), asmaer_sw_base(:,ilat,ilon,:),&
                       tauaer_lw_base(:,ilat,ilon,:), fds_1d, fus_1d, htr_sw_1d,&
                       fdl_1d, ful_1d, htr_lw_1d, fd_1d_q, fu_1d_q, htr_1d, rad_1d_q, lw_1d, sw_1d)
               end if

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

              ! Aerosol — skipped if global skip_aerosol flag set (saves 1 + nspecies rad_driver calls per column)
              if (.not. skip_aerosol) then
               call rad_driver(nlayer, iaer, icld, co2ppmv_base, ch4ppmv, n2oppmv, ps_base(ilat,ilon)/100.0,&
                   ts_base(ilat,ilon), zenith_base(ilat,ilon), albedo_sw_base(ilat,ilon), albedo_lw_base(ilat,ilon),&
                   plev, t_base(:,ilat,ilon), q_base(:,ilat,ilon), o3_base(:,ilat,ilon), cldfrac_base(:,ilat,ilon), cldlwc_base(:,ilat,ilon),&
                   cldiwc_base(:,ilat,ilon), tauaer_sw_warm(:,ilat,ilon,:), ssaaer_sw_warm(:,ilat,ilon,:), asmaer_sw_warm(:,ilat,ilon,:),&
                   tauaer_lw_warm(:,ilat,ilon,:), fds_1d, fus_1d, htr_sw_1d,&
                   fdl_1d, ful_1d, htr_lw_1d, fd_1d_aerosol, fu_1d_aerosol, htr_1d, rad_1d_aerosol, lw_1d, sw_1d)
              else
                  ! Identical to base (no perturbation): copy base diagnostics so frc_aerosol = 0
                  rad_1d_aerosol = rad_1d_base
                  fd_1d_aerosol  = fd_1d_base
                  fu_1d_aerosol  = fu_1d_base
              end if

               ! Per-species aerosol perturbation (Phase 3) — same skip guard.
               if (skip_aerosol) then
                   frc_spc = 0.0
                   ! fall through to forcing computation; per-species frcs are zero
                   goto 7777
               end if
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
7777           continue   ! per-species skip target

               ! ============== Lu/Cai full-state perturbation ==============
               ! All variables at perturbed values INCLUDING T_warm, ts_warm.
               ! Used to compute ΔQ_dyn = (rad_1d_full - rad_1d_base) for closure.
               ! Aerosol always uses warm optics here (aerosol skip only affects the
               ! standalone aerosol partial; the bulk warm state is part of the full).
               call rad_driver(nlayer, iaer, icld, co2ppmv_warm, ch4ppmv, n2oppmv, ps_warm(ilat,ilon)/100.0,&
                   ts_warm(ilat,ilon), zenith_warm(ilat,ilon), albedo_sw_warm(ilat,ilon), albedo_lw_warm(ilat,ilon),&
                   plev, t_warm(:,ilat,ilon), q_warm(:,ilat,ilon), o3_warm(:,ilat,ilon), cldfrac_warm(:,ilat,ilon), cldlwc_warm(:,ilat,ilon),&
                   cldiwc_warm(:,ilat,ilon), tauaer_sw_warm(:,ilat,ilon,:), ssaaer_sw_warm(:,ilat,ilon,:), asmaer_sw_warm(:,ilat,ilon,:),&
                   tauaer_lw_warm(:,ilat,ilon,:), fds_1d_full, fus_1d_full, htr_sw_1d,&
                   fdl_1d_full, ful_1d_full, htr_lw_1d, fd_1d_full, fu_1d_full, htr_1d, rad_1d_full, lw_1d_full, sw_1d_full)

               ! Compute forcing (W/m2)
               frc_warm(1:nlayer)    = rad_1d_warm(1:nlayer) - rad_1d_base(1:nlayer)
               frc_warm(nlayer+1)    = fd_1d_warm(nlayer+1)-fu_1d_warm(nlayer+1)-(fd_1d_base(nlayer+1)-fu_1d_base(nlayer+1))
               frc_co2(1:nlayer)     = rad_1d_co2(1:nlayer) - rad_1d_base(1:nlayer)
               frc_co2(nlayer+1)     = fd_1d_co2(nlayer+1)-fu_1d_co2(nlayer+1)-(fd_1d_base(nlayer+1)-fu_1d_base(nlayer+1))
               ! frc_q: standard CFRAM subtracts the base radiation; Manabe
               ! feedback subtracts the warm-T q-base reference (computed above).
               if (use_q_feedback) then
                   frc_q(1:nlayer) = rad_1d_q(1:nlayer) - rad_1d_q_ref(1:nlayer)
                   frc_q(nlayer+1) = fd_1d_q(nlayer+1)-fu_1d_q(nlayer+1)-(fd_1d_q_ref(nlayer+1)-fu_1d_q_ref(nlayer+1))
               else
                   frc_q(1:nlayer) = rad_1d_q(1:nlayer) - rad_1d_base(1:nlayer)
                   frc_q(nlayer+1) = fd_1d_q(nlayer+1)-fu_1d_q(nlayer+1)-(fd_1d_base(nlayer+1)-fu_1d_base(nlayer+1))
               end if
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

               ! Lu/Cai full-state forcing: rad_1d_full uses T_warm + all-perturbed.
               ! Python computes  dT_dry = -drdt_inv @ (-frc_full) = drdt_inv @ frc_full
               ! which captures ΔQ_dyn from full-state energy balance (Lu/Cai Eq.4).
               frc_full(1:nlayer) = rad_1d_full(1:nlayer) - rad_1d_base(1:nlayer)
               frc_full(nlayer+1) = fd_1d_full(nlayer+1)-fu_1d_full(nlayer+1)-(fd_1d_base(nlayer+1)-fu_1d_base(nlayer+1))

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
                  frc_full_output(itime,ilayer,ilat,ilon)    = frc_full(ilayer)
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
                  frc_full_output(itime,nlayer+1:nlev,ilat,ilon)    = -999.0
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
               frc_full_output(itime,nlev+1,ilat,ilon)    = frc_full(nlayer+1)
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
    call write_out_3d(frc_full_output,   "frc_full.dat    ",      ntime,nlev,nlat,nlon)
    call write_out_3d(frc_solar_output,  "frc_solar.dat   ",      ntime,nlev,nlat,nlon)
    ! Per-species aerosol forcing (Phase 3)
    call write_out_3d(frc_bc_output,     "frc_bc.dat      ",      ntime,nlev,nlat,nlon)
    call write_out_3d(frc_ocphi_output,  "frc_ocphi.dat   ",      ntime,nlev,nlat,nlon)
    call write_out_3d(frc_ocpho_output,  "frc_ocpho.dat   ",      ntime,nlev,nlat,nlon)
    call write_out_3d(frc_sulf_output,   "frc_sulf.dat    ",      ntime,nlev,nlat,nlon)
    call write_out_3d(frc_ss_output,     "frc_ss.dat      ",      ntime,nlev,nlat,nlon)
    call write_out_3d(frc_dust_output,   "frc_dust.dat    ",      ntime,nlev,nlat,nlon)

 End Program CFRAM
