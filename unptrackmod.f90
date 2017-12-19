!============================================================================**
!*---------**                 UNPTRACKMOD.F90                     **---------**
!*===========================================================================**
! Version 2.0                                                       June 2017 *
!-----------------------------------------------------------------------------*
!-----------------------------------------------------------------------------*
!
!=============================================================================*

      MODULE UNPTRACKMOD
      
      SAVE
      
!-| Grid node and elements
      integer :: gr_nn_sw, gr_ne_sw, gr_nghmax_sw
      integer :: gr_nlayers_sw = 1
      integer :: gr_jpcell_sw
      integer, allocatable :: gr_ntve_sw(:)
      integer, allocatable :: gr_nodes_sw(:,:)
      integer, allocatable :: gr_nbve_sw(:,:)
      real, allocatable :: gr_x_sw(:), gr_y_sw(:)
      real, allocatable :: gr_xc_sw(:), gr_yc_sw(:)
      real, allocatable :: gr_z_sw(:), gr_dz_sw(:)
      real, allocatable :: gr_depth_sw(:), gr_eledep_sw(:)
      real, allocatable :: gr_area_sw(:), gr_vol_sw(:,:)
      
!-|Time interval between velocity fields from hydrodynamic model (s)
      real :: tm_dt_sw

!-|Time
      integer :: tm_ks_sw, dd, hh, mm, ss

!-|Velocity array dimensions
      real, allocatable :: u(:,:)
      real, allocatable :: v(:,:)
      real, allocatable :: w(:,:)
      real, allocatable :: tem(:,:)
      
!-|Density and buoyancy arrays
      real :: rho0 = 1025.
      real, allocatable :: rho(:,:)

!-|Surface elevation
      real, allocatable :: eta(:)

!-|Wind-driven surface velocity field
      real, allocatable :: usw(:)
      real, allocatable :: vsw(:)

!-|Wind forcing data
      real :: winddata(8760,10)
      real :: xwind(5), ywind(5)
      integer :: nwinterval, nwsite

!-|Horizontal and vertical eddy diffusion
      real :: kv0 = 0.001
      real :: kh0 = 0.1
      real, allocatable :: tu_kh_sw(:,:)
      real, allocatable :: tu_kv_sw(:,:)

!-|Particle numbers
      integer :: NPMAX
      
!-|Particle density
      real, allocatable :: pden(:,:), psum(:,:)

!-|Time step of lagrangian model (s)
      real :: p_dt

!-|Duration of simulation
      real :: lsim

!-|Output frequency
      real :: fout

!-|Particle location and characteristics
      integer :: tr_n_sw = 0
      integer :: tr_init_sw = 0
      integer, allocatable :: tr_stage_sw(:)
      integer, allocatable :: tr_igo_sw(:)
      integer, allocatable :: tr_element_sw(:)
      real, allocatable :: tr_x_sw(:)
      real, allocatable :: tr_y_sw(:)
      real, allocatable :: tr_z_sw(:)
      real, allocatable :: tr_s_sw(:)
      real, allocatable :: tr_age_sw(:)
      real, allocatable :: tr_pfac_sw(:)
      real, allocatable :: tr_prho_sw(:)
      real, allocatable :: tr_pvol_sw(:)

!-|Particle sources
      integer :: tr_nsource_sw
      real, allocatable :: tr_x0_sw(:), tr_y0_sw(:), tr_z0_sw(:)
      real, allocatable :: tr_xrange_sw(:), tr_yrange_sw(:), tr_zrange_sw(:)
      real, allocatable :: tr_mass_sw(:), tr_volume_sw(:), tr_density_sw(:)
      real, allocatable :: tr_start_sw(:), tr_stop_sw(:)
      real, allocatable :: tr_wsettle_sw(:)

!-|Particle counting
      integer :: tr_jpcell_sw

!-|Variables for the tracking routine
      real :: tpassive, tmobile          ! Life cycle events
      real :: cmort, pconv, pfac, tdecay
      real :: wsink, wswim
      real :: wsettle = 0
      real :: hpmax
      real :: shrmin
      real, allocatable :: dzz_sv(:)
      integer :: np_steps
      integer :: JPCELL

!-|Switches
      integer :: alpha = 0                  ! Switch for growth
      integer :: beta = 0                   ! Switch for mortality
      integer :: wind_io = 0                ! Switch for direct wind forcing
      integer :: utype = 0                  ! Type of velocity data
      integer :: advscheme = 0              ! Advection scheme (rk4 or euler)
      integer :: buoyancy = 0               ! Switch to include buoyancy forcing on particles
      integer :: resus = 0                  ! Switch to include resuspension of particles

      END
      
