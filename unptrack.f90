!============================================================================**
!*---------**                    UNPTRACK.F90                     **---------**
!!===========================================================================**
! Version 1.0                                                       June 2017 *
!-----------------------------------------------------------------------------*
!                                                                             *
! Lagrangian (particle-tracking) model to simulate the transport pathways     *
! of pelagic biota or chemical contaminants. The model runs offline (i.e.     *
! velocity fields from a 3D hydrodynamic model are input) and uses a simple   *
! estimate of turbulent diffusion calculated from the current velocity.       *
!                                                                             *
! This version uses velocity fields from unstructured (flexible) mesh         *
! hydrodynamic models e.g. FVCOM, RiCOM.                                      *
!                                                                             *
!-----------------------------------------------------------------------------*

!=============================================================================*
      program unptrack
!=============================================================================*
!--External variables
      use unptrackmod
      
      implicit none

!--Description
!    The core of the model, which calls all necessary subroutines.

!--Initialise model parameters
      call unptrack_init

!--Run model
      call unptrack_run

!--Finish
      end

!=============================================================================*
      subroutine unptrack_init
!=============================================================================*
!--External variables
      use unptrackmod

      implicit none

!--Description
!   Initializes the ptrack module.

!--Local variable
      real :: mtotal = 0
      real :: H = 0
      integer :: i, j, k
      integer :: iun
      integer :: nwrec
      integer :: stat
      integer :: starttime
      integer :: ierr, myproc
      integer :: ilen, iloc1, iloc2
      integer :: inode, nn, ne
      character(80) :: dataline,nextline,fname,udata,adv_scheme

!...................................

      write(6,*) ' Initialising model'

      iun = 20
      open(iun,file='inputs/grid-input.txt',status='unknown')
      do
        read(iun,'(a)',IOSTAT=stat) nextline
        if(stat.gt.0) then
          write(*,*) 'ERROR reading file'
          ierr = 71
          call terminate(ierr,myproc)
        elseif(stat.lt.0) then
          exit
        endif
        ilen = len_trim(nextline)
!        write(*,*) 'ilen=',ilen
        if(ilen.eq.0) then !blank line
          cycle
        elseif(nextline(1:1).eq.'#') then !comment line
!          write(*,*) nextline
          cycle
        elseif(nextline(1:3).eq.'END') then !end of data
          exit
        else  !parse for keyword
          iloc1 = index(nextline,'=')
          iloc2 = index(nextline,' ')
          !          write(*,*) 'iloc=',ilen
          if(nextline(1:iloc1-1).eq.'NN') then
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) gr_nn_sw
            if(stat.ne.0) then
              write(*,*) 'ERROR reading NN parameter gr_nn_sw'
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' gr_nn_sw = ', &
                                     gr_nn_sw
            endif
          elseif(nextline(1:iloc1-1).eq.'NE') then
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) gr_ne_sw
            if(stat.ne.0) then
              write(*,*) 'ERROR reading NE parameter gr_ne_sw'
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' gr_ne_sw = ', &
                                     gr_ne_sw
            endif
          elseif(nextline(1:iloc1-1).eq.'NGHMAX') then
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) gr_nghmax_sw
            if(stat.ne.0) then
              write(*,*) 'ERROR reading NGHMAX parameter gr_nghmax_sw'
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1), &
                                     ' gr_nghmax_sw = ', gr_nghmax_sw
            endif
          elseif(nextline(1:iloc1-1).eq.'BATHY') then
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) H
            if(stat.ne.0) then
              write(*,*) 'ERROR reading BATHY parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' H = ', &
                                     H
            endif
          elseif(nextline(1:iloc1-1).eq.'NPMAX') then
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) NPMAX
            if(stat.ne.0) then
              write(*,*) 'ERROR reading NPMAX parameter'
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' = ', &
                                     NPMAX
            endif
          elseif(nextline(1:iloc1-1).eq.'JPCELL') then
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) tr_jpcell_sw
            if(stat.ne.0) then
              write(*,*) 'ERROR reading JPCELL parameter'
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' = ', &
                                     tr_jpcell_sw
            endif
          elseif(nextline(1:iloc1-1).eq.'NLAYERS') then
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) gr_nlayers_sw
            if(stat.ne.0) then
              write(*,*) 'ERROR reading NLAYERS parameter gr_nlayers_sw'
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' gr_nlayers_sw = ' &
                                    ,gr_nlayers_sw
              if (gr_nlayers_sw.gt.0) then
!               Read vertical level and layer depths
                allocate (gr_z_sw(gr_nlayers_sw), dzz_sv(gr_nlayers_sw+1))
                read(iun,*) (gr_z_sw(k), k = 1, gr_nlayers_sw)
                read(iun,*) (dzz_sv(k), k = 1, gr_nlayers_sw+1)
              endif
            endif
          endif
        endif
      enddo
      close(iun)
      
!     Allocate all required arrays
      call arrayalloc
    
!     Get the x and y grid coordinates, and element information 
      open(iun,file='inputs/grid-data.dat',status='unknown')
      read(iun,*) nn, ne
      write(*,*) 'Reading grid data: np = ',nn,' ne = ',ne
      if (nn .ne. gr_nn_sw .or. ne .ne. gr_ne_sw) then
        write(*,*) 'ERROR: Grid file does not match input file spec.'
        stop
      else
        do i = 1, gr_nn_sw
          read(iun,*) inode, gr_x_sw(i), gr_y_sw(i), gr_depth_sw(i)
        enddo
        do i = 1, gr_ne_sw
          read(iun,*) gr_nodes_sw(i,1),gr_nodes_sw(i,2),gr_nodes_sw(i,3) 
        enddo
      endif
      close(iun)
      
!     Get data on elements neighbouring each node
!      open(iun,file='inputs/grid-elements.dat',status='unknown')
!      read(iun,*) nn
!      write(*,*) 'Reading grid node neighbours: np = ',nn
!      if (nn .ne. gr_nn_sw) then
!        write(*,*) 'ERROR: Number of elements mis-match. Stopping'
!        stop
!      else
!        do i = 1, gr_nn_sw
!          read(iun,*) gr_ntve_sw(i),(gr_nbve_sw(i,k),k=1,gr_ntve_sw(i))
!          write(*,*) gr_ntve_sw(i),(gr_nbve_sw(i,k),k=1,gr_ntve_sw(i))
!        enddo
!      endif

!     Derive data on elements surrounding each node (needed for particle searching)
      gr_ntve_sw = 0
      gr_nbve_sw = 0
      do i = 1, gr_ne_sw
        do k = 1, 3
          gr_ntve_sw(gr_nodes_sw(i,k)) = gr_ntve_sw(gr_nodes_sw(i,k))+1
          gr_nbve_sw(gr_nodes_sw(i,k),gr_ntve_sw(gr_nodes_sw(i,k))) = i
        enddo
      enddo
      open(iun,file='Results/node-elements-out.dat',status='unknown')
      do i = 1, gr_nn_sw
        write(iun,'(12i6)') gr_ntve_sw(i),(gr_nbve_sw(i,k),k=1,gr_nghmax_sw)
      enddo
      close(iun)

!     Initialise time counter
      tm_ks_sw = 0

!     Load constant bathymetry into gr_depth_sw
      if (H .gt. 0) gr_depth_sw = H
!     Derive grid parameters from bathymetry
      write(*,*) 'Define grid parameters'
      call derive_grid_param

!     Initialise hydrodynamic variables
      write(*,*) 'Initialise variables'
      eta = 0.
      u = 0.
      v = 0.
      w = 0.
      tu_kv_sw = 0
      tu_kh_sw = 0
      winddata = 0.
      usw = 0.
      vsw = 0.

!     Open input data file to get parameters for particle tracking model
      iun = 20
      open(iun,file='inputs/run-input.dat',status='old')

!     *** read keyword inputs ***
      do
        read(iun,'(a)',IOSTAT=stat) nextline
        if(stat.gt.0) then
          write(*,*) 'ERROR reading file'
          ierr = 71
          call terminate(ierr,myproc)
        elseif(stat.lt.0) then
          exit
        endif
        ilen = len_trim(nextline)
!        write(*,*) 'ilen=',ilen
        if(ilen.eq.0) then !blank line
          cycle
        elseif(nextline(1:1).eq.'#') then !comment line
!          write(*,*) nextline
          cycle
        elseif(nextline(1:3).eq.'END') then !end of data
          exit
        else  !parse for keyword
          iloc1 = index(nextline,'=')
          iloc2 = index(nextline,' ')
          !          write(*,*) 'iloc=',ilen
          if(nextline(1:iloc1-1).eq.'DELTAT') then
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) p_dt
            if(stat.ne.0) then
              write(*,*) 'ERROR reading DELTAT parameter p_dt'
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' p_dt = ',p_dt
            endif
          elseif(nextline(1:iloc1-1).eq.'DT') then
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) tm_dt_sw
            if(stat.ne.0) then
              write(*,*) 'ERROR reading DT parameter tm_dt_sw'
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' tm_dt_sw = ', &
                                     tm_dt_sw
            endif
          elseif(nextline(1:iloc1-1).eq.'STARTDAY') then
!            write(*,*) 'keyword=',nextline(1:iloc-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) dd
            if(stat.ne.0) then
              write(*,*) 'ERROR reading STARTDAY parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' dd = ', dd
            endif
          elseif(nextline(1:iloc1-1).eq.'STARTTIME') then
!            write(*,*) 'keyword=',nextline(1:iloc-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) starttime
            hh = int(starttime/10000)
            mm = int((starttime - hh * 10000)/100)
            ss = int(starttime - hh * 10000 - mm * 100)
            if(stat.ne.0) then
              write(*,*) 'ERROR reading STARTTIME parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' hhmmss = ', &
                                      starttime
              write(*,*) hh,mm,ss
            endif
          elseif(nextline(1:iloc1-1).eq.'DURATION') then
!            write(*,*) 'keyword=',nextline(1:iloc1-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) lsim
            if(stat.ne.0) then
              write(*,*) 'ERROR reading DURATION parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' lsim = ',lsim
            endif
          elseif(nextline(1:iloc1-1).eq.'OUTPUTFREQ') then
!            write(*,*) 'keyword=',nextline(1:iloc1-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) fout
            if(stat.ne.0) then
              write(*,*) 'ERROR reading OUTPUTFREQ parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' fout = ',fout
            endif
          elseif(nextline(1:iloc1-1).eq.'VELOCITYDATA') then
!            write(*,*) 'keyword=',nextline(1:iloc1-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) udata
            if(stat.ne.0) then
              write(*,*) 'ERROR reading VELOCITYDATA parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' udata = ',&
                                       udata
              select case (udata)
                case ('obs')
                  utype = 0
                case ('mesh')
                  utype = 1
                case default
                  write(*,*) 'ERROR: Velocity data type incorrectly ', &
                             'specified. Exit.'
                  stop
              end select
            endif
          elseif(nextline(1:iloc1-1).eq.'ADV_SCHEME') then
!            write(*,*) 'keyword=',nextline(1:iloc1-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) adv_scheme
            if(stat.ne.0) then
              write(*,*) 'ERROR reading ADV_SCHEME parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' adv_scheme',&
                                     ' = ',adv_scheme
              select case (adv_scheme)
                case ('rk4')
                  advscheme = 0
                case ('euler')
                  advscheme = 1
                case default
                  write(*,*) 'ERROR: Advection scheme incorrectly ', &
                             'specified. Exit.'
                  stop
              end select
            endif
          elseif(nextline(1:iloc1-1).eq.'PASSIVESTAGE') then
!            write(*,*) 'keyword=',nextline(1:iloc1-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) tpassive
            if(stat.ne.0) then
              write(*,*) 'ERROR reading PASSIVESTAGE parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' tpassive = ', &
                                     tpassive
            endif
          elseif(nextline(1:iloc1-1).eq.'MOBILESTAGE') then
!            write(*,*) 'keyword=',nextline(1:iloc1-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) tmobile
            if(stat.ne.0) then
              write(*,*) 'ERROR reading MOBILESTAGE parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' tmobile = ', &
                                     tmobile
            endif
          elseif(nextline(1:iloc1-1).eq.'SWIMSPEEDUP') then
!            write(*,*) 'keyword=',nextline(1:iloc1-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) wswim
            if(stat.ne.0) then
              write(*,*) 'ERROR reading SWIMSPEEDUP parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' wswim = ', &
                                     wswim
            endif
          elseif(nextline(1:iloc1-1).eq.'SWIMSPEEDDOWN') then
!            write(*,*) 'keyword=',nextline(1:iloc1-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) wsink
            if(stat.ne.0) then
              write(*,*) 'ERROR reading SWIMSPEEDDOWN parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' wsink = ', &
                                     wsink
            endif
          elseif(nextline(1:iloc1-1).eq.'SETTLINGVELOCITY') then
!            write(*,*) 'keyword=',nextline(1:iloc1-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) wsettle
            if(stat.ne.0) then
              write(*,*) 'ERROR reading SETTLINGVELOCITY parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' wsettle = ', &
                                     wsettle
            endif
          elseif(nextline(1:iloc1-1).eq.'CELLGROWTH') then
!            write(*,*) 'keyword=',nextline(1:iloc1-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) alpha
            if(stat.ne.0) then
              write(*,*) 'ERROR reading CELLGROWTH parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' alpha = ', &
                                     alpha
            endif
          elseif(nextline(1:iloc1-1).eq.'CELLMORTALITY') then
!            write(*,*) 'keyword=',nextline(1:iloc1-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) beta
            if(stat.ne.0) then
              write(*,*) 'ERROR reading CELLMORTALITY parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' beta = ',beta
            endif
          elseif(nextline(1:iloc1-1).eq.'MORTALITYCONST') then
!            write(*,*) 'keyword=',nextline(1:iloc1-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) cmort
            if(stat.ne.0) then
              write(*,*) 'ERROR reading MORTALITYCONST parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' cmort = ', &
                                     cmort
            endif
          elseif(nextline(1:iloc1-1).eq.'MINSHEAR') then
!            write(*,*) 'keyword=',nextline(1:iloc1-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) shrmin
            if(stat.ne.0) then
              write(*,*) 'ERROR reading MINSHEAR parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' shrmin = ', &
                                     shrmin
            endif
          elseif(nextline(1:iloc1-1).eq.'MASSPERPARTICLE') then
!            write(*,*) 'keyword=',nextline(1:iloc1-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) pfac
            if(stat.ne.0) then
              write(*,*) 'ERROR reading MASSPERPARTICLE parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' pfac = ',pfac
            endif
          elseif(nextline(1:iloc1-1).eq.'CHLPERCELL') then
!            write(*,*) 'keyword=',nextline(1:iloc1-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) pconv
            if(stat.ne.0) then
              write(*,*) 'ERROR reading CHLPERCELL parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' pconv = ', &
                                     pconv
            endif
          elseif(nextline(1:iloc1-1).eq.'DEPTHLIMIT') then
!            write(*,*) 'keyword=',nextline(1:iloc1-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) hpmax
            if(stat.ne.0) then
              write(*,*) 'ERROR reading DEPTHLIMIT parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' hpmax = ', &
                                     hpmax
            endif
          elseif(nextline(1:iloc1-1).eq.'PINITIAL') then
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) tr_init_sw
            if(stat.ne.0) then
              write(*,*) 'ERROR reading PINITIAL parameter tr_init_sw'
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),  &
                            ' tr_init_sw = ', tr_init_sw
            endif
          elseif(nextline(1:iloc1-1).eq.'NSOURCE') then
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) tr_nsource_sw
            if(stat.ne.0) then
              write(*,*) 'ERROR reading NSOURCE parameter tr_nsource_sw'
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),  &
                            ' tr_nsource_sw = ',tr_nsource_sw
              if (tr_nsource_sw.gt.0) then
                allocate (tr_x0_sw(tr_nsource_sw))
                allocate (tr_y0_sw(tr_nsource_sw))
                allocate (tr_z0_sw(tr_nsource_sw))
                allocate (tr_xrange_sw(tr_nsource_sw)) 
                allocate (tr_yrange_sw(tr_nsource_sw))
                allocate (tr_zrange_sw(tr_nsource_sw))
                allocate (tr_start_sw(tr_nsource_sw))
                allocate (tr_stop_sw(tr_nsource_sw))
                allocate (tr_mass_sw(tr_nsource_sw))
!               Read particle source details
                do i = 1, tr_nsource_sw
                  read(iun,*) tr_x0_sw(i), tr_y0_sw(i), tr_z0_sw(i), &
                  tr_xrange_sw(i), tr_yrange_sw(i), tr_zrange_sw(i), &
                  tr_start_sw(i), tr_stop_sw(i), tr_mass_sw(i)
                  mtotal = mtotal + tr_mass_sw(i)
                enddo
                if (pfac .eq. 0.) pfac = mtotal / float(NPMAX)
              endif
            endif
          elseif(nextline(1:iloc1-1).eq.'VERTICALDIFF') then
!            write(*,*) 'keyword=',nextline(1:iloc1-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) kv0
            if(stat.ne.0) then
              write(*,*) 'ERROR reading VERTICALDIFF parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' kv0 = ',kv0
            endif
          elseif(nextline(1:iloc1-1).eq.'HORIZONTALDIFF') then
!            write(*,*) 'keyword=',nextline(1:iloc1-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) kh0
            if(stat.ne.0) then
              write(*,*) 'ERROR reading HORIZONTALDIFF parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' kh0 = ',kh0
            endif
          elseif(nextline(1:iloc1-1).eq.'WINDFORCING') then
!            write(*,*) 'keyword=',nextline(1:iloc1-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) wind_io
            if(stat.ne.0) then
              write(*,*) 'ERROR reading WINDFORCING parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' wind_io = ', &
                                     wind_io
            endif
          elseif(nextline(1:iloc1-1).eq.'BUOYANCY') then
!            write(*,*) 'keyword=',nextline(1:iloc1-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) buoyancy
            if(stat.ne.0) then
              write(*,*) 'ERROR reading BUOYANCY parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' buoyancy = ', &
                                     buoyancy
              if (buoyancy .ne. 0) then
                allocate (tr_volume_sw(tr_nsource_sw))
                allocate (tr_density_sw(tr_nsource_sw))
                read(iun,*) rho0
                do i = 1, tr_nsource_sw
                  read(iun,*) tr_volume_sw(i), tr_density_sw(i)
                enddo
              endif
            endif
          elseif(nextline(1:iloc1-1).eq.'RESUSPENSION') then
!            write(*,*) 'keyword=',nextline(1:iloc1-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) resus
            if(stat.ne.0) then
              write(*,*) 'ERROR reading RESUSPENSION parameter '
            else
              write(*,*) 'keyword: ',nextline(1:iloc1-1),' resus = ', &
                                     resus
            endif
          endif
        endif
      enddo
      close(iun)

!     Read wind forcing data if needed
      iun = 20
      if (wind_io .eq. 1) then
        open(iun,file='inputs/wind_velocity.dat')
        read(iun,*) nwsite
        read(iun,*) nwrec
        read(iun,*) nwinterval
        do i = 1, nwsite
          read(iun,*) xwind(i), ywind(i)
        enddo
        do i = 1, nwrec
          read(iun,*) (winddata(i,k), k = 1, nwsite*2)
        enddo
      endif
      close(iun)
      
!     Set up density fields if needed for buoyancy calculations
      if (buoyancy .ne. 0) then
        allocate (tr_prho_sw(NPMAX), tr_pvol_sw(NPMAX))
        allocate (rho(tr_jpcell_sw,gr_ne_sw))
        tr_prho_sw = 0.
        tr_pvol_sw = 0.
        rho = rho0
      endif
          
!     Number of particle tracking model time steps for each
!     data record from hydrodynamic model.
      np_steps = tm_dt_sw / p_dt
      if (mod(tm_dt_sw,p_dt) .ne. 0) then
        write(*,*) 'Particle timestep not compatible with data timestep'
        write(*,*) 'Particle time step = ',p_dt,' s, data time step = ', &
                    tm_dt_sw,' s'
        do while (mod(tm_dt_sw, p_dt) .ne. 0)
          p_dt = p_dt - 1
          np_steps = tm_dt_sw / p_dt
        enddo
        write(*,*) 'Changing particle time step to ',p_dt,' seconds'
      endif

!     Convert particle stage durations to seconds
      tpassive = tpassive * 24. * 3600.
      tmobile = tmobile * 24. * 3600.

!     Initialise particle parameters
      tr_n_sw = 0
      tr_x_sw = 0.
      tr_y_sw = 0.
      tr_z_sw = 0.
      tr_s_sw = 0.
      tr_age_sw = 0.
      tr_stage_sw = 0
      tr_pfac_sw = 0.
      tr_element_sw = 0

!     Initialise diffusion coefficients
      tu_kh_sw = kh0
      tu_kv_sw = kv0
!     Set stratification-modified diffusion if required
      if (kv0 .lt. 0.0 .and. buoyancy .ne. 0) call calc_diff
      if (kh0 .gt. 0. .and. utype .eq. 0) then
        write(*,*) 'Constant diffusion with observed velocity field specified.'
        write(*,*) '      Okubo diffusion will be applied.'
      endif
      
!     Initialise particle locations
      call p_start

!     Specify any particle sources
      call p_source

!     Get initial particle densities
      call p_density

!     Output initial particle locations
      call ptrack_out

!     Output grid coordinates
      open(iun,file='Results/grid_coordinates.dat',status='unknown')
      write(*,*) 'Writing grid coordinates'
      write(iun,*) gr_nn_sw, gr_ne_sw, gr_nlayers_sw
      do i = 1, gr_nn_sw
        write(iun,*) gr_x_sw(i), gr_y_sw(i),gr_depth_sw(i)
      enddo
      do k = 1, gr_ne_sw
        write(iun,*) gr_nodes_sw(k,1),gr_nodes_sw(k,2),gr_nodes_sw(k,3)
      enddo
      close(iun)

      end

!=============================================================================*
      subroutine unptrack_run
!=============================================================================*
!--External variables
      use unptrackmod

      implicit none

!--Description
!   Initializes the ptrack module.

!--Local variable
      real :: pstep, hour
      integer :: iper, istep, iday, idiff, ip, k, npl, npd
      integer :: nperiod, dstep, dd1
      integer :: nwread, npout
      integer :: iun
      character*12 :: ufile, vfile, tfile
!...................................

      write(6,*) ' Starting simulation'
      write(6,*) ' Time =   0'

!     Output frequency in time steps
      npout = fout / p_dt

!     Determine the number of hydrodynamic time steps are required
      nperiod = ceiling(lsim * 3600 / tm_dt_sw)

!     Note start day of simulation
      dd1 = dd

!     Set up main loop
      do iper = 1, nperiod

!       Initialise secondary loop counter to zero.
        istep = 0
        pstep = 0.0

        if (utype.eq.0) then
!         Load observed velocity data
!         Specify input file name
          iun = 22
          if (iper .eq. 1) &
            open(iun,file='inputs/flow-obs.dat',status='unknown')
          call load_obs(iun,iper,pstep)
          if (iper .eq. nperiod) close(iun)
        elseif (utype.eq.1) then
          iun = 22
          if (iper .eq. 1) &
!           open(iun,file='inputs/flow-regular.bin',status='old',
!               form='unformatted',convert='BIG_ENDIAN')
            open(iun,file='inputs/flow-mesh.dat',status='old')
          call load_mesh(iun,iper,pstep)
          if (iper .eq. nperiod) close(iun)
        endif

!       Sensitivity test: T = T + 10%
!        tem = tem * 1.05

!       Calculate horizontal and vertical eddy diffusivities
        if (kv0 .lt. 0.0 .and. buoyancy .ne. 0) call calc_diff

!       Adjust np_steps during first day, if necessary, to allow
!       for delayed start (e.g. hh = 15) on Day 1
        dstep = 0
        if (iper .eq. 1) dstep = (hh * 3600 + mm * 60 + ss) / p_dt

!       Set up secondary loop
        do istep = 1, np_steps - dstep

!         Update time step counter
          tm_ks_sw = tm_ks_sw + 1

!         Update time
          ss = ss + p_dt
          if (ss .ge. 60) then
            mm = mm + ss/60
            ss = 0
            if (mm .ge. 60) then
              hh = hh + mm/60
              mm = 0
              if (hh .ge. 24) then
                dd = dd + hh/24
                hh = 0
              end if
            end if
          end if

!         write time
          hour = (tm_ks_sw * p_dt) / 3600.
          if (mod(hour,1.0) .eq. 0) then
            npl = 0
            npd = 0
            do ip = 1, tr_n_sw
              if (tr_stage_sw(ip) .eq. 1 .or. tr_stage_sw(ip) .eq. 2) &
                                                    npl = npl + 1
              if (tr_stage_sw(ip) .eq. 3) npd = npd + 1
            enddo
            write(6,'(2x,7hTime = ,i5,29h hours; live particle count =, &
                    i7,22h; dead particle count:,i7)') int(hour),npl,npd
          endif

!         Calculate fraction of primary time step
          pstep = float(istep)/float(np_steps)

          if (utype.eq.0) then
!           Load observed velocity data
!           Specify input file name
            iun = 22
            call load_obs(iun,iper,pstep)
          elseif (utype .eq. 1) then
            iun = 22
            call load_mesh(iun,iper,pstep)
          endif

!         Calculate wind-driven velocity fields
          if (wind_io .eq. 1) then
            nwread = nwinterval * 3600 / p_dt
    	    if (mod(tm_ks_sw-1, nwread) .eq. 0) call wind_forcing
          end if

!          write(6,*) ''

!         Release particles from any sources
          do ip = 1, tr_n_sw
            idiff = tm_ks_sw - tr_igo_sw(ip)
            if (idiff .eq. 0) then
              tr_age_sw(ip) = p_dt
              tr_stage_sw(ip) = 1
              k = tr_element_sw(ip) 
              tr_s_sw(ip) = (tr_z_sw(ip)-eta(k)) / (gr_eledep_sw(k)+eta(k))
            elseif (idiff .gt. 0) then
!             Establish stage and behavioural characteristics
!             of sea lice larvae
              if (tr_stage_sw(ip) .lt. 3) call life_cycle(ip)
            endif
          enddo

!         Recalculate particle densities and tracer concentration
          call p_density

!         Get water density if dynamic tracer
          if (buoyancy .ne. 0) call wdensity

!         Move particles
          if (advscheme .eq. 0) then
            call ptrack_move_rk4
          else
            call ptrack_move_euler
          endif

!         Write results
          if (mod(tm_ks_sw,npout) .eq. 0) then
            call p_density
            call ptrack_out
          endif

!         Check for completion
          if (hour .ge. lsim) then
            write(*,*) 'Simulation completed. Stopping.'
            close(21)
            stop
          endif
          
        end do

      end do

      end

!=============================================================================*
      subroutine ptrack_out
!=============================================================================*
!--External variables
      use unptrackmod

      implicit none

!--Description
!     Outputs particle locations and "concentrations"

!--Local variable
      real :: ptime
      real :: zout
      integer :: j, k, ip
      integer :: iun
      integer :: nplive
      character*28 :: filename

!...................................

!     Output particle results
      ptime = tm_ks_sw * p_dt
      write(filename(1:15),'(a15)') 'Results/ptrack-'
      write(filename(16:24),'(i9)') int(ptime) + 100000000
      write(filename(16:16),'(a)') '0'
      write(filename(25:28),'(a4)') '.dat'
      open(53,file=filename,status='unknown')
      do ip = 1, tr_n_sw
        write(53,'(6f12.4,i5)') tr_x_sw(ip), tr_y_sw(ip), tr_s_sw(ip), &
	     tr_z_sw(ip),tr_pfac_sw(ip),tr_age_sw(ip)/3600.,tr_stage_sw(ip)
      enddo
      close(53)

!     Output particle densities in units of kg m-3
      write(filename(10:14),'(a5)') 'dnsty'
      open(53,file=filename,status='unknown')
      do k = 1, gr_ne_sw
        do j = 1, tr_jpcell_sw
          zout = float(-j) + 0.5
          if (pden(j,k) .gt. 0.0) write(53,*) &
                k, gr_xc_sw(k), gr_yc_sw(k), zout, pden(j,k)
        enddo
      enddo
      close(53)

!     Output estimated surface layer chlorophyll-a concentration
!      filename(10:14) = 'chl_a'
!      open(53,file=filename,status='unknown')
!      j = 1
!      do k = 1, gr_ne_sw
!        write(53,'(192f10.3)') (pconv * pden(j,k), k = 1, gr_ne_sw)
!        write(53,'(192f10.3)') (usw(k), k = 1, gr_ne_sw)
!      end do
!      do i = 1, gr_nn_sw
!        write(53,'(192f10.3)') (vsw(k), k = 1, gr_ne_sw)
!      end do
!      close(53)

!     Output particle mass deposited on the seabed
      write(filename(10:14),'(a5)') 'depos'
      write(filename(16:24),'(i9)') int(ptime) + 100000000
      filename(16:16) = '0'
      open(53,file=filename,status='unknown')
      write(53,*) '"x","y","ug/m2","ug/kg"'
      j = tr_jpcell_sw + 1
      do k = 1, gr_ne_sw
          write(53,'(2(f9.1,a),f13.1,a,f13.1)') gr_xc_sw(k),',', &
              gr_yc_sw(k),',',pden(j,k)*1e9,',',pden(j,k)*1e9/120.8
      enddo
      close(53)

!     output model parameter settings
      if (tm_ks_sw .le. 1) then
        iun = 21
        open(iun,file='Results/model_parameters.log',status='unknown')
        write(iun,'(a,i8)')     'N                : ', tr_n_sw
        write(iun,'(a,i8)')     'Npmax            : ', NPMAX
        write(iun,'(a,4i5)')    'NN, NE, JLO      : ', gr_nn_sw, &
                                                gr_ne_sw, gr_nlayers_sw
        write(iun,'(a,2f8.1)')  'tm_dt_sw         : ', tm_dt_sw
        write(iun,'(a,10f6.1)') 'dz               : ', gr_dz_sw
        write(iun,'(a,10f6.1)') 'dzz_sv           : ', dzz_sv
        write(iun,'(a,10f6.1)') 'z                : ', gr_z_sw
        write(iun,'(a,f8.4)')     'Time step (s)    : ', p_dt
        write(iun,'(a,i4)')     'Start day        : ', dd
        write(iun,'(a,3i4)')     'Start time       : ', hh, mm, ss
        write(iun,'(a,f12.2)')     'Run duration (h) : ', lsim
        write(iun,'(a,f8.2)')     'Output freq (s)  : ', fout
        write(iun,'(a,f10.1)')     'Passive time (h) : ', tpassive/3600.
        write(iun,'(a,f10.1)')     'Mobile time (h)  : ', tmobile/3600.
        write(iun,'(a,f8.4)')     'Swim speed (m/h) : ', wswim
        write(iun,'(a,f8.4)')     'Sink speed (m/h) : ', wsink
        write(iun,'(a,f8.2)')     'Max depth (m)    : ', hpmax
        write(iun,'(a,i8)')     'Growth on/off    : ', alpha
        write(iun,'(a,i8)')     'Mortality on/off : ', beta
        write(iun,'(a,e11.4)')     'Mortality const  : ', cmort
        write(iun,'(a,f8.4)')     'Minimum shear    : ', shrmin
        write(iun,'(a,e11.4)')     'Mass/particle   : ', pfac
        write(iun,'(a,e11.4)')     'Chlorophyll/cell : ', pconv
        write(iun,'(a,f8.4)')     'Horizontal diffusivity: ', kh0
        write(iun,'(a,f8.5)')     'Vertical diffusivity : ', kv0
        write(iun,'(a,i8)')     'Wind on/off      : ',wind_io
        write(iun,'(a)')        '  Time (h)  No. particles'
      else
        iun = 21
        nplive = 0
        do ip = 1, tr_n_sw
          if (tr_stage_sw(ip) .eq. 1 .or. tr_stage_sw(ip) .eq. 2) &
                                                    nplive = nplive + 1
        enddo
        write(iun,'(f10.4,2i12)') ptime/3600., nplive, tr_n_sw
      endif

      end

!=============================================================================*
      subroutine ptrack_move_euler
!=============================================================================*
!--External variables
      use unptrackmod
!      use IFPORT
      implicit none

!--Description
!   This subroutine advects (and possibly diffuses) the tracer
!   over the current time step.
!
!   The subroutine solves d(x,y,z) / dt = (u,v,w) in the Arakawa
!   C-grid using an Euler forward method.

!--Local variables
      integer :: idiff, ierr, ip, i, k, ihost
      integer :: ifound = 0
      real :: up, vp, wp
      real :: x_new, y_new, z_new, s_new
      real :: x_mid, y_mid, z_mid, s_mid
      real :: kx, ky, kz
      real :: rwx, rwy, rwz
      real :: z, z0, zc, zfac
      real*4 :: rnd

!--Functions called
      real :: field_value
      real :: vertmig
      logical :: isintriangle
!...................................

!     Initialise advection and diffusion coefficients
      up = 0.
      vp = 0.
      wp = 0.
      kx = 0.
      ky = 0.
      kz = 0.

!     wind forcing parameters (from Elliott, 1986)
      z0 = 0.1
      zc = 30
          
      do ip = 1, tr_n_sw
        idiff = tm_ks_sw - tr_igo_sw(ip)
        if (idiff .lt. 0 .or. tr_stage_sw(ip) .eq. 3) cycle

!       Derive local velocity field
        k = tr_element_sw(ip) 
        up = field_value(u,2,2,tr_x_sw(ip),tr_y_sw(ip),tr_s_sw(ip),k)
        vp = field_value(v,2,2,tr_x_sw(ip),tr_y_sw(ip),tr_s_sw(ip),k)
        wp = field_value(w,2,2,tr_x_sw(ip),tr_y_sw(ip),tr_s_sw(ip),k) &
                   + vertmig(ip, tr_stage_sw(ip))
                   
        if (wind_io .eq. 1) then
          z = eta(k) - tr_z_sw(ip)
          if (z .lt. z0) z = z0
          zfac = (1. - (log10(z/z0) / log10(zc/z0)))
          up = up + usw(k) * zfac
          vp = vp + usw(k) * zfac
        endif

        x_new = tr_x_sw(ip) + p_dt * up
        y_new = tr_y_sw(ip) + p_dt * vp
        s_new = tr_s_sw(ip) + p_dt * wp
!       update position due to advection
        if (resus .ne. 0 .and. tr_stage_sw(ip) .eq. 5) then
!         Call resuspension routine
!          call bedtransport (ip)
        else
          ihost = tr_element_sw(ip)
          call drycheck(x_new, y_new, s_new, ihost, ierr)
!          if (ierr .lt. 2) then
            x_mid = (tr_x_sw(ip) + x_new) / 2.
            y_mid = (tr_y_sw(ip) + y_new) / 2.
            s_mid = (tr_s_sw(ip) + s_new) / 2.
            tr_x_sw(ip) = x_new
            tr_y_sw(ip) = y_new
            tr_s_sw(ip) = s_new
            if (ierr .lt. 2) tr_element_sw(ip) = ihost
!          else
!            call drycheck(tr_x_sw(ip), tr_y_sw(ip), s_new, ihost, ierr)
!            if (ierr .lt. 2) then
!              x_mid = tr_x_sw(ip)
!              y_mid = tr_y_sw(ip)
!              s_mid = (tr_s_sw(ip) + s_new) / 2.
!              tr_s_sw(ip) = s_new
!            else
!              x_mid = tr_x_sw(ip)
!              y_mid = tr_y_sw(ip)
!              s_mid = tr_s_sw(ip)
!            endif
!          endif
        endif
        if (ierr .gt. 0) then
          tr_stage_sw(ip) = 3
          cycle
        endif

!       Add "diffusive" random walk.

!       Get local value of horizontal diffusivities
        if (kh0 .ne. 0) then
          kx = kh0
          if (kh0 .lt. 0.0) &
              kx = field_value(tu_kh_sw, 2, 2, x_mid, y_mid, s_mid, ip)
          ky = kx
!         Calculate displacement due to diffusion
          call random_number(rnd)
          rnd = 2. * rnd - 1.
          rwx = rnd * sqrt(6. * kx * p_dt)
          call random_number(rnd)
          rnd = 2. * rnd - 1.
          rwy = rnd * sqrt(6. * ky * p_dt)
!         Update position due to diffusion
          x_new = tr_x_sw(ip) + rwx
          y_new = tr_y_sw(ip) + rwy
        end if

!       Get local gradient of vertical diffusivity at current position
        if (kv0 .ne. 0) then
          kz = kv0
          if (kv0 .lt. 0.0) &
              kz = field_value(tu_kv_sw, 2, 2, x_mid, y_mid, s_mid, ip)
          k = tr_element_sw(ip) 
          z = eta(k) + gr_eledep_sw(k)
!         Calculate displacement due to diffusion
          call random_number(rnd)
          rnd = 2. * rnd - 1.
          rwz = rnd * sqrt(6. * kz * p_dt)
!         Update position due to diffusion
          s_new = tr_s_sw(ip) + (rwz / z)
        endif

!       Check updated position
        if (kh0 .ne. 0.0 .or. kv0 .ne. 0.0) then
          ihost = tr_element_sw(ip)
          call drycheck(x_new, y_new, s_new, ihost, ierr)
!        if (ierr .lt. 2) then
            tr_x_sw(ip) = x_new
            tr_y_sw(ip) = y_new
            tr_s_sw(ip) = s_new
            if (ierr .lt. 2) tr_element_sw(ip) = ihost
!          else
!            call drycheck(tr_x_sw(ip), tr_y_sw(ip), s_new, &
!                                       tr_element_sw(ip), ierr)
!            if (ierr .lt. 2) tr_s_sw(ip) = s_new
!          endif         
        endif

!       For particles on seabed or outside domain, no more movement        
        if (ierr .gt. 0) tr_stage_sw(ip) = 3

!       Update vertical position in z coordinates
        k = tr_element_sw(ip)
        tr_z_sw(ip) = (eta(k) + gr_eledep_sw(k)) * tr_s_sw(ip) + eta(k)
                          
      end do

      end

!-|Private functions

!=============================================================================*
      subroutine ptrack_move_rk4
!=============================================================================*
!--External variables

      use unptrackmod

      implicit none

!--Description
!   This subroutine advects (and possibly diffuses) the tracer
!   over the current time step.
!
!   The subroutine solves d(x,y,z) / dt = (u,v,w) in the Arakawa
!   C-grid using a fourth order Runge-Kutta method and tri-linear
!   serendipity shape functions.

!--Local variables
      integer :: idiff, ierr, ip, i, k, ihost
      integer :: ifound = 0
      real :: up, vp, wp
      real :: ke1(3), ke2(3), ke3(3), ke4(3)
      real :: x_new, y_new, z_new, s_new
      real :: x_mid, y_mid, z_mid, s_mid
      real :: kx, ky, kz
      real :: rwx, rwy, rwz
      real :: z, z0, zc, zfac
      real :: t
      real*4 :: rnd

!--Functions called
      real :: field_value
      real :: vertmig
      logical :: isintriangle
!...................................

!     Initialise advection and diffusion coefficients
      up = 0.
      vp = 0.
      wp = 0.
      kx = 0.
      ky = 0.
      kz = 0.

!     wind forcing parameters (from Elliott, 1986)
      z0 = 0.1
      zc = 30
          
      do ip = 1, tr_n_sw
        idiff = tm_ks_sw - tr_igo_sw(ip)
        if (idiff .eq. 0) then
          tr_age_sw(ip) = p_dt
          tr_stage_sw(ip) = 1
          k = tr_element_sw(ip) 
          tr_s_sw(ip) = (tr_z_sw(ip) - eta(k)) / &
                            (gr_eledep_sw(k) + eta(k))          
        elseif (idiff .gt. 0) then
          if (tr_stage_sw(ip) .lt. 3) call life_cycle(ip)
          if (tr_stage_sw(ip) .eq. 3) cycle
        else
          cycle
        endif

!       Calculate terms for RK4 advection
!       Term ke1
        k = tr_element_sw(ip)
        up = field_value(u, 2, 2, tr_x_sw(ip), tr_y_sw(ip), &
                          tr_s_sw(ip), k)
        vp = field_value(v, 2, 2, tr_x_sw(ip), tr_y_sw(ip), &
                          tr_s_sw(ip), k)
        wp = field_value(w, 1, 1, tr_x_sw(ip), tr_y_sw(ip), &
                     tr_s_sw(ip), k) + vertmig(ip, tr_stage_sw(ip))
        if (wind_io .eq. 1) then
          z = eta(k) - tr_z_sw(ip)
          if (z .lt. z0) z = z0
          zfac = (1. - (log10(z/z0) / log10(zc/z0)))
          up = up + usw(k) * zfac
          vp = vp + usw(k) * zfac
        end if
        ke1(1) = p_dt * up
        ke1(2) = p_dt * vp
        ke1(3) = p_dt * wp

!       Term ke2
        x_new = tr_x_sw(ip) + ke1(1) / 2
        y_new = tr_y_sw(ip) + ke1(2) / 2
        s_new = tr_s_sw(ip) + ke1(3) / 2
        if(isintriangle(k,x_new,y_new) .eqv. .true.) then
          ifound = 1
        else
          call QuickSearch(x_new,y_new,k,ifound)
          if (ifound .eq. 0) call FullSearch(x_new,y_new,k,ifound)
        endif
        if (ifound .eq. 1) then
          up = field_value(u, 2, 2, x_new, y_new, s_new, k)
          vp = field_value(v, 2, 2, x_new, y_new, s_new, k)
          wp = field_value(w, 1, 1, x_new, y_new, s_new, k) &
                     + vertmig(ip, tr_stage_sw(ip))
        else
          up = 0.0
          vp = 0.0
          wp = 0.0
        endif
        if (wind_io .eq. 1) then
          up = up + usw(k) * zfac
          vp = vp + usw(k) * zfac
        endif
        ke2(1) = p_dt * up
        ke2(2) = p_dt * vp
        ke2(3) = p_dt * wp

!       Term ke3
        x_new = tr_x_sw(ip) + ke2(1) / 2
        y_new = tr_y_sw(ip) + ke2(2) / 2
        s_new = tr_s_sw(ip) + ke2(3) / 2
        k = tr_element_sw(ip)
        if(isintriangle(k,x_new,y_new) .eqv. .true.) then
          ifound = 1
        else
          call QuickSearch(x_new,y_new,k,ifound)
          if (ifound .eq. 0) call FullSearch(x_new,y_new,k,ifound)
        endif
        if (ifound .eq. 1) then
          up = field_value(u, 2, 2, x_new, y_new, s_new, k)
          vp = field_value(v, 2, 2, x_new, y_new, s_new, k)
          wp = field_value(w, 1, 1, x_new, y_new, s_new, k) &
                     + vertmig(ip, tr_stage_sw(ip))
        else
          up = 0.0
          vp = 0.0
          wp = 0.0
        endif
        if (wind_io .eq. 1) then
          up = up + usw(k) * zfac
          vp = vp + usw(k) * zfac
        endif
        ke3(1) = p_dt * up
        ke3(2) = p_dt * vp
        ke3(3) = p_dt * wp

!       Term ke4
        x_new = tr_x_sw(ip) + ke3(1)
        y_new = tr_y_sw(ip) + ke3(2)
        s_new = tr_s_sw(ip) + ke3(3)
        k = tr_element_sw(ip)
        if(isintriangle(k,x_new,y_new) .eqv. .true.) then
          ifound = 1
        else
          call QuickSearch(x_new,y_new,k,ifound)
          if (ifound .eq. 0) call FullSearch(x_new,y_new,k,ifound)
        endif
        if (ifound .eq. 1) then
          up = field_value(u, 2, 2, x_new, y_new, s_new, k)
          vp = field_value(v, 2, 2, x_new, y_new, s_new, k)
          wp = field_value(w, 1, 1, x_new, y_new, s_new, k) &
                     + vertmig(ip, tr_stage_sw(ip))
        else
          up = 0.0
          vp = 0.0
          wp = 0.0
        endif
        if (wind_io .eq. 1) then
          up = up + usw(k) * zfac
          vp = vp + usw(k) * zfac
        end if
        ke4(1) = p_dt * up
        ke4(2) = p_dt * vp
        ke4(3) = p_dt * wp

!       Calculate new particle position
        x_new = tr_x_sw(ip) + (ke1(1) + 2 * ke2(1) &
                                + 2 * ke3(1) + ke4(1)) / 6
        y_new = tr_y_sw(ip) + (ke1(2) + 2 * ke2(2) &
                                + 2 * ke3(2) + ke4(2)) / 6
        s_new = tr_s_sw(ip) + (ke1(3) + 2 * ke2(3) &
                                + 2 * ke3(3) + ke4(3)) / 6

!       Check new position and update if OK
        if (resus .ne. 0 .and. tr_stage_sw(ip) .eq. 5) then
!         Call resuspension routine
!         call bedtransport (ip)
        else
          ihost = tr_element_sw(ip)
          call drycheck(x_new, y_new, s_new, ihost, ierr)
!          if (ierr .lt. 2) then
            x_mid = (tr_x_sw(ip) + x_new) / 2.
            y_mid = (tr_y_sw(ip) + y_new) / 2.
            s_mid = (tr_s_sw(ip) + s_new) / 2.
            tr_x_sw(ip) = x_new
            tr_y_sw(ip) = y_new
            tr_s_sw(ip) = s_new
            if (ierr .lt. 2) tr_element_sw(ip) = ihost
!          else
!            call drycheck(tr_x_sw(ip), tr_y_sw(ip), s_new, &
!                                           ihost, ierr)
!            if (ierr .lt. 2) then
!              x_mid = tr_x_sw(ip)
!              y_mid = tr_y_sw(ip)
!              s_mid = (tr_s_sw(ip) + s_new) / 2.
!              tr_s_sw(ip) = s_new
!            else
!              x_mid = tr_x_sw(ip)
!              y_mid = tr_y_sw(ip)
!              s_mid = tr_s_sw(ip)
!           endif
!            tr_element_sw(ip) = ihost
!          endif
        endif
        if (ierr .gt. 0) then
          tr_stage_sw(ip) = 3
          cycle
        endif

!       Add "diffusive" random walk.

!       Get local value of horizontal diffusivities
        if (kh0 .ne. 0) then
          kx = kh0
          if (kh0 .lt. 0.) then
            kx = field_value(tu_kh_sw,2,2,x_mid,y_mid,s_mid,ip)
          elseif (kh0 .gt. 0. .and. utype .eq. 0) then
!           Apply Okubo (1971) formula for ocean diffusion
            t = tr_age_sw(ip)
            kx = max(2.7e-7 * t**1.34, 0.1)
          endif
          ky = kx
          call random_number(rnd)
          rnd = 2. * rnd - 1.
          rwx = rnd * sqrt(6. * kx * p_dt)
          call random_number(rnd)
          rnd = 2. * rnd - 1.
          rwy = rnd * sqrt(6. * ky * p_dt)
          x_new = tr_x_sw(ip) + rwx
          y_new = tr_y_sw(ip) + rwy
        endif

!       Get local gradient of vertical diffusivity at current position
        if (kv0 .ne. 0) then
          kz = kv0
!         Constant vertical diffusivity gradient
          if (kv0 .lt. 0.) &
                  kz = field_value(tu_kv_sw,2,2,x_mid,y_mid,s_mid,ip)
          k = tr_element_sw(ip) 
          z = eta(k) + gr_eledep_sw(k)
          call random_number(rnd)
          rnd = 2. * rnd - 1.
          rwz = rnd * sqrt(6. * kz * p_dt)
          s_new = tr_s_sw(ip) + (rwz / z)
        endif

!       update position due to diffusion
        if (kh0 .ne. 0.0 .or. kv0 .ne. 0.0) then
          ihost = tr_element_sw(ip)
          call drycheck(x_new, y_new, s_new, ihost, ierr)
          if (ierr .lt. 2) then
            tr_x_sw(ip) = x_new
            tr_y_sw(ip) = y_new
            tr_s_sw(ip) = s_new
          else
            call drycheck(tr_x_sw(ip), tr_y_sw(ip), s_new, &
                                           ihost, ierr)
            if (ierr .lt. 2) tr_s_sw(ip) = s_new
          endif
          tr_element_sw(ip) = ihost
        endif
            
!       For particles on seabed, no more movement        
        if (ierr .gt. 0) tr_stage_sw(ip) = 3

!       Update vertical position in z coordinates
        k = tr_element_sw(ip)
        tr_z_sw(ip) = (eta(k)+gr_eledep_sw(k))*tr_s_sw(ip) + eta(k)

      enddo

      end

!-|Private functions

!=============================================================================*
      real function field_value(field,hor_index,ver_index,xp,yp,zp,k)
                                
!=============================================================================*
!-- Description
!
!   This subroutine interpolates a field value at a given point
!   (xp,yp,zp) where zp is in meters,
!   negative down and zero at the surface.
!   This version is for an unstructured mesh, with horizontal 
!   velocities given at mid-layer depth at the nodes, and the 
!   vertical velocity given at the cell centre at sigma levels.
!   The field value is calculated using a barycentric interpolation.
!   Note that in this version, zp is in sigma coordinates.
!
!   Meaning of index:
!      hor_index = 1 : for w, sigmat
!      hor_index = 2 : for u, v
!
!      ver_index = 1 : for w
!      ver_index = 2 : for sigmat, u, v

!--External variables
      use unptrackmod

      implicit none

!--Input arguments
      integer, intent(in) :: hor_index
      integer, intent(in) :: ver_index
      integer, intent(in) :: k
      real, intent(in)    :: xp, yp, zp
      real, intent(in) :: field(gr_nlayers_sw,gr_nn_sw)

!--Local variables
      real :: A1 = 0
      real :: A2 = 0
      real :: A3 = 0
      real :: f1(2) = 0
      real :: f2(2) = 0
      real :: f3(2) = 0
      real :: xn(3), yn(3), Un(3)
      real :: dsig = 0
      integer :: j, i, jm1

!...................................

      field_value = 0.0
      Un = 0.0

!     - Get the coordinates of the nodes making up element k
      do i = 1, 3
        xn(i) = gr_x_sw(gr_nodes_sw(k,i))
        yn(i) = gr_y_sw(gr_nodes_sw(k,i))
      enddo

!     - First check if below bottom
      if (zp .le. -1.0) return
      if (gr_eledep_sw(k) .eq. 0) return
      
!     - Set the horizontal and vertical weighting for v-point.
      if (hor_index .eq. 2) then
        f1(1) = xn(1) - xp
        f1(2) = yn(1) - yp
        f2(1) = xn(2) - xp
        f2(2) = yn(2) - yp
        f3(1) = xn(3) - xp
        f3(2) = yn(3) - yp
        
        A1 = abs(f2(1) * f3(2) - f2(2) * f3(1))
        A2 = abs(f3(1) * f1(2) - f3(2) * f1(1)) 
        A3 = abs(f1(1) * f2(2) - f1(2) * f2(1))
      endif
      
!     - Interpolate the node velocities to the same sigma depth as zp
      if (gr_nlayers_sw .eq. 1) then
        do i = 1, 3
          Un(i) = field(1,gr_nodes_sw(k,i))
        enddo
      elseif (ver_index .eq. 2) then
        if (zp .gt. 0.0 .or. zp .lt. -1.0) then
                       write(*,*) 'ERROR: particle not in water column'
                       write(*,*) xp, yp, zp
                       stop
        endif
        do i = 1, 3
          if (zp .ge. gr_z_sw(1)) then
            Un(i) = field(1,gr_nodes_sw(k,i))
          elseif (zp .le. gr_z_sw(gr_nlayers_sw)) then
            Un(i) = field(gr_nlayers_sw,gr_nodes_sw(k,i))
          else
            j = 1
            do while (gr_z_sw(j) .gt. zp) 
              j = j + 1
            enddo
            jm1 = j - 1
            dsig = (zp - gr_z_sw(jm1)) / gr_dz_sw(jm1)
            Un(i) = (1 - dsig) * field(jm1,gr_nodes_sw(k,i)) + &
                         dsig * field(j,gr_nodes_sw(k,i))
          endif
        enddo
      endif
      
      field_value = (A1*Un(1)+A2*Un(2)+A3*Un(3)) / (2.*gr_area_sw(k))
      end

!=============================================================================*
      real function vertmig(ip, istage)
!=============================================================================*
!--External variables
      use unptrackmod

      implicit none

!--Input arguments
      integer :: ip, istage

!--Description
!  Can be used to add to tracer some vertical motion with respect
!  to ambient water.

!--Local variables
      real :: sal, rho_amb
      real :: speed_up, speed_dn
      real :: z, wb, h
      integer :: i, j, k, kn, nn
      integer :: iney, nele
!...................................

      vertmig = 0.
      
!     k is element number containing particle
      k = tr_element_sw(ip)
      h = 1. / (gr_eledep_sw(k) + eta(k))

!     z is depth of particle in m
      z = (tr_s_sw(ip) / h) + eta(k) 

      if (istage .eq. 2) then

!       specify vertical migration speeds (m/s)
        speed_dn = -wsink / 3600.
        speed_up =  wswim / 3600.

!       Convert swimming speeds into sigma space
        speed_dn = speed_dn * h
        speed_up = speed_up * h
        
!       Apply sinking speed at night, and
!       upward swimming speed during daylight
        vertmig = speed_dn
        if (hh .ge. 6 .and. hh .le. 18) vertmig = speed_up

!       depth avoidance
        if (z .lt. -hpmax) vertmig = speed_up
        
      end if

!     Add settling velocity
      vertmig = vertmig - (tr_wsettle_sw(ip) * h)
      
!     If buoyancy is included in simulations, calculate vertical buoyancy 
      if (buoyancy .ne. 0) then
        j = int(-z) + 1
        wb = 0.0
        if (j .gt. 1 .and. j .le. tr_jpcell_sw) then
          rho_amb = 0.0
          nele = 0
          do i=1, 3
            nn = gr_nodes_sw(k,i)
            do kn=1, gr_ntve_sw(nn)
              iney = gr_nbve_sw(nn,kn) 
              if (iney .ne. k) then
                rho_amb = rho_amb + rho(j,iney)
                nele = nele + 1
              endif
            enddo
          enddo
          if (nele .gt. 0) then
            rho_amb = rho_amb / float(nele)
            wb = (rho_amb - rho(j,k)) * 9.81 * p_dt / &
                           (2.0 * rho_amb)
!            write(*,*) 'Vertmig: ',ip, j, k, nele, rho(j,k), rho_amb,wb &
!                 ,gr_area_sw(k)
            if (wb .lt. 0.0) wb = 0.0
!           Convert to sigma space
            wb = wb * h
          endif
        endif
        vertmig = vertmig + wb
      endif
      
      end
      
!=============================================================================*
      subroutine drycheck(x, y, z, ihost, ierr)
!=============================================================================*
!--External variables
      use unptrackmod

      implicit none

!--Input arguments
      real, intent(inout) :: x, y, z
      integer, intent(inout) :: ihost

!--Output argument
      integer, intent(out) :: ierr

!--Description
!  Checks for status of particle position (x,y,z) and returns an "error" code.
!  ierr = 0 : particle remains within the water column.
!  ierr = 1 : particle is settled on the seabed.
!  ierr = 2 : particle is grounded on dry land.
!  ierr = 3 : particle cannot be found on the grid - assumed exported from the domain.

!--Local variables
      integer :: i,k
      integer :: ifound
      logical :: isintriangle      
!...................................

!     Initialise.
      ierr = 0
      ifound = 0

!     Find element containing (x,y) coordinate
      if(isintriangle(ihost,x,y) .eqv. .false.) then
        call QuickSearch(x,y,ihost,ifound)
        if (ifound .eq. 0) call FullSearch(x,y,ihost,ifound)
        if (ifound .eq. 0) then
          ierr = 3
          return
        endif
      endif
      k = ihost

      if (gr_eledep_sw(k) .le. 0) then
        ierr = 2
        return
      end if
      
      if (wsettle .ne. 0. .and. z .le. -1.0) then
        z = -1.0
        if (resus .eq. 0) ierr = 1
      elseif (z .le. -1.0) then
        z = -(2. + z)
      endif
      
      if (z .gt. 0.0) z = 0.0

      end

!=============================================================================*
      subroutine life_cycle(ip)
!=============================================================================*
!--External variables
      use unptrackmod

      implicit none

!--Input arguments
      integer :: ip
      
!--Description
!   .. ==>

!--Functions called
      real :: field_value
!      real :: field_grad

!--Local variables
      integer :: i, j, k
      real(4) :: rnd
      real :: t, mu, gamma
      real :: p_growth, p_mortality
      real :: ushr, vshr, shr
      real :: x, y, z

!     particle location
      x = tr_x_sw(ip)
      y = tr_y_sw(ip)
      z = tr_z_sw(ip)

!     Update particle age
      tr_age_sw(ip) = tr_age_sw(ip) + p_dt

!     check age and stage of particle
!     1. change from passive to mobile stage
      if (tr_stage_sw(ip) .eq. 1 .and. tr_age_sw(ip) .gt. tpassive) &
                      tr_stage_sw(ip) = 2
!     2. mortality after n days
      if (tr_age_sw(ip) .gt. tpassive+tmobile) tr_stage_sw(ip) = 3

!     Allow for mortality of cells
      gamma = 0.
      if (beta .eq. 1) then
        k = tr_element_sw(ip) 
        j = int(-z) + 1
        if (z .ge. 0.) j = 1
!        ushr = field_grad(u, 2, x, y, z)
!        vshr = field_grad(v, 2, x, y, z)
        shr = sqrt(ushr * ushr + vshr * vshr)
        shr = amax1(shr, shrmin)
        gamma = cmort * shr * (pden(j,k) / 1000.)
!
!-------Do mortality by removing particles.------------------------
!
!	p_mortality = 1. - exp(-gamma * p_dt)
!        rnd = ran
! 	if (rnd .lt. p_mortality) then
!	  tr_stage_sw(ip) = 3
!	  return
!	end if
!------------------------------------------------------------------
      end if

!     Allow for particle growth if switched on.
!     Growth equation from Gentien et al (2007)
      mu = 0.
      if (alpha .eq. 1) then
        t = field_value(tem, 1, 2, x, y, z, ip)
!       calculate temperature-dependent growth in divisions per day
        mu = (0.0025*t*t*t - 0.15*t*t + 2.8775*t - 17.25)
!       set minimum growth rate (see Vanhoutte-Brunier et al, 2008)
        mu = amax1(mu, 0.05)
!       convert units to divisions per second
        mu = mu / 86400

!------ Do growth by splitting cells and creating new particles.-------
!(float(ihr)) * 3600 / p_dt + 1
!        p_growth = 1. - exp(-mu * p_dt)
!        rnd = ran
!        if (rnd .lt. p_growth) then
!          i = tr_n_sw + 1
!          write(6,*) i
!          do while (i .le. NPMAX)
!            if (tr_stage_sw(i) .eq. 0 .or. tr_stage_sw(i) .eq. 3) exit
!              i = i + 1
!          end do
!          if (i .le. NPMAX) then
!            tr_x_sw(i) = tr_x_sw(ip)
!            tr_y_sw(i) = tr_y_sw(ip)
!            tr_z_sw(i) = tr_z_sw(ip)
!            tr_stage_sw(i) = 1
!            tr_age_sw(i) = 0
!            tr_igo_sw(i) = tm_ks_sw + 1
!            tr_pfac_sw(ip) = pfac
!            tr_n_sw = i
!          end if
!        end if
!
      end if

!
!----- Do growth and mortality by modifying particle properties.----------
!
      tr_pfac_sw(ip) = tr_pfac_sw(ip) * (1.0 + p_dt * (mu - gamma))
!
!-------------------------------------------------------------------------
!
      end

!=============================================================================*
      subroutine p_start
!=============================================================================*
!--External variables

      use unptrackmod

      implicit none
      
!--Description
!   Sets up 2D field of initial particle locations

!--Local variables
      real :: x, y, z, pcon

      integer :: i, j, k
      integer :: ip, np, npoints
      integer :: iun
      integer :: ifound

!...................................

!     Initialise particle counter and other variables
      ip = 0
      pcon = 0.
      np = 0
      npoints = 0
      iun = 20
      ifound = 0
      
!     Either assign initial value or load initial concentrations from file
      if (tr_init_sw .gt. 0.0) then
        do k = 1, gr_ne_sw
          do j = 1, gr_nlayers_sw
            do np = 1, tr_init_sw
              ip = ip + 1
              tr_x_sw(ip) = (gr_x_sw(gr_nodes_sw(k,1)) + &
                             gr_x_sw(gr_nodes_sw(k,2)) + &
                             gr_x_sw(gr_nodes_sw(k,3)))/3.
              tr_y_sw(ip) = (gr_y_sw(gr_nodes_sw(k,1)) + &
                             gr_y_sw(gr_nodes_sw(k,2)) + &
                             gr_y_sw(gr_nodes_sw(k,3)))/3.
              tr_z_sw(ip) = gr_z_sw(j)
              tr_age_sw(ip) = 0.
              tr_stage_sw(ip) = 0
              tr_igo_sw(ip) = 1
              tr_pfac_sw(ip) = pfac
              tr_element_sw(ip) = k
              tr_s_sw(ip) = tr_z_sw(ip) / gr_eledep_sw(k)
            enddo
          enddo
        enddo
      elseif (tr_init_sw .lt. 0.0) then
!       Open and read input file
        open(iun,file='inputs/initial_positions.dat',status = 'old')
!       read initial concentrations and convert to particles
        ip = 0
        read(iun,*) npoints
        do j = 1, npoints
          read (iun,*) x, y, z, pcon
          np = int(pcon / pfac)
          do i = 1, np
            ip = ip + 1
            tr_x_sw(ip) = x
            tr_y_sw(ip) = y
            tr_z_sw(ip) = z
            tr_age_sw(ip) = 0.
            tr_stage_sw(ip) = 0
            tr_igo_sw(ip) = 1
            tr_pfac_sw(ip) = pfac
            call FullSearch(x,y,tr_element_sw(ip),ifound)
            if (ifound .eq. 0) then
              write(*,*) 'Cannot find particle location in the grid.', &
                         'Stopping.'
              stop
            endif
            tr_s_sw(ip) = tr_z_sw(ip) / gr_eledep_sw(tr_element_sw(ip))
          enddo
        enddo
      endif

      tr_n_sw = ip
      if (tr_n_sw .gt. NPMAX) then
        write(6,'(2x,1h )')
        write(6,'(2x,31hERROR: Too many particles (N = ,i7)') tr_n_sw
        write(6,'(2x,43h       Decrease the value of variable pfac )')
        write(6,'(2x,37h       or increase the size of NPMAX.)')
        write(6,'(2x,1h )')
        STOP
      else
        write(6,'(40h Initial Number of particles released = ,i7)') ip
      endif

      close(iun)

      end

!=============================================================================*
      subroutine p_source
!=============================================================================*
!--External variables

      use unptrackmod

      implicit none
      
!--Description
!   Sets up 2D field of initial particle locations

!--Local variables
      integer :: i, k, ip, ippr
      integer :: nppr = 0
      integer :: elementno, ifound
      real*4 :: rnd
      real :: Trelease
      real :: tvol = 0
      real*4 :: x1, x2, y1
      real :: diam(NPMAX)
      
      real, parameter :: pi = 3.1415927
      real, parameter :: g = 9.81
      real, parameter :: kv = 1.83e-6
      real, parameter :: sixth = 1./6.
      real, parameter :: sigma = 1./3.

!--Functions called

!...................................

!     initialise random number generator
      call random_seed()

!     Initialise variables
      ip = tr_n_sw
      ifound = 0
      elementno = 0
      diam = 0.
      
!     Number of particles released per release (based on NPMAX)
!      if (tr_nsource_sw .gt. 0) &
!                nppr = int((NPMAX - tr_n_sw) / tr_nsource_sw)
      
!     Loop through releases
      do i = 1, tr_nsource_sw

!       Calculate number of particles for this release
        nppr = max(int(tr_mass_sw(i) / pfac), 1)

!       Duration (in hours) of this release
        Trelease = tr_stop_sw(i) - tr_start_sw(i)

!       Find element number of source
        call FullSearch(tr_x0_sw(i),tr_y0_sw(i),elementno,ifound)
        if (ifound .eq. 0) then
          write(*,*) 'Cannot find particle source on the grid.', &
                     'Stopping.'
          stop
        else
          write(*,*) 'Source ',i,': Element Number ',elementno
        endif

!       Define release positions including random component over given radius        
        do ippr = 1, nppr
            ip = ip + 1
            call random_number(rnd)
            rnd = rnd * Trelease
            tr_igo_sw(ip) = (tr_start_sw(i) + rnd) * 3600 / p_dt + 1
            call random_number(rnd)
            rnd = (rnd - 0.5) * 2.0
            tr_x_sw(ip) = tr_x0_sw(i) + tr_xrange_sw(i) * rnd
            call random_number(rnd)
            rnd = (rnd - 0.5) * 2.0
            tr_y_sw(ip) = tr_y0_sw(i) + tr_yrange_sw(i) * rnd
            call random_number(rnd)
            rnd = (rnd - 0.5) * 2.0
            tr_z_sw(ip) = tr_z0_sw(i) + tr_zrange_sw(i) * rnd
            tr_pfac_sw(ip) = tr_mass_sw(i) / float(nppr)
            if (buoyancy .ne. 0) then
              tr_pvol_sw(ip) = tr_volume_sw(i) / float(nppr)
              tr_prho_sw(ip) = tr_density_sw(i)
            endif
            tr_age_sw(ip) = 0.
            tr_stage_sw(ip) = 0
            tr_element_sw(ip) = elementno
            tr_s_sw(ip) = tr_z_sw(ip) / gr_eledep_sw(tr_element_sw(ip))
          enddo
      enddo

      tr_n_sw = ip
      write(*,*) 'No. of particles to be released = ', tr_n_sw

!     Add random 10% component to settling velocity if used
      if (wsettle .gt. 0.) then
        do ip = 1, tr_n_sw
          call random_number(rnd)
          tr_wsettle_sw(ip) = (0.9 + (2. * rnd / 10.)) * wsettle
        enddo
      elseif (wsettle .lt. 0.) then
!        Apply normal distribution of settling velocity and mass to particles.
!        This distribution specifically for a particulate tracer release 
!        experiment by Partrac during May 2017.
         ip = 0
         do i = 1, tr_nsource_sw
           do while (ip .lt. i*nppr)
             call random_number(rnd)
             x1 = 2. * rnd - 1.
             call random_number(rnd)
             x2 = rnd
             y1 = exp(-x1*x1/(2.0*sigma*sigma))
             if (x2 .lt. y1) then
               ip = ip + 1
               diam(ip) = 10.**(2.34 + x1 * 2.75 * 0.163)
               diam(ip) = diam(ip) * 1e-6
               tr_wsettle_sw(ip) = (10.36 * kv / diam(ip)) * (sqrt(1 + &
                   0.156 * (3000./1025. - 1) * 9.81 * diam(ip)**3 / &
                   (16 * kv * kv)) - 1)
               tr_pfac_sw(ip) = sixth * pi * diam(ip)**3
               tvol = tvol + tr_pfac_sw(ip)
             endif
           enddo
           if (tvol .gt. 0.) tr_pfac_sw = tr_pfac_sw * tr_mass_sw(i) / &
                                           tvol
         enddo
      endif

!     Write settling velocities to filename
      open(20,file='Results/ParticleCharacteristics.dat',status='unknown')
      do ip = 1, tr_n_sw
        write(20,'(i7,3f10.6)') ip, tr_wsettle_sw(ip), &
                        diam(ip), tr_pfac_sw(ip)
      enddo
      close(20)
      
      end

!=============================================================================*
      subroutine derive_grid_param
!=============================================================================*
!--External variables
      use unptrackmod

      implicit none
      
!--Description
!   Derives grid parameters from bathymetry and layer thicknesses

!--Local variables
      integer :: i, j, k
      real :: xn(3), yn(3)
      real :: f1(2), f2(2)
      real :: areamin, areamax, areamean

!...................................
      if (gr_nlayers_sw .eq. 1) then
        gr_z_sw(1) = -0.5
        dzz_sv(1) = 0.0
        dzz_sv(2) = -1.0
      endif
      
!     Derive layer thicknesses
      gr_dz_sw(1) = -dzz_sv(2)
      do j = 2, gr_nlayers_sw
        gr_dz_sw(j) = dzz_sv(j) - dzz_sv(j+1)
      end do
      
!     Calculate water depth at each element
      do k = 1, gr_ne_sw
        gr_eledep_sw(k) = 0.
        do i = 1, 3
          gr_eledep_sw(k) = gr_eledep_sw(k) + &
                      gr_depth_sw(gr_nodes_sw(k,i))
        enddo
        gr_eledep_sw(k) = gr_eledep_sw(k) / 3.
      enddo

!     Get coordinates of element centroids
      areamean = 0.0
      areamax = 0.0
      areamin = 100000.0
      do k = 1, gr_ne_sw
        do i = 1, 3
          xn(i) = gr_x_sw(gr_nodes_sw(k,i))
          yn(i) = gr_y_sw(gr_nodes_sw(k,i))
        enddo
        gr_xc_sw(k) = (xn(1) + xn(2) + xn(3)) / 3.
        gr_yc_sw(k) = (yn(1) + yn(2) + yn(3)) / 3.
      
!       Calculate surface area of elements
        f1(1) = xn(1) - xn(2)
        f1(2) = yn(1) - yn(2)
        f2(1) = xn(1) - xn(3)
        f2(2) = yn(1) - yn(3)
        gr_area_sw(k) = abs(f1(1) * f2(2) - f1(2) * f2(1)) / 2.0
        areamean = areamean + gr_area_sw(k)
        if (gr_area_sw(k) .gt. areamax) areamax = gr_area_sw(k)
        if (gr_area_sw(k) .lt. areamin) areamin = gr_area_sw(k)
      enddo
      areamean = areamean / gr_ne_sw
      
      write(*,*) 'Mean element area = ',areamean,' m2'
      write(*,*) 'Max element area = ',areamax,' m2'
      write(*,*) 'Min element area = ',areamin,' m2'

      end

!=============================================================================*
      subroutine calc_diff
!=============================================================================*
!--External variables
      use unptrackmod

      implicit none
      
!--Description
!   Sets up 3D eddy diffusion fields: tu_kh_sw, tu_kv_sw

!--Local variables
      integer :: i, j, k
      real :: N, N2

!...................................

!     Include damping of vertical diffusion by stratification.
!     Use simple algorithm based on Gargett & Holloway (1984), where
!     Kv = kv0. N^-1
!     where N is the buoyancy frequency.

      do i = 1, gr_nn_sw
        do j = 2, gr_nlayers_sw
          N2 = 9.81 * (rho(j,i) - rho(j-1,i)) / rho0
          if (N2 .le. 0.0) N2 = 1.0e-4
          N = 1. / sqrt(N2)
          tu_kv_sw(j,i) = -kv0 * N
!          if (tu_kv_sw(j,i) .lt. 0.0009) write(*,*) j,rho(j-1,i), &
!          rho(j,i), tu_kv_sw(j,i)
        enddo
      enddo
      
      end

!=============================================================================*
      subroutine p_density
!=============================================================================*
!--External variables
      use unptrackmod

      implicit none
      
!--Description
!   Calculates particle density throughout the model domain.
!   Units are cell numbers per litre.

!--Local variable
      integer :: i, j, k, ip
!...................................

!     Initialise density and other variables to zero
      pden = 0.
      psum = 0
      j = 0

!     Loop through particles, summing each grid cell over 1m thick layers
      do ip = 1, tr_n_sw
        k = tr_element_sw(ip)
        if (tr_stage_sw(ip) .gt. 0 .and. tr_stage_sw(ip) .lt. 3) then
          j = int(-tr_z_sw(ip)) + 1
          if (j .le. tr_jpcell_sw) psum(j,k) = psum(j,k) + & 
                                               tr_pfac_sw(ip)
        endif
        
!       Record particles settled on the seabed
        j = tr_jpcell_sw + 1
        if (tr_z_sw(ip) .le. -gr_eledep_sw(k)) psum(j,k) = &
                    psum(j,k) + tr_pfac_sw(ip)
      enddo

!     Concentration in mass per cubic metre (assuming 1m thick layers)
      do k = 1, gr_ne_sw
        do j = 1, tr_jpcell_sw
          if (psum(j,k) .gt. 0) pden(j,k) = psum(j,k) / &
                                    gr_area_sw(k)
        enddo
          
!       Record settled particles in mass per unit area          
        j = tr_jpcell_sw + 1
        pden(j,k) = psum(j,k) / gr_area_sw(k)
      enddo

      end

!=============================================================================*
      subroutine wdensity
!=============================================================================*
!--External variables
      use unptrackmod

      implicit none
      
!--Description
!   Calculates water density throughout the model domain in cases for dynamic
!   particles i.e. where particles have a density and influence water density.

!--Local variable
      integer :: i, j, k, ip
      integer :: jmin, kmin
      real :: rhomin, prhomin
!...................................

!     Initialise particle counter to zero
      psum = 0.0

!     Loop through particles, summing each grid cell
      do ip = 1, tr_n_sw
        if (tr_stage_sw(ip) .gt. 0 .and. tr_stage_sw(ip) .lt. 3) then
          k = tr_element_sw(ip)
          j = int(-tr_z_sw(ip)) + 1
          if (tr_z_sw(ip) .ge. 0.) j = 1
          psum(j,k) = psum(j,k) + tr_pvol_sw(ip) * &
                          (rho0 - tr_prho_sw(ip))
        endif
      enddo

!     Calculate water density in cell
      rhomin = rho0
      do k = 1, gr_ne_sw
        do j = 1, tr_jpcell_sw
          rho(j,k) = rho0 - (psum(j,k) / gr_area_sw(k))
!          if (rho(j,k) .lt. rho0) &
!           write(*,*) 'Water density', j, k, rho(j,k)
          if (rho(j,k) .lt. rhomin) then
            rhomin = rho(j,k) 
            jmin = j
            kmin = k
          endif
        enddo
      enddo
      
!     Check for sensible outcome
      if (rhomin .lt. minval(tr_density_sw)) then 
        write(*,*) 'ERROR: Unrealistic contaminant source.'
        write(*,*) 'Time step ', tm_ks_sw,', Density = ', rhomin
        write(*,*) 'j,k = ', jmin, kmin
      endif
         
      end

!=============================================================================*
      subroutine wind_forcing
!=============================================================================*
!--External variables
      use unptrackmod

      implicit none
      
!--Description
!   Calculates the wind-driven velocity at the particle location

!--Local variable
      integer :: i, j, k
      integer :: isite
      integer :: irec
      integer :: iuw, ivw
      real :: dist, dist2

!...................................

!     identify wind record number for current time step
      irec = int(((dd - 1) * 24 + hh) / nwinterval + 1)

!     Work through grid, finding the closest wind data site

      do i = 1, gr_nn_sw
!         find the closest wind data site
          isite = 1
          dist = (gr_x_sw(i) - xwind(isite))**2 + &
                                 (gr_y_sw(i) - ywind(isite))**2
          do j = 2, nwsite
            dist2 = (gr_x_sw(i) - xwind(j))**2 + &
                                     (gr_y_sw(i) - ywind(j))**2
            if (dist2 .lt. dist) then
              dist = dist2
              isite = j
            endif
          enddo

!         calculate wind-driven surface velocities at each grid cell
          iuw = (isite * 2) - 1
          ivw = isite * 2
          usw(i) = 0.5 * sqrt(1.2 / 1025.) * winddata(irec,iuw)
          vsw(i) = 0.5 * sqrt(1.2 / 1025.) * winddata(irec,ivw)
      enddo

!      write(6,'(a,2f10.2)') 'Wind velocity (m/s) : ',winddata(irec,iuw), &
!                                                     winddata(irec,ivw)

      end

!=============================================================================*
      subroutine load_obs(iun,iper,pstep)
!=============================================================================*
!--External variables
      use unptrackmod

      implicit none

!--Input arguments
      integer :: iun, iper
      real :: pstep
      
!--Description
!   Loads observed velocity data. Format includes a header, and comma-separated
!   data.

      integer :: i, j, k
      integer :: irec
      integer :: stat
      integer :: iloc1, iloc2, ilen
      integer :: ierr, myproc
      real :: uin(3), vin(3)
      real :: h_obs
      real, allocatable, save :: z_obs(:)
      real, allocatable, save :: u_obs1(:), v_obs1(:)
      real, allocatable, save :: u_obs2(:), v_obs2(:)
      character(80) :: nextline

      if (iper.eq.1 .and. pstep .eq. 0.) then
        allocate (u_obs1(gr_nlayers_sw), u_obs2(gr_nlayers_sw), &
                  v_obs1(gr_nlayers_sw), v_obs2(gr_nlayers_sw))
        allocate (z_obs(3))
        do
          read(iun,'(a)',IOSTAT=stat) nextline
          iloc1 = index(nextline,'=')
          iloc2 = index(nextline,' ')
          if(stat.gt.0) then
            write(*,*) 'ERROR reading file'
            ierr = 71
            call terminate(ierr,myproc)
          elseif(stat.lt.0) then
            exit
          elseif(nextline(1:iloc1-1).eq.'Flowmetry.siteDepth') then
!            write(*,*) 'keyword=',nextline(1:iloc1-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) h_obs
            h_obs = h_obs
            if(stat.ne.0) then
              write(*,*) 'ERROR reading Flowmetry.siteDepth parameter'
              stop
            endif
          elseif(nextline(1:iloc1-1).eq.'Flowmetry.meterDepths') then
!            write(*,*) 'keyword=',nextline(1:iloc1-1)
            read(nextline(iloc1+1:iloc2-1),*,IOSTAT=stat) (z_obs(j), &
                                          j=1,3)
            z_obs = -z_obs / h_obs
            if(stat.ne.0) then
              write(*,*) 'ERROR reading Flowmetry.meterDepths parameter'
              stop
            endif
          endif
          ilen = len_trim(nextline)
!         write(*,*) 'ilen=',ilen
          if(ilen.eq.0) then !blank line
            cycle
          elseif(nextline(1:6).eq.'#start') then !start of data
            read(iun,*) irec,uin(1),vin(1),uin(2),vin(2),uin(3),vin(3)
            do j = 1, gr_nlayers_sw
              if (gr_z_sw(j) .ge. z_obs(1)) then
                    u_obs2(j) = uin(1)
                    v_obs2(j) = vin(1)
              elseif (gr_z_sw(j) .le. z_obs(3)) then
                    u_obs2(j) = uin(3)
                    v_obs2(j) = vin(3)
              elseif (gr_z_sw(j) .lt. z_obs(1) .and. &
                      gr_z_sw(j) .ge. z_obs(2)) then
                    u_obs2(j) = uin(1) + (uin(2) - uin(1)) * &
                        (gr_z_sw(j) - z_obs(1)) / (z_obs(2) - z_obs(1))
                    v_obs2(j) = vin(1) + (vin(2) - vin(1)) * &
                        (gr_z_sw(j) - z_obs(1)) / (z_obs(2) - z_obs(1))
              elseif (gr_z_sw(j) .lt. z_obs(2) .and. & 
                      gr_z_sw(j) .ge. z_obs(3)) then
                    u_obs2(j) = uin(2) + (uin(3) - uin(2)) * &
                        (gr_z_sw(j) - z_obs(2)) / (z_obs(3) - z_obs(2)) 
                    v_obs2(j) = vin(2) + (vin(3) - vin(2)) * &
                        (gr_z_sw(j) - z_obs(2)) / (z_obs(3) - z_obs(2))
              endif
            enddo
            exit
          else
            cycle
          endif
        enddo
      endif

      if (pstep .eq. 0.) then
        u_obs1 = u_obs2
        v_obs1 = v_obs2
        read(iun,*) irec, uin(1), vin(1), uin(2), vin(2),uin(3), vin(3)
!        write(*,*) irec, uin(1), vin(1), uin(2), vin(2),uin(3), vin(3)
        do j = 1, gr_nlayers_sw
          if (gr_z_sw(j) .ge. z_obs(1)) then
            u_obs2(j) = uin(1)
            v_obs2(j) = vin(1)
          elseif (gr_z_sw(j) .le. z_obs(3)) then
            u_obs2(j) = uin(3)
            v_obs2(j) = vin(3)
          elseif (gr_z_sw(j) .lt. z_obs(1) .and. &
                  gr_z_sw(j) .ge. z_obs(2)) then
                u_obs2(j) = uin(1) + (uin(2) - uin(1)) * &
                        (gr_z_sw(j) - z_obs(1)) / (z_obs(2) - z_obs(1))
                v_obs2(j) = vin(1) + (vin(2) - vin(1)) * &
                        (gr_z_sw(j) - z_obs(1)) / (z_obs(2) - z_obs(1))
          elseif (gr_z_sw(j) .lt. z_obs(2) .and. &
                  gr_z_sw(j) .ge. z_obs(3)) then
            u_obs2(j) = uin(2) + (uin(3) - uin(2)) * &
                        (gr_z_sw(j) - z_obs(2)) / (z_obs(3) - z_obs(2))
            v_obs2(j) = vin(2) + (vin(3) - vin(2)) * &
                        (gr_z_sw(j) - z_obs(2)) / (z_obs(3) - z_obs(2))
          endif
        enddo
      endif
!     Interpolate observed velocity onto model grid
      do i = 1, gr_nn_sw
        do j = 1, gr_nlayers_sw
          u(j,i) = u_obs1(j) + pstep * (u_obs2(j) - u_obs1(j))
          v(j,i) = v_obs1(j) + pstep * (v_obs2(j) - v_obs1(j))
        enddo
      enddo

      end

!=============================================================================*
      subroutine load_mesh(iun,iper,pstep)
!=============================================================================*
!--External variables
      use unptrackmod

      implicit none

!--Input arguments
      integer :: iun, iper
      real  :: pstep
      
!--Description
!     Loads velocity data from regular grid model.

!--Local variables
      integer :: i, j, k, nn
      integer :: stat
      integer, save :: ntimesteps
      integer(8) :: time
      real,save :: ufactor
      real, allocatable, save :: uin1(:,:), vin1(:,:), win1(:,:)
      real, allocatable, save :: uin2(:,:), vin2(:,:), win2(:,:)

      if (iper .eq. 1 .and. pstep .eq. 0.) then
        allocate (uin1(gr_nlayers_sw,gr_nn_sw), &
                  vin1(gr_nlayers_sw,gr_nn_sw), &
                  win1(gr_nlayers_sw,gr_nn_sw), &
                  uin2(gr_nlayers_sw,gr_nn_sw), &
                  vin2(gr_nlayers_sw,gr_nn_sw), &
                  win2(gr_nlayers_sw,gr_nn_sw))

!       Read first set of 3D velocity data
        read(iun,*) ntimesteps, ufactor
        read(iun,*) time
        do i=1, gr_nn_sw
          do j=1, gr_nlayers_sw
            read(iun,*) nn, uin2(j,i), vin2(j,i), win2(j,i)
          enddo
        enddo
!       Apply velocity scaling
        uin2 = uin2 * ufactor
        vin2 = vin2 * ufactor
        win2 = win2 * ufactor
      endif

      if (iper .lt. ntimesteps .and. pstep .eq. 0.) then
!       Update previous values (e.g. uin1) with new values (e.g. uin2)
        uin1 = uin2
        vin1 = vin2
        win1 = win2
!       Read in new 3D velocity data
        read(iun,*) time
        do i=1, gr_nn_sw
          do j=1, gr_nlayers_sw
            read(iun,*) nn, uin2(j,i), vin2(j,i), win2(j,i)
          enddo
        enddo
!       Apply velocity scaling
        uin2 = uin2 * ufactor
        vin2 = vin2 * ufactor
        win2 = win2 * ufactor
      elseif (iper .ge. ntimesteps .and. pstep .eq. 0.) then
!       Maintain steady flow field
        uin2 = uin1
        vin2 = vin1
        win2 = win1
      endif

!     Interpolate observed velocity to current particle time step
      do i = 1, gr_nn_sw
        do j = 1, gr_nlayers_sw
          u(j,i) = uin1(j,i) + pstep * (uin2(j,i) - uin1(1,i))
          v(j,i) = vin1(j,i) + pstep * (vin2(j,i) - vin1(1,i))
          w(j,i) = win1(j,i) + pstep * (win2(j,i) - win1(1,i))
          
!         Convert vertical velocity to sigma space
          w(j,i) = w(j,i) / (eta(i) + gr_depth_sw(i))
        enddo
      enddo

      end

!=============================================================================*
      logical function isintriangle(ihost,xp,yp) 
!==============================================================================|
!     Determine if point (xp,yp) is in triangle defined by nodes (xt(3),yt(3))    |
!     using algorithm used for scene rendering in computer graphics               |
!==============================================================================|
!------------------------------------------------------------------------------|

      use  unptrackmod
      implicit none
!------------------------------------------------------------------------------!
      real, intent(in) :: xp,yp
      integer, intent(in) :: ihost
      real :: xt(3), yt(3)
      real :: f1,f2,f3
      integer :: i
      
!------------------------------------------------------------------------------|

      isintriangle = .false.  

!     Find element containing (x,y) coordinate
      do i = 1, 3
        xt(i) = gr_x_sw(gr_nodes_sw(ihost,i))
        yt(i) = gr_y_sw(gr_nodes_sw(ihost,i))
      enddo

!     Determine presence or absence of (xp, yp)
      f1 = (yp-yt(1))*(xt(2)-xt(1)) - (xp-xt(1))*(yt(2)-yt(1))
      f2 = (yp-yt(3))*(xt(1)-xt(3)) - (xp-xt(3))*(yt(1)-yt(3))
      f3 = (yp-yt(2))*(xt(3)-xt(2)) - (xp-xt(2))*(yt(3)-yt(2))
      if(f1*f3 >= 0.0 .and. f3*f2 >= 0.0) isintriangle = .true.

      return
      end function isintriangle
      

!=============================================================================*
      subroutine QuickSearch(xp,yp,ihost,ifound)
!==============================================================================|
!     Determine which element a particle resides in by searching       |
!     neighboring elements.  Updates host component of lagrangian particle     |  
!     type and updates logical array "elem_found" flagging whether the host    |
!     has been found							       |
!==============================================================================|

      use unptrackmod
      implicit none
!------------------------------------------------------------------------------!
      real, intent(in)    :: xp,yp
      integer, intent(inout) :: ihost,ifound
!------------------------------------------------------------------------------!
      integer j, k, iney, ncheck
      logical :: isintriangle
!==============================================================================|

      ifound = 0
      if(isintriangle(ihost,xp,yp)) then      !!particle remains in element
        ifound = 1 
      else                                          !!check neighbors           
        outer: do j=1,3
          ncheck = gr_nodes_sw(ihost,j)
          do k=1,gr_ntve_sw(ncheck)
            iney = gr_nbve_sw(ncheck,k) 
            if (isintriangle(iney,xp,yp)) then
              ifound = 1 
              ihost = iney
              exit outer
            endif
          enddo
        enddo outer
      endif

      return
      end subroutine

!==============================================================================|
      subroutine FullSearch(xp,yp,ihost,ifound)
!==============================================================================|
!     Find home element for points (x,y)                                       |
!  search nearest element to progressively further elements. updates lagrangian| 
!  component "host" and marks lagrangian component "ifound" with 1 if          |
!  found.  returns logical variable "all_found" if all lagrangian variables    |
!  have a known host element.  the host element may have been found prior to   |
!  entry in this routine.                                                      |
!==============================================================================|

      use unptrackmod
      implicit none
!------------------------------------------------------------------------------!
      real, intent(in)    :: xp,yp 
      integer, intent(inout) :: ihost,ifound
!------------------------------------------------------------------------------|
      integer :: i, k
      real :: xn(3), yn(3)
      logical :: isintriangle
      
      if (tm_ks_sw .gt. 0) write(*,*) 'Conducting Full Search'

      ifound = 0
      do k = 1, gr_ne_sw
        if (isintriangle(k,xp,yp)) then
          ihost = k
          ifound = 1
          return
        endif
      enddo
        
      end subroutine

!==============================================================================|
      subroutine arrayalloc
!=============================================================================*
!--External variables
      use unptrackmod

      implicit none
      
!--Description
!   Reads 3D velocity fields e.g. u, v, t from file

!--Local variables
      integer :: NN, NE, JLO, NMAX

      NN = gr_nn_sw
      NE = gr_ne_sw
      JLO = gr_nlayers_sw
      NMAX = gr_nghmax_sw

!-|Grid coordinates and bathymetry at nodes
      allocate (gr_x_sw(NN), gr_y_sw(NN))
      allocate (gr_depth_sw(NN))
      allocate (gr_dz_sw(JLO))
      allocate (gr_ntve_sw(NN), gr_nbve_sw(NN,NMAX))
      
!-|Grid element information    
      allocate (gr_xc_sw(NE), gr_yc_sw(NE))
      allocate (gr_nodes_sw(NE,3))
      allocate (gr_eledep_sw(NE))
      allocate (gr_area_sw(NE))
      allocate (gr_vol_sw(JLO,NE))

!-|Velocity array dimensions
      allocate (u(JLO,NN), v(JLO,NN), w(JLO,NN))
      allocate (tem(JLO,NN))

!-|Surface elevation
      allocate (eta(NE))

!-|Wind-driven surface velocity field
      allocate (usw(NN), vsw(NN)) 

!-|Horizontal and vertical eddy diffusion
      allocate (tu_kh_sw(JLO,NN), tu_kv_sw(JLO,NE))

!-|Particle density
      allocate (pden(tr_jpcell_sw+1,NE), psum(tr_jpcell_sw+1,NE))

!-|Particle location and characteristics
      allocate (tr_x_sw(NPMAX), tr_y_sw(NPMAX))
      allocate (tr_s_sw(NPMAX), tr_z_sw(NPMAX))
      allocate (tr_age_sw(NPMAX), tr_pfac_sw(NPMAX), &
                tr_stage_sw(NPMAX), tr_igo_sw(NPMAX))
      allocate (tr_element_sw(NPMAX), tr_wsettle_sw(NPMAX))

      end
!=============================================================================*
      subroutine terminate(ierr,iproc)
!=============================================================================*

      implicit none

      integer :: iproc,ierr
      integer, parameter :: lst=6

      write(lst,'(2i5)') ierr,iproc
      stop
      
      end subroutine

!=============================================================================*

