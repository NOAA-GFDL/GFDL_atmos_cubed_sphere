#define NC_ERR_STOP(status) \
    if (status /= nf90_noerr) write(0,*) "file: ", __FILE__, " line: ", __LINE__, trim(nf90_strerror(status)); \
    if (status /= nf90_noerr) call ESMF_Finalize(endflag=ESMF_END_ABORT)

module module_get_cubed_sphere_inc

  use netcdf
  use esmf
  use mpp_mod,            only: mpp_pe, mpp_get_current_pelist
  use tracer_manager_mod, only: get_tracer_index, get_number_tracers
  use field_manager_mod,  only: MODEL_ATMOS
  use CCPP_data,          only: GFS_control
  use time_manager_mod,   only: time_type, get_time, set_time
  use fv_arrays_mod,      only: fv_atmos_type
  use mpi

  implicit none
  type iau_internal_data_type
    real,allocatable :: ua_inc(:,:,:)
    real,allocatable :: va_inc(:,:,:)
    real,allocatable :: temp_inc(:,:,:)
    real,allocatable :: delp_inc(:,:,:)
    real,allocatable :: delz_inc(:,:,:)
    real,allocatable :: tracer_inc(:,:,:,:)
  end type iau_internal_data_type

  type netcdf_data
    character(len=NF90_MAX_NAME) :: varname
    integer                      :: varid
    integer                      :: dimid
    integer                      :: dimsize
  end type netcdf_data

  public read_netcdf, calc_iau_tendency, read_netcdf_inc, iau_internal_data_type

  logical :: par

  contains

!----------------------------------------------------------------------------------------
  logical function compvals(a, b)
    real, intent(in)     :: a
    real, intent(in)     :: b

    if(abs(a - b) < 1.0e-22) then
       compvals = .true.
    else
       compvals = .false.
    endif
    write(6,*) 'vals are ',a,b
    write(6,*) 'diff is ',abs(a - b),compvals

  end function
  logical function is_tracer(varname,tracer_idx)
    character(*), intent(in)      ::     varname
    integer, intent(out)          ::    tracer_idx 
    select case (varname)
      case("liq_wat")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'liq_wat')
      case("ice_wat")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'ice_wat')
      case("rainwat")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'rainwat')
      case("snowwat")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'snowwat')
      case("spfh")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'sphum')
      case("sphum")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'sphum')
      case("spo3")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'spo3')
      case("o3mr")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'o3mr')
      case("graupel")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'graupel')
      case("dust1")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'dust1')
      case("dust2")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'dust2')
      case("dust3")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'dust3')
      case("dust4")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'dust4')
      case("dust5")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'dust5')
      case("seas1")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'seas1')
      case("seas2")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'seas2')
      case("seas3")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'seas3')
      case("seas4")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'seas4')
      case("seas5")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'seas5')
      case("so4")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'so4')
      case("bc1")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'bc1')
      case("bc2")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'bc2')
      case("oc1")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'oc1')
      case("oc2")
        is_tracer = .true.
        tracer_idx = get_tracer_index(MODEL_ATMOS, 'oc2')
      case default
        tracer_idx = get_tracer_index(MODEL_ATMOS, trim(varname))
        write(6,*) 'tracer index is ',tracer_idx,' and name is ',trim(varname)
        is_tracer = .false.
      end select
  end function
  subroutine read_netcdf_inc(filename, increment_data, Atm, mygrid, &
                          testing, im_ret, jm_ret, pf_ret, tileCount, tests_passed,rc)
    character(*), intent(in)                         :: filename
    type(iau_internal_data_type), intent(inout)      :: increment_data
    type (fv_atmos_type), allocatable, intent(inout) :: Atm(:)
    integer, intent(in)                              :: mygrid
    logical, intent(in)                              :: testing
    integer, optional, intent(out)                   :: im_ret
    integer, optional, intent(out)                   :: jm_ret
    integer, optional, intent(out)                   :: pf_ret
    integer, optional, intent(out)                   :: tileCount
    logical, optional,intent(out)                    :: tests_passed
    integer, optional,intent(out)                    :: rc

! local variables
    integer :: i,j
    integer :: mpi_comm
    integer :: mype

    real(ESMF_KIND_R8), dimension(:,:,:,:,:), pointer :: array_r8_3d_tiled
    real(ESMF_KIND_R8)  :: testval
    integer :: ncerr
    integer :: ncid
    integer :: im_dimid, jm_dimid, tile_dimid, pfull_dimid, phalf_dimid, time_dimid 
    integer :: tm,ph
    integer :: im_varid, jm_varid, tile_varid, lon_varid, lat_varid, timeiso_varid
    integer :: liqwat_varid, sphum_varid, o3mr_varid, icewat_varid
    integer :: ugrd_varid, vgrd_varid, dpres_varid, tmp_varid 

    integer :: par_access, nvar, xtype, ndims, nAtts 
    integer :: sphum_idx, ntcw, ntiw, ntoz, ntrw, ntsw, ntgl, ntdu1, ntdu2, ntdu3, ntdu4, ntdu5
    integer :: ntss1, ntss2, ntss3, ntss4,ntss5, ntsu, ntbcb, ntbcl,ntocb, ntocl
    integer :: isc, iec, jsc, jec, TC
    integer :: im, jm, pf
    integer :: mytile, num_tracers, idx_val, varid_val
    integer, dimension(:), allocatable :: varids
    integer, dimension(:), allocatable :: tracer_idx
    character(len=NF90_MAX_NAME), dimension(:), allocatable :: varnames
    character(len=NF90_MAX_NAME) :: varname
    type(netcdf_data), allocatable    :: incvars(:)    
!
    TC = 6
    write(6,*) 'in read_netcdf_inc, testing is',testing
    if(present(tileCount)) TC = tileCount
    if(present(im_ret)) im = im_ret
    if(present(jm_ret)) jm = jm_ret
    if(present(pf_ret)) pf = pf_ret
    if(present(tests_passed)) tests_passed = .true.
    if(testing) then
      testval = 0.0
      mytile = 1
      mype = 0
      mpi_comm = 0
    else 
      testval = 1.0
      mytile = Atm(mygrid)%tile_of_mosaic
      mype = mpp_pe()
      call mpp_get_current_pelist(Atm(mygrid)%pelist, commID=mpi_comm)
    endif
    par = .false.

    if (par) then
       ! not implemented yet
       ncerr = nf90_open(trim(filename),&
               mode=nf90_nowrite, &
               comm=mpi_comm, info = MPI_INFO_NULL, ncid=ncid); NC_ERR_STOP(ncerr)
    else
       ncerr = nf90_open(trim(filename),&
               mode=nf90_nowrite, &
               ncid=ncid); NC_ERR_STOP(ncerr)
    end if
    ncerr = nf90_inquire(ncid, nvariables = nvar); NC_ERR_STOP(ncerr)
    write(6,*) 'nvars is ',nvar
    allocate(incvars(nvar))
    allocate(varids(nvar))
    allocate(tracer_idx(nvar)) ! we don't need it to be this large, but don't yet know how many tracers there are
    ncerr = nf90_inq_varids(ncid, nvar, varids); NC_ERR_STOP(ncerr)
    num_tracers = 0
    do i=1,nvar
      ncerr = nf90_inquire_variable(ncid, varids(i), name=varname, xtype=xtype, ndims=ndims, nAtts=nAtts)
      NC_ERR_STOP(ncerr)
      incvars(i)%varname = varname
      if(.not.(testing)) then
        idx_val = get_tracer_index(MODEL_ATMOS, trim(incvars(i)%varname))
        if((idx_val > 0 ).and.(idx_val <= nvar)) then
!       if( is_tracer(incvars(i)%varname,idx_val)) then 
           num_tracers = num_tracers + 1
           ! store tracer_idx by tracer number so we can iterate through 1,num_tracers later
           tracer_idx(num_tracers) = idx_val
           ! store/overwrite varids by tracer number so we can iterate through 1,num_tracers later
           varids(num_tracers) = varids(i)
        endif 
      else
        num_tracers=4
      endif
      incvars(i)%varid = varids(i)
      enddo
    write(6,*) "Number of tracers is ",num_tracers
    !get dimensions of fields in the file
    varname = "grid_yt"
    ncerr = nf90_inq_dimid(ncid, trim(varname), jm_dimid) ; NC_ERR_STOP(ncerr)
    ncerr = nf90_inquire_dimension(ncid,jm_dimid,len=jm); NC_ERR_STOP(ncerr)
    varname = "grid_xt"
    ncerr = nf90_inq_dimid(ncid, trim(varname), im_dimid) ; NC_ERR_STOP(ncerr)
    ncerr = nf90_inquire_dimension(ncid,im_dimid,len=im); NC_ERR_STOP(ncerr)
    varname = "time"
    ncerr = nf90_inq_dimid(ncid, trim(varname), time_dimid) ; NC_ERR_STOP(ncerr)
    ncerr = nf90_inquire_dimension(ncid,time_dimid,len=tm); NC_ERR_STOP(ncerr)
    varname = "tile"
    ncerr = nf90_inq_dimid(ncid, trim(varname), tile_dimid) ; NC_ERR_STOP(ncerr)
    ncerr = nf90_inquire_dimension(ncid,tile_dimid,len=TC); NC_ERR_STOP(ncerr)
    varname = "pfull"
    ncerr = nf90_inq_dimid(ncid, trim(varname), pfull_dimid) ; NC_ERR_STOP(ncerr)
    ncerr = nf90_inquire_dimension(ncid,pfull_dimid,len=pf); NC_ERR_STOP(ncerr)
    varname = "phalf"
    ncerr = nf90_inq_dimid(ncid, trim(varname), phalf_dimid) ; NC_ERR_STOP(ncerr)
    ncerr = nf90_inquire_dimension(ncid,phalf_dimid,len=ph); NC_ERR_STOP(ncerr)

    ! return sizes of arrays, if requested
    if(present(im_ret)) im_ret = im
    if(present(jm_ret)) jm_ret = jm
    if(present(pf_ret)) pf_ret = pf
    if(present(tileCount)) tileCount = TC

    if(.not. allocated(increment_data%ua_inc)) then
      ! Allocate space in increment 
      allocate(increment_data%ua_inc(im,jm,pf))
      allocate(increment_data%va_inc(im,jm,pf))
      allocate(increment_data%temp_inc(im,jm,pf))
      allocate(increment_data%delp_inc(im,jm,pf))
      allocate(increment_data%delz_inc(im,jm,pf))
      allocate(increment_data%tracer_inc(im,jm,pf,num_tracers)) 
    endif
    if(testing) then
      ! assign dummy indices
      varname = "ice_wat"
      ncerr = nf90_inq_varid(ncid,trim(varname),icewat_varid); NC_ERR_STOP(ncerr); NC_ERR_STOP(ncerr)
      varname = "liq_wat"
      ncerr = nf90_inq_varid(ncid,trim(varname),liqwat_varid); NC_ERR_STOP(ncerr); NC_ERR_STOP(ncerr)
      varname = "spfh"
      ncerr = nf90_inq_varid(ncid,trim(varname),sphum_varid); NC_ERR_STOP(ncerr); NC_ERR_STOP(ncerr)
      varname = "o3mr"
      ncerr = nf90_inq_varid(ncid,trim(varname),o3mr_varid); NC_ERR_STOP(ncerr)
      varids(1) = sphum_varid ! spfh
      varids(2) = icewat_varid ! ice_wat
      varids(3) = liqwat_varid ! liq_wat
      varids(4) = o3mr_varid ! o3mr
      tracer_idx(1) = 1
      tracer_idx(2) = 2
      tracer_idx(3) = 3
      tracer_idx(4) = 4
      tracer_idx(5) = 5
      write(6,*) 'test varid is ',varids(1),' tracer_idx is ',tracer_idx(1)
      write(6,*) 'test varid is ',varids(2),' tracer_idx is ',tracer_idx(2)
      write(6,*) 'test varid is ',varids(3),' tracer_idx is ',tracer_idx(3)
      write(6,*) 'test varid is ',varids(4),' tracer_idx is ',tracer_idx(4)
      ! for non-testing, indices were returned from check for is_tracer
    endif
    ! allocate temporary array to hold variables
    ! TODO read only what we need instead of the whole field
    allocate(array_r8_3d_tiled(im,jm,pf,TC,tm))
    isc = GFS_control%isc
    iec = GFS_control%isc+GFS_control%nx-1
    jsc = GFS_control%jsc
    jec = GFS_control%jsc+GFS_control%ny-1


    do i = 1,nvar
      select case (incvars(i)%varname)
        case("ugrd")
          ncerr = nf90_get_var(ncid, incvars(i)%varid, array_r8_3d_tiled) ; NC_ERR_STOP(ncerr)
          increment_data%ua_inc(isc:iec,jsc:jec,:) = array_r8_3d_tiled(isc:iec,jsc:jec,:,mytile,1)
        case("ua")
          ncerr = nf90_get_var(ncid, incvars(i)%varid, array_r8_3d_tiled) ; NC_ERR_STOP(ncerr)
          increment_data%ua_inc(isc:iec,jsc:jec,:) = array_r8_3d_tiled(isc:iec,jsc:jec,:,mytile,1)
        case("vgrd")
          ncerr = nf90_get_var(ncid, incvars(i)%varid, array_r8_3d_tiled) ; NC_ERR_STOP(ncerr)
          increment_data%va_inc(isc:iec,jsc:jec,:) = array_r8_3d_tiled(isc:iec,jsc:jec,:,mytile,1)
        case("va")
          ncerr = nf90_get_var(ncid, incvars(i)%varid, array_r8_3d_tiled) ; NC_ERR_STOP(ncerr)
          increment_data%va_inc(isc:iec,jsc:jec,:) = array_r8_3d_tiled(isc:iec,jsc:jec,:,mytile,1)
        case("tmp")
          ncerr = nf90_get_var(ncid, incvars(i)%varid, array_r8_3d_tiled) ; NC_ERR_STOP(ncerr)
          increment_data%temp_inc(isc:iec,jsc:jec,:) = array_r8_3d_tiled(isc:iec,jsc:jec,:,mytile,1)
        case("T")
          ncerr = nf90_get_var(ncid, incvars(i)%varid, array_r8_3d_tiled) ; NC_ERR_STOP(ncerr)
          increment_data%temp_inc(isc:iec,jsc:jec,:) = array_r8_3d_tiled(isc:iec,jsc:jec,:,mytile,1)
        case("dpres")
          ncerr = nf90_get_var(ncid, incvars(i)%varid, array_r8_3d_tiled) ; NC_ERR_STOP(ncerr)
          increment_data%delp_inc(isc:iec,jsc:jec,:) = array_r8_3d_tiled(isc:iec,jsc:jec,:,mytile,1)
        case("DELP")
          ncerr = nf90_get_var(ncid, incvars(i)%varid, array_r8_3d_tiled) ; NC_ERR_STOP(ncerr)
          increment_data%delp_inc(isc:iec,jsc:jec,:) = array_r8_3d_tiled(isc:iec,jsc:jec,:,mytile,1)
        case("delz")
          ncerr = nf90_get_var(ncid, incvars(i)%varid, array_r8_3d_tiled) ; NC_ERR_STOP(ncerr)
          increment_data%delz_inc(isc:iec,jsc:jec,:) = array_r8_3d_tiled(isc:iec,jsc:jec,:,mytile,1)
      end select
    enddo
    do i = 1,num_tracers
      write(6,*) 'varid is ',varids(i),' tracer_idx is ',tracer_idx(i)
      if((tracer_idx(i) > 0) .and. (tracer_idx(i) <= num_tracers)) then
        ncerr = nf90_get_var(ncid, varids(i), array_r8_3d_tiled) ; NC_ERR_STOP(ncerr)
        increment_data%tracer_inc(isc:iec,jsc:jec,:,tracer_idx(i)) = array_r8_3d_tiled(isc:iec,jsc:jec,:,mytile,1)
      endif
    enddo

    deallocate(array_r8_3d_tiled)
    deallocate(incvars)
    deallocate(varids)
    deallocate(tracer_idx)
  
  end subroutine read_netcdf_inc
  subroutine read_netcdf(filename, Atm, mygrid, &
                          use_parallel_netcdf, &
                          testing,tests_passed,rc)
!
    character(*), intent(in)                         :: filename
    type (fv_atmos_type), allocatable, intent(inout) :: Atm(:)
    integer,          intent(in)                     :: mygrid
    logical, intent(in)                              :: use_parallel_netcdf
    logical, intent(in)                              :: testing
    logical, optional,intent(out)                    :: tests_passed
    integer, optional,intent(out)                    :: rc
!
!** local vars
    type(iau_internal_data_type)               :: increment_data
    integer :: i,j,k
    integer :: im, jm
    integer :: mpi_comm
    integer :: mype

    real(ESMF_KIND_R8), dimension(:,:,:,:,:), pointer :: array_r8_3d_tiled
    real(ESMF_KIND_R8)  :: testval

    integer :: ncerr
    integer :: ncid
    integer :: im_dimid, jm_dimid, tile_dimid, pfull_dimid, phalf_dimid, time_dimid 
    integer :: tm,pf,ph, tileCount
    integer :: im_varid, jm_varid, tile_varid, lon_varid, lat_varid, timeiso_varid
    integer :: liqwat_varid, sphum_varid, o3mr_varid, icewat_varid
    integer :: ugrd_varid, vgrd_varid, dpres_varid, tmp_varid 

    integer :: par_access, nvar, xtype, ndims, nAtts 
    integer :: sphum_idx, ntcw, ntiw, ntoz
    integer :: isc, iec, jsc, jec
    integer :: mytile
    integer, dimension(:), allocatable :: varids
    character(len=NF90_MAX_NAME) :: varname

    if(testing) then
      testval = 0.0
      mytile = 1
      tileCount = 6
      mype = 0
      mpi_comm = 0
      allocate(increment_data%ua_inc(96,96,127))
      allocate(increment_data%va_inc(96,96,127))
      allocate(increment_data%temp_inc(96,96,127))
      allocate(increment_data%delp_inc(96,96,127))
      allocate(increment_data%delz_inc(96,96,127))
      allocate(increment_data%tracer_inc(96,96,127,6))
    else 
      testval = 1.0
      mytile = Atm(mygrid)%tile_of_mosaic
      mype = mpp_pe()
      call mpp_get_current_pelist(Atm(mygrid)%pelist, commID=mpi_comm)
    endif
    par = .false.
    call read_netcdf_inc(filename, increment_data,Atm, mygrid, testing,im_ret=im, jm_ret=jm, pf_ret=pf, tileCount=tileCount, tests_passed=tests_passed, rc=rc)
    if(testing) then
      ! allocate 6 tiles for Atm
      write(6,*) "im, jm, etc are",im,jm,pf
      allocate(Atm(6))
      ! assign dummy indices
      sphum_idx   = 1
      ntiw = 2
      ntcw = 3
      ntoz    = 4
      ! Allocate space in Atm for testing 
      do i=1,tileCount
        allocate(Atm(i)%u(im,jm,pf))
        allocate(Atm(i)%v(im,jm,pf))
        allocate(Atm(i)%pt(im,jm,pf))
        allocate(Atm(i)%delp(im,jm,pf))
        allocate(Atm(i)%q(im,jm,pf,4)) ! currently, just need specific humdity, liq_wat, ice_wat, and o3mr
        Atm(i)%u(:,:,:) = 0.0
        Atm(i)%v(:,:,:) = 0.0
        Atm(i)%pt(:,:,:) = 0.0
        Atm(i)%delp(:,:,:) = 0.0
        Atm(i)%q(:,:,:,:) = 0.0
      enddo
    else
      ! get the tracer index to update variables in Atm(t)%tracer_inc
      ! this in not available when on the write component
      sphum_idx   = get_tracer_index(MODEL_ATMOS, 'sphum')
      ntiw = get_tracer_index(MODEL_ATMOS, 'ice_wat')
      ntcw = get_tracer_index(MODEL_ATMOS, 'liq_wat')
      ntoz    = get_tracer_index(MODEL_ATMOS, 'o3mr')
    endif
    isc = GFS_control%isc
    iec = GFS_control%isc+GFS_control%nx-1
    jsc = GFS_control%jsc
    jec = GFS_control%jsc+GFS_control%ny-1
    !Update u
    Atm(mygrid)%u(isc:iec,jsc:jec,:) = Atm(mygrid)%u(isc:iec,jsc:jec,:) + &
      testval * increment_data%ua_inc(isc:iec,jsc:jec,:)
    
    !Update v
    Atm(mygrid)%v(isc:iec,jsc:jec,:) = Atm(mygrid)%v(isc:iec,jsc:jec,:) + &
      testval * increment_data%va_inc(isc:iec,jsc:jec,:)

    !Update potential temp
    Atm(mygrid)%pt(isc:iec,jsc:jec,:) = Atm(mygrid)%pt(isc:iec,jsc:jec,:) + &
      testval * increment_data%temp_inc(isc:iec,jsc:jec,:)

    !Update delp
    Atm(mygrid)%delp(isc:iec,jsc:jec,:) = Atm(mygrid)%delp(isc:iec,jsc:jec,:) + &
      testval * increment_data%delp_inc(isc:iec,jsc:jec,:)

    !Update sphum
    Atm(mygrid)%q(isc:iec,jsc:jec,:,sphum_idx) = Atm(mygrid)%q(isc:iec,jsc:jec,:,sphum_idx) + &
      testval * increment_data%tracer_inc(isc:iec,jsc:jec,:,sphum_idx)

    !Update ice_wat
    Atm(mygrid)%q(isc:iec,jsc:jec,:,ntiw) = Atm(mygrid)%q(isc:iec,jsc:jec,:,ntiw) + &
      testval * increment_data%tracer_inc(isc:iec,jsc:jec,:,ntiw)

    !Update liq_wat
    Atm(mygrid)%q(isc:iec,jsc:jec,:,ntcw) = Atm(mygrid)%q(isc:iec,jsc:jec,:,ntcw) + &
      testval * increment_data%tracer_inc(isc:iec,jsc:jec,:,ntcw)
    !Update o3mr 
    Atm(mygrid)%q(isc:iec,jsc:jec,:,ntoz) = Atm(mygrid)%q(isc:iec,jsc:jec,:,ntoz) + &
      testval * increment_data%tracer_inc(isc:iec,jsc:jec,:,ntoz)

    if(testing) then
        tests_passed = tests_passed .and. compvals(increment_data%ua_inc(isc,jsc,1) ,-1.8169951587765354E-007)
        tests_passed = tests_passed .and. compvals(increment_data%va_inc(isc,jsc,1) , 3.0927015537418612E-007) 
        tests_passed = tests_passed .and. compvals(increment_data%delp_inc(isc,jsc,1), -1.8873791418627661E-015)
        tests_passed = tests_passed .and. compvals(increment_data%temp_inc(isc,jsc,1) , -6.8499218741635559E-008)
        tests_passed = tests_passed .and. compvals(increment_data%tracer_inc(isc,jsc,1,ntiw) , -4.4666691960233954E-019)
        tests_passed = tests_passed .and. compvals(increment_data%tracer_inc(isc,jsc,1,ntcw) , -7.3514147181582930E-022)
        tests_passed = tests_passed .and. compvals(increment_data%tracer_inc(isc,jsc,1,ntoz) , -1.2690141736645521E-017)
        tests_passed = tests_passed .and. compvals(increment_data%tracer_inc(isc,jsc,1,sphum_idx) , -2.9569591390493574E-006)

      do i=1,tileCount
        deallocate(Atm(i)%u)
        deallocate(Atm(i)%v)
        deallocate(Atm(i)%pt)
        deallocate(Atm(i)%delp)
        deallocate(Atm(i)%q) ! currently, just need specific humdity, liq_wat, ice_wat, and o3mr
      enddo
      deallocate(Atm) ! currently, just need specific humdity, liq_wat, ice_wat, and o3mr
    endif 
      
  end subroutine read_netcdf
  subroutine calc_iau_tendency(filenames, &
                          Atm, mygrid, &
                          tendency, increment, &
                          testing, rc)
!
    character(len=128), dimension(3), intent(in) :: filenames
    type (fv_atmos_type), allocatable, intent(inout) :: Atm(:)
    integer,                             intent(in) :: mygrid
    type(iau_internal_data_type), intent(inout)     :: tendency
    type(iau_internal_data_type), intent(out),optional :: increment
    logical, intent(in), optional                   :: testing
    integer, intent(out),optional                   :: rc
!
!** local vars
    integer :: mpi_comm
    integer :: mype
    type(iau_internal_data_type)               :: increment_data(3)
    integer :: i,j,k
    integer :: im, jm

    real(ESMF_KIND_R8), dimension(:,:,:,:,:), pointer :: array_r8_3d_tiled
    real(ESMF_KIND_R8)  :: testval

    integer :: ncerr
    integer :: ncid
    integer :: im_dimid, jm_dimid, tile_dimid, pfull_dimid, phalf_dimid, time_dimid 
    integer :: tm,pf,ph, tileCount
    integer :: im_varid, jm_varid, tile_varid, lon_varid, lat_varid, timeiso_varid
    integer :: liqwat_varid, sphum_varid, o3mr_varid, icewat_varid
    integer :: ugrd_varid, vgrd_varid, dpres_varid, tmp_varid 

    integer :: par_access, nvar, xtype, ndims, nAtts 
    integer :: sphum_idx, ntcw, ntiw, ntoz, ntrw, ntsw, ntgl, ntdu1, ntdu2, ntdu3, ntdu4, ntdu5
    integer :: ntss1, ntss2, ntss3, ntss4,ntss5, ntsu, ntbcb, ntbcl,ntocb, ntocl
    integer :: isc, iec, jsc, jec
    integer :: mytile
    real :: dt_atmos
    integer, dimension(:), allocatable :: varids
    character(len=NF90_MAX_NAME) :: varname
    logical do_testing

    do_testing = .false.
    dt_atmos = 250.
    if(present(testing)) do_testing=testing
    if(do_testing) then
      testval = 0.0
      mytile = 1
      mype = 0
      mpi_comm = 0
    else 
      testval = 1.0
      mytile = Atm(mygrid)%tile_of_mosaic
      mype = mpp_pe()
      call mpp_get_current_pelist(Atm(mygrid)%pelist, commID=mpi_comm)
    endif
    par = .false.
    ! read in 3? increments
    do i=1,3
      call read_netcdf_inc(filenames(i), increment_data(i),Atm,mygrid,do_testing,im_ret=im, jm_ret=jm, pf_ret=pf, tileCount=tileCount, rc=rc)
    enddo
    ! allocate space for tendencies
    allocate(tendency%ua_inc(im,jm,pf))
    allocate(tendency%va_inc(im,jm,pf))
    allocate(tendency%temp_inc(im,jm,pf))
    allocate(tendency%delp_inc(im,jm,pf))
    allocate(tendency%tracer_inc(im,jm,pf,size(increment_data(1)%tracer_inc,4))) 
    ! Calculate the tendencies by subtracting the first increment from the second and dividing by dt_atmos
    if(present(increment)) increment = increment_data(2)
    tendency%ua_inc(:,:,:) = (increment_data(2)%ua_inc(:,:,:) - increment_data(1)%ua_inc(:,:,:))/dt_atmos
    tendency%va_inc(:,:,:) = (increment_data(2)%va_inc(:,:,:) - increment_data(1)%va_inc(:,:,:))/dt_atmos
    tendency%temp_inc(:,:,:) = (increment_data(2)%temp_inc(:,:,:) - increment_data(1)%temp_inc(:,:,:))/dt_atmos
    tendency%delp_inc(:,:,:) = (increment_data(2)%delp_inc(:,:,:) - increment_data(1)%delp_inc(:,:,:))/dt_atmos
    tendency%tracer_inc(:,:,:,:) = (increment_data(2)%tracer_inc(:,:,:,:) - increment_data(1)%tracer_inc(:,:,:,:))/dt_atmos

  end subroutine calc_iau_tendency
  
  
!----------------------------------------------------------------------------------------
end module module_get_cubed_sphere_inc
