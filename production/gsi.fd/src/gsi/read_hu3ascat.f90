subroutine read_hu3ascat(nread,ndata,nodata,infile,obstype,lunout,gstime,twind,sis,&
     prsl_full,nobs)
!$$$  subprogram documentation block
! subroutine based primarily on read_sfcwnd.f90: documentation block for that
! routine follows the history log.
! program history log:
!   2022-05-12: Hoover - HU-3 (AI-derived) ASCAT processed as ASCAT winds
!                      - Optimal uob, vob determined from max likelihood of
!                        ambiguities
!                .      .    .                                       .
! subprogram:    read_sfcwnd                    read scatterometer winds
!   prgmmr: Li Bi                               date: 2012-08-20
!
! abstract:  This routine reads OSCAT scatterometer winds from dump.        
!            it also has options to thin the data by using conventional 
!            thinning programs 
!            When running the gsi in regional mode, the code only
!            retains those observations that fall within the regional
!            domain
!
! program history log:
!   2012-08-20 Li Bi      
!   2014-04-15 Su -  new error table
!   2015-02-23  Rancic/Thomas - add thin4d to time window logical
!   2015-03-23  Su      -fix array size with maximum message and subset number from fixed number to
!                        dynamic allocated array
!   2015-10-01  guo     - consolidate use of ob location (in deg)
!   2016-03-15  Su      - modified the code so that the program won't stop when
!                         no subtype is found in non linear qc error table and b table
!   2020-05-04  wu   - no rotate_wind for fv3_regional
!
!   input argument list:
!     ithin    - flag to thin data
!     rmesh    - thinning mesh size (km)
!     gstime   - analysis time in minutes from reference date
!     infile   - unit from which to read BUFR data
!     lunout   - unit to which to write data for further processing
!     obstype  - observation type to process
!     twind    - input group time window (hours)
!     sis      - satellite/instrument/sensor indicator
!
!   output argument list:
!     nread    - number of satellite winds read 
!     ndata    - number of satellite winds retained for further processing
!     nodata   - number of satellite winds retained for further processing
!     nobs     - array of observations on each subdomain for each processor
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
! Changes relative to read_sfcwnd
! 
!$$$
  use kinds, only: r_kind,r_double,i_kind,r_single
  use gridmod, only: diagnostic_reg,regional,nlon,nlat,nsig,&
       tll2xy,txy2ll,rotate_wind_ll2xy,rotate_wind_xy2ll,&
       rlats,rlons,fv3_regional
  use qcmod, only: errormod,noiqc,njqc

  use convthin, only: make3grids,map3grids,del3grids,use_all
  use constants, only: deg2rad,zero,rad2deg,one_tenth,&
        tiny_r_kind,huge_r_kind,r60inv,one_tenth,&
        one,two,three,four,five,half,quarter,r60inv,r10,r100,r2000
  use converr,only: etabl
  use converr_uv,only: etabl_uv,isuble_uv,maxsub_uv
  use convb_uv,only: btabl_uv
  use obsmod, only: ran01dom,bmiss
  use convinfo, only: nconvtype, &
       icuse,ictype,icsubtype,ioctype, &
       ithin_conv,rmesh_conv,pmesh_conv
  use gsi_4dvar, only: l4dvar,l4densvar,iwinbgn,winlen,time_4dvar,thin4d
  use deter_sfc_mod, only: deter_sfc_type,deter_sfc2
  use mpimod, only: npe
  implicit none

! Declare passed variables
  character(len=*)                      ,intent(in   ) :: infile,obstype
  character(len=20)                     ,intent(in   ) :: sis
  integer(i_kind)                       ,intent(in   ) :: lunout
  integer(i_kind)                       ,intent(inout) :: nread,ndata,nodata
  integer(i_kind),dimension(npe)        ,intent(inout) :: nobs
  real(r_kind)                          ,intent(in   ) :: twind
  real(r_kind),dimension(nlat,nlon,nsig),intent(in   ) :: prsl_full

! Declare local parameters

  real(r_kind),parameter:: r6= 6.0_r_kind
  real(r_kind),parameter:: r90= 90.0_r_kind
  real(r_kind),parameter:: r110= 110.0_r_kind
  real(r_kind),parameter:: r360 = 360.0_r_kind
  real(r_kind),parameter:: r421=421.0_r_kind
  real(r_kind),parameter:: r1200= 1200.0_r_kind
  
  

! Declare local variables
  logical outside
  logical luse,ithinp
  logical,allocatable,dimension(:,:):: lmsg     ! set true when convinfo entry id found in a message

  character(70) obstr,hdrtr,wndstr
  character(8) subset
  character(8) c_prvstg,c_sprvstg

  integer(i_kind) ireadmg,ireadsb,iuse,mxtb,nmsgmax
  integer(i_kind) i,maxobs,idomsfc,nsattype,j,ncount
  integer(i_kind) nc,nx,isflg,itx,nchanl
  integer(i_kind) ntb,ntmatch,ncx,ncsave,ntread
  integer(i_kind) kk,klon1,klat1,klonp1,klatp1
  integer(i_kind) nmind,lunin,idate,ilat,ilon,iret,k
  integer(i_kind) nreal,ithin,iout,ntmp,icount,iiout,ii
  integer(i_kind) itype,iosub,ixsub,isubsub,iobsub 
  integer(i_kind) lim_qm
  integer(i_kind) nlevp         ! vertical level for thinning
  integer(i_kind) pflag
  integer(i_kind) ntest,nvtest
  integer(i_kind) kl,k1,k2
  integer(i_kind) nmsg                ! message index
  integer(i_kind) qc1,ierr
  
  
 
  integer(i_kind),dimension(nconvtype) :: ntxall 
  integer(i_kind),dimension(nconvtype+1) :: ntx  
  
  integer(i_kind),dimension(5):: idate5 
  integer(i_kind),allocatable,dimension(:):: isort,iloc,nrep
  integer(i_kind),allocatable,dimension(:,:)::tab

! integer(i_kind) itypex,lcount,iflag,m
  integer(i_kind),dimension(1):: opt_idx
  integer(i_kind) itypey
  real(r_kind) toff,t4dv
  real(r_kind) rmesh,ediff,usage,tdiff
  real(r_kind) u0,v0,uob,vob,dx,dy,dx1,dy1,w00,w10,w01,w11
  real(r_kind) dlnpob,ppb,ppb2,qifn,qify,ee,var_jb
  real(r_kind) woe,dlat,dlon,dlat_earth,dlon_earth,oelev
  real(r_kind) dlat_earth_deg,dlon_earth_deg
  real(r_kind) cdist,disterr,disterrmax,rlon00,rlat00
  real(r_kind) vdisterrmax,u00,v00,uob1,vob1
  real(r_kind) del,werrmin,obserr,ppb1,wjbmin
  real(r_kind) tsavg,ff10,sfcr,sstime,gstime,zz
  real(r_kind) crit1,timedif,xmesh,pmesh
  real(r_kind),dimension(nsig):: presl
  real(r_kind) uob_1,uob_2,uob_3,uob_4,vob_1,vob_2,vob_3,vob_4
  real(r_kind),dimension(4):: uob_vec
  real(r_kind),dimension(4):: vob_vec
  real(r_kind) uob_rng, vob_rng ! "extreme" u- and v-wind ranges
  real(r_kind) woe_est ! estimated woe value from "extreme" u/v ranges
  real(r_kind) err_spd, err_dir ! wind speed and direction error

  real(r_double),dimension(8):: hdrdat
  real(r_double),dimension(2):: satqc
  real(r_double),dimension(5):: obsdat
  real(r_double),dimension(3,4):: wnddat
  real(r_double),dimension(1,1):: r_prvstg,r_sprvstg
  real(r_kind),allocatable,dimension(:):: presl_thin
  real(r_kind),allocatable,dimension(:,:):: cdata_all,cdata_out

! equivalence to handle character names
  equivalence(r_prvstg(1,1),c_prvstg)
  equivalence(r_sprvstg(1,1),c_sprvstg)


!******** Modify below from the bufrtable: 
  data hdrtr /'SAID CLATH CLONH YEAR MNTH DAYS HOUR MINU'/ 
  data obstr/'DBID EESSM MWS10 MWD10 WVCQ'/  
  
  
  data ithin / -9 /
  data lunin / 11 /
  data rmesh / -99.999_r_kind /

!**************************************************************************

write(6,*) ' READ_HU3ASCAT: entering routine'

! Set lower limits for observation errors
! ** keep this way for now
! nreal keep the dimension of cdata_all 
  werrmin=one
  nsattype=0
  nreal=24
  if (noiqc) then
     lim_qm=8
  else
     lim_qm=4
  endif

! ** read convtype from convinfo file 
! ** only read in ASCAT 290 for now ** 
  ntread=1
  ntmatch=0
  ntx(ntread)=0
  ntxall=0
  do nc=1,nconvtype
     if(trim(ioctype(nc)) == 'uv' .and. ictype(nc) == 290) then
        ntmatch=ntmatch+1
        ntxall(ntmatch)=nc
        ithin=ithin_conv(nc)
        if(ithin > 0)then
           ntread=ntread+1
           ntx(ntread)=nc
        end if
     end if
  end do
  if(ntmatch == 0)then
     write(6,*) ' READ_HU3ASCAT: no matching obstype found in convinfo ',obstype
     return
  end if
      
!!  go through the satedump to find out how many subset to process
!** Open and read data from bufr data file

  open(lunin,file=trim(infile),form='unformatted')
  call openbf(lunin,'IN',lunin)
  call datelen(10)

!! get message and subset counts

  call getcount_bufr(infile,nmsgmax,mxtb)

  allocate(lmsg(nmsgmax,ntread),tab(mxtb,3),nrep(nmsgmax))

  lmsg = .false.
  maxobs=0
  tab=0
  nmsg=0
  nrep=0
  ntb =0
  ncount=0
  msg_report: do while (ireadmg(lunin,subset,idate) == 0)

!    Time offset
     if(nmsg == 0) call time_4dvar(idate,toff)
     nmsg=nmsg+1
     if (nmsg>nmsgmax) then
        write(6,*)'READ_HU3ASCAT: messages exceed maximum ',nmsgmax
        call stop2(49)
     endif
     loop_report: do while (ireadsb(lunin) == 0)
        ntb = ntb+1
        maxobs=maxobs+1
        nrep(nmsg)=nrep(nmsg)+1
        if (ntb>mxtb) then
           write(6,*)'READ_HU3ASCAT: reports exceed maximum ',mxtb   
           call stop2(49)
        endif
!** Extract sat ID information from BUFR and assign type 
!   This part will extract information from the bufrtable
!   iobsub is the prepbufr subtype, still part of the convinfo file
!********* iobsub= 3 or 4 for existing ASCAT winds
!********* we will use iobsub=6 to distinguish HU-3 ASCAT winds
!********* NOTE: Max subtype=7 for uv obs

        call ufbint(lunin,hdrdat,8,1,iret,hdrtr)
!       determine the satellite wind type 
!       BTH: Spoofing with HU3ASCAT: subset=NC012122, hdrdat(1)=5.0                                       
        iobsub=6
        itype=-1
        if(trim(subset) == 'NC012122') then
           if( hdrdat(1) < r6 ) then           !    HU3ASCAT data == 5.0
             itype=290
           endif
        endif
!  Match ob to proper convinfo type
        ncsave=0
        matchloop:do ncx=1,ntmatch
           nc=ntxall(ncx)
           if (itype /= ictype(nc)) cycle matchloop
        write(6,*) 'BTH: nc=',nc
!  Find convtype which match ob type and subtype
           if(icsubtype(nc) == iobsub) then
              ncsave=nc
              exit matchloop
           else
!  Find convtype which match ob type and subtype group (isubtype == ?*)
!       where ? specifies the group and icsubtype = ?6)
              ixsub=icsubtype(nc)/10
              iosub=iobsub/10
              isubsub=icsubtype(nc)-ixsub*10
              if(ixsub == iosub .and. isubsub == 0) then
                 ncsave=nc
!  Find convtype which match ob type and subtype is all remaining
!       (icsubtype(nc) = 6)
              else if (ncsave == 0 .and. icsubtype(nc) == 6) then
                 ncsave=nc
              end if
           end if
        end do matchloop
!  Save information for next read
        if(ncsave /= 0) then
           maxobs=maxobs+1
           nx=1
           if(ithin_conv(ncsave) > 0)then
              do ii=2,ntread
                 if(ntx(ii) == ncsave)nx=ii
              end do
           end if
           tab(ntb,1)=ncsave
           tab(ntb,2)=nx
           tab(ntb,3)=1
           lmsg(nmsg,nx) = .true.
        end if
     enddo loop_report
  enddo msg_report

! Loop over convinfo file entries; operate on matches
  write(6,*) 'BTH: maxobs final=',maxobs
  allocate(cdata_all(nreal,maxobs),isort(maxobs))
  isort = 0
  cdata_all=zero
  nread=0
  ntest=0
  nvtest=0
  nchanl=0
  ilon=2
  ilat=3

!!  read satellite winds one type a time
!   same as in the read_prepbufr.f90 file

  loop_convinfo: do nx=1,ntread 
     use_all = .true.
     ithin=0
     if(nx >1) then
        nc=ntx(nx)
        ithin=ithin_conv(nc)
        if (ithin > 0 ) then
           rmesh=rmesh_conv(nc)
           pmesh=pmesh_conv(nc)
           use_all = .false.
           if(pmesh > zero) then
              pflag=1
              nlevp=r1200/pmesh
           else
              pflag=0
              nlevp=nsig
           endif
           xmesh=rmesh

           call make3grids(xmesh,nlevp)

           if (.not.use_all) then
              allocate(presl_thin(nlevp))
              if (pflag==1) then
                 do k=1,nlevp
                    presl_thin(k)=(r1200-(k-1)*pmesh)*one_tenth
                 enddo
              endif
           endif

           write(6,*)'READ_HU3ASCAT: ictype(nc),rmesh,pflag,nlevp,pmesh,nc ',&
                   ioctype(nc),ictype(nc),rmesh,pflag,nlevp,pmesh,nc
        endif
     endif

     call closbf(lunin)
     close(lunin)
     open(lunin,file=trim(infile),form='unformatted')
     call openbf(lunin,'IN',lunin)
     call datelen(10)

! Big loop over BUFR file

     ntb = 0
     nmsg = 0
     loop_msg:  do while(ireadmg(lunin,subset,idate) == 0)
        nmsg = nmsg+1
        if(.not.lmsg(nmsg,nx)) then
           ntb=ntb+nrep(nmsg)
           write(6,*) 'BTH: nmsg no useable reports, skip ahead',nmsg
           cycle loop_msg ! no useable reports this mesage, skip ahead report count
        end if

        loop_readsb: do while(ireadsb(lunin) == 0)
!          use msg lookup table to decide which messages to skip
!          use report id lookup table to only process matching reports
           ntb = ntb+1
           nc=tab(ntb,1)
           if(nc <= 0 .or. tab(ntb,2) /= nx) cycle loop_readsb

           hdrdat=bmiss
           obsdat=bmiss
           wnddat=bmiss
           satqc=bmiss
           iobsub=6
           itype=-1
           uob=bmiss
           vob=bmiss
           ppb=bmiss
           ppb1=bmiss
           ppb2=bmiss
           uob1=bmiss
           vob1=bmiss
           ee=r110
           qifn=r110
           qify=r110
           uob_1=bmiss
           vob_1=bmiss
           uob_2=bmiss
           vob_2=bmiss
           uob_3=bmiss
           vob_3=bmiss
           uob_4=bmiss
           vob_4=bmiss
           uob_rng=bmiss
           vob_rng=bmiss
           woe_est=bmiss
           err_spd=bmiss
           err_dir=bmiss
           qc1=bmiss
 
! Extract type, date, and location information
           call ufbint(lunin,hdrdat,8,1,iret,hdrtr) 
           call ufbint(lunin,obsdat,5,1,iret,obstr)
           !call ufbrep(lunin,wnddat,3,4,iret,wndstr)


!** potential can reject bad cell, etc, place holder for now
!   reject the data with bad quality mark from SDM
!** cycle loop means skip to the next record     

           if(hdrdat(1) >= r6) cycle loop_readsb
!       Compare relative obs time with window.  If obs 
!       falls outside of window, don't use this obs
           idate5(1) = hdrdat(4)     ! year
           idate5(2) = hdrdat(5)     ! month
           idate5(3) = hdrdat(6)     ! day
           idate5(4) = hdrdat(7)     ! hours
           idate5(5) = hdrdat(8)     ! minutes
           call w3fs21(idate5,nmind)
           t4dv = real((nmind-iwinbgn),r_kind)*r60inv
           if (l4dvar.or.l4densvar) then
              if (t4dv<zero .OR. t4dv>winlen) cycle loop_readsb 
           else
              sstime = real(nmind,r_kind) 
              tdiff=(sstime-gstime)*r60inv
              if (abs(tdiff)>twind) cycle loop_readsb 
           endif


!       determine the satellite wind type as in prepbufr
!       290: ASCAT winds                              
           iosub=6
           if(abs(hdrdat(2)) >r90 ) cycle loop_readsb 
           if(hdrdat(3) <zero) hdrdat(3)=hdrdat(3)+r360
           if(hdrdat(3) == r360) hdrdat(3)=hdrdat(3)-r360
           if(hdrdat(3) >r360) cycle loop_readsb 
           !if(abs(obsdat(2)) >= 100) cycle loop_readsb
             ! applied to wind direction, move to lower bloc
             ! when optimal wind vector has been selected
           if(obsdat(1) >=1) cycle loop_readsb

           if(trim(subset) == 'NC012122') then    ! HU3ASCAT wind
              if( hdrdat(1) < r6) then          
                   itype=290
              endif
           endif

           nread=nread+2
           dlon_earth_deg=hdrdat(3)
           dlat_earth_deg=hdrdat(2)
           dlon_earth=hdrdat(3)*deg2rad
           dlat_earth=hdrdat(2)*deg2rad
                              
!       If regional, map obs lat,lon to rotated grid.
           if(regional)then
              call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)
              if(diagnostic_reg) then
                 call txy2ll(dlon,dlat,rlon00,rlat00)
                 ntest=ntest+1
                 cdist=sin(dlat_earth)*sin(rlat00)+cos(dlat_earth)*cos(rlat00)* &
                       (sin(dlon_earth)*sin(rlon00)+cos(dlon_earth)*cos(rlon00))
                 cdist=max(-one,min(cdist,one))
                 disterr=acos(cdist)*rad2deg
                 disterrmax=max(disterrmax,disterr)
              end if
              if(outside) cycle loop_readsb ! check to see if outside regional 
           else
              dlon=dlon_earth
              dlat=dlat_earth
              call grdcrd1(dlat,rlats,nlat,1)
              call grdcrd1(dlon,rlons,nlon,1)
           endif

!     If ASCAT data, determine primary surface type.  If not open sea,
!     skip this observation.  This check must be done before thinning.
!     isflg    - surface flag 0:sea 1:land 2:sea ice 3:snow 4:mixed

           if (itype==290) then                              
              call deter_sfc_type(dlat_earth,dlon_earth,t4dv,isflg,tsavg)
              if (isflg /= 0) cycle loop_readsb
              if (tsavg <= 273.0_r_kind) cycle loop_readsb
           endif

       
!!    convert from wind direction and speed to u,v component
           !uob=-obsdat(2)*sin(obsdat(1)*deg2rad)
           !vob=-obsdat(2)*cos(obsdat(1)*deg2rad)
               ! no optimal wind vector in hu3_ascat BUFR file,
               ! we will select one from the available ambiguities
               ! based on assigned likelihood, below

           err_spd=obsdat(1) ! DBID: contains wind direction error (deg)
           err_dir=obsdat(2) ! EESM: contains wind speed error (m/s)
           ! Convert "optimal" 10m speed/dir to u-, v-wind
           uob=-obsdat(3)*sin(obsdat(4)*deg2rad)
           vob=-obsdat(3)*cos(obsdat(4)*deg2rad)
           ! Pull quality-flag
           qc1=obsdat(5)
           ! Ob errors are provided in speed/direction space. To model these
           ! errors in u/v component space, we will use the following logic:
           ! 1. Assume the wind speed is a uniform distribution between
           !    MWS10-0.5*err_spd and MWS10+0.5*err_spd
           ! 2. Assume the wind direction is a uniform distribution between
           !    MWD10-0.5*err_dir and MWD10+0.5*err_dir
           ! 3. Therefore, we can use trigonometric identities:
           !     3a. sin(a +- b) = sin(a)cos(b) +- cos(a)sin(b)
           !     3b. cos(a +- b) = cos(a)cos(b) -+ sin(a)sin(b)
           !    To define the range of u/v components by four possible values
           !    based on the error bounds as:
           ! 4. uob_1 = -(obsdat(3)+0.5*err_spd)*(
           !                                      sin(obsdat(4)*deg2rad)*cos(err_dir*deg2rad)+ 
           !                                      cos(obsdat(4)*deg2rad)*sin(err_dir*deg2rad)
           !                                     )
           !    uob_2 = -(obsdat(3)+0.5*err_spd)*(
           !                                      sin(obsdat(4)*deg2rad)*cos(err_dir*deg2rad)-
           !                                      cos(obsdat(4)*deg2rad)*sin(err_dir*deg2rad)
           !                                     )
           !    uob_3 = -(obsdat(3)-0.5*err_spd)*(
           !                                      sin(obsdat(4)*deg2rad)*cos(err_dir*deg2rad)+
           !                                      cos(obsdat(4)*deg2rad)*sin(err_dir*deg2rad)
           !                                     )
           !    uob_4 = -(obsdat(3)-0.5*err_spd)*(
           !                                      sin(obsdat(4)*deg2rad)*cos(err_dir*deg2rad)-
           !                                      cos(obsdat(4)*deg2rad)*sin(err_dir*deg2rad)
           !                                     )                                      
           !    vob_1 = -(obsdat(3)+0.5*err_spd)*(
           !                                      cos(obsdat(4)*deg2rad)*cos(err_dir*deg2rad)- 
           !                                      sin(obsdat(4)*deg2rad)*sin(err_dir*deg2rad)
           !                                     )
           !    vob_2 = -(obsdat(3)+0.5*err_spd)*(
           !                                      cos(obsdat(4)*deg2rad)*cos(err_dir*deg2rad)+
           !                                      sin(obsdat(4)*deg2rad)*sin(err_dir*deg2rad)
           !                                     )
           !    vob_3 = -(obsdat(3)-0.5*err_spd)*(
           !                                      cos(obsdat(4)*deg2rad)*cos(err_dir*deg2rad)-
           !                                      sin(obsdat(4)*deg2rad)*sin(err_dir*deg2rad)
           !                                     )
           !    vob_4 = -(obsdat(3)-0.5*err_spd)*(
           !                                      cos(obsdat(4)*deg2rad)*cos(err_dir*deg2rad)+
           !                                      sin(obsdat(4)*deg2rad)*sin(err_dir*deg2rad)
           !                                     )
           ! 5. The error-bounds on u and v are then computed as:
           !    uob_err = 0.5*(
           !                   max(uob_1,uob_2,uob_3,uob_4)-
           !                   min(uob_1,uob_2,uob_3,uob_4)
           !                  )
           !    vob_err = 0.5*(
           !                   max(vob_1,vob_2,vob_3,vob_4)-
           !                   min(vob_1,vob_2,vob_3,vob_4)
           !                  )
           !    This presupposes that the range of possible u/v values given the
           !    speed/direction errors is centered on the observed u/v, which is
           !    NOT necessarily a good assumption, but it's usually very close.
           ! 6. Since woe is taken as a single value to express error in either
           !    the u- or v-component, we calculate woe as a simple average of
           !    uob_err and vob_error: woe = 0.5*(uob_err+vob_err)
           ! 7. If either spd_err or dir_err is zero, then insufficient
           !    information is provided for this calculation and woe will revert to
           !    the error table.
           ! Compute 4 potential extremes for u,v
           uob_1=-(obsdat(3)+0.5*err_spd)*(sin(obsdat(4)*deg2rad)*cos(err_dir*deg2rad)+cos(obsdat(4)*deg2rad)*sin(err_dir*deg2rad))
           uob_2=-(obsdat(3)+0.5*err_spd)*(sin(obsdat(4)*deg2rad)*cos(err_dir*deg2rad)-cos(obsdat(4)*deg2rad)*sin(err_dir*deg2rad))
           uob_3=-(obsdat(3)-0.5*err_spd)*(sin(obsdat(4)*deg2rad)*cos(err_dir*deg2rad)+cos(obsdat(4)*deg2rad)*sin(err_dir*deg2rad))
           uob_4=-(obsdat(3)-0.5*err_spd)*(sin(obsdat(4)*deg2rad)*cos(err_dir*deg2rad)-cos(obsdat(4)*deg2rad)*sin(err_dir*deg2rad))
           vob_1=-(obsdat(3)+0.5*err_spd)*(cos(obsdat(4)*deg2rad)*cos(err_dir*deg2rad)-sin(obsdat(4)*deg2rad)*sin(err_dir*deg2rad))
           vob_2=-(obsdat(3)+0.5*err_spd)*(cos(obsdat(4)*deg2rad)*cos(err_dir*deg2rad)+sin(obsdat(4)*deg2rad)*sin(err_dir*deg2rad))
           vob_3=-(obsdat(3)-0.5*err_spd)*(cos(obsdat(4)*deg2rad)*cos(err_dir*deg2rad)-sin(obsdat(4)*deg2rad)*sin(err_dir*deg2rad))
           vob_4=-(obsdat(3)-0.5*err_spd)*(cos(obsdat(4)*deg2rad)*cos(err_dir*deg2rad)+sin(obsdat(4)*deg2rad)*sin(err_dir*deg2rad))
           ! load uob, vob values into vectors
           uob_vec(1)=uob_1
           uob_vec(2)=uob_2
           uob_vec(3)=uob_3
           uob_vec(4)=uob_4
           vob_vec(1)=vob_1
           vob_vec(2)=vob_2
           vob_vec(3)=vob_3
           vob_vec(4)=vob_4
           ! compute uob_rng and vob_rng based on min/max values in uob_vec,
           ! vob_vec
           uob_rng=maxval(uob_vec)-minval(uob_vec)
           vob_rng=maxval(vob_vec)-minval(vob_vec)
           ! compute estimated woe value based no half of the average of uob_rng
           ! and vob_rng
           woe_est=0.25*(uob_rng+vob_rng)
           write(6,*) 'BTH: selected wind uob,vob,woe_est=',uob,vob,woe_est           
!!  Get observation error from PREPBUFR observation error table
!   only need read the 4th column for type 290 from the right
 
           ppb=max(zero,min(ppb,r2000))
           itypey=itype
           if(njqc) then
              ierr=0
              do i =1,maxsub_uv
                 if( icsubtype(nc) == isuble_uv(itypey,i) ) then
                    ierr=i+1
                    exit
                 else if( i == maxsub_uv .and. icsubtype(nc) /= isuble_uv(itypey,i)) then
                    ncount=ncount+1
                    do j=1,maxsub_uv
                       if(isuble_uv(itypey,j) ==0 ) then
                          ierr=j+1
                          exit
                       endif
                    enddo
                    if (ncount ==1) then
                       write(6,*) 'READ_HU3ASCAT,WARNING!! cannot find subtyep in the error table,&
                                   itype,iobsub=',itypey,icsubtype
                       write(6,*) 'read error table at colomn subtype as 0,error table column=',j
                    endif
                 endif
              enddo
              if(ppb>=etabl_uv(itypey,1,1)) k1=1
              do kl=1,32
                 if(ppb>=etabl_uv(itypey,kl+1,1).and.ppb<=etabl_uv(itypey,kl,1)) k1=kl
              end do
              if(ppb<=etabl_uv(itypey,33,1)) k1=5
              k2=k1+1
              ediff = etabl_uv(itypey,k2,1)-etabl_uv(itypey,k1,1)
              if (abs(ediff) > tiny_r_kind) then
                 del = (ppb-etabl_uv(itypey,k1,1))/ediff
              else
                 del = huge_r_kind
              endif
              del=max(zero,min(del,one))
              obserr=(one-del)*etabl_uv(itypey,k1,ierr)+del*etabl_uv(itypey,k2,ierr)
              obserr=max(obserr,werrmin)
! get non linear qc parameter from b table
              var_jb=(one-del)*btabl_uv(itypey,k1,ierr)+del*btabl_uv(itypey,k2,ierr)
              var_jb=max(var_jb,wjbmin)
              if (var_jb >= 10.0_r_kind) var_jb=zero
          else
             if(ppb>=etabl(itype,1,1)) k1=1
              do kl=1,32
                 if(ppb>=etabl(itype,kl+1,1).and.ppb<=etabl(itype,kl,1)) k1=kl
              end do
              if(ppb<=etabl(itype,33,1)) k1=5
              k2=k1+1
              ediff = etabl(itype,k2,1)-etabl(itype,k1,1)
              if (abs(ediff) > tiny_r_kind) then
                 del = (ppb-etabl(itype,k1,1))/ediff
              else
                 del = huge_r_kind
              endif
              del=max(zero,min(del,one))
              obserr=(one-del)*etabl(itype,k1,4)+del*etabl(itype,k2,4)
              obserr=max(obserr,werrmin)
           endif

!         Set usage variable
           usage = 0 
           iuse=icuse(nc)
           if(iuse <= 0)usage=r100

! Get information from surface file necessary for conventional data here
! This is different from the previous sfc_type call
           call deter_sfc2(dlat_earth,dlon_earth,t4dv,idomsfc,tsavg,ff10,sfcr,zz)
 
!!    process the thining procedure
                
           ithin=ithin_conv(nc)
           ithinp = ithin > 0 .and. pflag /= 0
           if(ithinp   )then
!          Interpolate guess pressure profile to observation location
              klon1= int(dlon);  klat1= int(dlat)
              dx   = dlon-klon1; dy   = dlat-klat1
              dx1  = one-dx;     dy1  = one-dy
              w00=dx1*dy1; w10=dx1*dy; w01=dx*dy1; w11=dx*dy

              klat1=min(max(1,klat1),nlat); klon1=min(max(0,klon1),nlon)
              if (klon1==0) klon1=nlon
              klatp1=min(nlat,klat1+1); klonp1=klon1+1
              if (klonp1==nlon+1) klonp1=1
              do kk=1,nsig
                 presl(kk)=w00*prsl_full(klat1 ,klon1 ,kk) +  &
                           w10*prsl_full(klatp1,klon1 ,kk) + &
                           w01*prsl_full(klat1 ,klonp1,kk) + &
                           w11*prsl_full(klatp1,klonp1,kk)
              end do
           end if

! The following two lines are new in SATWND:
           dlnpob=log(one_tenth*ppb)  ! ln(pressure in cb)
           ppb=one_tenth*ppb         ! from mb to cb

 !         Special block for data thinning - if requested
           if (ithin > 0 .and. iuse >=0) then
              ntmp=ndata  ! counting moved to map3gridS

 !         Set data quality index for thinning
              if (thin4d) then
                 timedif = zero
              else
                 timedif=abs(t4dv-toff)
              endif

              crit1 = timedif/r6+half

              if (pflag==0) then
                 do kk=1,nsig
                    presl_thin(kk)=presl(kk)
                 end do
              endif
 
              call map3grids(-1,pflag,presl_thin,nlevp,dlat_earth,dlon_earth,&
                              ppb,crit1,ndata,iout,ntb,iiout,luse,.false.,.false.)

              if (.not. luse) cycle loop_readsb
              if(iiout > 0) isort(iiout)=0
              if (ndata > ntmp) then
                 nodata=nodata+2
              endif
              isort(ntb)=iout

           else
              ndata=ndata+1
              nodata=nodata+2
              iout=ndata
              isort(ntb)=iout
           endif
           ! if woe_est is non-zero and non-bmiss, use it for woe in place of
           ! obserr, otherwise use obserr
           if(woe_est > 0. .and. woe_est < r110)then 
              woe=woe_est
           else
              woe=obserr
           endif
           oelev=r10

           if(regional .and. .not. fv3_regional)then
              u0=uob
              v0=vob
              call rotate_wind_ll2xy(u0,v0,uob,vob,dlon_earth,dlon,dlat)
              if(diagnostic_reg) then
                 call rotate_wind_xy2ll(uob,vob,u00,v00,dlon_earth,dlon,dlat)
                 nvtest=nvtest+1
                 disterr=sqrt((u0-u00)**2+(v0-v00)**2)
                 vdisterrmax=max(vdisterrmax,disterr)
              end if
           endif
           write(6,*) 'BTH: Getting to cdata_all write-out'
! Output result to array, need restruct from here. 


           cdata_all(1,iout)=woe                  ! wind error
           cdata_all(2,iout)=dlon                 ! grid relative longitude
           cdata_all(3,iout)=dlat                 ! grid relative latitude
           cdata_all(4,iout)=dlnpob               ! ln(pressure in cb)
           cdata_all(5,iout)=oelev                 ! index of height         
           cdata_all(6,iout)=uob                  ! u obs
           cdata_all(7,iout)=vob                  ! v obs 
           cdata_all(8,iout)=ndata                ! station id 
           cdata_all(9,iout)=t4dv                 ! time
           cdata_all(10,iout)=nc                  ! index of type in convinfo file
           cdata_all(11,iout)=0                   ! index of station elevation
           cdata_all(12,iout)=qc1                 ! index of quality mark
           cdata_all(13,iout)=obserr              ! original obs error
           cdata_all(14,iout)=usage               ! usage parameter
           cdata_all(15,iout)=idomsfc             ! dominate surface type
           cdata_all(16,iout)=tsavg               ! skin temperature
           cdata_all(17,iout)=ff10                ! 10 meter wind factor
           cdata_all(18,iout)=sfcr                ! surface roughness
           cdata_all(19,iout)=dlon_earth_deg      ! earth relative longitude (degrees)
           cdata_all(20,iout)=dlat_earth_deg      ! earth relative latitude (degrees)
           cdata_all(21,iout)=zz                  ! terrain height at ob location
           cdata_all(22,iout)=r_prvstg(1,1)       ! provider name
           cdata_all(23,iout)=r_sprvstg(1,1)      ! subprovider name
           cdata_all(24,iout)=var_jb              ! non linear qc parameter

        enddo  loop_readsb

     enddo loop_msg

!    Close unit to bufr file
!    Deallocate arrays used for thinning data
     if (.not.use_all) then
        deallocate(presl_thin)
        call del3grids
     endif
! Normal exit

  enddo loop_convinfo! loops over convinfo entry matches
  call closbf(lunin)
  deallocate(lmsg,nrep,tab)
 

  ! Write header record and data to output file for further processing
  allocate(iloc(ndata))
  icount=0
  do i=1,maxobs
     if(isort(i) > 0)then
        icount=icount+1
        iloc(icount)=isort(i)
     end if
  end do
  if(ndata /= icount)then
     write(6,*) ' READ_HU3ASCAT: mix up in read_hu3ascat ,ndata,icount ',ndata,icount
     call stop2(49)
  end if

  allocate(cdata_out(nreal,ndata))
  do i=1,ndata
     itx=iloc(i)
     do k=1,nreal
        cdata_out(k,i)=cdata_all(k,itx)
     end do
  end do
  deallocate(iloc,isort,cdata_all)
!  deallocate(etabl)
  
  call count_obs(ndata,nreal,ilat,ilon,cdata_out,nobs)
  write(lunout) obstype,sis,nreal,nchanl,ilat,ilon
  write(lunout) cdata_out

  deallocate(cdata_out)
900 continue
  if(diagnostic_reg .and. ntest>0) write(6,*)'READ_HU3ASCAT:  ',&
       'ntest,disterrmax=',ntest,disterrmax
  if(diagnostic_reg .and. nvtest>0) write(6,*)'READ_HU3ASCAT:  ',&
       'nvtest,vdisterrmax=',ntest,vdisterrmax

  if (ndata == 0) then
     write(6,*)'READ_HU3ASCAT:  closbf(',lunin,') no data'
  endif
  
  write(6,*) 'READ_HU3ASCAT,nread,ndata,nreal,nodata=',nread,ndata,nreal,nodata

  close(lunin)

! End of routine
  return



end subroutine read_hu3ascat
