
    !run using 
    !gfortran meshgen.f90 flow_ini.f90 util_solve.f90 util_general.f90 && ./a.out
    module global_variables
    integer ::Nx,Ny,readmesh,Tstep,Tsave,giter_ADE,giter_PPE,probei,probej
    real*8    ::Lx,Ly,dummy,dt,t,tEnd,Re,ADS_tol,PPE_tol,ADS_error,PPE_error,coff_p,coff_u,uerr,verr,su,sv,m_in,m_out,m_corr
    real*8    ::omega_pp,omega_ad,res_u,res_v,res_p,probex,probey,probe_dummy,xcntr,ycntr,angle,pi
    real*8 :: f_drag, f_lift,f_drag_left,f_drag_right,f_lift_top,f_lift_bottom
    real*8, dimension(:), allocatable :: x,y,xc,yc,dx,dy,lam_e,lam_w,lam_n,lam_s
    integer, dimension(:,:), allocatable :: iblank,ghost,IP,IM,JP,JM
    real*8, dimension(:,:), allocatable ::adc1,adc2,adc3,adc4,adc5,adc6,adc7,adc8,adc9,adc10,ppc1,ppc2,ppc3,ppc4,ppc5
    real*8, dimension(:,:), allocatable :: un,vn,unp1,vnp1,uf,vf,ufnp1,vfnp1,pold,pnp1,source_p,ufs,vfs,us,vs,us_old,vs_old
    real*8, dimension(:,:), allocatable :: uen,uwn,usn,unn,ven,vwn,vsn,vnn,uenm1,uwnm1,usnm1,unnm1,venm1,vwnm1,vsnm1,vnnm1
    real*8, dimension(:,:), allocatable :: unm1,vnm1,ufnm1,vfnm1,source_u,source_v, pfe,pfw,pfn,pfs,probe_dis
    character(len=64) :: fname
    CHARACTER *500 :: cmd,dir
    LOGICAL :: L_EXISTS
    end module global_variables

    

program project

    use global_variables
    open(3,file='log.dat')
    write(3,*) 'Starting simulation'
    write(3,*) 'Reading input file'
    !read input
    open(11,file='input.dat')
    read(11,*) 
    read(11,*) Lx,Ly,Nx,Ny,xcntr,ycntr
    read(11,*) 
    read(11,*) readmesh
    read(11,*)
    read(11,*) dt,tEnd,Tsave
    read(11,*)
    read(11,*) Re, omega_ad,omega_pp
    read(11,*)
    read(11,*) ADS_tol,PPE_tol
    read(11,*)
    read(11,*) probex,probey

   
  PI=4.D0*DATAN(1.D0)
   t=0
   Tstep=0

   Nx=Nx+2
   Ny=Ny+2 !for ghost cells
   open(31,file='drag_lift.dat')
   write(31,*) 't','drag','lift'
   open(32,file='ADE_residual.dat')
   open (33,file='PPE_residual.dat')
    open (34,file='velprobe.dat')
    open (35,file='pressprobe.dat')
    open (36,file='surfacepress.dat')

    !generate mesh files
    call meshgen()

  

    call flow_ini()
    call set_coeff()
    
    
    !call flow_save()

    do while (t<tEnd)
        write(3,*) 'Time step',t
        call ADS()

        call PPE()

        call velocity_update()
        write(34,*) t,unp1(probei,probej),vnp1(probei,probej)
        write(35,*) t, pnp1(probei,probej)
      
        call drag_lift()
       

        if (mod(Tstep,Tsave).eq.0) call flow_save()



        unm1(:,:)=un(:,:)
        vnm1(:,:)=vn(:,:)
        un(:,:)=unp1(:,:)
        vn(:,:)=vnp1(:,:)
        ufnm1(:,:)=uf(:,:)
        vfnm1(:,:)=vf(:,:)    !is uf at current time?
        uf(:,:)=ufnp1(:,:)
        vf(:,:)=vfnp1(:,:)
        
        
        
    

        t=t+dt
        write(*,*) 'time is',t
        Tstep=Tstep+1


    enddo
    call surface_pressure()


    close(3)
    close(31)  !drag_lift.dat file
    close(32)  !ADE_residual.dat file
    close(33)  !PPE_residual.dat file
    close(34)  !probe.dat file
    close(35)  !presssurf.dat file
end program 


subroutine flow_save()
    use global_variables
    write(3,*) 'Writing Tecplot output'
    WRITE(dir, '("./Qfiles/")')
    !    INQUIRE(DIRECTORY= trim(dir), EXIST=L_EXISTS)      !manually make dir if code gives error
    !        IF(.NOT.L_EXISTS) THEN
    !             WRITE(cmd, '(" mkdir -p   ",A)') trim(dir)
    !             call system(trim(cmd))
    !        END IF                                          !!!!!
    WRITE(fname,"(A,'q.',I7.7,'.dat')") trim(dir),Tstep
    open(15, file=fname, status='unknown')
    WRITE(15,*) 'TITLE = "Post Processing Tecplot"'
    WRITE(15,*) 'VARIABLES = "X", "Y", "iblank", "u", "v", "pressure"'
    WRITE(15,*) 'ZONE T="BIG ZONE", I=',nx-2,', J=',ny-2,', DATAPACKING=POINT'

    DO j=2,ny-1
     DO i = 2,nx-1   !dont save ghost cells
         WRITE(15,*) xc(i), yc(j), iblank(i,j), unp1(i,j), vnp1(i,j), Pnp1(i,j)
     END DO
    END DO
    close(15)
    end subroutine
