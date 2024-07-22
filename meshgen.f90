
subroutine meshgen()

    use global_variables

    write(3,*) 'entering meshgen'

    allocate(x(Nx+1),y(Ny+1))  !ghost cells included


    if (readmesh.eq.0) then
        do i=1,Nx+1
             x(i)=Lx*(i-2)/(Nx-2)
        enddo
        !x(1)=x(2)-Lx/(Nx-2)

        do j=1,Ny+1
             y(j)=Ly*(j-2)/(Ny-2)
        enddo    

        write(3,*) 'writing mesh files'
        open(12,file='xgrid.dat')
        do i=2,Nx
          write(12,*) i-1, x(i)
        enddo
        open(13,file='ygrid.dat')
        do j=2,Ny
             write(13,*) j-1, y(j)  !dont write ghhost cell
         enddo
         close(12)
         close(13)
    elseif (readmesh.eq.1) then  !make sure this is correct later
             open(12,file='xgrid.dat')
               do i=2,Nx
                    read(12,*) dummy, x(i)    !reading only interior cells
               enddo
             open(13,file='ygrid.dat')
               do j=2,Ny
                    read(13,*) dummmy, y(j)
               enddo
            close(12)
            close(13)
            x(1)=x(2)-(x(3)-x(2))
            y(1)=y(2)-(y(3)-y(2))
            x(Nx+1)=x(Nx)+(x(Nx)-x(Nx-1))
            y(Ny+1)=y(Ny)+(y(Ny)-y(Ny-1))
    endif

     !cell centers
    write(3,*) 'Computing cell centers'
    allocate(xc(Nx),yc(Ny))
    do i=1,Nx
        xc(i)=0.5*(x(i)+x(i+1))
    enddo
    do i=1,Ny
        yc(i)=0.5*(y(i)+y(i+1))
    enddo

    !compute cell size
    write(3,*) 'Computing grid size'
    allocate(dx(Nx),dy(Ny))
    do i=1,Nx
        dx(i)=x(i+1)-x(i)
    enddo
    do i=1,Ny
        dy(i)=y(i+1)-y(i)
    enddo
     dx(1)=0
     dx(nx)=0
     dy(1)=0
     dy(ny)=0 !outside cells are not used


     !circle geometry
    write(3,*) 'Finding iBlank'
    allocate(iblank(Nx,Ny))  !ghost cells
    iblank(:,:)=1
    do i=1,Nx
        do j=1,Ny
            if( (xc(i)-xcntr)**2+(yc(j)-ycntr)**2 <= 0.5**2 ) iblank(i,j)=0    !should be 0 for inside boundary  
        enddo
    enddo

    write(3,*) 'Finding iBlank boundary'
    allocate(ghost(Nx,Ny))
    ghost(:,:)=0
    do i=2,Nx-1     !note that boundary cant be on the edge
        do j=2,Ny-1
            if (iblank(i,j).eq.0) then     !only need cell inside boundary
                if (abs(iBlank(i,j)-iblank(i+1,j)).eq.1.or.&
                    abs(iBlank(i,j)-iblank(i-1,j)).eq.1.or.&
                    abs(iBlank(i,j)-iblank(i,j+1)).eq.1.or.&
                    abs(iBlank(i,j)-iblank(i,j-1)).eq.1) then
                        ghost(i,j)=1
                endif
            endif
        enddo
    enddo

    allocate(IP(Nx,Ny),JP(Nx,Ny),IM(Nx,Ny),JM(Nx,Ny))

    IP(:,:)=0
    JP(:,:)=0
    IM(:,:)=0
    JM(:,:)=0


    do i=2,Nx-1
        do j=2,Ny-1    !note that boundary cant be on the edge
            if (ghost(i,j).eq.0.and.ghost(i+1,j).eq.1) then
                IP(i,j)=1
            endif
            if (ghost(i,j).eq.0.and.ghost(i-1,j).eq.1) then
                IM(i,j)=1
            endif
            if (ghost(i,j).eq.0.and.ghost(i,j+1).eq.1) then
                JP(i,j)=1
            endif
            if (ghost(i,j).eq.0.and.ghost(i,j-1).eq.1) then
                JM(i,j)=1
            endif
        enddo
    enddo
    !write(*,*) 'ip'
    !do j=1,ny
    !    write(*,*) ip(:,j)
     ! enddo
      !write(*,*) 'im'
     ! do j=1,ny
      !  write(*,*) im(:,j)
      !enddo
      !write(*,*) 'jp'
      !do j=1,ny
      !  write(*,*) jp(:,j)
      !enddo
      !write(*,*) 'jm'
      !do j=1,ny
      !  write(*,*) jm(:,j)
      !enddo

    write(3,*) 'exiting meshgen'
allocate(lam_e(Nx),lam_w(Nx),lam_n(Ny),lam_s(Ny))
lam_e(:)=0
lam_w(:)=0
lam_n(:)=0
lam_s(:)=0

do i=2,nx-1
    lam_e(i)=(dx(i+1)/(dx(i)+dx(i+1)))
    lam_w(i)=(dx(i-1)/(dx(i)+dx(i-1)))
enddo

do j=2,ny-1
    lam_n(j)=(dy(j+1)/(dy(j)+dy(j+1)))
    lam_s(j)=(dy(j-1)/(dy(j)+dy(j-1)))
enddo

allocate(probe_dis(Nx,Ny))
probe_dis(:,:)=0

do i =2,nx-1
    do j=2,ny-1
        probe_dis(i,j)=sqrt((xc(i)-probex)**2+(yc(j)-probey)**2)
    enddo
enddo
probe_dummy=10
do j=2,ny-1
    do i=2,nx-1
        if (probe_dis(i,j).le.probe_dummy) then
            probe_dummy=probe_dis(i,j)
            probei=i
            probej=j
        endif
    enddo
enddo
write(34,*) 'probe_i', probei, 'probe_j', probej
write(35,*) 'probe_i', probei, 'probe_j', probej


    

end subroutine