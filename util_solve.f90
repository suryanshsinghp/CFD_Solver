
subroutine set_coeff()
  use global_variables
  allocate(adc1(nx,ny),adc2(nx,ny),adc3(nx,ny),adc4(nx,ny),adc5(nx,ny))
  allocate(adc6(nx,ny),adc7(nx,ny),adc8(nx,ny),adc9(nx,ny),adc10(nx,ny))
  allocate(ppc1(nx,ny),ppc2(nx,ny),ppc3(nx,ny),ppc4(nx,ny),ppc5(nx,ny))
  adc1(:,:)=0
  adc2(:,:)=0
  adc3(:,:)=0
  adc4(:,:)=0
  adc5(:,:)=0
  adc6(:,:)=0
  adc7(:,:)=0
  adc8(:,:)=0
  adc9(:,:)=0
  adc10(:,:)=0
  ppc1(:,:)=0
  ppc2(:,:)=0
  ppc3(:,:)=0
  ppc4(:,:)=0
  ppc5(:,:)=0



  do i=2,nx-1
    do j=2,ny-1
      adc1(i,j)=(2/(dx(i)*(dx(i)+dx(i+1))))
      adc2(i,j)=(2/(dx(i)*(dx(i)+dx(i-1))))
      adc4(i,j)=(2/(dy(j)*(dy(j)+dy(j-1))))
      adc5(i,j)=(2/(dy(j)*(dy(j)+dy(j+1))))
      !adc3(i,j)=1-( adc1(i,j)*(1-ip(i,j))+adc2(i,j)*(1-im(i,j))&  
       !           +adc4(i,j)*(1-jm(i,j))+adc5(i,j)*(1-jp(i,j)) +2*(dt/(2*Re))*ip(i,j)/dx(i)**2+2*(dt/(2*Re))*im(i,j)/dx(i)**2 &
       !                                                        +2*(dt/(2*Re))*jp(i,j)/dy(j)**2+2*(dt/(2*Re))*jm(i,j)/dy(j)**2 )
      adc3(i,j)=-((1-ip(i,j))*adc1(i,j)+(1-IM(i,j))*adc2(i,j)+(1-JM(i,j))*adc4(i,j)+ &
                 +(1-JP(i,j))*adc5(i,j)+&
                 2*(((IP(i,j)+IM(i,j))/(dx(i)**2))+((JP(i,j)+JM(i,j))/(dy(j)**2))))
      adc3(i,j)=1-((dt/(2*Re))*adc3(i,j))
      ppc3(i,j)=-((1-ip(i,j))*adc1(i,j)+(1-im(i,j))*adc2(i,j)+(1-jm(i,j))*adc4(i,j)+(1-jp(i,j))*adc5(i,j))
  
      adc1(i,j)=(-dt/(2*Re))*(adc1(i,j)*(1-ip(i,j)))*(iblank(i,j))   
      adc2(i,j)=(-dt/(2*Re))*(adc2(i,j)*(1-im(i,j)))*(iblank(i,j))  
      adc4(i,j)=(-dt/(2*Re))*(adc4(i,j)*(1-jm(i,j)))*(iblank(i,j))  
      adc5(i,j)=(-dt/(2*Re))*(adc5(i,j)*(1-jp(i,j)))*(iblank(i,j))  
    enddo
  enddo

  !do j=1,ny
    !write(*,*) adc3(:,j)
  !enddo
  
  do i=2,nx-1
    do j=2,ny-1
      adc6(i,j)=-adc1(i,j)
      adc7(i,j)=-adc2(i,j)
      adc9(i,j)=-adc4(i,j)
      adc10(i,j)=-adc5(i,j)
      adc8(i,j)=2-adc3(i,j)
    enddo
  enddo

  ppc1(:,:)=(adc1(:,:)*2*Re)/dt    !check this value if PPE does not converge
  ppc2(:,:)=(adc2(:,:)*2*Re)/dt
  ppc4(:,:)=(adc4(:,:)*2*Re)/dt
  ppc5(:,:)=(adc5(:,:)*2*Re)/dt
  !do i=2,nx-1
  !  do j=2,ny-1
  !    ppc3(i,j)=-ppc1(i,j)-ppc2(i,j)-ppc4(i,j)-ppc5(i,j)
   !   ppc3(i,j)=(ppc3(i,j)*2*Re)/(dt)
   ! enddo
  !enddo

  !do j=1,ny
  !  write(*,*) ppc5(:,j)
  !enddo

end subroutine

subroutine ADS()
    use global_variables
  !interpolate to get face velocities at n and n-1
    do i=2,nx-1
      do j=2,ny-1
        uen(i,j)=lam_e(i)*un(i,j)+(1-lam_e(i))*un(i+1,j)
        uwn(i,j)=lam_w(i)*un(i,j)+(1-lam_w(i))*un(i-1,j)
        unn(i,j)=lam_n(j)*un(i,j)+(1-lam_n(j))*un(i,j+1)
        usn(i,j)=lam_s(j)*un(i,j)+(1-lam_s(j))*un(i,j-1)
        uen(i,j)=uen(i,j)*(1-ip(i,j))
        uwn(i,j)=uwn(i,j)*(1-im(i,j))
        unn(i,j)=unn(i,j)*(1-jp(i,j))
        usn(i,j)=usn(i,j)*(1-jm(i,j))

        ven(i,j)=lam_e(i)*vn(i,j)+(1-lam_e(i))*vn(i+1,j)
        vwn(i,j)=lam_w(i)*vn(i,j)+(1-lam_w(i))*vn(i-1,j)
        vnn(i,j)=lam_n(j)*vn(i,j)+(1-lam_n(j))*vn(i,j+1)
        vsn(i,j)=lam_s(j)*vn(i,j)+(1-lam_s(j))*vn(i,j-1)
        ven(i,j)=ven(i,j)*(1-ip(i,j))
        vwn(i,j)=vwn(i,j)*(1-im(i,j))
        vnn(i,j)=vnn(i,j)*(1-jp(i,j))
        vsn(i,j)=vsn(i,j)*(1-jm(i,j))

        uenm1(i,j)=lam_e(i)*unm1(i,j)+(1-lam_e(i))*unm1(i+1,j)
        uwnm1(i,j)=lam_w(i)*unm1(i,j)+(1-lam_w(i))*unm1(i-1,j)
        unnm1(i,j)=lam_n(j)*unm1(i,j)+(1-lam_n(j))*unm1(i,j+1)
        usnm1(i,j)=lam_s(j)*unm1(i,j)+(1-lam_s(j))*unm1(i,j-1)
        uenm1(i,j)=uenm1(i,j)*(1-ip(i,j))
        uwnm1(i,j)=uwnm1(i,j)*(1-im(i,j))
        unnm1(i,j)=unnm1(i,j)*(1-jp(i,j))
        usnm1(i,j)=usnm1(i,j)*(1-jm(i,j))

        venm1(i,j)=lam_e(i)*vnm1(i,j)+(1-lam_e(i))*vnm1(i+1,j)
        vwnm1(i,j)=lam_w(i)*vnm1(i,j)+(1-lam_w(i))*vnm1(i-1,j)
        vnnm1(i,j)=lam_n(j)*vnm1(i,j)+(1-lam_n(j))*vnm1(i,j+1)
        vsnm1(i,j)=lam_s(j)*vnm1(i,j)+(1-lam_s(j))*vnm1(i,j-1)
        venm1(i,j)=venm1(i,j)*(1-ip(i,j))
        vwnm1(i,j)=vwnm1(i,j)*(1-im(i,j))
        vnnm1(i,j)=vnnm1(i,j)*(1-jp(i,j))
        vsnm1(i,j)=vsnm1(i,j)*(1-jm(i,j))
        
        
      enddo
    enddo
    !do j=1,ny
    !write(*,*) venm1(:,j)
    !enddo 

    !do j=1,ny
    !write(*,*) uen(:,j)
    !enddo 

    !do j=1,ny
    !write(*,*) ven(:,j)
    !enddo 

    
    do i=2,nx-1
      do j=2,ny-1
        source_u(i,j)= adc6(i,j)*un(i+1,j)+adc7(i,j)*un(i-1,j)+adc8(i,j)*un(i,j)+adc9(i,j)*un(i,j-1)+adc10(i,j)*un(i,j+1) &
                          -dt*( (3/(2*dx(i)))*(uen(i,j)*uf(i,j-1)-uwn(i,j)*uf(i-1,j-1)) &
                              + (3/(2*dy(j)))*(unn(i,j)*vf(i-1,j)-usn(i,j)*vf(i-1,j-1))  &
                              - (1/(2*dx(i)))*(uenm1(i,j)*ufnm1(i,j-1)-uwnm1(i,j)*ufnm1(i-1,j-1)) &
                              - (1/(2*dy(j)))*(unnm1(i,j)*vfnm1(i-1,j)-usnm1(i,j)*vfnm1(i-1,j-1)) )
        source_v(i,j)= adc6(i,j)*vn(i+1,j)+adc7(i,j)*vn(i-1,j)+adc8(i,j)*vn(i,j)+adc9(i,j)*vn(i,j-1)+adc10(i,j)*vn(i,j+1) &
                          -dt*( (3/(2*dx(i)))*(ven(i,j)*uf(i,j-1)-vwn(i,j)*uf(i-1,j-1)) &
                              + (3/(2*dy(j)))*(vnn(i,j)*vf(i-1,j)-vsn(i,j)*vf(i-1,j-1))  &
                              - (1/(2*dx(i)))*(venm1(i,j)*ufnm1(i,j-1)-vwnm1(i,j)*ufnm1(i-1,j-1)) &
                              - (1/(2*dy(j)))*(vnnm1(i,j)*vfnm1(i-1,j)-vsnm1(i,j)*vfnm1(i-1,j-1)) )        
        source_u(i,j)=source_u(i,j)*iBlank(i,j)
        source_v(i,j)=source_v(i,j)*iBlank(i,j)

      enddo
    enddo
    !do j=1,ny
    !write(*,*) source_u(:,j)
    !enddo 

    !do j=1,ny
    !write(*,*) source_v(:,j)
    !enddo 

 !gauss sidel
    ADS_error=1
    us(:,:)=un(:,:)
    vs(:,:)=vn(:,:)   !initialization
    us_old(:,:)=us(:,:)
      vs_old(:,:)=vs(:,:)   
    do while (ADS_error.gt.ADS_tol)
      
      do j=2,ny-1
        do i=2,nx-1
          us(i,j)=omega_ad*( -adc1(i,j)*us(i+1,j)-adc2(i,j)*us(i-1,j) &
                          -adc4(i,j)*us(i,j-1)-adc5(i,j)*us(i,j+1)+source_u(i,j))
          us(i,j)=us(i,j)/adc3(i,j)
          us(i,j)= (1-omega_ad)*us_old(i,j)+us(i,j)
          vs(i,j)=omega_ad*( -adc1(i,j)*vs(i+1,j)-adc2(i,j)*vs(i-1,j) &
                          -adc4(i,j)*vs(i,j-1)-adc5(i,j)*vs(i,j+1)+source_v(i,j))
          vs(i,j)=vs(i,j)/adc3(i,j)
          vs(i,j)= (1-omega_ad)*vs_old(i,j)+vs(i,j)
        enddo
      enddo


      !BC   
      us(nx,:)=us(nx-1,:)
      us(:,ny)=us(:,ny-1)
      us(:,1)=us(:,2)
      vs(nx,:)=vs(nx-1,:)
      us_old(:,:)=us(:,:)
      vs_old(:,:)=vs(:,:)   
      res_x=0
      res_y=0
      do i=2,nx-1
        do j=2,ny-1
          res_x=res_x+abs(source_u(i,j)-adc1(i,j)*us(i+1,j)-adc2(i,j)*us(i-1,j)-adc3(i,j)*us(i,j)-&
                           adc4(i,j)*us(i,j-1)-adc5(i,j)*us(i,j+1))
          res_y=res_y+abs(source_v(i,j)-adc1(i,j)*vs(i+1,j)-adc2(i,j)*vs(i-1,j)-adc3(i,j)*vs(i,j)-&
                          adc4(i,j)*vs(i,j-1)-adc5(i,j)*vs(i,j+1))
        enddo
      enddo
      res_x=res_x/((nx-2)*(ny-2))
      res_y=res_y/((nx-2)*(ny-2))
      ADS_error=max(res_x,res_y)
      !write(*,*) 'ads error',ADS_error
      giter_ADE=giter_ADE+1
      write(32,*) giter_ADE,ADS_error
    enddo
 

    !do j=1,ny
     ! write(*,*) us(:,j)
      !enddo 
    !correct mass
    m_in=0
    m_out=0
    do j=2,ny-1
      m_in=m_in+us(1,j)*dy(j)
      m_out=m_out+us(nx-1,j)*dy(j)
    enddo
    m_corr=(m_out-m_in)/Ly
    us(nx-1,:)=us(nx-1,:)-m_corr
    us(nx,:)=us(nx-1,:)   !du/dx=0
  !do j=1,ny
    ! write(*,*) us(:,j)
      !enddo 

     !do j=1,ny
    ! write(*,*) vs(:,j)
      !enddo 
  !face star velocity
    do i=2,nx-1
      do j=2,ny-1
        ufs(i,j-1)=(1-ip(i,j))*lam_e(i)*us(i,j)+(1-im(i,j))*(1-lam_e(i))*us(i+1,j)
        vfs(i-1,j)=(1-jp(i,j))*lam_n(j)*vs(i,j)+(1-jm(i,j))*(1-lam_n(j))*vs(i,j+1)
      enddo
    enddo

     !do j=1,ny
    ! write(*,*) ufs(:,j)
      !enddo 

     !do j=1,ny
    ! write(*,*) vfs(:,j)
      !enddo 

end subroutine


subroutine PPE()
    use global_variables
    do i=2,nx-1
      do j=2,ny-1
        !source_p(i,j)=(ufs(i,j-1)-ufs(i-1,j-1))/dx(i) + (vfs(i-1,j)-vfs(i-1,j-1))/dy(j)  !!!!!!!!!!!!!!!!!!!!CHECK THIS!!!!!!!
        source_p(i,j)=(((1-IP(i,j))*(lam_e(i)*us(i,j)+(1-lam_e(i))*us(i+1,j))) &
                      -((1-IM(i,j))*(lam_w(i)*us(i,j)+(1-lam_w(i))*us(i-1,j))))/(dx(i))
        source_p(i,j)=source_p(i,j)+(((1-JP(i,j))*(lam_n(j)*vs(i,j)+(1-lam_n(j))*vs(i,j+1))) &
                      -((1-JM(i,j))*(lam_s(j)*vs(i,j)+(1-lam_s(j))*vs(i,j-1))))/(dy(j))
        source_p(i,j)=source_p(i,j)/dt
        source_p(i,j)=source_p(i,j)*iblank(i,j)
      enddo
   enddo
   !do j=1,ny
     !write(*,*) source_p(:,j)
   !enddo 
   PPE_error=1
   pold(:,:)=pnp1(:,:)
   do while (PPE_error.gt.PPE_tol)    !!!!!!!!!
  
   
    do i=2,nx-1
      do j=2,ny-1
        pnp1(i,j)=omega_pp*(ppc1(i,j)*pnp1(i+1,j)+ppc2(i,j)*pnp1(i-1,j) +&
                             ppc4(i,j)*pnp1(i,j-1)+ppc5(i,j)*pnp1(i,j+1)+source_p(i,j))
        pnp1(i,j)=pnp1(i,j)/ppc3(i,j)
        pnp1(i,j)=pnp1(i,j)+(1-omega_pp)*pold(i,j)
      enddo
    enddo

    pnp1(1,:)=pnp1(2,:)
    pnp1(:,1)=pnp1(:,2)
    pnp1(nx,:)=pnp1(nx-1,:)
    pnp1(:,ny)=pnp1(:,ny-1)
    pold(:,:)=pnp1(:,:)
    PPE_error=0
    do i=2,nx-1
      do j=2,ny-1
        PPE_error=PPE_error+abs(ppc3(i,j)*pnp1(i,j)-ppc1(i,j)*pnp1(i+1,j)-ppc2(i,j)*pnp1(i-1,j)&
                                                   -ppc4(i,j)*pnp1(i,j-1)-ppc5(i,j)*pnp1(i,j+1)-source_p(i,j))
      enddo
    enddo
    PPE_error=PPE_error/((nx-2)*(ny-2))
    giter_PPE=giter_PPE+1
      write(33,*) giter_PPE,PPE_error
    !write(*,*) 'PPE_error',PPE_error
    !write(*,*) 'ppc1'
   ! do j=1,ny
      !write(*,*) ppc1(:,j)
    !enddo 
   
   

   enddo !!!!!!!!!!!!
   
    
     

end subroutine


subroutine velocity_update()
    use global_variables

    do i=2,nx-1
      do j=2,ny-1
        pfe(i,j)=lam_e(i)*pnp1(i,j)+(1-lam_e(i))*pnp1(i+1,j)      
        pfw(i,j)=lam_w(i)*pnp1(i,j)+(1-lam_w(i))*pnp1(i-1,j)       
        pfn(i,j)=lam_n(j)*pnp1(i,j)+(1-lam_n(j))*pnp1(i,j+1)      
        pfs(i,j)=lam_s(j)*pnp1(i,j)+(1-lam_s(j))*pnp1(i,j-1)
        pfe(i,j)=(1-ip(i,j))*pfe(i,j)+ip(i,j)*pnp1(i,j)
        pfw(i,j)=(1-im(i,j))*pfw(i,j)+im(i,j)*pnp1(i,j)
        pfn(i,j)=(1-jp(i,j))*pfn(i,j)+jp(i,j)*pnp1(i,j)
        pfs(i,j)=(1-jm(i,j))*pfs(i,j)+jm(i,j)*pnp1(i,j)
      enddo
    enddo
    
    do i=2,nx-1
      do j=2,ny-1
        unp1(i,j)=us(i,j)-(dt/dx(i))*(pfe(i,j)-pfw(i,j))
        vnp1(i,j)=vs(i,j)-(dt/dy(j))*(pfn(i,j)-pfs(i,j))
        unp1(i,j)=unp1(i,j)*iblank(i,j)
        vnp1(i,j)=vnp1(i,j)*iblank(i,j)
      enddo
    enddo
    unp1(nx,:)=unp1(nx-1,:)
    unp1(:,1)=unp1(:,2)
    unp1(:,ny)=unp1(:,ny-1)
    vnp1(nx,:)=vnp1(nx-1,:)
 
    
    do i=1,nx-1
      do j=2,ny-1
        ufnp1(i,j-1)=ufs(i,j-1)-(((1-im(i+1,j))*pnp1(i+1,j)-(1-ip(i,j))*pnp1(i,j))/(0.5*(dx(i)+dx(i+1))))*dt
      enddo
    enddo

    do i=2,nx-1
      do j=2,ny-2
        vfnp1(i-1,j)=vfs(i-1,j)-(((1-jm(i,j+1))*pnp1(i,j+1)-(1-jp(i,j))*pnp1(i,j))/(0.5*(dy(j)+dy(j+1))))*dt
      enddo
    enddo



end subroutine

subroutine drag_lift()
  use global_variables
  f_drag_left=0
  f_drag_right=0
  f_lift_top=0
  f_lift_bottom=0
  f_drag=0
  f_lift=0
  do i=2,nx-1
    do j=2,ny-1
      f_drag_left=f_drag_left+(pnp1(i,j)*(1-iblank(i+1,j)))*dy(j)
      f_drag_right=f_drag_right+(pnp1(i,j)*(1-iblank(i-1,j)))*dy(j)
      f_drag=2*(f_drag_left-f_drag_right)   !drag coeff factor=2
      f_lift_top=f_lift_top+(pnp1(i,j)*(1-iblank(i,j-1)))*dx(i)
      f_lift_bottom=f_lift_bottom+(pnp1(i,j)*(1-iblank(i,j+1)))*dx(i)
      f_lift=2*(f_lift_bottom-f_lift_top)   !lift coeff factor=2
    enddo
  enddo
  write(31,*) t,f_drag,f_lift!,f_drag_left,f_drag_right,f_lift_top,f_lift_bottom


End subroutine

subroutine surface_pressure()

  use global_variables
  do i=2,nx-1
    do j=2,ny-1
      f_drag_left=(pnp1(i,j)*(1-iblank(i+1,j)))*dy(j)
      f_drag_right=(pnp1(i,j)*(1-iblank(i-1,j)))*dy(j)
      angle=asin(abs((yc(j)-ycntr)/(sqrt((xc(i)-xcntr)**2+(yc(j)-ycntr)**2))))
      if (xc(i).lt.xcntr.and.yc(j).lt.ycntr) then
        angle=pi+angle
      elseif (xc(i).gt.xcntr.and.yc(j).lt.ycntr) then
        angle=2*pi-angle
      elseif (xc(i).lt.xcntr.and.yc(j).gt.ycntr) then
        angle=pi-angle
      elseif (xc(i).gt.xcntr.and.yc(j).gt.ycntr) then
        angle=angle
      end if
      if (f_drag_left.ne.0) then
        write(36,*) angle,f_drag_left
      endif
      if (f_drag_right.ne.0) then
        write(36,*) angle,f_drag_right
      endif

    enddo
  enddo



  end subroutine  