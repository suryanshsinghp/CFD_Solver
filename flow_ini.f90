subroutine flow_ini()
    use global_variables


    allocate(un(nx,ny),vn(nx,ny),unp1(nx,ny),vnp1(nx,ny),unm1(nx,ny),vnm1(nx,ny))
    allocate(us(nx,ny),vs(nx,ny),us_old(nx,ny),vs_old(nx,ny))
    allocate(uf(nx-1,ny-2),vf(nx-2,ny-1),ufnp1(nx-1,ny-2),vfnp1(nx-2,ny-1))
    allocate(ufs(nx-1,ny-2),vfs(nx-2,ny-1),ufnm1(nx-1,ny-2),vfnm1(nx-2,ny-1))
    allocate(pold(nx,ny),pnp1(nx,ny))
    allocate(source_p(nx,ny))
    allocate(uen(nx,ny),uwn(nx,ny),unn(nx,ny),usn(nx,ny))
    allocate(ven(nx,ny),vwn(nx,ny),vnn(nx,ny),vsn(nx,ny))
    allocate(uenm1(nx,ny),uwnm1(nx,ny),unnm1(nx,ny),usnm1(nx,ny))
    allocate(venm1(nx,ny),vwnm1(nx,ny),vnnm1(nx,ny),vsnm1(nx,ny))
    allocate(source_u(nx,ny),source_v(nx,ny))
    allocate(pfe(nx,ny),pfw(nx,ny),pfn(nx,ny),pfs(nx,ny))


    un(:,:)=0
    vn(:,:)=0
    unp1(:,:)=0
    vnp1(:,:)=0
    uf(:,:)=0
    vf(:,:)=0
    ufnp1(:,:)=0
    vfnp1(:,:)=0
    ufs(:,:)=0
    vfs(:,:)=0
    pold(:,:)=0
    pnp1(:,:)=0
    source_p(:,:)=0
    us(:,:)=0
    vs(:,:)=0
    us_old(:,:)=0
    vs_old(:,:)=0
    ufnm1(:,:)=0
    vfnm1(:,:)=0
    unm1(:,:)=0
    vnm1(:,:)=0
    pfe(:,:)=0
    pfw(:,:)=0
    pfn(:,:)=0
    pfs(:,:)=0



    !BC

    uf(1,:)=1
    ufs(1,:)=1    !check this
    ufnm1(1,:)=1
    ufnp1(1,:)=1
    un(1,:)=1     !dx=0
    us(1,:)=1
    unm1(1,:)=1
    unp1(1,:)=1

    uen(:,:)=0
    uwn(:,:)=0
    unn(:,:)=0
    usn(:,:)=0
    ven(:,:)=0
    vwn(:,:)=0
    vnn(:,:)=0
    vsn(:,:)=0
    uenm1(:,:)=0
    uwnm1(:,:)=0
    unnm1(:,:)=0
    usnm1(:,:)=0
    venm1(:,:)=0
    vwnm1(:,:)=0
    vnnm1(:,:)=0
    vsnm1(:,:)=0

    source_u(:,:)=0
    source_v(:,:)=0

    giter_ADE=0
    giter_PPE=0   !global iteration counter






end subroutine 