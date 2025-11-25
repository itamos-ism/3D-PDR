!Written by Zhengping Zhu
module m_Ray_box
    use m_Mesh
    use healpix_types
    use maincode_module,only:pdr,nrays,maxpoints
    implicit none
    private

    type, public :: box
        real(RK) :: min(3), max(3)
    end type box

    type :: slab
        real(RK) :: origin(3),  dir(3),  dir_inv(3) 
    end type slab

    type, public :: HEALPix_ray
        integer :: eval
        real(RK) :: length 
        real(RK) :: origin(3), angle(2)
    end type HEALPix_ray

    integer,public :: levels
    integer,public,allocatable :: maxpoints_ray(:)
    type(box),public :: box1
    type(HEALPix_ray),public :: ray
    real(kind=dp),public :: thfpix, phfpix, contribution
    real(kind=dp),public :: corner_min(3), corner_max(3)

#ifdef RAYTHEIA_MO
    public:: raytheia, raytheia_maxpoints
#else
    public:: raytheia_table
#endif

contains
    logical function isfinite(x)
            real(RK), intent(in) :: x
            isfinite = (abs(x) < huge(x) .and. x == x) 
    end function isfinite

    subroutine intersections(ray_xyz, box_in, length)
    implicit none
    type(slab), intent(in) :: ray_xyz
    type(box), intent(in) :: box_in
    real(RK), intent(inout) :: length
    integer :: i, d
    real(RK) :: tmin, tmax, t1, t2
    real(RK) :: Pmin(3), Pmax(3)
    real(RK) :: distance, diff
    real(RK) :: temp

    ! 初始化 tmin 和 tmax
    tmin = 0.0
    tmax = huge(0.0)  ! 设置 tmax 为无穷大

    ! 遍历 x, y, z 三个维度
    do d = 1, 3
        if (isfinite(ray_xyz%dir_inv(d))) then
            ! 计算 t1 和 t2
            t1 = (box_in%min(d) - ray_xyz%origin(d)) * ray_xyz%dir_inv(d)
            t2 = (box_in%max(d) - ray_xyz%origin(d)) * ray_xyz%dir_inv(d)

            ! 确保 t1 是较小值，t2 是较大值
            if (t1 > t2) then
                temp = t1
                t1 = t2
                t2 = temp
            end if

            ! 更新 tmin 和 tmax
            tmin = max(tmin, t1)
            tmax = min(tmax, t2)
        else if (ray_xyz%origin(d) < box_in%min(d) .or. ray_xyz%origin(d) > box_in%max(d)) then
            ! 射线与某个维度的边界盒平行且射线起点不在盒子内
            tmin = huge(0.0)
            exit
        end if
    end do

    ! 判断射线是否与盒子相交
    if (tmin <= tmax) then
        ! 计算交点坐标
        do d = 1, 3
            Pmin(d) = ray_xyz%origin(d) + tmin / ray_xyz%dir_inv(d)
            Pmax(d) = ray_xyz%origin(d) + tmax / ray_xyz%dir_inv(d)
        end do

        ! 计算欧几里得距离
        distance = 0.0
        do d = 1, 3
            diff = Pmax(d) - Pmin(d)
            distance = distance + diff**2
        end do
        length = sqrt(distance)  ! 射线穿过盒子的欧几里得距离
    else
        length = 0.0  ! 如果不相交，长度为 0
    end if

    end subroutine intersections

#ifdef RAYTHEIA_MO
    recursive subroutine raytheia(ray, parent, level, ipix, epray, projected, plength)
        type(HEALPix_ray), intent(in) :: ray
        type(slab) :: ray_xyz
        type(box), intent(in) :: parent
        type(box) :: children(0:7)
        integer, intent(in) :: level, ipix
        integer :: epray(0:nrays-1)
        integer :: projected(0:nrays-1,0:maxpoints)
        real(RK) :: plength(0:nrays-1,0:maxpoints)
        integer :: i, j, k, d, m, parent_index, ir, node_count, start_index, cI, id
        real(RK) :: mid(3), extent(3), center(3)
        logical :: intersect
        real(RK) :: x, y, z, r, xnode, ynode, znode
!        real(RK) :: Aij,c,nu,pc2cm,f,g_i,g_j,thfpix,phfpix,length
        real(RK) :: thfpix,phfpix,length

        if (level > 0) then
            ray_xyz%origin = ray%origin
            thfpix = ray%angle(1)
            phfpix = ray%angle(2)
            ray_xyz%dir(1) = sin(thfpix)*cos(phfpix)
            ray_xyz%dir(2) = sin(thfpix)*sin(phfpix)
            ray_xyz%dir(3) = cos(thfpix)
            do m = 1, 3
                if (abs(ray_xyz%dir(m)) > 1.0e-6) then
                    ray_xyz%dir_inv(m) = 1.0 / ray_xyz%dir(m)
                else
                    ray_xyz%dir_inv(m) = huge(0.D0)  ! 避免除以零
                end if
            end do
            length = 0.0
        
            call intersections(ray_xyz, parent, length)

            intersect = .false.
            if (length > 0.0) then
                intersect = .true.
            endif

            if(intersect) then
                ! generate eight child boxes
                start_index = 0
                do i = 0, 7
                    do d = 1, 3
                        mid(d) = (parent%min(d) + parent%max(d)) / 2.0
                        if (btest(i, d-1)) then
                            children(start_index)%min(d) = mid(d)
                            children(start_index)%max(d) = parent%max(d)
                        else
                            children(start_index)%min(d) = parent%min(d)
                            children(start_index)%max(d) = mid(d)
                        endif
                    enddo
                    start_index = start_index + 1
                enddo

                do i = 0, 7
                    call raytheia(ray, children(i), level-1, ipix, epray, projected, plength)
                enddo
            endif
        endif

        ! leaf nodes calculations
        if(level == 0) then

            I = nint(parent%max(1)/dx)
            J = nint(parent%max(2)/dy)
            K = nint(parent%max(3)/dz)

#ifdef XYZ
            cI = I + (J-1)*nxc + (K-1)*(nxc*nyc)
#else
            cI = K + (J-1)*nzc + (I-1)*(nzc*nyc)
#endif

            ! penetration length
            ray_xyz%origin = ray%origin
            thfpix = ray%angle(1)
            phfpix = ray%angle(2)
            ray_xyz%dir(1) = sin(thfpix)*cos(phfpix)
            ray_xyz%dir(2) = sin(thfpix)*sin(phfpix)
            ray_xyz%dir(3) = cos(thfpix)
            do m = 1, 3
                if (abs(ray_xyz%dir(m)) > 1.0e-6) then
                    ray_xyz%dir_inv(m) = 1.0 / ray_xyz%dir(m)
                else
                    ray_xyz%dir_inv(m) = huge(0.D0)  ! 避免除以零
                end if
            end do
            length = 0.0
        
            call intersections(ray_xyz, parent, length)
            if(length.ne.0.D0) then
                epray(ipix) = epray(ipix) + 1
                id = epray(ipix)
                projected(ipix,id) = cI
                plength(ipix,id) = length
            endif

        endif

    end subroutine raytheia

#else

    recursive subroutine raytheia_table(ray, parent, level, source, ipix)
        type(HEALPix_ray), intent(in) :: ray
        type(slab) :: ray_xyz
        type(box), intent(in) :: parent
        type(box) :: children(0:7)
        integer, intent(in) :: level, source, ipix
        integer :: i, j, k, d, m, parent_index, ir, node_count, start_index, cI, id
        real(RK) :: mid(3), extent(3), center(3)
        logical :: intersect
        real(RK) :: x, y, z, r, xnode, ynode, znode
!        real(RK) :: Aij,c,nu,pc2cm,f,g_i,g_j,thfpix,phfpix,length
        real(RK) :: thfpix,phfpix,length

        if (level > 0) then
            ray_xyz%origin = ray%origin
            thfpix = ray%angle(1)
            phfpix = ray%angle(2)
            ray_xyz%dir(1) = sin(thfpix)*cos(phfpix)
            ray_xyz%dir(2) = sin(thfpix)*sin(phfpix)
            ray_xyz%dir(3) = cos(thfpix)
            do m = 1, 3
                if (abs(ray_xyz%dir(m)) > 1.0e-6) then
                    ray_xyz%dir_inv(m) = 1.0 / ray_xyz%dir(m)
                else
                    ray_xyz%dir_inv(m) = huge(0.D0)  ! 避免除以零
                end if
            end do
            length = 0.0
        
            call intersections(ray_xyz, parent, length)

            intersect = .false.
            if (length > 0.0) then
                intersect = .true.
            endif

            if(intersect) then
                ! generate eight child boxes
                start_index = 0
                do i = 0, 7
                    do d = 1, 3
                        mid(d) = (parent%min(d) + parent%max(d)) / 2.0
                        if (btest(i, d-1)) then
                            children(start_index)%min(d) = mid(d)
                            children(start_index)%max(d) = parent%max(d)
                        else
                            children(start_index)%min(d) = parent%min(d)
                            children(start_index)%max(d) = mid(d)
                        endif
                    enddo
                    start_index = start_index + 1
                enddo

                do i = 0, 7
                    call raytheia_table(ray, children(i), level-1, source, ipix)
                enddo
            endif
        endif

        ! leaf nodes calculations
        if(level == 0) then

            I = nint(parent%max(1)/dx)
            J = nint(parent%max(2)/dy)
            K = nint(parent%max(3)/dz)
            
#ifdef XYZ
            cI = I + (J-1)*nxc + (K-1)*(nxc*nyc)
#else
            cI = K + (J-1)*nzc + (I-1)*(nzc*nyc)
#endif

            ! penetrate length
            ray_xyz%origin = ray%origin
            thfpix = ray%angle(1)
            phfpix = ray%angle(2)
            ray_xyz%dir(1) = sin(thfpix)*cos(phfpix)
            ray_xyz%dir(2) = sin(thfpix)*sin(phfpix)
            ray_xyz%dir(3) = cos(thfpix)
            do m = 1, 3
                if (abs(ray_xyz%dir(m)) > 1.0e-6) then
                    ray_xyz%dir_inv(m) = 1.0 / ray_xyz%dir(m)
                else
                    ray_xyz%dir_inv(m) = huge(0.D0)  ! 避免除以零
                end if
            end do
            length = 0.0
        
            call intersections(ray_xyz, parent, length)
            if(length.ne.0.D0) then
                pdr(source)%epray(ipix) = pdr(source)%epray(ipix) + 1
                id = pdr(source)%epray(ipix)
                pdr(source)%projected(ipix,id) = cI
                pdr(source)%length(ipix,id) = length
            endif
        endif

    end subroutine raytheia_table
#endif

    recursive subroutine raytheia_maxpoints(ray, parent, level, ipix, epray)
        type(HEALPix_ray), intent(in) :: ray
        type(slab) :: ray_xyz
        type(box), intent(in) :: parent
        type(box) :: children(0:7)
        integer, intent(in) :: level, ipix
        integer :: epray(0:nrays-1)
        integer :: projected(0:nrays-1,0:maxpoints)
        real(RK) :: plength(0:nrays-1,0:maxpoints)
        integer :: i, j, k, d, m, parent_index, ir, node_count, start_index, cI, id
        real(RK) :: mid(3), extent(3), center(3)
        logical :: intersect
        real(RK) :: x, y, z, r, xnode, ynode, znode
!        real(RK) :: Aij,c,nu,pc2cm,f,g_i,g_j,thfpix,phfpix,length
        real(RK) :: thfpix,phfpix,length

        if (level > 0) then
            ray_xyz%origin = ray%origin
            thfpix = ray%angle(1)
            phfpix = ray%angle(2)
            ray_xyz%dir(1) = sin(thfpix)*cos(phfpix)
            ray_xyz%dir(2) = sin(thfpix)*sin(phfpix)
            ray_xyz%dir(3) = cos(thfpix)
            do m = 1, 3
                if (abs(ray_xyz%dir(m)) > 1.0e-6) then
                    ray_xyz%dir_inv(m) = 1.0 / ray_xyz%dir(m)
                else
                    ray_xyz%dir_inv(m) = huge(0.D0)  ! 避免除以零
                end if
            end do
            length = 0.0
        
            call intersections(ray_xyz, parent, length)

            intersect = .false.
            if (length > 0.0) then
                intersect = .true.
            endif

            if(intersect) then
                ! generate eight child boxes
                start_index = 0
                do i = 0, 7
                    do d = 1, 3
                        mid(d) = (parent%min(d) + parent%max(d)) / 2.0
                        if (btest(i, d-1)) then
                            children(start_index)%min(d) = mid(d)
                            children(start_index)%max(d) = parent%max(d)
                        else
                            children(start_index)%min(d) = parent%min(d)
                            children(start_index)%max(d) = mid(d)
                        endif
                    enddo
                    start_index = start_index + 1
                enddo

                do i = 0, 7
                    call raytheia_maxpoints(ray, children(i), level-1, ipix, epray)
                enddo
            endif
        endif

        ! leaf nodes calculations
        if(level == 0) then

            I = nint(parent%max(1)/dx)
            J = nint(parent%max(2)/dy)
            K = nint(parent%max(3)/dz)

#ifdef XYZ
            cI = I + (J-1)*nxc + (K-1)*(nxc*nyc)
#else
            cI = K + (J-1)*nzc + (I-1)*(nzc*nyc)
#endif

            ! penetration length
            ray_xyz%origin = ray%origin
            thfpix = ray%angle(1)
            phfpix = ray%angle(2)
            ray_xyz%dir(1) = sin(thfpix)*cos(phfpix)
            ray_xyz%dir(2) = sin(thfpix)*sin(phfpix)
            ray_xyz%dir(3) = cos(thfpix)
            do m = 1, 3
                if (abs(ray_xyz%dir(m)) > 1.0e-6) then
                    ray_xyz%dir_inv(m) = 1.0 / ray_xyz%dir(m)
                else
                    ray_xyz%dir_inv(m) = huge(0.D0)  ! 避免除以零
                end if
            end do
            length = 0.0
        
            call intersections(ray_xyz, parent, length)
            if(length.ne.0.D0) then
                epray(ipix) = epray(ipix) + 1
            endif

        endif

    end subroutine raytheia_maxpoints

end module m_Ray_box
