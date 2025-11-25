module m_Ray_box
    use m_Mesh
    implicit none
    private
    real(RK),parameter,public::PC = 3.08568025D+18

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

    public:: octree

contains
    logical function isfinite(x)
            real(RK), intent(in) :: x
            isfinite = (abs(x) < huge(x) .and. x == x) 
    end function isfinite

    subroutine intersections(ray, box_in, length)
    implicit none
    type(slab), intent(in) :: ray
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
        if (isfinite(ray%dir_inv(d))) then
            ! 计算 t1 和 t2
            t1 = (box_in%min(d) - ray%origin(d)) * ray%dir_inv(d)
            t2 = (box_in%max(d) - ray%origin(d)) * ray%dir_inv(d)

            ! 确保 t1 是较小值，t2 是较大值
            if (t1 > t2) then
                temp = t1
                t1 = t2
                t2 = temp
            end if

            ! 更新 tmin 和 tmax
            tmin = max(tmin, t1)
            tmax = min(tmax, t2)
        else if (ray%origin(d) < box_in%min(d) .or. ray%origin(d) > box_in%max(d)) then
            ! 射线与某个维度的边界盒平行且射线起点不在盒子内
            tmin = huge(0.0)
            exit
        end if
    end do

    ! 判断射线是否与盒子相交
    if (tmin <= tmax) then
        ! 计算交点坐标
        do d = 1, 3
            Pmin(d) = ray%origin(d) + tmin / ray%dir_inv(d)
            Pmax(d) = ray%origin(d) + tmax / ray%dir_inv(d)
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

    recursive subroutine octree(ray, parent, level, contribution,contribution_1D,density_3D,density_1D)
        type(HEALPix_ray), intent(in) :: ray
        type(slab) :: ray_xyz
        type(box), intent(in) :: parent
        type(box) :: children(0:7)
        integer, intent(in) :: level
        real(RK), intent(inout) :: contribution,contribution_1D
        integer :: i, j, k, d, m, parent_index, ir, node_count, start_index, p
        real(RK) :: mid(3), extent(3), center(3)
        logical :: intersect
        real(RK) :: x, y, z, r, xnode, ynode, znode
        real(RK) :: Aij,c,nu,pc2cm,f,g_i,g_j,thfpix,phfpix,length
        real(RK) :: density_3D(:,:,:)
        real(RK) :: density_1D(:)

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
                    call octree(ray, children(i), level-1, contribution,contribution_1D,density_3D,density_1D)
                enddo
            endif
        endif

        ! leaf nodes calculations
        if(level == 0) then

            I = nint(parent%max(1)/dx)
            J = nint(parent%max(2)/dy)
            K = nint(parent%max(3)/dz)
            p = I + (J-1)*nxc + (K-1)*(nxc*nyc)

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
            if(I.eq.33.and.J.eq.1.and.K.eq.32) print*,ray_xyz%origin, ray_xyz%dir,parent%min,parent%max, length
            if(length.ne.0) print*, I,J,K,length
            contribution = contribution + density_3D(I,J,K)*length*pc
            contribution_1D = contribution_1D + density_1D(p)*length*pc
        endif

    end subroutine octree

end module m_Ray_box
