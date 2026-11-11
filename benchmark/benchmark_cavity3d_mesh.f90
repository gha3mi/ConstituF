! SPDX-License-Identifier: Apache-2.0

module benchmark_cavity3d_mesh
    !! Conforming cubed-sphere hexahedral shell: a cube containing a spherical cavity.
    !! The six structured face patches share nodes along every patch edge.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf, only: wp
    implicit none (type, external)
    private

    integer, parameter, public :: n_side=10, n_layers=12
    real(wp), parameter, public :: half_width=0.030_wp, cavity_radius=0.012_wp
    integer, parameter :: nodes_per_shell=6*n_side*n_side+2
    integer, parameter :: faces_per_shell=6*n_side*n_side
    integer, parameter :: hex_corners=8

    type, public :: cavity_mesh_t
        integer :: n_nodes=0, n_cells=0, n_faces=0
        real(wp), allocatable :: nodes(:,:), centers(:,:), volumes(:)
        real(wp), allocatable :: face_centers(:,:), face_areas(:)
        integer, allocatable :: cells(:,:), faces(:,:)
    end type cavity_mesh_t

    public :: build_cavity_mesh, radial_coordinate, unit_normal

contains

    subroutine build_cavity_mesh(mesh)
        type(cavity_mesh_t), intent(out) :: mesh
        integer :: lookup(0:n_side,0:n_side,0:n_side)
        integer :: i, j, k, face, p, n, corner, layer
        integer :: ix(4), iy(4), iz(4), v(4)
        real(wp) :: radius, t
        real(wp) :: cube(3), x(3,4), cross1(3), cross2(3)

        lookup=0
        n=0
        do k=0,n_side
            do j=0,n_side
                surface_scan: do i=0,n_side
                    if (min(i,j,k)>0.and.max(i,j,k)<n_side) cycle surface_scan
                    n=n+1
                    lookup(i,j,k)=n
                end do surface_scan
            end do
        end do
        if (n/=nodes_per_shell) error stop "Invalid cubed-sphere surface indexing"

        mesh%n_nodes=(n_layers+1)*nodes_per_shell
        mesh%n_faces=faces_per_shell
        mesh%n_cells=n_layers*faces_per_shell
        allocate(mesh%nodes(3,mesh%n_nodes),mesh%centers(3,mesh%n_cells))
        allocate(mesh%cells(8,mesh%n_cells),mesh%volumes(mesh%n_cells))
        allocate(mesh%faces(4,mesh%n_faces),mesh%face_centers(3,mesh%n_faces))
        allocate(mesh%face_areas(mesh%n_faces))

        do layer=0,n_layers
            t=real(layer,wp)/real(n_layers,wp)
            do k=0,n_side
                do j=0,n_side
                    node_scan: do i=0,n_side
                        n=lookup(i,j,k)
                        if (n==0) cycle node_scan
                        cube=2.0_wp*real([i,j,k],wp)/real(n_side,wp)-1.0_wp
                        radius=sqrt(sum(cube*cube))
                        mesh%nodes(:,layer*nodes_per_shell+n)= &
                            (1.0_wp-t)*cavity_radius*cube/radius+t*half_width*cube
                    end do node_scan
                end do
            end do
        end do

        p=0
        do layer=1,n_layers
            do face=1,6
                do j=0,n_side-1
                    facet_scan: do i=0,n_side-1
                        call surface_corner(face,i,j,ix(1),iy(1),iz(1))
                        call surface_corner(face,i+1,j,ix(2),iy(2),iz(2))
                        call surface_corner(face,i+1,j+1,ix(3),iy(3),iz(3))
                        call surface_corner(face,i,j+1,ix(4),iy(4),iz(4))
                        do corner=1,4
                            v(corner)=lookup(ix(corner),iy(corner),iz(corner))
                        end do
                        p=p+1
                        mesh%cells(1:4,p)=v+(layer-1)*nodes_per_shell
                        mesh%cells(5:8,p)=v+layer*nodes_per_shell
                        mesh%centers(:,p)=sum(mesh%nodes(:,mesh%cells(:,p)),dim=2)/8.0_wp
                        mesh%volumes(p)=hexahedron_volume(mesh%nodes(:,mesh%cells(:,p)))
                        if (layer/=1) cycle facet_scan
                        n=p
                        mesh%faces(:,n)=v
                        x=mesh%nodes(:,v)
                        mesh%face_centers(:,n)=sum(x,dim=2)/4.0_wp
                        cross1=cross_product(x(:,2)-x(:,1),x(:,3)-x(:,1))
                        cross2=cross_product(x(:,3)-x(:,1),x(:,4)-x(:,1))
                        mesh%face_areas(n)=0.5_wp*(sqrt(sum(cross1**2))+sqrt(sum(cross2**2)))
                    end do facet_scan
                end do
            end do
        end do

        if (p/=mesh%n_cells.or.any(mesh%cells<=0).or.any(mesh%cells>mesh%n_nodes)) then
            error stop "Invalid cubed-sphere connectivity"
        end if
        if (.not.all(ieee_is_finite(mesh%nodes)).or..not.all(ieee_is_finite(mesh%volumes))) then
            error stop "Nonfinite cubed-sphere geometry"
        end if
        if (any(mesh%volumes<=0.0_wp).or.any(mesh%face_areas<=0.0_wp)) then
            error stop "Inverted or degenerate cavity element"
        end if
        ! Convergence check for the faceted inner boundary (not an exact sphere).
        if (abs(sum(mesh%volumes)-(2.0_wp*half_width)**3+ &
                4.0_wp*acos(-1.0_wp)*cavity_radius**3/3.0_wp)> &
                0.01_wp*(2.0_wp*half_width)**3) then
            error stop "Cavity volume is inconsistent with the analytical geometry"
        end if
    end subroutine build_cavity_mesh

    pure subroutine surface_corner(face,u,v,i,j,k)
        !! Each face uses a right-handed tangential frame relative to
        !! the cube-outward radial direction; VTK hexahedra have positive volume.
        integer, intent(in) :: face,u,v
        integer, intent(out) :: i,j,k

        select case(face)
        case(1)                    ! +x; e_y cross e_z = e_x
            i=n_side
            j=u
            k=v
        case(2)                    ! -x; e_z cross e_y = -e_x
            i=0
            j=v
            k=u
        case(3)                    ! +y; e_z cross e_x = e_y
            i=v
            j=n_side
            k=u
        case(4)                    ! -y; e_x cross e_z = -e_y
            i=u
            j=0
            k=v
        case(5)                    ! +z; e_x cross e_y = e_z
            i=u
            j=v
            k=n_side
        case(6)                    ! -z; e_y cross e_x = -e_z
            i=v
            j=u
            k=0
        case default
            i=0
            j=0
            k=0
        end select
    end subroutine surface_corner

    pure function cross_product(a,b) result(c)
        real(wp), intent(in) :: a(3),b(3)
        real(wp) :: c(3)

        c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
    end function cross_product

    pure real(wp) function hexahedron_volume(x) result(volume)
        !! \(2\times2\times2\) Gauss integration of the trilinear isoparametric Jacobian.
        real(wp), intent(in) :: x(3,hex_corners)
        integer, parameter :: sx(hex_corners)=[-1,1,1,-1,-1,1,1,-1]
        integer, parameter :: sy(hex_corners)=[-1,-1,1,1,-1,-1,1,1]
        integer, parameter :: sz(hex_corners)=[-1,-1,-1,-1,1,1,1,1]
        real(wp) :: xi, eta, zeta, g
        real(wp) :: jac(3,3)
        integer :: a,b,c,k

        volume=0.0_wp
        g=1.0_wp/sqrt(3.0_wp)
        do a=-1,1,2
            do b=-1,1,2
                do c=-1,1,2
                    xi=real(a,wp)*g
                    eta=real(b,wp)*g
                    zeta=real(c,wp)*g
                    jac=0.0_wp
                    do k=1,8
                        jac(:,1)=jac(:,1)+x(:,k)*real(sx(k),wp)* &
                            (1.0_wp+real(sy(k),wp)*eta)*(1.0_wp+real(sz(k),wp)*zeta)/8.0_wp
                        jac(:,2)=jac(:,2)+x(:,k)*real(sy(k),wp)* &
                            (1.0_wp+real(sx(k),wp)*xi)*(1.0_wp+real(sz(k),wp)*zeta)/8.0_wp
                        jac(:,3)=jac(:,3)+x(:,k)*real(sz(k),wp)* &
                            (1.0_wp+real(sx(k),wp)*xi)*(1.0_wp+real(sy(k),wp)*eta)/8.0_wp
                    end do
                    volume=volume+dot_product(jac(:,1),cross_product(jac(:,2),jac(:,3)))
                end do
            end do
        end do
    end function hexahedron_volume

    pure real(wp) function radial_coordinate(position) result(radius)
        real(wp), intent(in) :: position(3)

        radius=sqrt(sum(position**2))
    end function radial_coordinate

    pure function unit_normal(position) result(normal)
        real(wp), intent(in) :: position(3)
        real(wp) :: normal(3)

        normal=position/radial_coordinate(position)
    end function unit_normal

end module benchmark_cavity3d_mesh
