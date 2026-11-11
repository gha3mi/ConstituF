! SPDX-License-Identifier: Apache-2.0

module benchmark_annulus_mesh
    !! Consistent SI geometry shared by the transport and mechanical probes.
    use constituf, only: wp
    implicit none (type, external)
    private

    integer, parameter, public :: nradial = 24, nangle = 64
    real(wp), parameter, public :: pi = acos(-1.0_wp)
    real(wp), parameter, public :: r_inner = 0.010_wp, r_outer = 0.030_wp
    real(wp), parameter, public :: interface_radius = 0.020_wp
    real(wp), parameter, public :: dr = (r_outer-r_inner)/real(nradial,wp)

    public :: cell_radius, cell_angle, node_radius, node_angle

contains

    pure real(wp) function node_radius(i) result(value)
        integer, intent(in) :: i
        value = r_inner + real(i-1,wp)*dr
    end function node_radius

    pure real(wp) function cell_radius(i) result(value)
        integer, intent(in) :: i
        value = r_inner + (real(i,wp)-0.5_wp)*dr
    end function cell_radius

    pure real(wp) function node_angle(j) result(value)
        integer, intent(in) :: j
        value = 2.0_wp*pi*real(j-1,wp)/real(nangle,wp)
    end function node_angle

    pure real(wp) function cell_angle(j) result(value)
        integer, intent(in) :: j
        value = 2.0_wp*pi*(real(j,wp)-0.5_wp)/real(nangle,wp)
    end function cell_angle
end module benchmark_annulus_mesh
