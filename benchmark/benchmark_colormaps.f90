! SPDX-License-Identifier: Apache-2.0

module benchmark_colormaps
    !! Benchmark palette selection and physical plotting limits.
    !! ForColormap provides the RGB values for all scientific palettes.
    use benchmark_models, only: cases, ncases
    use constituf, only: wp
    use forcolormap, only: Colormap, color_wp => wp
    implicit none (type, external)
    private

    integer, parameter, public :: cm_vik=1
    integer, parameter, public :: cm_magma=2
    integer, parameter, public :: cm_viridis=3
    integer, parameter, public :: cm_batlow=4

    public :: Colormap, colormap_for, colormap_name, colormap_limits
    public :: init_colormap, colormap_color

contains

    pure integer function colormap_for(id) result(kind)
        integer, intent(in) :: id

        if (id<1.or.id>ncases) error stop "Unrecognized benchmark model in colormap selection"
        kind=cases(id)%colormap_kind
    end function colormap_for

    pure function colormap_name(kind) result(name)
        integer, intent(in) :: kind
        character(12) :: name

        select case(kind)
        case(cm_vik)
            name="vik"
        case(cm_magma)
            name="magma"
        case(cm_viridis)
            name="viridis"
        case(cm_batlow)
            name="batlow"
        case default
            error stop "Unrecognized scientific colormap"
        end select
    end function colormap_name

    pure logical function colormap_is_diverging(kind) result(diverging)
        integer, intent(in) :: kind

        diverging=kind==cm_vik
    end function colormap_is_diverging

    pure subroutine colormap_limits(kind,values,lo,hi)
        !! Diverging: zero maps to the neutral midpoint, with symmetric
        !! limits \((-M,M)\) for a positive scale \(M\). Sequential: use the actual linear range.
        !! Near-constant fields are broadened to avoid emphasizing roundoff.
        integer, intent(in) :: kind
        real(wp), intent(in) :: values(:)
        real(wp), intent(out) :: lo,hi
        real(wp) :: magnitude, midpoint, halfwidth, span

        lo=minval(values)
        hi=maxval(values)
        if (colormap_is_diverging(kind)) then
            magnitude=max(abs(lo),abs(hi))
            if (magnitude<=tiny(1.0_wp)) magnitude=1.0_wp
            lo=-magnitude
            hi=magnitude
        else
            magnitude=max(abs(lo),abs(hi))
            span=hi-lo
            if (span<=128.0_wp*epsilon(1.0_wp)*magnitude) then
                midpoint=lo+0.5_wp*span
                halfwidth=0.05_wp*magnitude
                if (halfwidth<=tiny(1.0_wp)) halfwidth=1.0_wp
                lo=midpoint-halfwidth
                hi=midpoint+halfwidth
            end if
        end if
    end subroutine colormap_limits

    pure subroutine init_colormap(kind,cmap)
        integer, intent(in) :: kind
        type(Colormap), intent(inout) :: cmap

        call cmap%set(trim(colormap_name(kind)),0.0_color_wp,1.0_color_wp)
    end subroutine init_colormap

    pure function colormap_color(cmap,fraction) result(hex)
        type(Colormap), intent(in) :: cmap
        real(wp), intent(in) :: fraction
        character(7) :: hex
        integer :: red
        integer :: green
        integer :: blue

        call cmap%compute_RGB(real(max(0.0_wp,min(1.0_wp,fraction)),color_wp),red,green,blue)
        write(hex,'("#",3(z2.2))') red,green,blue
    end function colormap_color

end module benchmark_colormaps
