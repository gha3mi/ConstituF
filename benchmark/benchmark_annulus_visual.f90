! SPDX-License-Identifier: Apache-2.0

module benchmark_annulus_visual
    !! Scientific response contours and analytical parameter heatmaps.
    !! Heatmap CSV data retain normalized values and resolution categories.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use benchmark_annulus_mesh, only: r_inner, r_outer, interface_radius, nradial, nangle, &
                                     node_radius, node_angle, cell_radius, cell_angle, pi
    use benchmark_colormaps, only: Colormap, cm_viridis, colormap_for, colormap_name, colormap_color, colormap_limits, init_colormap
    use benchmark_metrics, only: sensitivity_metric_t
    use benchmark_models, only: cases, names, interface_case, max_parameters
    use benchmark_paths, only: annulus_directory
    use constituf, only: wp
    implicit none (type, external)
    private

    real(wp), parameter :: center_x=315.0_wp, center_y=323.0_wp, plot_radius=248.0_wp
    integer, parameter :: plot_x=348, plot_width=800, plot_y=46, row_height=48
    integer, parameter :: canvas_width=1430, legend_x=1200
    integer, parameter :: undefined_cell=0, zero_cell=1, unresolved_cell=2, valid_cell=3
    public :: write_annulus_preview, write_annulus_heatmap

contains

    subroutine write_annulus_preview(path,id,values,label,units)
        character(*), intent(in) :: path,label,units
        integer, intent(in) :: id
        real(wp), intent(in) :: values(:)
        integer :: u, i, j, p, m, scheme
        real(wp) :: lo, hi, r0, r1, a0, a1, x0, x1, x2, x3, y0, y1, y2, y3
        real(wp) :: t, y, tick_value
        type(Colormap) :: cmap
        character(7) :: color
        character(24) :: display_label

        display_label=label
        if (index(label,"d(")/=1) then
        select case(trim(cases(id)%response_label))
        case("heat_flux_radial")
            display_label="q_r"
        case("molar_flux_radial")
            display_label="j_r"
        case("current_density_radial")
            display_label="J_r"
        case("electric_displacement_radial")
            display_label="D_r"
        case("darcy_velocity_radial")
            display_label="v_r"
        case default
        end select
        end if

        scheme=colormap_for(id)
        call init_colormap(scheme,cmap)
        call colormap_limits(scheme,values,lo,hi)

        open(newunit=u,file=path,status="replace",action="write")
        write(u,"(a)") '<?xml version="1.0" encoding="UTF-8"?>'
        write(u,"(a)") '<svg xmlns="http://www.w3.org/2000/svg" width="1080" height="810" viewBox="0 40 840 630">'
        write(u,"(a)") '<rect x="0" y="40" width="840" height="630" fill="#ffffff"/>'
        write(u,"(a)") "<desc>Constitutive response benchmark; annular radial response; "//trim(display_label)//" ["//trim(units)// &
            "]. Linear "//trim(colormap_name(scheme))//" color scale. See colorbar for physical values.</desc>"
        if (interface_case(id)) then
            write(u,"(a)") '<circle cx="315" cy="323" r="248" fill="#e7ecef"/>'
            write(u,"(a)") '<circle cx="315" cy="323" r="82.67" fill="white"/>'
            do j=1,nangle
                a0=node_angle(j)
                a1=node_angle(j+1)
                x0=svg_x(interface_radius,a0)
                y0=svg_y(interface_radius,a0)
                x1=svg_x(interface_radius,a1)
                y1=svg_y(interface_radius,a1)
                color=colormap_color(cmap,(values(j)-lo)/(hi-lo))
                write(u,"(a)") '<line x1="'//trim(num(x0))//'" y1="'//trim(num(y0))// &
                    '" x2="'//trim(num(x1))//'" y2="'//trim(num(y1))// &
                    '" stroke="'//color//'" stroke-width="13" stroke-linecap="round"/>'
            end do
        else
            do i=1,nradial
                r0=node_radius(i)
                r1=node_radius(i+1)
                do j=1,nangle
                    p=(i-1)*nangle+j
                    a0=node_angle(j)
                    a1=node_angle(j+1)
                    x0=svg_x(r0,a0)
                    y0=svg_y(r0,a0)
                    x1=svg_x(r1,a0)
                    y1=svg_y(r1,a0)
                    x2=svg_x(r1,a1)
                    y2=svg_y(r1,a1)
                    x3=svg_x(r0,a1)
                    y3=svg_y(r0,a1)
                    color=colormap_color(cmap,(values(p)-lo)/(hi-lo))
                    write(u,"(a)") '<polygon points="'//trim(num(x0))//","//trim(num(y0))// &
                        " "//trim(num(x1))//","//trim(num(y1))// &
                        " "//trim(num(x2))//","//trim(num(y2))// &
                        " "//trim(num(x3))//","//trim(num(y3))// &
                        '" fill="'//color//'" stroke="'//color//'" stroke-width="0.5"/>'
                end do
            end do
            write(u,"(a)") '<circle cx="315" cy="323" r="82.67" fill="white"/>'
        end if

        write(u,"(a)") '<circle cx="315" cy="323" r="248" fill="none" stroke="#344251" stroke-width="1.3"/>'
        write(u,"(a)") '<circle cx="315" cy="323" r="82.67" fill="none" stroke="#344251" stroke-width="1.3"/>'
        write(u,"(a)") '<line x1="67" y1="600" x2="563" y2="600" stroke="#344251"/>'
        do i=-2,2
            x0=center_x+plot_radius*real(i,wp)/2.0_wp
            write(u,"(a)") '<line x1="'//trim(num(x0))//'" y1="600" x2="'//trim(num(x0))// &
                '" y2="607" stroke="#344251"/>'
            write(u,"(a)") '<text x="'//trim(num(x0))// &
                '" y="629" text-anchor="middle" font-size="15" fill="#263645">'// &
                trim(integer_str(i*15))//"</text>"
        end do
        write(u,"(a)") '<text x="315" y="655" text-anchor="middle" font-size="15" fill="#263645">x [mm]</text>'
        write(u,"(a)") '<line x1="67" y1="75" x2="67" y2="571" stroke="#344251"/>'
        do i=-2,2
            y=center_y-plot_radius*real(i,wp)/2.0_wp
            write(u,"(a)") '<line x1="61" y1="'//trim(num(y))//'" x2="67" y2="'//trim(num(y))// &
                '" stroke="#344251"/>'
            write(u,"(a)") '<text x="54" y="'//trim(num(y+5.0_wp))// &
                '" text-anchor="end" font-size="15" fill="#263645">'//trim(integer_str(i*15))//"</text>"
        end do
        write(u,"(a)") '<text x="13" y="323" transform="rotate(-90 13 323)" text-anchor="middle" font-size="15" fill="#263645">y [mm]</text>'

        ! Match the contour normalization exactly; 256 short swatches have
        ! overlapping edges so there are no white seams in SVG viewers.
        do m=0,255
            t=real(m,wp)/255.0_wp
            color=colormap_color(cmap,1.0_wp-t)
            y=153.0_wp+1.5_wp*real(m,wp)
            write(u,"(a)") '<rect x="649" y="'//trim(num(y))// &
                '" width="27" height="1.55" fill="'//color//'"/>'
        end do
        write(u,"(a)") '<rect x="649" y="153" width="27" height="384" fill="none" stroke="#344251" stroke-width="0.7"/>'
        write(u,"(a)") '<text x="648" y="110" font-size="15" fill="#263645">'//trim(display_label)//"</text>"
        write(u,"(a)") '<text x="648" y="133" font-size="14" fill="#506174">['//trim(units)//"]</text>"
        do m=0,4
            t=real(m,wp)/4.0_wp
            y=153.0_wp+384.0_wp*t
            tick_value=hi-t*(hi-lo)
            write(u,"(a)") '<line x1="678" y1="'//trim(num(y))// &
                '" x2="686" y2="'//trim(num(y))//'" stroke="#344251" stroke-width="1.2"/>'
            write(u,"(a)") '<text x="692" y="'//trim(num(y+5.0_wp))// &
                '" font-size="14" font-variant-numeric="tabular-nums" fill="#263645">'// &
                trim(scientific_label(tick_value))//"</text>"
        end do
        write(u,"(a)") "</svg>"
        close(u)
    end subroutine write_annulus_preview

    pure real(wp) function svg_x(r,angle) result(value)
        real(wp), intent(in) :: r,angle
        value=center_x+plot_radius*r/r_outer*cos(angle)
    end function svg_x

    pure real(wp) function svg_y(r,angle) result(value)
        real(wp), intent(in) :: r,angle
        value=center_y-plot_radius*r/r_outer*sin(angle)
    end function svg_y

    function num(value) result(s)
        real(wp), intent(in) :: value
        character(32) :: s
        write(s,"(f12.3)") value
        s=adjustl(s)
    end function num

    function scientific_label(value) result(s)
        real(wp), intent(in) :: value
        character(32) :: s
        write(s,"(es12.3e2)") value
        s=adjustl(s)
    end function scientific_label

    function integer_str(value) result(s)
        integer, intent(in) :: value
        character(16) :: s
        write(s,"(i0)") value
    end function integer_str

    subroutine write_annulus_heatmap(id,parameters,units,metrics,fields)
        !! Response \(y\) shares units with its scaled derivative \(s_p\,\mathrm{d}y/\mathrm{d}p\),
        !! where \(p\) is a scalar parameter and \(s_p\) its reference scale.
        !! At each radial station, the normalized magnitude is
        !! \(\operatorname{RMS}_{\theta}(s_p\,\mathrm{d}y/\mathrm{d}p)/\operatorname{RMS}_{\Omega}(y)\).
        !! Cohesive interfaces use \(|s_p\,\mathrm{d}y/\mathrm{d}p|/\operatorname{RMS}_{\mathrm{interface}}(y)\).
        !! Heatmap colors show the base-10 logarithm of this magnitude.
        integer, intent(in) :: id
        character(*), intent(in) :: parameters(:),units(:)
        type(sensitivity_metric_t), intent(in) :: metrics(:)
        real(wp), intent(in) :: fields(:,:)
        real(wp) :: local_ratio(nangle,max_parameters)
        real(wp) :: signed_ratio(nangle,max_parameters)
        integer :: categories(nangle,max_parameters)
        integer :: n,nsample,j,i,first,last
        real(wp) :: reference,magnitude,lo,hi,maximum,minimum
        character(:), allocatable :: prefix

        n=size(parameters)
        nsample=nradial
        if (interface_case(id)) nsample=nangle
        if (n<=0.or.n>max_parameters) error stop "Invalid heatmap parameter count"
        if (size(units)/=n.or.size(metrics)/=n.or.size(fields,2)/=n) then
            error stop "Invalid per-material heatmap metadata"
        end if
        if (size(fields,1)/=nangle*nsample.and..not.interface_case(id)) then
            error stop "Incorrect bulk annulus sample shape"
        end if
        if (interface_case(id).and.size(fields,1)/=nangle) error stop "Incorrect interface sample shape"
        if (.not.all(ieee_is_finite(fields))) error stop "Nonfinite analytical heatmap input"

        local_ratio=0.0_wp
        signed_ratio=0.0_wp
        categories=undefined_cell
        reference=metrics(1)%response_rms
        if (reference>0.0_wp) then
            do j=1,n
                do i=1,nsample
                    if (interface_case(id)) then
                        local_ratio(i,j)=abs(fields(i,j))/reference
                        signed_ratio(i,j)=fields(i,j)/reference
                    else
                        first=(i-1)*nangle+1
                        last=i*nangle
                        magnitude=maxval(abs(fields(first:last,j)))
                        if (magnitude>0.0_wp) then
                            local_ratio(i,j)=(magnitude/reference)* &
                                sqrt(sum((fields(first:last,j)/magnitude)**2)/real(nangle,wp))
                            signed_ratio(i,j)=(magnitude/reference)* &
                                sum(fields(first:last,j)/magnitude)/real(nangle,wp)
                        end if
                    end if
                    if (.not.ieee_is_finite(local_ratio(i,j)).or. &
                        .not.ieee_is_finite(signed_ratio(i,j))) error stop "Invalid normalized heatmap cell"
                    categories(i,j)=valid_cell
                    if (local_ratio(i,j)<=0.0_wp) categories(i,j)=zero_cell
                    if (local_ratio(i,j)>0.0_wp.and. &
                        local_ratio(i,j)<=256.0_wp*epsilon(1.0_wp)) categories(i,j)=unresolved_cell
                end do
            end do
        end if

        minimum=huge(1.0_wp)
        maximum=0.0_wp
        do j=1,n
            resolved_bin_loop: do i=1,nsample
                if (categories(i,j)/=valid_cell) cycle resolved_bin_loop
                minimum=min(minimum,local_ratio(i,j))
                maximum=max(maximum,local_ratio(i,j))
            end do resolved_bin_loop
        end do
        lo=-1.0_wp
        hi=1.0_wp
        if (maximum>0.0_wp) then
            lo=real(floor(log10(minimum)),wp)
            hi=real(ceiling(log10(maximum)),wp)
            if (hi<=lo) hi=lo+1.0_wp
        end if

        prefix=annulus_directory//trim(names(id))//"_sensitivity"
        call export_heatmap_csv(prefix//".csv",id,parameters,units,metrics, &
                                local_ratio,signed_ratio,categories,nsample)
        call export_heatmap_svg(prefix//".svg",id,parameters,units,metrics, &
                                local_ratio,signed_ratio,categories,nsample,lo,hi)
    end subroutine write_annulus_heatmap

    subroutine export_heatmap_csv(filename,id,parameters,units,metrics,ratios,signed,categories,nsample)
        character(*), intent(in) :: filename
        character(*), intent(in) :: parameters(:),units(:)
        integer, intent(in) :: id,nsample
        type(sensitivity_metric_t), intent(in) :: metrics(:)
        real(wp), intent(in) :: ratios(:,:),signed(:,:)
        integer, intent(in) :: categories(:,:)
        real(wp) :: coordinate
        integer :: unit,j,i
        character(24) :: log_label
        character(12) :: coordinate_unit

        coordinate_unit="m"
        if (interface_case(id)) coordinate_unit="rad"
        open(newunit=unit,file=filename,status="replace",action="write")
        write(unit,"(a)") "parameter_index,parameter,parameter_units,parameter_value,parameter_scale,coordinate,coordinate_units,relative_rms,log10_relative_rms,signed_mean_relative,status"
        do j=1,size(parameters)
            do i=1,nsample
                coordinate=cell_radius(i)
                if (interface_case(id)) coordinate=cell_angle(i)
                log_label="NA"
                if (categories(i,j)==valid_cell.or.categories(i,j)==unresolved_cell) then
                    write(log_label,"(es21.12e3)") log10(ratios(i,j))
                end if
                write(unit,'(i0,",",a,",",a,",",es22.14e3,",",es22.14e3,",",es22.14e3,",",a,",",es22.14e3,",",a,",",es22.14e3,",",a)') &
                    j,trim(parameters(j)),trim(units(j)),metrics(j)%reference_value, &
                    metrics(j)%parameter_scale,coordinate,trim(coordinate_unit),ratios(i,j), &
                    trim(adjustl(log_label)),signed(i,j),trim(status_name(categories(i,j)))
            end do
        end do
        close(unit)
    end subroutine export_heatmap_csv

    subroutine export_heatmap_svg(filename,id,parameters,units,metrics,ratios,signed,categories,nsample,lo,hi)
        character(*), intent(in) :: filename
        character(*), intent(in) :: parameters(:),units(:)
        integer, intent(in) :: id,nsample
        type(sensitivity_metric_t), intent(in) :: metrics(:)
        real(wp), intent(in) :: ratios(:,:),signed(:,:)
        real(wp), intent(in) :: lo,hi
        integer, intent(in) :: categories(:,:)
        real(wp) :: x,y,dx,fraction,coordinate,tick_value
        integer :: unit,height,i,j,tick,m,legend_height,legend_bottom
        type(Colormap) :: cmap
        character(7) :: color
        character(32) :: exponent
        character(18) :: coord_label
        character(6) :: coordinate_name
        character(3) :: coordinate_unit

        coordinate_name="radius"
        coordinate_unit="m"
        if (interface_case(id)) then
            coordinate_name="angle"
            coordinate_unit="rad"
        end if
        legend_height=max(120,size(parameters)*row_height-6)
        legend_bottom=plot_y+legend_height
        height=max(plot_y+size(parameters)*row_height+90,legend_bottom+98)
        dx=real(plot_width,wp)/real(nsample,wp)
        call init_colormap(cm_viridis,cmap)
        open(newunit=unit,file=filename,status="replace",action="write")
        write(unit,"(a)") '<?xml version="1.0" encoding="UTF-8"?>'
        write(unit,"(a)") '<svg xmlns="http://www.w3.org/2000/svg" width="'//trim(integer_string(canvas_width))// &
            '" height="'//trim(integer_string(height))//'" viewBox="0 0 '//trim(integer_string(canvas_width))// &
            " "//trim(integer_string(height))//'">'
        write(unit,"(a)") '<rect width="100%" height="100%" fill="#ffffff"/>'
        write(unit,"(a)") "<desc>Analytical material-parameter sensitivity heatmap. Rows identify scalar parameters;"// &
            " columns identify radial positions or interface angles. Colors encode log10 normalized sensitivity RMS."// &
            " Cell tooltips provide numerical values, coordinate units and differentiation status.</desc>"

        do j=1,size(parameters)
            y=real(plot_y+(j-1)*row_height,wp)
            if (mod(j,2)==0) then
                write(unit,"(a)") '<rect x="15" y="'//trim(decimal(y))// &
                    '" width="1140" height="'//trim(integer_string(row_height))//'" fill="#F3F6F9"/>'
            end if
            write(unit,"(a)") '<text x="24" y="'//trim(decimal(y+29.0_wp))// &
                '" font-family="Georgia,serif" font-size="13" fill="#23374a">P'// &
                trim(integer_string(j))//"  "//trim(parameters(j))//" ["//trim(units(j))//"]</text>"
            do i=1,nsample
                x=real(plot_x,wp)+real(i-1,wp)*dx
                color="#E5B57B"
                select case(categories(i,j))
                case(zero_cell)
                    color="#D7DBE2"
                case(unresolved_cell)
                    color="#BAC5D1"
                case(valid_cell)
                    fraction=(log10(ratios(i,j))-lo)/(hi-lo)
                    color=colormap_color(cmap,fraction)
                case default
                end select
                coordinate=cell_radius(i)
                if (interface_case(id)) coordinate=cell_angle(i)
                write(unit,"(a)") '<g class="sensitivity-cell" data-parameter="'// &
                    trim(integer_string(j))//'" data-bin="'//trim(integer_string(i))//'">'
                write(unit,"(a)") "<title>"//trim(names(id))//" / "//trim(parameters(j))// &
                    " ["//trim(units(j))//"]; "//trim(coordinate_name)//"="//trim(scientific(coordinate))// &
                    " "//trim(coordinate_unit)//"; normalized RMS="//trim(scientific(ratios(i,j)))// &
                    "; normalized signed mean="//trim(scientific(signed(i,j)))// &
                    "; status="//trim(status_name(categories(i,j)))//"</title>"
                write(unit,"(a)") '<rect x="'//trim(decimal(x))//'" y="'//trim(decimal(y+3.0_wp))// &
                    '" width="'//trim(decimal(dx+0.02_wp))//'" height="'// &
                    trim(integer_string(row_height-6))//'" fill="'//color//'"/>'
                write(unit,"(a)") "</g>"
            end do
        end do

        do tick=0,4
            x=real(plot_x,wp)+real(plot_width,wp)*real(tick,wp)/4.0_wp
            coordinate=1000.0_wp*r_inner+(1000.0_wp*(r_outer-r_inner))*real(tick,wp)/4.0_wp
            if (interface_case(id)) coordinate=360.0_wp*real(tick,wp)/4.0_wp
            write(coord_label,"(f7.1)") coordinate
            write(unit,"(a)") '<line x1="'//trim(decimal(x))//'" y1="'// &
                trim(integer_string(plot_y+size(parameters)*row_height))// &
                '" x2="'//trim(decimal(x))//'" y2="'// &
                trim(integer_string(plot_y+size(parameters)*row_height+7))//'" stroke="#42566a"/>'
            write(unit,"(a)") '<text x="'//trim(decimal(x))//'" y="'// &
                trim(integer_string(plot_y+size(parameters)*row_height+28))// &
                '" text-anchor="middle" font-family="Georgia,serif" font-size="14" fill="#30475e">'// &
                trim(adjustl(coord_label))//"</text>"
        end do
        coord_label="radius [mm]"
        if (interface_case(id)) coord_label="angle [deg]"
        write(unit,"(a)") '<text x="'//trim(integer_string(plot_x+plot_width/2))//'" y="'// &
            trim(integer_string(plot_y+size(parameters)*row_height+65))// &
            '" text-anchor="middle" font-family="Georgia,serif" font-size="15" fill="#30475e">'// &
            trim(coord_label)//"</text>"

        write(unit,"(a)") '<text x="'//trim(integer_string(legend_x))// &
            '" y="'//trim(integer_string(plot_y-16))// &
            '" font-family="Georgia,serif" font-size="14" fill="#293b4c">log10 relative RMS</text>'
        do m=0,255
            fraction=1.0_wp-real(m,wp)/255.0_wp
            color=colormap_color(cmap,fraction)
            write(unit,"(a)") '<rect x="'//trim(integer_string(legend_x))//'" y="'// &
                trim(decimal(real(plot_y,wp)+real(legend_height,wp)*real(m,wp)/255.0_wp))// &
                '" width="27" height="'//trim(decimal(real(legend_height,wp)/255.0_wp+0.08_wp))// &
                '" fill="'//color//'"/>'
        end do
        write(unit,"(a)") '<rect x="'//trim(integer_string(legend_x))// &
            '" y="'//trim(integer_string(plot_y))//'" width="27" height="'// &
            trim(integer_string(legend_height))//'" fill="none" stroke="#344251" stroke-width="0.7"/>'
        do tick=0,4
            tick_value=hi-(hi-lo)*real(tick,wp)/4.0_wp
            write(exponent,"(f9.2)") tick_value
            write(unit,"(a)") '<text x="'//trim(integer_string(legend_x+36))//'" y="'// &
                trim(decimal(real(plot_y+5,wp)+real(legend_height,wp)*real(tick,wp)/4.0_wp))// &
                '" font-family="Georgia,serif" font-size="13" fill="#293b4c">'// &
                trim(adjustl(exponent))//"</text>"
        end do
        call draw_status_legend(unit,legend_x,legend_bottom+18,"#D7DBE2","exact zero")
        call draw_status_legend(unit,legend_x,legend_bottom+43,"#BAC5D1","below resolution")
        call draw_status_legend(unit,legend_x,legend_bottom+68,"#E5B57B","undefined")
        write(unit,"(a)") "</svg>"
        close(unit)
    end subroutine export_heatmap_svg

    subroutine draw_status_legend(unit,x,y,color,label)
        integer, intent(in) :: unit,x,y
        character(*), intent(in) :: color,label

        write(unit,"(a)") '<rect x="'//trim(integer_string(x))//'" y="'// &
            trim(integer_string(y))//'" width="19" height="16" fill="'//color//'"/>'
        write(unit,"(a)") '<text x="'//trim(integer_string(x+27))//'" y="'// &
            trim(integer_string(y+13))//'" font-family="Georgia,serif" font-size="12"'// &
            ' fill="#293b4c">'//trim(label)//"</text>"
    end subroutine draw_status_legend

    pure function status_name(classification) result(name)
        integer, intent(in) :: classification
        character(18) :: name

        name="undefined"
        select case(classification)
        case(zero_cell)
            name="exact_zero"
        case(unresolved_cell)
            name="below_resolution"
        case(valid_cell)
            name="ok"
        case default
            name="undefined"
        end select
    end function status_name

    function scientific(value) result(s)
        real(wp), intent(in) :: value
        character(24) :: s

        write(s,"(es13.5e2)") value
        s=adjustl(s)
    end function scientific

    function decimal(value) result(s)
        real(wp), intent(in) :: value
        character(24) :: s

        write(s,"(f12.3)") value
        s=adjustl(s)
    end function decimal

    function integer_string(value) result(s)
        integer, intent(in) :: value
        character(18) :: s

        write(s,"(i0)") value
    end function integer_string

end module benchmark_annulus_visual
