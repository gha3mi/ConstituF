! SPDX-License-Identifier: Apache-2.0

module benchmark_cavity3d_visual
    !! Dependency-free 3D projected response previews and parameter heatmaps.
    !! The VTK export contains the full uncut three-dimensional geometry.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use benchmark_cavity3d_mesh, only: cavity_mesh_t, n_layers, unit_normal, half_width, cavity_radius
    use benchmark_colormaps, only: Colormap, cm_viridis, colormap_for, colormap_color, colormap_limits, init_colormap
    use benchmark_metrics, only: sensitivity_metric_t
    use benchmark_models, only: names
    use benchmark_paths, only: cavity3d_directory
    use constituf, only: wp, lk
    implicit none (type, external)
    private

    integer, parameter :: n_hexa_faces=6
    integer, parameter :: local_face(4,n_hexa_faces)=reshape([ &
        1,4,3,2, 5,6,7,8, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8],[4,n_hexa_faces])
    integer, parameter :: exact_zero=1, below_resolution=2, undefined_ratio=3, resolved=4
    integer, parameter :: depth_buckets=120, hash_size=131071

    public :: write_cavity_preview, write_cavity_heatmap

contains

    subroutine write_cavity_preview(mesh,id,is_interface,values,label,units)
        !! SVG painters' projection of the outer surface, three cut planes
        !! and the spherical wall exposed by removing the positive octant.
        !! The removal is a view operation only; VTK retains all cells.
        type(cavity_mesh_t), intent(in) :: mesh
        integer, intent(in) :: id
        logical, intent(in) :: is_interface
        real(wp), intent(in) :: values(:)
        character(*), intent(in) :: label,units
        integer, allocatable :: vertices(:,:),owner(:),hits(:),hash_table(:),keys(:,:)
        integer :: max_faces,n_faces,cell,f,slot,existing,side,bucket,u,k,scheme,p,m
        integer :: nodes(4),sorted(4)
        real(wp) :: depth, min_depth, max_depth, lo, hi, fraction
        real(wp) :: center(3)
        real(wp) :: tick_value
        real(wp) :: screen(2)
        type(Colormap) :: cmap
        character(7) :: color
        character(:), allocatable :: path

        if (.not.all(ieee_is_finite(values))) error stop "Nonfinite cavity preview field"
        scheme=colormap_for(id)
        call init_colormap(scheme,cmap)
        call colormap_limits(scheme,values,lo,hi)
        max_faces=mesh%n_faces
        if (.not.is_interface) max_faces=n_hexa_faces*mesh%n_cells
        allocate(vertices(4,max_faces),owner(max_faces),hits(max_faces))
        allocate(keys(4,max_faces),hash_table(hash_size))
        hits=0
        hash_table=0
        n_faces=0

        if (is_interface) then
            do f=1,mesh%n_faces
                n_faces=n_faces+1
                vertices(:,n_faces)=mesh%faces(:,f)
                owner(n_faces)=f
                hits(n_faces)=1
            end do
        else
            visible_cells: do cell=1,mesh%n_cells
                ! Remove one octant to reveal the actual internal cavity.
                if (all(mesh%centers(:,cell)>0.0_wp)) cycle visible_cells
                do side=1,n_hexa_faces
                    nodes=mesh%cells(local_face(:,side),cell)
                    sorted=sorted_four(nodes)
                    slot=hash_slot(sorted)
                    face_lookup: do
                        existing=hash_table(slot)
                        if (existing==0) then
                            n_faces=n_faces+1
                            vertices(:,n_faces)=nodes
                            keys(:,n_faces)=sorted
                            owner(n_faces)=cell
                            hits(n_faces)=1
                            hash_table(slot)=n_faces
                            exit face_lookup
                        end if
                        if (all(keys(:,existing)==sorted)) then
                            hits(existing)=hits(existing)+1
                            exit face_lookup
                        end if
                        slot=mod(slot,hash_size)+1
                    end do face_lookup
                end do
            end do visible_cells
        end if

        ! Projection depth is a dot product with a fixed isometric view.
        min_depth=-2.0_wp*half_width
        max_depth=2.0_wp*half_width
        path=cavity3d_directory//trim(names(id))//"_response.svg"
        open(newunit=u,file=path,status="replace",action="write")
        write(u,"(a)") '<?xml version="1.0" encoding="UTF-8"?>'
        write(u,"(a)") '<svg xmlns="http://www.w3.org/2000/svg" width="1080" height="810" viewBox="0 20 1080 810">'
        write(u,"(a)") '<rect x="0" y="20" width="1080" height="810" fill="#ffffff"/>'
        write(u,"(a)") "<desc>Prescribed local response of "//trim(names(id))//"; "//trim(label)// &
            " ["//trim(units)//"]. Three-dimensional volume cutaway or cohesive interface.</desc>"
        write(u,"(a)") '<g stroke-width="0.4" stroke-linejoin="round">'
        do bucket=0,depth_buckets-1
            view_faces: do f=1,n_faces
                if (hits(f)/=1) cycle view_faces
                center=sum(mesh%nodes(:,vertices(:,f)),dim=2)/4.0_wp
                depth=(center(1)+center(2)+0.9_wp*center(3))/sqrt(2.81_wp)
                k=min(depth_buckets-1,max(0,int(depth_buckets*(depth-min_depth)/ &
                    (max_depth-min_depth))))
                if (k/=bucket) cycle view_faces
                fraction=(values(owner(f))-lo)/(hi-lo)
                color=colormap_color(cmap,fraction)
                write(u,"(a)",advance="no") '<polygon fill="'//color//'" stroke="'//color//'" points="'
                do p=1,4
                    screen=projection(mesh%nodes(:,vertices(p,f)),is_interface)
                    write(u,"(a)",advance="no") trim(decimal(screen(1)))//","//trim(decimal(screen(2)))//" "
                end do
                write(u,"(a)") '"/>'
            end do view_faces
        end do
        write(u,"(a)") "</g>"
        write(u,"(a)") '<text x="875" y="150" font-family="Georgia,serif" font-size="19" fill="#1e293b">'// &
            trim(label)//"</text>"
        write(u,"(a)") '<text x="875" y="179" font-family="sans-serif" font-size="15" fill="#475569">['// &
            trim(units)//"]</text>"
        do m=0,255
            fraction=1.0_wp-real(m,wp)/255.0_wp
            color=colormap_color(cmap,fraction)
            write(u,"(a)") '<rect x="881" y="'//trim(decimal(210.0_wp+1.65_wp*real(m,wp)))// &
                '" width="27" height="1.72" fill="'//color//'"/>'
        end do
        write(u,"(a)") '<rect x="881" y="210" width="27" height="421" fill="none" stroke="#334155"/>'
        do m=0,4
            tick_value=hi-(hi-lo)*real(m,wp)/4.0_wp
            write(u,"(a)") '<text x="925" y="'//trim(decimal(216.0_wp+420.0_wp*real(m,wp)/4.0_wp))// &
                '" font-size="15" font-family="sans-serif" fill="#334155">'//trim(scientific(tick_value))//"</text>"
        end do
        write(u,"(a)") "</svg>"
        close(u)
    end subroutine write_cavity_preview

    pure function sorted_four(nodes) result(sorted)
        integer, intent(in) :: nodes(4)
        integer :: i, j, temp
        integer :: sorted(4)

        sorted=nodes
        do i=2,4
            temp=sorted(i)
            j=i-1
            insertion_scan: do while(j>=1)
                if (sorted(j)<=temp) exit insertion_scan
                sorted(j+1)=sorted(j)
                j=j-1
            end do insertion_scan
            sorted(j+1)=temp
        end do
    end function sorted_four

    pure integer function hash_slot(sorted) result(slot)
        integer, intent(in) :: sorted(4)
        integer(lk) :: key

        key=int(sorted(1),lk)*73856093_lk+int(sorted(2),lk)*19349663_lk+ &
            int(sorted(3),lk)*83492791_lk+int(sorted(4),lk)*2654435761_lk
        slot=int(modulo(key,int(hash_size,lk)))+1
    end function hash_slot

    pure function projection(position,is_interface) result(screen)
        real(wp), intent(in) :: position(3)
        logical, intent(in) :: is_interface
        real(wp) :: screen(2),q(3)
        real(wp) :: view_scale

        ! Uniform display magnification for the inner interface; VTK stays in SI units.
        view_scale=half_width
        if (is_interface) view_scale=cavity_radius
        q=position/view_scale
        screen(1)=430.0_wp+235.0_wp*(q(2)-q(1))/sqrt(2.0_wp)
        screen(2)=420.0_wp+235.0_wp*(0.38_wp*(q(1)+q(2))-0.84_wp*q(3))
    end function projection

    subroutine write_cavity_heatmap(mesh,id,is_interface,parameters,metrics,scaled_fields)
        !! At each layer (bulk), or z-direction bin (interface), plot
        !! the weighted RMS of the scaled parameter derivative relative to
        !! the global weighted RMS of the response. All derivatives are analytic.
        type(cavity_mesh_t), intent(in) :: mesh
        integer, intent(in) :: id
        logical, intent(in) :: is_interface
        character(*), intent(in) :: parameters(:)
        type(sensitivity_metric_t), intent(in) :: metrics(:)
        real(wp), intent(in) :: scaled_fields(:,:)
        real(wp), allocatable :: ratio(:,:)
        integer, allocatable :: category(:,:)
        real(wp) :: reference, magnitude, lo, hi
        real(wp) :: weights(n_layers), sum_square(n_layers)
        real(wp) :: extent,coordinate,t,x,y,tick
        integer :: n,np,j,k,p,u,m,rows,legend_height,legend_bottom
        integer, parameter :: plot_y=46, row_height=48, plot_x=326, bin_width=65, legend_x=1151
        type(Colormap) :: cmap
        character(7) :: color
        character(:), allocatable :: base
        character(24) :: log_label
        character(:), allocatable :: axis_label

        np=size(scaled_fields,1)
        n=size(parameters)
        if (n<1.or.n/=size(metrics).or.size(scaled_fields,2)/=n) then
            error stop "Invalid cavity heatmap metadata"
        end if
        if (is_interface) then
            if (np/=mesh%n_faces) error stop "Invalid cavity interface heatmap size"
        else
            if (np/=mesh%n_cells) error stop "Invalid cavity volume heatmap size"
        end if
        if (.not.all(ieee_is_finite(scaled_fields))) error stop "Nonfinite cavity analytical fields"

        allocate(ratio(n_layers,n),category(n_layers,n))
        ratio=0.0_wp
        category=undefined_ratio
        parameter_rows: do j=1,n
            weights=0.0_wp
            sum_square=0.0_wp
            reference=metrics(j)%response_rms
            magnitude=maxval(abs(scaled_fields(:,j)))
            if (reference<=0.0_wp) cycle parameter_rows
            do p=1,np
                k=sample_bin(mesh,is_interface,p)
                if (is_interface) then
                    extent=mesh%face_areas(p)
                else
                    extent=mesh%volumes(p)
                end if
                weights(k)=weights(k)+extent
                if (magnitude>0.0_wp) then
                    sum_square(k)=sum_square(k)+extent*(scaled_fields(p,j)/magnitude)**2
                end if
            end do
            sample_bins: do k=1,n_layers
                if (weights(k)<=0.0_wp) cycle sample_bins
                category(k,j)=exact_zero
                if (magnitude<=0.0_wp) cycle sample_bins
                ratio(k,j)=(magnitude/reference)*sqrt(sum_square(k)/weights(k))
                if (ratio(k,j)>0.0_wp) category(k,j)=resolved
                if (ratio(k,j)>0.0_wp.and.ratio(k,j)<=256.0_wp*epsilon(1.0_wp)) then
                    category(k,j)=below_resolution
                end if
            end do sample_bins
        end do parameter_rows

        lo=huge(1.0_wp)
        hi=-huge(1.0_wp)
        value_rows: do j=1,n
            value_bins: do k=1,n_layers
                if (category(k,j)/=resolved) cycle value_bins
                lo=min(lo,log10(ratio(k,j)))
                hi=max(hi,log10(ratio(k,j)))
            end do value_bins
        end do value_rows
        if (hi<lo) then
            lo=-1.0_wp
            hi=1.0_wp
        else
            lo=real(floor(lo),wp)
            hi=real(ceiling(hi),wp)
            if (hi<=lo) hi=lo+1.0_wp
        end if

        base=cavity3d_directory//trim(names(id))//"_sensitivity"
        open(newunit=u,file=base//".csv",status="replace",action="write")
        write(u,"(a)") "parameter_index,parameter,parameter_units,bin,coordinate,coordinate_units,relative_rms,log10_relative_rms,status"
        do j=1,n
            do k=1,n_layers
                coordinate=(real(k,wp)-0.5_wp)/real(n_layers,wp)
                if (is_interface) coordinate=2.0_wp*coordinate-1.0_wp
                log_label="NA"
                if (ratio(k,j)>0.0_wp) write(log_label,"(es21.12e3)") log10(ratio(k,j))
                write(u,'(i0,",",a,",",a,",",i0,",",es22.14e3,",1,",es22.14e3,",",a,",",a)') &
                    j,trim(parameters(j)),trim(metrics(j)%units),k,coordinate,ratio(k,j), &
                    trim(adjustl(log_label)),trim(category_name(category(k,j)))
            end do
        end do
        close(u)

        call init_colormap(cm_viridis,cmap)
        axis_label="&#955; [-]"
        if (is_interface) axis_label="n_z [-]"
        legend_height=max(120,n*row_height-6)
        legend_bottom=plot_y+legend_height
        rows=max(plot_y+n*row_height+90,legend_bottom+98)
        open(newunit=u,file=base//".svg",status="replace",action="write")
        write(u,"(a)") '<?xml version="1.0" encoding="UTF-8"?>'
        write(u,"(a)") '<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="'// &
            trim(int_string(rows))//'" viewBox="0 0 1400 '//trim(int_string(rows))//'">'
        write(u,"(a)") '<rect width="1400" height="100%" fill="#ffffff"/>'
        write(u,"(a)") "<desc>Logarithmic ratio of the binned weighted RMS of each scaled analytical"// &
            " parameter derivative to the global response RMS. Colorbar and categories give numerical"// &
            " values and resolution status.</desc>"
        do j=1,n
            y=real(plot_y+(j-1)*row_height,wp)
            if (mod(j,2)==0) then
                write(u,"(a)") '<rect x="12" y="'//trim(decimal(y))// &
                    '" width="1115" height="'//trim(int_string(row_height))//'" fill="#f3f6f9"/>'
            end if
            write(u,"(a)") '<text x="18" y="'//trim(decimal(y+31.0_wp))// &
                '" font-family="sans-serif" font-size="15" fill="#25364b">P'//trim(int_string(j))// &
                " "//trim(parameters(j))//" ["//trim(metrics(j)%units)//"]</text>"
            do k=1,n_layers
                x=real(plot_x+(k-1)*bin_width,wp)
                select case(category(k,j))
                case(exact_zero)
                    color="#d7dbe2"
                case(below_resolution)
                    color="#bac5d1"
                case(undefined_ratio)
                    color="#e5b57b"
                case(resolved)
                    color=colormap_color(cmap,(log10(ratio(k,j))-lo)/(hi-lo))
                case default
                    error stop "Invalid cavity heatmap classification"
                end select
                write(u,"(a)") '<g class="sensitivity-cell" data-parameter="'// &
                    trim(int_string(j))//'" data-bin="'//trim(int_string(k))//'">'
                write(u,"(a)") "<title>Scaled relative RMS: "//trim(scientific(ratio(k,j)))// &
                    "; status: "//trim(category_name(category(k,j)))//"</title>"
                write(u,"(a)") '<rect x="'//trim(decimal(x))//'" y="'//trim(decimal(y+3.0_wp))// &
                    '" width="64.2" height="44" fill="'//color//'"/>'
                write(u,"(a)") "</g>"
            end do
        end do
        do m=0,4
            x=real(plot_x,wp)+real(n_layers*bin_width,wp)*real(m,wp)/4.0_wp
            coordinate=real(m,wp)/4.0_wp
            if (is_interface) coordinate=2.0_wp*coordinate-1.0_wp
            write(u,"(a)") '<line x1="'//trim(decimal(x))//'" y1="'// &
                trim(int_string(plot_y+n*row_height))//'" x2="'//trim(decimal(x))//'" y2="'// &
                trim(int_string(plot_y+n*row_height+6))//'" stroke="#334155"/>'
            write(u,"(a)") '<text x="'//trim(decimal(x))//'" y="'// &
                trim(int_string(plot_y+n*row_height+27))// &
                '" text-anchor="middle" font-family="sans-serif" font-size="13" fill="#334155">'// &
                trim(decimal(coordinate))//"</text>"
        end do
        write(u,"(a)") '<text x="'//trim(int_string(plot_x+n_layers*bin_width/2))//'" y="'// &
            trim(int_string(plot_y+n*row_height+64))// &
            '" font-family="sans-serif" text-anchor="middle" font-size="15" fill="#334155">'// &
            axis_label//"</text>"
        write(u,"(a)") '<text x="'//trim(int_string(legend_x))//'" y="'// &
            trim(int_string(plot_y-16))// &
            '" font-family="sans-serif" font-size="14" fill="#334155">log10 relative RMS</text>'
        do m=0,255
            t=1.0_wp-real(m,wp)/255.0_wp
            color=colormap_color(cmap,t)
            write(u,"(a)") '<rect x="'//trim(int_string(legend_x))//'" y="'// &
                trim(decimal(real(plot_y,wp)+real(legend_height,wp)*real(m,wp)/255.0_wp))// &
                '" width="24" height="'//trim(decimal(real(legend_height,wp)/255.0_wp+0.08_wp))// &
                '" fill="'//color//'"/>'
        end do
        write(u,"(a)") '<rect x="'//trim(int_string(legend_x))//'" y="'// &
            trim(int_string(plot_y))//'" width="24" height="'// &
            trim(int_string(legend_height))//'" fill="none" stroke="#334155"/>'
        do m=0,4
            tick=hi-(hi-lo)*real(m,wp)/4.0_wp
            y=real(plot_y+5,wp)+real(legend_height,wp)*real(m,wp)/4.0_wp
            write(u,"(a)") '<text x="'//trim(int_string(legend_x+34))//'" y="'//trim(decimal(y))// &
                '" font-family="sans-serif" font-size="13" fill="#334155">'//trim(decimal(tick))//"</text>"
        end do
        write(u,"(a)") '<rect x="'//trim(int_string(legend_x))//'" y="'// &
            trim(int_string(legend_bottom+18))//'" width="18" height="15" fill="#d7dbe2"/>'
        write(u,"(a)") '<text x="'//trim(int_string(legend_x+28))//'" y="'// &
            trim(int_string(legend_bottom+31))// &
            '" font-family="sans-serif" font-size="12" fill="#334155">exact zero</text>'
        write(u,"(a)") '<rect x="'//trim(int_string(legend_x))//'" y="'// &
            trim(int_string(legend_bottom+43))//'" width="18" height="15" fill="#bac5d1"/>'
        write(u,"(a)") '<text x="'//trim(int_string(legend_x+28))//'" y="'// &
            trim(int_string(legend_bottom+56))// &
            '" font-family="sans-serif" font-size="12" fill="#334155">below resolution</text>'
        write(u,"(a)") '<rect x="'//trim(int_string(legend_x))//'" y="'// &
            trim(int_string(legend_bottom+68))//'" width="18" height="15" fill="#e5b57b"/>'
        write(u,"(a)") '<text x="'//trim(int_string(legend_x+28))//'" y="'// &
            trim(int_string(legend_bottom+81))// &
            '" font-family="sans-serif" font-size="12" fill="#334155">undefined</text>'
        write(u,"(a)") "</svg>"
        close(u)
    end subroutine write_cavity_heatmap

    pure integer function sample_bin(mesh,is_interface,p) result(k)
        type(cavity_mesh_t), intent(in) :: mesh
        logical, intent(in) :: is_interface
        integer, intent(in) :: p
        real(wp) :: z
        real(wp) :: normal(3)

        if (is_interface) then
            normal=unit_normal(mesh%face_centers(:,p))
            z=normal(3)
            k=min(n_layers,max(1,1+int((z+1.0_wp)*real(n_layers,wp)/2.0_wp)))
        else
            k=(p-1)/mesh%n_faces+1
        end if
    end function sample_bin

    pure function category_name(category) result(name)
        integer, intent(in) :: category
        character(18) :: name

        select case(category)
        case(exact_zero)
            name="exact_zero"
        case(below_resolution)
            name="below_resolution"
        case(undefined_ratio)
            name="undefined"
        case(resolved)
            name="ok"
        case default
            name="undefined"
        end select
    end function category_name

    function scientific(value) result(s)
        real(wp), intent(in) :: value
        character(24) :: s

        write(s,"(es13.5e2)") value
        s=adjustl(s)
    end function scientific

    function decimal(value) result(s)
        real(wp), intent(in) :: value
        character(24) :: s

        write(s,"(f13.3)") value
        s=adjustl(s)
    end function decimal

    function int_string(value) result(s)
        integer, intent(in) :: value
        character(18) :: s

        write(s,"(i0)") value
    end function int_string

end module benchmark_cavity3d_visual
