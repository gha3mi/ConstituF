! SPDX-License-Identifier: Apache-2.0

module benchmark_annulus_export
    !! Numeric VTK/CSV outputs and generated per-material FORD result pages.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use benchmark_annulus_fields, only: annulus_displacement
    use benchmark_annulus_mesh, only: interface_radius, nradial, nangle, node_radius, node_angle, cell_radius, cell_angle
    use benchmark_annulus_visual, only: write_annulus_preview, write_annulus_heatmap
    use benchmark_metrics, only: sensitivity_metric_t
    use benchmark_models, only: names, interface_case, latex_response, latex_units, latex_number
    use benchmark_paths, only: annulus_directory
    use constituf, only: wp
    implicit none (type, external)
    private

    public :: write_annulus_response, write_annulus_sensitivity

contains

    subroutine write_annulus_response(id, values, label, units)
        integer, intent(in) :: id
        real(wp), intent(in) :: values(:)
        character(*), intent(in) :: label, units
        character(:), allocatable :: base

        if (.not. all(ieee_is_finite(values))) error stop "Annulus field contains nonfinite data"
        if (interface_case(id)) then
            if (size(values) /= nangle) error stop "Incorrect interface output size"
        else
            if (size(values) /= nangle*nradial) error stop "Incorrect annulus output size"
        end if

        base = annulus_directory//trim(names(id))//"_response"
        call write_field_csv(base//".csv",id,values)
        call write_field_vtk(base//".vtk",id,values,label)
        call write_annulus_preview(base//".svg",id,values,label,units)
    end subroutine write_annulus_response


    subroutine write_field_csv(path,id,values)
        character(*), intent(in) :: path
        integer, intent(in) :: id
        real(wp), intent(in) :: values(:)
        integer :: u, p, i, j
        real(wp) :: r

        open(newunit=u,file=path,status="replace",action="write")
        write(u,"(a)") "radius_m,angle_rad,response"
        do p=1,size(values)
            j=mod(p-1,nangle)+1
            i=(p-1)/nangle+1
            r=cell_radius(i)
            if (interface_case(id)) r=interface_radius
            write(u,'(es24.16e3,",",es24.16e3,",",es24.16e3)') r,cell_angle(j),values(p)
        end do
        close(u)
    end subroutine write_field_csv

    subroutine write_field_vtk(path,id,values,label)
        character(*), intent(in) :: path,label
        integer, intent(in) :: id
        real(wp), intent(in) :: values(:)
        integer :: u, i, j, m, ncell, nnode, id0, id1, id2, id3
        real(wp) :: r, angle, displacement

        open(newunit=u,file=path,status="replace",action="write")
        write(u,"(a)") "# vtk DataFile Version 3.0"
        write(u,"(a)") "ConstituF annulus (SI units, see benchmark page)"
        write(u,"(a)") "ASCII"
        if (interface_case(id)) then
            write(u,"(a)") "DATASET POLYDATA"
            write(u,"(a,i0,a)") "POINTS ",nangle," double"
            do j=1,nangle
                angle=node_angle(j)
                write(u,"(3(es24.16e3,1x))") interface_radius*cos(angle), &
                    interface_radius*sin(angle),0.0_wp
            end do
            write(u,"(a,i0,1x,i0)") "LINES ",nangle,3*nangle
            do j=1,nangle
                write(u,"(a,i0,1x,i0)") "2 ",j-1,mod(j,nangle)
            end do
            ncell=nangle
        else
            nnode=(nradial+1)*nangle
            ncell=nradial*nangle
            write(u,"(a)") "DATASET UNSTRUCTURED_GRID"
            write(u,"(a,i0,a)") "POINTS ",nnode," double"
            do i=1,nradial+1
                r=node_radius(i)
                do j=1,nangle
                    angle=node_angle(j)
                    write(u,"(3(es24.16e3,1x))") r*cos(angle),r*sin(angle),0.0_wp
                end do
            end do
            write(u,"(a,i0,1x,i0)") "CELLS ",ncell,5*ncell
            do i=1,nradial
                do j=1,nangle
                    m=mod(j,nangle)
                    id0=(i-1)*nangle+j-1
                    id1=i*nangle+j-1
                    id2=i*nangle+m
                    id3=(i-1)*nangle+m
                    write(u,"(a,4(i0,1x))") "4 ",id0,id1,id2,id3
                end do
            end do
            write(u,"(a,i0)") "CELL_TYPES ",ncell
            do m=1,ncell
                write(u,"(a)") "9"
            end do
            write(u,"(a,i0)") "POINT_DATA ",nnode
            write(u,"(a)") "VECTORS displacement_m double"
            do i=1,nradial+1
                r=node_radius(i)
                displacement=annulus_displacement(id,r)
                do j=1,nangle
                    angle=node_angle(j)
                    write(u,"(3(es24.16e3,1x))") displacement*cos(angle), &
                        displacement*sin(angle),0.0_wp
                end do
            end do
        end if
        write(u,"(a,i0)") "CELL_DATA ",ncell
        write(u,"(a)") "SCALARS "//trim(label)//" double 1"
        write(u,"(a)") "LOOKUP_TABLE default"
        do m=1,ncell
            write(u,"(es24.16e3)") values(m)
        end do
        close(u)
    end subroutine write_field_vtk

    subroutine write_annulus_sensitivity(id,parameter_names,parameter_units,values, &
                                        metrics,response,scaled_fields,response_name,response_units,domain,steps)
        !! For response \(y\) and scalar parameter \(p_j\), `scaled_fields(point,j)` stores
        !! \(s_j\,\mathrm{d}y/\mathrm{d}p_j\), including derivatives of committed history.
        !! The reference scale \(s_j\) has the same units as \(p_j\).
        !! No finite differences, artificial spatial PDE solve, or global
        !! structural design-sensitivity claim is made.
        integer, intent(in) :: id,steps
        character(*), intent(in) :: parameter_names(:),parameter_units(:)
        real(wp), intent(in) :: values(:),response(:),scaled_fields(:,:)
        type(sensitivity_metric_t), intent(in) :: metrics(:)
        character(*), intent(in) :: response_name,response_units,domain
        integer :: n,np
        character(:), allocatable :: base

        n=size(parameter_names)
        np=size(response)
        if (n<1.or.size(parameter_units)/=n.or.size(values)/=n.or.size(metrics)/=n) then
            error stop "Material benchmark parameter metadata mismatch"
        end if
        if (size(scaled_fields,1)/=np.or.size(scaled_fields,2)/=n) then
            error stop "Material benchmark sensitivity field mismatch"
        end if
        if (.not.all(ieee_is_finite(scaled_fields)).or..not.all(ieee_is_finite(response))) then
            error stop "Material benchmark contains nonfinite results"
        end if
        if (interface_case(id)) then
            if (np/=nangle) error stop "Incorrect interface benchmark point count"
        else
            if (np/=nangle*nradial) error stop "Incorrect bulk benchmark point count"
        end if

        base=annulus_directory//trim(names(id))
        call write_annulus_heatmap(id,parameter_names,parameter_units,metrics,scaled_fields)
        call write_material_csv(base//"_sensitivity_fields.csv",id,parameter_names,response,scaled_fields)
        call write_material_page(id,parameter_names,parameter_units,values,metrics, &
                                 response_name,response_units,domain,steps)
    end subroutine write_annulus_sensitivity

    subroutine write_material_csv(filename,id,parameters,response,scaled_fields)
        !! Exact scaled derivatives, without a second response column.
        !! The companion _sensitivity.csv holds normalized heatmap values.
        character(*), intent(in) :: filename
        character(*), intent(in) :: parameters(:)
        integer, intent(in) :: id
        real(wp), intent(in) :: response(:),scaled_fields(:,:)
        integer :: unit,i,j,k
        real(wp) :: radius

        open(newunit=unit,file=filename,status="replace",action="write")
        write(unit,"(a)",advance="no") "radius_m,angle_rad"
        do j=1,size(parameters)
            write(unit,"(a)",advance="no") ",P"//trim(integer_string(j))//"_"//trim(parameters(j))
        end do
        write(unit,"(a)") ""
        do k=1,size(response)
            i=(k-1)/nangle+1
            radius=cell_radius(i)
            if (interface_case(id)) radius=interface_radius
            write(unit,'(es24.16e3,",",es24.16e3)',advance="no") &
                radius,cell_angle(mod(k-1,nangle)+1)
            do j=1,size(parameters)
                write(unit,'(",",es24.16e3)',advance="no") scaled_fields(k,j)
            end do
            write(unit,"(a)") ""
        end do
        close(unit)
    end subroutine write_material_csv

    subroutine write_material_page(id,parameters,units,values,metrics,response,response_units,domain,steps)
        !! Summarize each case in the generated FORD reference.
        integer, intent(in) :: id,steps
        character(*), intent(in) :: response,response_units,domain
        character(*), intent(in) :: parameters(:),units(:)
        real(wp), intent(in) :: values(:)
        type(sensitivity_metric_t), intent(in) :: metrics(:)
        integer :: unit,j
        character(:), allocatable :: filename,label,ratio
        character(18) :: status_label

        filename="docs/pages/benchmark-annulus-"//trim(names(id))//".md"
        label=material_label(names(id))
        open(newunit=unit,file=filename,status="replace",action="write")
        write(unit,"(a)") "---"
        write(unit,"(a)") "title: "//label//" - constitutive response and analytical sensitivity"
        write(unit,"(a)") "---"
        write(unit,"(a)") ""
        write(unit,"(a)") "Two analyses of the same prescribed annular loading: the constitutive response"
        write(unit,"(a)") "and the analytical derivatives returned by ConstituF for every declared parameter."
        write(unit,"(a)") ""
        write(unit,"(a)") "Response component: "//latex_response(response)//" ("//latex_units(response_units)//")."
        write(unit,"(a)") "Domain: "//trim(domain)//"; accepted loading steps: "//trim(integer_string(steps))//"."
        write(unit,"(a)") "Differentiation holds the prescribed input history fixed and includes committed material history, where present."
        write(unit,"(a)") ""
        write(unit,"(a)") "## 1. Constitutive response"
        write(unit,"(a)") ""
        write(unit,"(a)") "The response benchmark evaluates the material using ConstituF and plots"
        write(unit,"(a)") "the selected physical response at each prescribed sampling point."
        write(unit,"(a)") ""
        write(unit,"(a)") "![Response contour for "//label// &
            "](|media|/benchmarks/annulus/"//trim(names(id))//"_response.svg)"
        write(unit,"(a)") ""
        write(unit,"(a)") "[Constitutive response samples (CSV)](|media|/benchmarks/annulus/"// &
            trim(names(id))//"_response.csv) | [Response field (VTK)](|media|/benchmarks/annulus/"// &
            trim(names(id))//"_response.vtk)"
        write(unit,"(a)") ""
        write(unit,"(a)") "## 2. Analytical parameter sensitivities"
        write(unit,"(a)") ""
        write(unit,"(a)") "The sensitivity benchmark calls ConstituF's analytical derivative interface."
        write(unit,"(a)") "Rows correspond to declared scalar parameters; columns are radial positions"
        write(unit,"(a)") "or interface angles. Color represents the base-10 logarithm of the"
        write(unit,"(a)") "scaled derivative RMS divided by the spatial response RMS."
        write(unit,"(a)") "Zero, below-resolution and undefined values have separate categories."
        write(unit,"(a)") "See [sensitivity normalization](|page|/benchmark-annulus.html) for definitions."
        write(unit,"(a)") ""
        write(unit,"(a)") "![Parameter sensitivity heatmap for "//label// &
            "](|media|/benchmarks/annulus/"//trim(names(id))//"_sensitivity.svg)"
        write(unit,"(a)") ""
        write(unit,"(a)") "## Reference parameters"
        write(unit,"(a)") ""
        write(unit,"(a)") "| Index | Parameter | Unit | Value | Scale | Relative RMS | Status |"
        write(unit,"(a)") "| --- | --- | --- | ---: | ---: | ---: | --- |"
        do j=1,size(parameters)
            ratio="undefined"
            if (metrics(j)%relative_defined) ratio=latex_number(metrics(j)%relative_rms)
            status_label="ok"
            if (.not.metrics(j)%relative_defined) then
                status_label="undefined"
            else if (metrics(j)%sensitivity_rms<=0.0_wp) then
                status_label="exact zero"
            else if (metrics(j)%below_resolution) then
                status_label="below resolution"
            end if
            write(unit,"(a)") "| `P"//trim(integer_string(j))//"` | `"//trim(parameters(j))//"` | "// &
                latex_units(units(j))//" | "//latex_number(values(j))//" | "// &
                latex_number(metrics(j)%parameter_scale)//" | "//ratio//" | "//trim(status_label)//" |"
        end do
        write(unit,"(a)") ""
        write(unit,"(a)") "[Signed spatial derivatives (CSV)](|media|/benchmarks/annulus/"// &
            trim(names(id))//"_sensitivity_fields.csv)"
        write(unit,"(a)") ""
        write(unit,"(a)") "[Sensitivity heatmap data (CSV)](|media|/benchmarks/annulus/"// &
            trim(names(id))//"_sensitivity.csv)"
        write(unit,"(a)") ""
        write(unit,"(a)") "[Benchmark formulation](|page|/benchmark-annulus.html) | "// &
            "[Material results](|page|/benchmark-annulus-results.html)"
        close(unit)
    end subroutine write_material_page

    function material_label(slug) result(label)
        !! Human-readable heading from the stable underscore-separated case ID.
        character(*), intent(in) :: slug
        character(:), allocatable :: label
        integer :: i, first

        label=trim(slug)
        do i=1,len(label)
            if (label(i:i)=="_") label(i:i)=" "
        end do
        if (len(label)==0) return
        first=iachar(label(1:1))
        if (first>=iachar("a").and.first<=iachar("z")) label(1:1)=achar(first-32)
    end function material_label

    function integer_string(value) result(s)
        integer, intent(in) :: value
        character(18) :: s

        write(s,"(i0)") value
    end function integer_string

end module benchmark_annulus_export
