! SPDX-License-Identifier: Apache-2.0

module benchmark_cavity3d_export
    !! Numeric VTK/CSV outputs and generated per-material FORD result pages.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use benchmark_cavity3d_fields, only: cavity_displacement
    use benchmark_cavity3d_mesh, only: cavity_mesh_t, n_layers
    use benchmark_cavity3d_visual, only: write_cavity_preview, write_cavity_heatmap
    use benchmark_metrics, only: sensitivity_metric_t
    use benchmark_models, only: cases, names, load_mandel, load_F, load_thermoelastic, load_piezo, load_biot, &
                                latex_response, latex_units, latex_number
    use benchmark_paths, only: cavity3d_directory
    use constituf, only: wp
    implicit none (type, external)
    private

    public :: write_cavity_response, write_cavity_sensitivity

contains

    subroutine write_cavity_response(mesh,id,is_interface,field,label,units,amplitude)
        type(cavity_mesh_t), intent(in) :: mesh
        integer, intent(in) :: id
        logical, intent(in) :: is_interface
        real(wp), intent(in) :: amplitude
        real(wp), intent(in) :: field(:)
        character(*), intent(in) :: label,units
        integer :: p,u
        real(wp) :: x(3)
        character(:), allocatable :: base

        if (is_interface) then
            if (size(field)/=mesh%n_faces) error stop "Invalid cavity interface response size"
        else
            if (size(field)/=mesh%n_cells) error stop "Invalid cavity bulk response size"
        end if
        if (.not.all(ieee_is_finite(field))) error stop "Nonfinite cavity response"
        base=cavity3d_directory//trim(names(id))//"_response"

        open(newunit=u,file=base//".csv",status="replace",action="write")
        if (is_interface) then
            write(u,"(a)") "x_m,y_m,z_m,area_m2,response"
        else
            write(u,"(a)") "x_m,y_m,z_m,volume_m3,response"
        end if
        do p=1,size(field)
            if (is_interface) then
                x=mesh%face_centers(:,p)
                write(u,'(es24.16e3,",",es24.16e3,",",es24.16e3,",",es24.16e3,",",es24.16e3)') x,mesh%face_areas(p),field(p)
            else
                x=mesh%centers(:,p)
                write(u,'(es24.16e3,",",es24.16e3,",",es24.16e3,",",es24.16e3,",",es24.16e3)') x,mesh%volumes(p),field(p)
            end if
        end do
        close(u)

        call write_vtk(base//".vtk",mesh,id,is_interface,field,label,units,amplitude)
        call write_cavity_preview(mesh,id,is_interface,field,label,units)
    end subroutine write_cavity_response

    subroutine write_vtk(path,mesh,id,is_interface,field,label,units,amplitude)
        character(*), intent(in) :: path,label,units
        type(cavity_mesh_t), intent(in) :: mesh
        integer, intent(in) :: id
        logical, intent(in) :: is_interface
        real(wp), intent(in) :: amplitude
        real(wp), intent(in) :: field(:)
        integer :: u,p,j,nnode,ncell,ncorner,cell_type
        real(wp) :: displacement(3)
        logical :: mechanical

        nnode=mesh%n_nodes
        ncell=mesh%n_cells
        ncorner=8
        cell_type=12
        if (is_interface) then
            nnode=mesh%n_nodes/(n_layers+1)
            ncell=mesh%n_faces
            ncorner=4
            cell_type=9
        end if
        mechanical=any(cases(id)%input_kind==[load_mandel,load_F,load_thermoelastic,load_piezo,load_biot])
        open(newunit=u,file=path,status="replace",action="write")
        write(u,"(a)") "# vtk DataFile Version 3.0"
        write(u,"(a)") "ConstituF 3D spherical cavity; "//trim(label)//" ["//trim(units)//"]"
        write(u,"(a)") "ASCII"
        write(u,"(a)") "DATASET UNSTRUCTURED_GRID"
        write(u,"(a,i0,a)") "POINTS ",nnode," double"
        do p=1,nnode
            write(u,"(3(es24.16e3,1x))") mesh%nodes(:,p)
        end do
        write(u,"(a,i0,1x,i0)") "CELLS ",ncell,(ncorner+1)*ncell
        do p=1,ncell
            if (is_interface) then
                write(u,"(i0,1x,4(i0,1x))") 4,mesh%faces(:,p)-1
            else
                write(u,"(i0,1x,8(i0,1x))") 8,mesh%cells(:,p)-1
            end if
        end do
        write(u,"(a,i0)") "CELL_TYPES ",ncell
        do p=1,ncell
            write(u,"(i0)") cell_type
        end do
        if (mechanical) then
            write(u,"(a,i0)") "POINT_DATA ",nnode
            write(u,"(a)") "VECTORS prescribed_displacement_m double"
            do p=1,nnode
                displacement=cavity_displacement(mesh%nodes(:,p),amplitude)
                write(u,"(3(es24.16e3,1x))") displacement
            end do
        end if
        write(u,"(a,i0)") "CELL_DATA ",ncell
        write(u,"(a)") "SCALARS "//trim(label)//" double 1"
        write(u,"(a)") "LOOKUP_TABLE default"
        do j=1,ncell
            write(u,"(es24.16e3)") field(j)
        end do
        close(u)
    end subroutine write_vtk

    subroutine write_cavity_sensitivity(mesh,id,is_interface,parameter_names,metrics,scaled_fields)
        type(cavity_mesh_t), intent(in) :: mesh
        integer, intent(in) :: id
        logical, intent(in) :: is_interface
        character(*), intent(in) :: parameter_names(:)
        type(sensitivity_metric_t), intent(in) :: metrics(:)
        real(wp), intent(in) :: scaled_fields(:,:)
        integer :: u,p,j
        real(wp) :: x(3)
        character(:), allocatable :: path

        if (size(parameter_names)/=size(metrics).or.size(scaled_fields,2)/=size(metrics)) then
            error stop "Invalid cavity sensitivity parameters"
        end if
        path=cavity3d_directory//trim(names(id))//"_sensitivity_fields.csv"
        open(newunit=u,file=path,status="replace",action="write")
        write(u,"(a)",advance="no") "x_m,y_m,z_m"
        do j=1,size(metrics)
            write(u,"(a)",advance="no") ",P"//trim(int_string(j))//"_"//trim(parameter_names(j))
        end do
        write(u,"(a)") ""
        do p=1,size(scaled_fields,1)
            if (is_interface) then
                x=mesh%face_centers(:,p)
            else
                x=mesh%centers(:,p)
            end if
            write(u,'(es24.16e3,",",es24.16e3,",",es24.16e3)',advance="no") x
            do j=1,size(metrics)
                write(u,'(",",es24.16e3)',advance="no") scaled_fields(p,j)
            end do
            write(u,"(a)") ""
        end do
        close(u)

        call write_cavity_heatmap(mesh,id,is_interface,parameter_names,metrics,scaled_fields)
        call write_material_page(id,is_interface,parameter_names,metrics)
    end subroutine write_cavity_sensitivity

    subroutine write_material_page(id,is_interface,parameters,metrics)
        integer, intent(in) :: id
        logical, intent(in) :: is_interface
        character(*), intent(in) :: parameters(:)
        type(sensitivity_metric_t), intent(in) :: metrics(:)
        integer :: u,j
        character(:), allocatable :: slug,ratio

        slug=trim(names(id))
        open(newunit=u,file="docs/pages/benchmark-cavity3d-"//slug//".md",status="replace",action="write")
        write(u,"(a)") "---"
        write(u,"(a)") "title: "//slug//" - 3D cavity response and sensitivity"
        write(u,"(a)") "---"
        write(u,"(a)") ""
        write(u,"(a)") "These results sample the constitutive law at prescribed fields in the"
        if (is_interface) then
            write(u,"(a)") "spherical cavity interface (surface quadrilaterals)."
        else
            write(u,"(a)") "cubed-sphere hexahedral volume."
        end if
        write(u,"(a)") "They are not a boundary-value solution."
        write(u,"(a)") "See the [3D cavity formulation](|page|/benchmark-cavity3d.html)."
        write(u,"(a)") ""
        write(u,"(a)") "Response component: "//latex_response(cases(id)%response_label)// &
            " ("//latex_units(cases(id)%response_units)//")."
        write(u,"(a)") ""
        write(u,"(a)") "## Constitutive response"
        write(u,"(a)") ""
        if (is_interface) then
            write(u,"(a)") "![Spherical-interface response](|media|/benchmarks/cavity3d/"//slug//"_response.svg)"
        else
            write(u,"(a)") "![Three-dimensional response cutaway](|media|/benchmarks/cavity3d/"//slug//"_response.svg)"
        end if
        write(u,"(a)") ""
        if (is_interface) then
            write(u,"(a)") "[Interface VTK field](|media|/benchmarks/cavity3d/"//slug//"_response.vtk) | "// &
                "[Response samples (CSV)](|media|/benchmarks/cavity3d/"//slug//"_response.csv)"
        else
            write(u,"(a)") "[Volume VTK field](|media|/benchmarks/cavity3d/"//slug//"_response.vtk) | "// &
                "[Response samples (CSV)](|media|/benchmarks/cavity3d/"//slug//"_response.csv)"
        end if
        write(u,"(a)") ""
        write(u,"(a)") "## Analytical parameter sensitivity"
        write(u,"(a)") ""
        write(u,"(a)") "The rows correspond to the declared scalar parameters. Each cell contains a"
        write(u,"(a)") "weighted bin RMS of the scaled derivative divided by the global response RMS."
        write(u,"(a)") "The derivative includes committed material history, when present, and holds"
        write(u,"(a)") "the complete prescribed input history fixed."
        write(u,"(a)") ""
        write(u,"(a)") "![Per-parameter analytical sensitivity heatmap](|media|/benchmarks/cavity3d/"// &
            slug//"_sensitivity.svg)"
        write(u,"(a)") ""
        write(u,"(a)") "| Parameter | Unit | Nominal value | Reference scale | Relative RMS |"
        write(u,"(a)") "| --- | --- | ---: | ---: | ---: |"
        do j=1,size(parameters)
            ratio="undefined"
            if (metrics(j)%relative_defined) ratio=latex_number(metrics(j)%relative_rms)
            write(u,"(a)") "| `"//trim(parameters(j))//"` | "//latex_units(metrics(j)%units)//" | "// &
                latex_number(metrics(j)%reference_value)//" | "// &
                latex_number(metrics(j)%parameter_scale)//" | "//ratio//" |"
        end do
        write(u,"(a)") ""
        write(u,"(a)") "[Heatmap data (CSV)](|media|/benchmarks/cavity3d/"//slug//"_sensitivity.csv) | "// &
            "[Signed scaled derivative fields (CSV)](|media|/benchmarks/cavity3d/"// &
            slug//"_sensitivity_fields.csv)"
        write(u,"(a)") ""
        write(u,"(a)") "[All 3D cavity material results](|page|/benchmark-cavity3d-results.html)"
        close(u)
    end subroutine write_material_page

    function int_string(value) result(s)
        integer, intent(in) :: value
        character(18) :: s

        write(s,"(i0)") value
    end function int_string

end module benchmark_cavity3d_export
