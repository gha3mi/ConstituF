! SPDX-License-Identifier: Apache-2.0

module benchmark_cavity3d_runner
    !! Material-point response and analytical sensitivity on a common 3D mesh.
    !! Inputs are prescribed; derivatives are computed only by ConstituF.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use benchmark_cavity3d_export, only: write_cavity_response, write_cavity_sensitivity
    use benchmark_cavity3d_fields, only: supports_cavity3d, fill_cavity_inputs, project_cavity_response
    use benchmark_cavity3d_mesh, only: cavity_mesh_t, build_cavity_mesh
    use benchmark_colormaps, only: colormap_for, colormap_name, colormap_limits
    use benchmark_history, only: path
    use benchmark_metrics, only: sensitivity_metric_t, calculate_metric, parameter_scale
    use benchmark_models, only: ncases, names, max_parameters, cases, initialize_model, &
                              reference_parameters, interface_case, history_case
    use benchmark_pages, only: write_benchmark_index
    use benchmark_paths, only: cavity3d_directory
    use constituf, only: wp, material_t, cf_ok
    implicit none (type, external)
    private

    public :: run_cavity3d

contains

    subroutine run_cavity3d(mode,model_filter)
        !! Modes: all, response, sensitivity. An optional exact material slug
        !! selects only that model, without overwriting aggregate CSV reports.
        character(*), intent(in) :: mode
        character(*), intent(in), optional :: model_filter
        type(cavity_mesh_t) :: mesh
        type(material_t) :: material
        type(sensitivity_metric_t), allocatable :: metrics(:)
        real(wp), allocatable :: positions(:,:),input(:,:),output(:,:),derivatives(:,:,:)
        real(wp), allocatable :: field(:),derivative_field(:),scaled_fields(:,:),weights(:)
        real(wp) :: amplitude, scale, lo, hi
        real(wp) :: nominal(max_parameters)
        character(48) :: parameter_names(max_parameters)
        character(24) :: parameter_units(max_parameters)
        character(:), allocatable :: pname
        character(20) :: domain,status_label
        character(24) :: log_label
        integer :: id,npoint,ninput,noutput,step,nsteps,j,n,code,report_response,report_sens
        integer :: first,last,selected_id,completed,scheme
        logical :: response_mode,sensitivity_mode,filtered,found,is_interface

        response_mode=mode=="all".or.mode=="response"
        sensitivity_mode=mode=="all".or.mode=="sensitivity"
        if (.not.response_mode.and..not.sensitivity_mode) then
            error stop "Expected cavity3d mode: all, response, sensitivity"
        end if

        call build_cavity_mesh(mesh)
        filtered=present(model_filter)
        selected_id=0
        if (filtered) then
            do id=1,ncases
                if (trim(model_filter)==trim(names(id))) selected_id=id
            end do
            if (selected_id==0) error stop "Unknown 3D cavity material"
            if (.not.supports_cavity3d(selected_id)) error stop "Selected model is not a three-dimensional law"
        end if

        first=1
        last=ncases
        if (filtered) then
            first=selected_id
            last=selected_id
        end if

        if (response_mode.and..not.filtered) then
            open(newunit=report_response,file=cavity3d_directory//"response_manifest.csv", &
                status="replace",action="write")
            write(report_response,"(a)") &
                "model,domain,quantity,units,npoints,minimum,maximum,colormap,scale_min,scale_max"
        end if
        if (sensitivity_mode.and..not.filtered) then
            open(newunit=report_sens,file=cavity3d_directory//"sensitivity_metrics.csv", &
                status="replace",action="write")
            write(report_sens,"(a)") &
                "model,parameter_index,parameter,parameter_units,parameter_value,reference_scale,scale_basis,response,response_units,domain,steps,npoints,response_rms,scaled_sensitivity_rms,relative_rms,log10_relative_rms,sensitivity_linf,signed_mean,sign_coherence,status"
        end if

        completed=0
        material_loop: do id=first,last
            if (.not.supports_cavity3d(id)) cycle material_loop
            is_interface=interface_case(id)
            npoint=mesh%n_cells
            if (is_interface) npoint=mesh%n_faces
            allocate(positions(3,npoint),field(npoint),weights(npoint))
            if (is_interface) then
                positions=mesh%face_centers
                weights=mesh%face_areas
                domain="interface"
            else
                positions=mesh%centers
                weights=mesh%volumes
                domain="volume"
            end if

            if (sensitivity_mode) then
                call initialize_model(id,material,1)
                call reference_parameters(id,nominal,n)
                if (n/=material%parameter_count().or.n<1.or.n>max_parameters) then
                    error stop "Cavity parameter metadata mismatch"
                end if
                parameter_names=""
                parameter_units=cases(id)%parameter_units
                do j=1,n
                    call material%parameter_info(j,pname,found)
                    if (.not.found.or.len_trim(pname)==0.or.len_trim(parameter_units(j))==0) then
                        error stop "Missing cavity analytical parameter metadata"
                    end if
                    parameter_names(j)=pname
                end do
                call initialize_model(id,material,npoint,selected_parameters=parameter_names(:n))
                if (material%sensitivity_size()/=n) error stop "Incomplete cavity derivative selection"
            else
                call initialize_model(id,material,npoint)
            end if

            ninput=material%input_size()
            noutput=material%output_size()
            allocate(input(ninput,npoint),output(noutput,npoint))
            if (sensitivity_mode) then
                allocate(derivatives(noutput,n,npoint),derivative_field(npoint))
                allocate(scaled_fields(npoint,n),metrics(n))
            end if

            nsteps=1
            if (history_case(id)) nsteps=size(path)
            do step=1,nsteps
                amplitude=1.0_wp
                if (history_case(id)) amplitude=path(step)
                call fill_cavity_inputs(id,positions,amplitude,input)
                if (sensitivity_mode) then
                    call material%eval(input,output,status=code,dt=0.1_wp,dsdp=derivatives)
                else
                    call material%eval(input,output,status=code,dt=0.1_wp)
                end if
                if (code/=cf_ok) then
                    print "(a,1x,a,2(1x,i0))", "3D cavity material evaluation failed:",trim(names(id)),step,code
                    error stop 1
                end if
                call material%commit(code)
                if (code/=cf_ok) error stop "3D cavity history commit failed"
            end do
            if (.not.all(ieee_is_finite(output))) error stop "Nonfinite 3D cavity response"
            call project_cavity_response(id,positions,output,field)

            if (response_mode) then
                call write_cavity_response(mesh,id,is_interface,field, &
                    trim(cases(id)%response_label),trim(cases(id)%response_units),amplitude)
                if (.not.filtered) then
                    scheme=colormap_for(id)
                    call colormap_limits(scheme,field,lo,hi)
                    write(report_response, &
                        '(a,",",a,",",a,",",a,",",i0,",",es22.14e3,",",es22.14e3,",",a,",",es22.14e3,",",es22.14e3)') &
                        trim(names(id)),trim(domain),trim(cases(id)%response_label), &
                        trim(cases(id)%response_units),npoint,minval(field),maxval(field), &
                        trim(colormap_name(scheme)),lo,hi
                end if
            end if

            if (sensitivity_mode) then
                if (.not.all(ieee_is_finite(derivatives))) error stop "Nonfinite 3D cavity derivatives"
                do j=1,n
                    call project_cavity_response(id,positions,derivatives(:,j,:),derivative_field)
                    scale=parameter_scale(j,nominal(:n),parameter_units(:n))
                    if (scale<=0.0_wp) error stop "No positive reference scale for cavity parameter"
                    metrics(j)%parameter=parameter_names(j)
                    metrics(j)%units=parameter_units(j)
                    metrics(j)%reference_value=nominal(j)
                    metrics(j)%scale_source="abs_value"
                    if (abs(nominal(j))<=tiny(1.0_wp)) metrics(j)%scale_source="same_unit"
                    call calculate_metric(metrics(j),field,derivative_field,weights,scale)
                    scaled_fields(:,j)=scale*derivative_field

                    if (.not.filtered) then
                        status_label="ok"
                        log_label="NA"
                        if (.not.metrics(j)%relative_defined) then
                            status_label="zero_response"
                        else if (metrics(j)%relative_rms<=0.0_wp) then
                            status_label="zero_sensitivity"
                        else
                            if (metrics(j)%below_resolution) status_label="roundoff_level"
                            write(log_label,"(es21.12e3)") log10(metrics(j)%relative_rms)
                        end if
                        write(report_sens, &
                            '(a,",",i0,",",a,",",a,",",es22.14e3,",",es22.14e3,",",a,",",a,",",a,",",a,",",i0,",",i0,",",es22.14e3,",",es22.14e3,",",es22.14e3,",",a,",",es22.14e3,",",es22.14e3,",",es22.14e3,",",a)') &
                            trim(names(id)),j,trim(parameter_names(j)),trim(parameter_units(j)), &
                            nominal(j),scale,trim(metrics(j)%scale_source), &
                            trim(cases(id)%response_label),trim(cases(id)%response_units),trim(domain), &
                            nsteps,npoint,metrics(j)%response_rms,metrics(j)%sensitivity_rms, &
                            metrics(j)%relative_rms,trim(adjustl(log_label)), &
                            metrics(j)%sensitivity_linf,metrics(j)%signed_mean, &
                            metrics(j)%sign_coherence,trim(status_label)
                    end if
                end do
                call write_cavity_sensitivity(mesh,id,is_interface,parameter_names(:n),metrics,scaled_fields)
                deallocate(derivatives,derivative_field,scaled_fields,metrics)
            end if

            completed=completed+1
            print "(a,1x,a,1x,i0)", "3D cavity "//trim(domain),trim(names(id)),npoint
            deallocate(positions,field,weights,input,output)
        end do material_loop
        if (response_mode.and..not.filtered) close(report_response)
        if (sensitivity_mode.and..not.filtered) then
            close(report_sens)
            call write_benchmark_index("cavity3d","Spherical-cavity",[(supports_cavity3d(id),id=1,ncases)])
        end if
        print "(a,i0,a)", "Three-dimensional cavity benchmark: ",completed," material laws evaluated."
    end subroutine run_cavity3d

end module benchmark_cavity3d_runner
