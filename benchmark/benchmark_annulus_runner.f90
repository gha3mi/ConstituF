! SPDX-License-Identifier: Apache-2.0

module benchmark_annulus_runner
    !! Material-point response and analytical sensitivity on the annulus.
    !! Inputs are prescribed; derivatives are computed only by ConstituF.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use benchmark_annulus_export, only: write_annulus_response, write_annulus_sensitivity
    use benchmark_annulus_fields, only: fill_annulus_inputs, project_annulus_response
    use benchmark_annulus_mesh, only: nradial, nangle, cell_radius
    use benchmark_colormaps, only: colormap_for, colormap_name, colormap_limits
    use benchmark_history, only: npath, path
    use benchmark_metrics, only: sensitivity_metric_t, calculate_metric, parameter_scale, verify_metric_contract
    use benchmark_models, only: ncases, names, cases, initialize_model, response_kind, interface_case, &
                               history_case, max_parameters, reference_parameters
    use benchmark_pages, only: write_benchmark_index
    use benchmark_paths, only: annulus_directory
    use constituf, only: wp, material_t, cf_ok
    implicit none (type, external)
    private

    public :: run_annulus

contains

    subroutine run_annulus(mode,model_filter)
        !! Modes: all, response, sensitivity. An optional exact material slug
        !! selects only that model, without overwriting aggregate CSV reports.
        character(*), intent(in) :: mode
        character(*), intent(in), optional :: model_filter
        logical :: response_mode,sensitivity_mode

        response_mode=mode=="all".or.mode=="response"
        sensitivity_mode=mode=="all".or.mode=="sensitivity"
        if (.not.response_mode.and..not.sensitivity_mode) then
            error stop "Expected annulus mode: all, response, sensitivity"
        end if

        if (response_mode) call run_response(model_filter)
        if (sensitivity_mode) call run_sensitivity(model_filter)
    end subroutine run_annulus

    subroutine run_response(model_filter)
        !! An optional exact material slug selects only that model without
        !! overwriting the aggregate response manifest.
        character(*), intent(in), optional :: model_filter
        type(material_t) :: material
        real(wp), allocatable :: x(:,:), y(:,:), scalar(:)
        integer :: id, np, ni, no, step, steps, status, manifest, kind, scheme
        integer :: selected_id, first, last
        real(wp) :: scale_min, scale_max
        character(:), allocatable :: label, units
        character(12) :: domain
        logical :: single_model

        selected_id=0
        single_model=present(model_filter)
        if (single_model) then
            do id=1,ncases
                if (trim(model_filter)==trim(names(id))) selected_id=id
            end do
            if (selected_id==0) error stop "Unknown response material: use a declared model slug"
        end if
        first=1
        last=ncases
        if (single_model) then
            first=selected_id
            last=selected_id
        else
            open(newunit=manifest,file=annulus_directory//"response_manifest.csv", &
                 status="replace",action="write")
            write(manifest,"(a)") "model,domain,quantity,units,npoints,minimum,maximum,colormap,scale_min,scale_max"
        end if

        do id=first,last
            np=nradial*nangle
            if (interface_case(id)) np=nangle
            call initialize_model(id,material,np)
            ni=material%input_size()
            no=material%output_size()
            if (ni<1.or.no<1) error stop "Invalid constitutive field size"
            allocate(x(ni,np),y(no,np),scalar(np))

            steps=1
            if (history_case(id)) steps=size(path)
            do step=1,steps
                if (steps==1) then
                    call fill_annulus_inputs(id,x,1.0_wp)
                else
                    call fill_annulus_inputs(id,x,path(step))
                end if
                call material%eval(x,y,status=status,dt=0.1_wp)
                if (status/=cf_ok) then
                    print "(a,1x,a,1x,i0,1x,i0)", &
                        "Annulus local evaluation failed:",trim(names(id)),step,status
                    error stop 1
                end if
                call material%commit(status)
                if (status/=cf_ok) error stop "Annulus history commit failed"
            end do

            if (.not.all(ieee_is_finite(y))) error stop "Nonfinite constitutive response"
            call response_kind(id,kind,label,units)
            call project_annulus_response(id,kind,y,scalar)
            call write_annulus_response(id,scalar,trim(label),trim(units))
            domain="bulk"
            if (interface_case(id)) domain="interface"
            if (.not.single_model) then
                scheme=colormap_for(id)
                call colormap_limits(scheme,scalar,scale_min,scale_max)
                write(manifest,'(a,",",a,",",a,",",a,",",i0,",",es18.10e3,",",es18.10e3,",",a,",",es18.10e3,",",es18.10e3)') &
                    trim(names(id)),trim(domain),trim(label),trim(units),np, &
                    minval(scalar),maxval(scalar),trim(colormap_name(scheme)),scale_min,scale_max
            end if
            print "(a,1x,a,1x,i0)", "response "//trim(domain),trim(names(id)),np
            deallocate(x,y,scalar)
        end do
        if (single_model) then
            print "(a,1x,a)", "Standalone response material benchmark written:",trim(names(selected_id))
            return
        end if
        close(manifest)
        print "(a,i0,a)", "Constitutive response benchmark: ",ncases," materials; SVG/CSV/VTK written."
    end subroutine run_response

    subroutine run_sensitivity(model_filter)
        !! Optional exact model slug runs one complete, independent material
        !! benchmark without overwriting the full-atlas metrics or manifest.
        character(*), intent(in), optional :: model_filter
        !! Sensitivities differentiate *complete committed material history*
        !! at parameter-independent imposed fields and time increments.
        !! The scalar heatmap uses cell-area weighted response and derivative
        !! RMS norms, never elementwise division by a possibly zero response.
        type(material_t) :: material
        type(sensitivity_metric_t), allocatable :: metrics(:)
        real(wp), allocatable :: x(:,:),y(:,:)
        real(wp), allocatable :: grad(:,:,:)
        real(wp), allocatable :: field(:),dfield(:),weights(:)
        real(wp), allocatable :: per_parameter_fields(:,:)
        real(wp) :: values(max_parameters)
        real(wp) :: scale
        character(48) :: param_names(max_parameters)
        character(24) :: param_units(max_parameters)
        character(:), allocatable :: response_name,response_unit,parameter_name
        character(20) :: domain
        character(32) :: log_label,status_label
        integer :: id,np,nin,nout,status,istep,steps,kind,report
        integer :: n,j,k,derivative_columns,selected_id,first,last
        logical :: found,single_model

        selected_id=0
        single_model=present(model_filter)
        if (single_model) then
            do id=1,ncases
                if (trim(model_filter)==trim(names(id))) selected_id=id
            end do
            if (selected_id==0) error stop "Unknown sensitivity material: use a declared model slug"
        end if
        first=1
        last=ncases
        if (single_model) then
            first=selected_id
            last=selected_id
        else
            open(newunit=report,file=annulus_directory//"sensitivity_metrics.csv",status="replace",action="write")
            write(report,"(a)") "model,parameter_index,parameter,parameter_units,parameter_value,reference_scale,scale_basis,response,response_units,domain,steps,npoints,response_rms,scaled_sensitivity_rms,relative_rms,log10_relative_rms,sensitivity_linf,signed_mean,sign_coherence,status"
        end if

        call verify_metric_contract()
        derivative_columns=0
        do id=first,last
            ! First initialize a cheap one-point handle to interrogate the
            ! model's declared parameter names. The benchmark registry owns
            ! the SI labels for its chosen reference values.
            call initialize_model(id,material,1)
            call reference_parameters(id,values,n)
            if (n/=material%parameter_count()) then
                print "(a,1x,a,2(1x,i0))","Declared/reference parameter count mismatch:", &
                    trim(names(id)),material%parameter_count(),n
                error stop 1
            end if
            if (n<1.or.n>max_parameters) error stop "Unsupported parameter dimension"
            allocate(metrics(n))
            param_names=""
            param_units=cases(id)%parameter_units
            do j=1,n
                call material%parameter_info(j,parameter_name,found)
                if (.not.found.or.len_trim(parameter_name)==0.or.len_trim(param_units(j))==0) then
                    error stop "Missing analytical parameter metadata"
                end if
                param_names(j)=parameter_name
            end do

            np=nradial*nangle
            if (interface_case(id)) np=nangle
            call initialize_model(id,material,np,selected_parameters=param_names(:n))
            if (material%sensitivity_size()/=n) error stop "Incomplete parameter selection"
            nin=material%input_size()
            nout=material%output_size()
            allocate(x(nin,np),y(nout,np),grad(nout,n,np),field(np),dfield(np),weights(np))
            allocate(per_parameter_fields(np,n))

            do k=1,np
                weights(k)=1.0_wp
                if (.not.interface_case(id)) weights(k)=cell_radius((k-1)/nangle+1)
            end do
            steps=1
            if (history_case(id)) steps=npath
            do istep=1,steps
                if (steps==1) then
                    call fill_annulus_inputs(id,x,1.0_wp)
                else
                    call fill_annulus_inputs(id,x,path(istep))
                end if
                call material%eval(x,y,status=status,dt=0.1_wp,dsdp=grad)
                if (status/=cf_ok) then
                    print "(a,1x,a,2(1x,i0))", "All-parameter analytical evaluation:",trim(names(id)),istep,status
                    error stop 1
                end if
                call material%commit(status)
                if (status/=cf_ok) error stop "All-parameter derivative history commit"
            end do

            if (.not.all(ieee_is_finite(y)).or..not.all(ieee_is_finite(grad))) then
                error stop "Nonfinite constitutive response or analytical derivative"
            end if
            call response_kind(id,kind,response_name,response_unit)
            call project_annulus_response(id,kind,y,field)
            domain="bulk"
            if (interface_case(id)) domain="interface"

            do j=1,n
                call project_annulus_response(id,kind,grad(:,j,:),dfield)
                scale=parameter_scale(j,values(:n),param_units(:n))
                if (scale<=0.0_wp) error stop "No physical scaling reference for zero parameter"
                metrics(j)%parameter=param_names(j)
                metrics(j)%units=param_units(j)
                metrics(j)%reference_value=values(j)
                metrics(j)%scale_source="abs_value"
                if (abs(values(j))<=tiny(1.0_wp)) metrics(j)%scale_source="same_unit"
                call calculate_metric(metrics(j),field,dfield,weights,scale)
                per_parameter_fields(:,j)=scale*dfield
                derivative_columns=derivative_columns+1
                log_label="NA"
                status_label="ok"
                if (.not.metrics(j)%relative_defined) then
                    status_label="zero_response"
                else if (metrics(j)%relative_rms<=0.0_wp) then
                    status_label="zero_sensitivity"
                else
                    if (metrics(j)%below_resolution) status_label="roundoff_level"
                    write(log_label,"(es21.12e3)") log10(metrics(j)%relative_rms)
                end if
                if (.not.single_model) then
                    write(report,'(a,",",i0,",",a,",",a,",",es22.14e3,",",es22.14e3,",",a,",",a,",",a,",",a,",",i0,",",i0,",",es22.14e3,",",es22.14e3,",",es22.14e3,",",a,",",es22.14e3,",",es22.14e3,",",es22.14e3,",",a)') &
                        trim(names(id)),j,trim(param_names(j)),trim(param_units(j)),values(j),scale, &
                        trim(metrics(j)%scale_source),trim(response_name),trim(response_unit),trim(domain), &
                        steps,np,metrics(j)%response_rms,metrics(j)%sensitivity_rms, &
                        metrics(j)%relative_rms,trim(adjustl(log_label)),metrics(j)%sensitivity_linf, &
                        metrics(j)%signed_mean,metrics(j)%sign_coherence,trim(status_label)
                end if

            end do
            call write_annulus_sensitivity(id,param_names(:n),param_units(:n),values(:n),metrics, &
                                          field,per_parameter_fields,response_name,response_unit,trim(domain),steps)
            print "(a,1x,a,1x,i0)", "analytical parameters",trim(names(id)),n
            deallocate(x,y,grad,field,dfield,weights,per_parameter_fields,metrics)
        end do
        if (single_model) then
            print "(a,1x,a)", "Standalone analytical material benchmark written:",trim(names(selected_id))
            return
        end if
        close(report)
        call write_benchmark_index("annulus","Annulus")

        print "(a,i0,a)", "Analytical parameter sensitivity benchmark: ",derivative_columns," derivative fields evaluated; no FDM."
    end subroutine run_sensitivity

end module benchmark_annulus_runner
