! SPDX-License-Identifier: Apache-2.0

program test_axisymmetric_j2
    !! Test convention: stress and stress-like material parameters in Pa.
    !! Strain and deformation gradients are dimensionless. Explicit rescaling
    !! cases state their factors.
    !!
    !! Axisymmetric J2 full-history stress, tangent and sensitivity checks.
    !! Reference: the independent three-dimensional J2 material handle,
    !! evaluated with explicit suppressed r-theta/z-theta shear components.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
    use constituf, only: wp, material_t, verify_tangent, cf_ok, cf_invalid_argument, cf_evaluation_failed, cf_no_trial
    use constituf_axisymmetric_j2, only: axisymmetric_j2_t
    use constituf_j2_plasticity, only: j2_plasticity_t
    implicit none (type, external)

    integer, parameter :: naxis = 4, nmandel = 6, nstate = 7, npoints = 65
    integer, parameter :: nparameters = 4, nsteps = 6
    integer, parameter :: mapping(naxis) = [1,2,3,6]
    integer, parameter :: selected_ids(nparameters) = [2,4,1,3]
    real(wp), parameter :: base(nparameters) = [150000.0_wp,75000.0_wp,350.0_wp,2000.0_wp]
    real(wp), parameter :: factors(nsteps) = [0.6_wp,1.0_wp,1.5_wp,0.3_wp,-0.4_wp,1.7_wp]
    real(wp), parameter :: comparison_tolerance = 512.0_wp*epsilon(1.0_wp)
    character(len=18), parameter :: names(nparameters) = [character(len=18) :: &
        "shear_modulus", "hardening_modulus", "bulk_modulus", "yield_stress"]
    type(material_t) :: reduced, reference, invalid
    type(material_t) :: plus(nparameters,2), minus(nparameters,2)
    real(wp) :: input(naxis,npoints), output(naxis,npoints), tangent(naxis,naxis,npoints)
    real(wp) :: reference_input(nmandel,npoints), reference_output(nmandel,npoints)
    real(wp) :: reference_tangent(nmandel,nmandel,npoints)
    real(wp) :: dsdp(naxis,nparameters,npoints), reference_dsdp(nmandel,nparameters,npoints)
    real(wp) :: state_reduced(nstate,npoints), state_reference(nstate,npoints), before(nstate,npoints)
    real(wp) :: state_gradient_reduced(nstate,nparameters,npoints), state_gradient_reference(nstate,nparameters,npoints)
    real(wp) :: before_gradient(nstate,nparameters,npoints)
    real(wp) :: positive_parameters(nparameters), negative_parameters(nparameters)
    real(wp) :: perturbed_plus(naxis,npoints), perturbed_minus(naxis,npoints)
    real(wp) :: parameter_steps(nparameters,2), numerical(naxis,npoints,2)
    real(wp) :: derivative_tolerance(naxis,npoints)
    real(wp) :: step_size, error, reference_component, roundoff_tolerance
    integer :: status, step, k, i, j, point, refinement

    call reduced%init(axisymmetric_j2_t(base(1),base(2),base(3),base(4)),npoints,status, &
        sensitivity_parameters=names)
    if (status /= cf_ok) error stop "Axisymmetric J2 initialization"
    if (.not. reduced%field_matches("input","strain.small","custom.axisymmetric_mandel", &
        "small_strain","spatial",[3,3])) then
        error stop "Axisymmetric J2 has incorrect physical field convention"
    end if
    if (reduced%input_size() /= naxis .or. reduced%output_size() /= naxis .or. reduced%state_size() /= nstate) then
        error stop "Axisymmetric J2 declares incorrect packed field dimensions"
    end if

    call reference%init(j2_plasticity_t(base(1),base(2),base(3),base(4)),npoints,status, &
        sensitivity_parameters=names)
    if (status /= cf_ok) error stop "3D J2 reference initialization"

    do k = 1,nparameters
        parameter_steps(k,1) = epsilon(1.0_wp)**(1.0_wp/3.0_wp)*base(selected_ids(k))
        parameter_steps(k,2) = 0.5_wp*parameter_steps(k,1)
        do refinement = 1,2
            positive_parameters = base
            negative_parameters = base
            positive_parameters(selected_ids(k)) = positive_parameters(selected_ids(k))+parameter_steps(k,refinement)
            negative_parameters(selected_ids(k)) = negative_parameters(selected_ids(k))-parameter_steps(k,refinement)
            call plus(k,refinement)%init(axisymmetric_j2_t(positive_parameters(1),positive_parameters(2), &
                positive_parameters(3),positive_parameters(4)),npoints,status)
            if (status /= cf_ok) error stop "Positive axisymmetric parameter perturbation"
            call minus(k,refinement)%init(axisymmetric_j2_t(negative_parameters(1),negative_parameters(2), &
                negative_parameters(3),negative_parameters(4)),npoints,status)
            if (status /= cf_ok) error stop "Negative axisymmetric parameter perturbation"
        end do
    end do

    step_loop: do step = 1,nsteps
        input(:,1) = factors(step)*[0.006_wp,-0.003_wp,-0.003_wp,0.0010_wp]
        input(:,2) = factors(step)*[-0.004_wp,0.002_wp,0.002_wp,-0.0009_wp]
        do point = 3,npoints
            input(:,point) = input(:,1+mod(point-1,2))
        end do
        ! Point 65 has distinct history on the far side of the 64-point chunk boundary.
        input(:,npoints) = 1.25_wp*input(:,npoints)
        reference_input = 0.0_wp
        reference_input(mapping,:) = input

        call reduced%eval(input,output,tangent,status,dsdp=dsdp)
        if (status /= cf_ok) error stop "Axisymmetric J2 evaluated trial"
        call reference%eval(reference_input,reference_output,reference_tangent,status,dsdp=reference_dsdp)
        if (status /= cf_ok) error stop "3D J2 evaluated reference trial"
        if (.not. all(ieee_is_finite(output)) .or. .not. all(ieee_is_finite(reference_output)) .or. &
            .not. all(ieee_is_finite(tangent)) .or. .not. all(ieee_is_finite(reference_tangent))) then
            error stop "Nonfinite axisymmetric J2 response comparison"
        end if
        if (maxval(abs(output-reference_output(mapping,:))) > &
            comparison_tolerance*max(1.0_wp,maxval(abs(reference_output)))) then
            error stop "Axisymmetric J2 stress differs from 3D constitutive law"
        end if
        do j = 1,naxis
            do i = 1,naxis
                reference_component = maxval(abs(tangent(i,j,:)-reference_tangent(mapping(i),mapping(j),:)))
                if (reference_component > comparison_tolerance* &
                    max(1.0_wp,maxval(abs(reference_tangent(mapping(i),mapping(j),:))))) then
                    error stop "Axisymmetric J2 algorithmic tangent differs from 3D restriction"
                end if
            end do
        end do
        if (.not. all(ieee_is_finite(dsdp)) .or. .not. all(ieee_is_finite(reference_dsdp))) then
            error stop "Nonfinite axisymmetric J2 sensitivity comparison"
        end if
        do k = 1,nparameters
            if (maxval(abs(dsdp(:,k,:)-reference_dsdp(mapping,k,:))) > &
                comparison_tolerance*max(1.0_wp,maxval(abs(reference_dsdp(mapping,k,:))))) then
                error stop "Axisymmetric J2 analytical sensitivity does not equal full J2"
            end if
            ! Preserve both perturbed accepted histories for a central-difference truncation estimate.
            do refinement = 1,2
                call plus(k,refinement)%eval(input,perturbed_plus,status=status)
                if (status /= cf_ok) error stop "Positive full-history axisymmetric evaluation"
                call minus(k,refinement)%eval(input,perturbed_minus,status=status)
                if (status /= cf_ok) error stop "Negative full-history axisymmetric evaluation"
                numerical(:,:,refinement) = (perturbed_plus-perturbed_minus)/(2.0_wp*parameter_steps(k,refinement))
            end do
            derivative_tolerance = abs(numerical(:,:,2)-numerical(:,:,1))/3.0_wp + &
                64.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(perturbed_plus),abs(perturbed_minus))/parameter_steps(k,2)
            if (.not. all(ieee_is_finite(numerical))) error stop "Nonfinite axisymmetric parameter secant"
            if (any(abs(numerical(:,:,2)-dsdp(:,k,:)) > derivative_tolerance)) then
                error stop "Axisymmetric J2 sensitivity differs from independently perturbed full histories"
            end if
        end do

        call reduced%commit(status)
        if (status /= cf_ok) error stop "Axisymmetric J2 accept"
        call reference%commit(status)
        if (status /= cf_ok) error stop "3D J2 accept"
        do k = 1,nparameters
            do refinement = 1,2
                call plus(k,refinement)%commit(status)
                if (status /= cf_ok) error stop "Perturbed axisymmetric J2 accept"
                call minus(k,refinement)%commit(status)
                if (status /= cf_ok) error stop "Perturbed axisymmetric J2 accept"
            end do
        end do
        call reduced%export_committed(state_reduced,status,sensitivity_state=state_gradient_reduced)
        if (status /= cf_ok) error stop "Axisymmetric J2 accepted history export"
        call reference%export_committed(state_reference,status,sensitivity_state=state_gradient_reference)
        if (status /= cf_ok) error stop "3D J2 accepted history export"
        if (.not. all(ieee_is_finite(state_reduced)) .or. .not. all(ieee_is_finite(state_reference)) .or. &
            .not. all(ieee_is_finite(state_gradient_reduced)) .or. &
            .not. all(ieee_is_finite(state_gradient_reference))) then
            error stop "Nonfinite axisymmetric J2 history comparison"
        end if
        if (maxval(abs(state_reduced-state_reference)) > comparison_tolerance) then
            error stop "Axisymmetric J2 accepted history differs from original 3D J2"
        end if
        if (maxval(abs(state_gradient_reduced-state_gradient_reference)) > comparison_tolerance) then
            error stop "Axisymmetric J2 accepted derivative history differs from original 3D J2"
        end if
    end do step_loop

    input(:,1) = [0.014_wp,-0.006_wp,-0.008_wp,0.0023_wp]
    input(:,2) = [-0.010_wp,0.005_wp,0.005_wp,-0.0018_wp]
    do point = 3,npoints
        input(:,point) = input(:,1+mod(point-1,2))
    end do
    input(:,npoints) = 1.25_wp*input(:,npoints)
    call reduced%eval(input,output,status=status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Axisymmetric tangent reference"
    if (.not. all(ieee_is_finite(output))) error stop "Nonfinite axisymmetric tangent reference"
    step_size = maxval(abs(input))*epsilon(1.0_wp)**(1.0_wp/3.0_wp)
    roundoff_tolerance = 64.0_wp*epsilon(1.0_wp)*max(1.0_wp,maxval(abs(output)))/step_size
    call verify_tangent(reduced,input,error,status,step=[(step_size,i=1,naxis)], &
        absolute_tolerance=roundoff_tolerance)
    if (status /= cf_ok .or. .not. ieee_is_finite(error) .or. error > 1.0_wp) then
        error stop "History-bearing axisymmetric J2 tangent verification"
    end if
    call reduced%export_committed(before,status,sensitivity_state=before_gradient)
    if (status /= cf_ok .or. .not. all(abs(before-state_reduced) <= 0.0_wp) .or. &
        .not. all(abs(before_gradient-state_gradient_reduced) <= 0.0_wp)) then
        error stop "Tangent check overwrote axisymmetric J2 committed history"
    end if

    call reduced%eval(input,output,status=status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Axisymmetric J2 rejected trial setup"
    call reduced%revert()
    call reduced%export_committed(before,status,sensitivity_state=before_gradient)
    if (status /= cf_ok .or. .not. all(abs(before-state_reduced) <= 0.0_wp) .or. &
        .not. all(abs(before_gradient-state_gradient_reduced) <= 0.0_wp)) then
        error stop "Axisymmetric J2 revert changed accepted history"
    end if

    input(2,npoints) = ieee_value(0.0_wp,ieee_quiet_nan)
    call reduced%eval(input,output,status=status,dsdp=dsdp)
    if (status /= cf_evaluation_failed) error stop "Nonfinite axisymmetric J2 input accepted"
    if (reduced%last_failure_point() /= npoints) error stop "Axisymmetric J2 chunk failure offset"
    call reduced%commit(status)
    if (status /= cf_no_trial) error stop "Invalid axisymmetric J2 trial became committed"

    call invalid%init(axisymmetric_j2_t(base(1),base(2),0.0_wp,base(4)),npoints,status)
    if (status /= cf_invalid_argument) error stop "Axisymmetric J2 accepted zero yield stress"
    call invalid%init(axisymmetric_j2_t(base(1),base(2),base(3),-base(4)),npoints,status)
    if (status /= cf_invalid_argument) error stop "Axisymmetric J2 accepted negative hardening"
    call test_physical_hoop_and_shear()
    call test_exact_shear()

contains

    subroutine test_physical_hoop_and_shear()
        type(material_t) :: material
        real(wp) :: driving(naxis,2), stress(naxis,2), jacobian(naxis,naxis,2)
        real(wp) :: gradient(naxis,2,2), history(nstate,2)
        real(wp) :: history_gradient(nstate,2,2), expected_history_gradient(nstate,2,2)
        real(wp) :: expected_stress(naxis,2), expected_gradient(naxis,2,2), expected_history(nstate,2)
        real(wp) :: increment
        integer :: result, ipoint, component, parameter

        call material%init(axisymmetric_j2_t(100.0_wp,40.0_wp,5.0_wp,0.0_wp),2,result, &
            sensitivity_parameters=[character(len=12) :: "bulk_modulus","yield_stress"])
        if (result /= cf_ok) error stop "Physical axisymmetric reference setup"
        driving = 0.0_wp
        driving(3,1) = 0.01_wp
        ! Engineering rz shear 0.2 is packed as gamma/sqrt(2).
        driving(4,2) = 0.2_wp/sqrt(2.0_wp)
        call material%eval(driving,stress,jacobian,result,dsdp=gradient)
        if (result /= cf_ok) error stop "Physical axisymmetric reference evaluation"
        expected_stress = 0.0_wp
        expected_stress(:,1) = [11.0_wp/15.0_wp,11.0_wp/15.0_wp,23.0_wp/15.0_wp,0.0_wp]
        expected_stress(4,2) = 5.0_wp*sqrt(2.0_wp/3.0_wp)
        expected_gradient = 0.0_wp
        expected_gradient(1:3,1,1) = 0.01_wp
        expected_gradient(4,2,2) = sqrt(2.0_wp/3.0_wp)
        do ipoint = 1,2
            do component = 1,naxis
                call check(stress(component,ipoint),expected_stress(component,ipoint))
                do parameter = 1,2
                    call check(gradient(component,parameter,ipoint),expected_gradient(component,parameter,ipoint))
                end do
            end do
        end do
        call check(jacobian(4,4,2),0.0_wp)
        call material%commit(result)
        if (result /= cf_ok) error stop "Physical axisymmetric reference acceptance"
        call material%export_committed(history,result,sensitivity_state=history_gradient)
        if (result /= cf_ok) error stop "Physical axisymmetric history export"
        increment = (8.0_wp*sqrt(3.0_wp)-5.0_wp)/120.0_wp
        expected_history = 0.0_wp
        expected_history(6,2) = increment*sqrt(1.5_wp)
        expected_history(7,2) = increment
        expected_history_gradient = 0.0_wp
        expected_history_gradient(6,2,2) = -sqrt(1.5_wp)/120.0_wp
        expected_history_gradient(7,2,2) = -1.0_wp/120.0_wp
        do ipoint = 1,2
            do component = 1,nstate
                call check(history(component,ipoint),expected_history(component,ipoint))
                do parameter = 1,2
                    call check(history_gradient(component,parameter,ipoint), &
                        expected_history_gradient(component,parameter,ipoint))
                end do
            end do
        end do
    end subroutine test_physical_hoop_and_shear

    subroutine test_exact_shear()
        type(material_t) :: material
        integer, parameter :: count = 2
        real(wp) :: strain(naxis,count), stress(naxis,count), jacobian(naxis,naxis,count)
        real(wp) :: derivative(naxis,1,count), history(nstate,count), history_derivative(nstate,1,count)
        real(wp) :: direction
        integer :: result, point, component

        ! q_trial=12, delta_gamma=3/65 and q_new=48/13 follow directly from the radial return.
        ! Reverse the second point to check the shear sign and hardening response.
        call material%init(axisymmetric_j2_t(100.0_wp,60.0_wp,3.0_wp,15.0_wp),count,result, &
            sensitivity_parameters=[character(len=12) :: "yield_stress"])
        if (result /= cf_ok) error stop "Exact axisymmetric shear setup"
        strain = 0.0_wp
        strain(4,:) = sqrt(2.0_wp/3.0_wp)/10.0_wp
        strain(4,count) = -strain(4,count)
        call material%eval(strain,stress,jacobian,result,dsdp=derivative)
        if (result /= cf_ok) error stop "Exact axisymmetric shear evaluation"
        call material%commit(result)
        if (result /= cf_ok) error stop "Exact axisymmetric shear acceptance"
        call material%export_committed(history,result,sensitivity_state=history_derivative)
        if (result /= cf_ok) error stop "Exact axisymmetric shear history"
        do point = 1,count
            direction = sign(1.0_wp,strain(4,point))
            do component = 1,3
                call check(stress(component,point),0.0_wp)
                call check(derivative(component,1,point),0.0_wp)
            end do
            call check(stress(4,point),direction*(48.0_wp/13.0_wp)*sqrt(2.0_wp/3.0_wp))
            call check(jacobian(4,4,point),120.0_wp/13.0_wp)
            call check(derivative(4,1,point),direction*(12.0_wp/13.0_wp)*sqrt(2.0_wp/3.0_wp))
            do component = 1,5
                call check(history(component,point),0.0_wp)
                call check(history_derivative(component,1,point),0.0_wp)
            end do
            call check(history(6,point),direction*(3.0_wp/65.0_wp)*sqrt(3.0_wp/2.0_wp))
            call check(history(7,point),3.0_wp/65.0_wp)
            call check(history_derivative(6,1,point),-direction*sqrt(3.0_wp/2.0_wp)/195.0_wp)
            call check(history_derivative(7,1,point),-1.0_wp/195.0_wp)
        end do
    end subroutine test_exact_shear

    subroutine check(actual, expected)
        real(wp), intent(in) :: actual, expected

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
            error stop "Nonfinite axisymmetric physical reference"
        end if
        if (abs(actual-expected) > comparison_tolerance*max(1.0_wp,abs(expected))) then
            print *, "Expected",expected,"got",actual
            error stop "Axisymmetric physical reference mismatch"
        end if
    end subroutine check

end program test_axisymmetric_j2
