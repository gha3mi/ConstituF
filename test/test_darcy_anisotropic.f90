! SPDX-License-Identifier: Apache-2.0

program test_darcy_anisotropic
    !! Test convention: pressure in Pa, length in m and time in s.
    !! Permeability is m2, viscosity Pa*s and Darcy velocity m/s.
    !!
    !! Check Darcy tensor reference, dissipation, frame covariance and derivatives.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
    use constituf, only: wp, material_t, cf_ok, cf_invalid_argument, &
                        cf_invalid_shape, cf_evaluation_failed, cf_no_trial, cf_unsupported
    use constituf_darcy_anisotropic, only: darcy_anisotropic_t
    use constituf_darcy_isotropic, only: darcy_isotropic_t
    implicit none (type, external)

    integer, parameter :: ndim = 3, nparameters = 7
    character(len=18), parameter :: parameter_names(nparameters) = [character(len=18) :: &
        "dynamic_viscosity","permeability_23","permeability_13","permeability_12", &
        "permeability_33","permeability_22","permeability_11"]
    real(wp), parameter :: reference_sensitivities(ndim,nparameters) = reshape( &
        [11.0_wp/8.0_wp,-53.0_wp/16.0_wp,-0.25_wp, &
         0.0_wp,-0.25_wp,1.0_wp, -0.25_wp,0.0_wp,-0.5_wp, &
         1.0_wp,-0.5_wp,0.0_wp, 0.0_wp,0.0_wp,-0.25_wp, &
         0.0_wp,1.0_wp,0.0_wp, -0.5_wp,0.0_wp,0.0_wp],[ndim,nparameters])

    call test_reference_and_schema()
    call test_tangent_batch()
    call test_rotation_and_isotropic_limit()
    call test_sensitivities()
    call test_common_parameter_scale(sqrt(tiny(1.0_wp)))
    call test_common_parameter_scale(sqrt(huge(1.0_wp)))
    call test_cancelling_large_products()
    call test_cancellation_with_small_remainder()
    call test_sensitivity_range()
    call test_divided_mobility_underflow()
    call test_subnormal_product_sum()
    call test_failures()

contains

    function model_from_values(values) result(model)
        real(wp), intent(in) :: values(nparameters)
        type(darcy_anisotropic_t) :: model

        model = darcy_anisotropic_t(values(1),values(2),values(3),values(4), &
                                   values(5),values(6),values(7))
    end function model_from_values

    subroutine expect(actual, wanted)
        integer, intent(in) :: actual, wanted

        if (actual /= wanted) then
            print *, "Status", actual, "expected", wanted
            error stop "Anisotropic Darcy status mismatch"
        end if
    end subroutine expect

    subroutine close(actual, wanted)
        real(wp), intent(in) :: actual, wanted

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(wanted)) then
            error stop "Nonfinite anisotropic Darcy comparison"
        end if
        if (abs(actual-wanted) > 256.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(wanted))) then
            print *, "Actual", actual, "expected", wanted
            error stop "Anisotropic Darcy reference mismatch"
        end if
    end subroutine close

    subroutine permeability_reference(tensor)
        real(wp), intent(out) :: tensor(ndim,ndim)

        tensor(:,1) = [10.0_wp,2.0_wp,-1.0_wp]
        tensor(:,2) = [2.0_wp,8.0_wp,1.5_wp]
        tensor(:,3) = [-1.0_wp,1.5_wp,6.0_wp]
    end subroutine permeability_reference

    subroutine test_reference_and_schema()
        type(material_t) :: material
        real(wp) :: tensor(ndim,ndim)
        real(wp) :: gradient(ndim,3), velocity(ndim,3), tangent(ndim,ndim,3)
        integer :: status, ipoint, i, j

        call material%init(darcy_anisotropic_t(10.0_wp,8.0_wp,6.0_wp, &
                           2.0_wp,-1.0_wp,1.5_wp,2.0_wp),3,status)
        call expect(status,cf_ok)
        if (material%input_size() /= ndim .or. material%output_size() /= ndim) then
            error stop "Anisotropic Darcy field sizes"
        end if
        if (material%state_storage_bytes() /= 0 .or. material%state_size() /= 0) then
            error stop "Anisotropic Darcy must be stateless"
        end if
        if (.not. material%field_matches("input", "pressure.gradient", &
            "vector", "pressure_gradient", "spatial", [ndim])) then
            error stop "Anisotropic Darcy input metadata"
        end if
        if (.not. material%field_matches("output", "fluid.darcy_velocity", &
            "vector", "specific_discharge", "spatial", [ndim])) then
            error stop "Anisotropic Darcy output metadata"
        end if
        if (material%field_matches("output", "fluid.darcy_velocity", &
            "vector", "specific_discharge", "reference", [ndim])) then
            error stop "Anisotropic Darcy wrong frame accepted"
        end if

        call permeability_reference(tensor)
        gradient(:,1) = [1.0_wp,-2.0_wp,0.5_wp]
        gradient(:,2) = [-3.0_wp,4.0_wp,1.25_wp]
        gradient(:,3) = 0.0_wp
        call material%eval(gradient,velocity,tangent,status)
        call expect(status,cf_ok)
        call close(velocity(1,1),-2.75_wp)
        call close(velocity(2,1),6.625_wp)
        call close(velocity(3,1),0.5_wp)
        call close(-dot_product(gradient(:,1),velocity(:,1)),63.0_wp/4.0_wp)
        do ipoint = 1, 3
            do i = 1, ndim
                call close(velocity(i,ipoint), &
                           -dot_product(tensor(i,:),gradient(:,ipoint))/2.0_wp)
                do j = 1, ndim
                    call close(tangent(i,j,ipoint),-tensor(i,j)/2.0_wp)
                end do
            end do
            if (-dot_product(velocity(:,ipoint),gradient(:,ipoint)) < 0.0_wp) then
                error stop "Anisotropic Darcy dissipation is negative"
            end if
        end do
        call material%commit(status)
        call expect(status,cf_ok)
        call material%eval(gradient,velocity,status=status)
        call expect(status,cf_ok)
        call material%revert()
        call material%commit(status)
        call expect(status,cf_no_trial)
    end subroutine test_reference_and_schema

    subroutine test_tangent_batch()
        type(material_t) :: batch, single
        real(wp) :: gradient(ndim,2), velocity(ndim,2), tangent(ndim,ndim,2)
        real(wp) :: changed(ndim,2), plus(ndim,2), minus(ndim,2), alone(ndim,1)
        real(wp) :: h, derivative, error_bound
        integer :: status, ipoint, i, j

        call batch%init(darcy_anisotropic_t(10.0_wp,8.0_wp,6.0_wp, &
                        2.0_wp,-1.0_wp,1.5_wp,2.0_wp),2,status)
        call expect(status,cf_ok)
        call single%init(darcy_anisotropic_t(10.0_wp,8.0_wp,6.0_wp, &
                         2.0_wp,-1.0_wp,1.5_wp,2.0_wp),1,status)
        call expect(status,cf_ok)
        gradient(:,1) = [1.3_wp,-0.9_wp,2.1_wp]
        gradient(:,2) = [-1.2_wp,0.75_wp,0.4_wp]
        call batch%eval(gradient,velocity,tangent,status)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(tangent))) error stop "Nonfinite anisotropic Darcy tangent"
        call single%eval(gradient(:,1:1),alone,status=status)
        call expect(status,cf_ok)
        do i = 1, ndim
            call close(velocity(i,1),alone(i,1))
        end do

        do ipoint = 1, 2
            do j = 1, ndim
                ! The gradient enters affinely, so this finite step has no
                ! central-difference truncation and avoids cancellation.
                h = 0.25_wp*maxval(abs(gradient(:,ipoint)))
                changed = gradient
                changed(j,ipoint) = gradient(j,ipoint)+h
                call batch%eval(changed,plus,status=status)
                call expect(status,cf_ok)
                changed(j,ipoint) = gradient(j,ipoint)-h
                call batch%eval(changed,minus,status=status)
                call expect(status,cf_ok)
                if (.not. all(ieee_is_finite(plus)) .or. .not. all(ieee_is_finite(minus))) then
                    error stop "Nonfinite anisotropic Darcy tangent secant response"
                end if
                error_bound = 32.0_wp*epsilon(1.0_wp)* &
                    maxval(abs(plus(:,ipoint))+abs(minus(:,ipoint)))/(2.0_wp*h)
                if (.not. ieee_is_finite(error_bound)) error stop "Nonfinite anisotropic Darcy tangent bound"
                do i = 1, ndim
                    derivative = (plus(i,ipoint)-minus(i,ipoint))/(2.0_wp*h)
                    if (.not. ieee_is_finite(derivative) .or. &
                        abs(derivative-tangent(i,j,ipoint)) > error_bound) then
                        error stop "Anisotropic Darcy tangent finite difference"
                    end if
                end do
            end do
        end do
    end subroutine test_tangent_batch

    subroutine test_rotation_and_isotropic_limit()
        type(material_t) :: material, rotated, isotropic, diagonal
        real(wp) :: tensor(ndim,ndim), rotated_tensor(ndim,ndim)
        real(wp) :: rotation(ndim,ndim), gradient(ndim,1), rotated_gradient(ndim,1)
        real(wp) :: velocity(ndim,1), rotated_velocity(ndim,1)
        real(wp) :: iso_velocity(ndim,1), diag_velocity(ndim,1)
        integer :: status, i

        call permeability_reference(tensor)
        rotation(:,1) = [0.0_wp,1.0_wp,0.0_wp]
        rotation(:,2) = [-1.0_wp,0.0_wp,0.0_wp]
        rotation(:,3) = [0.0_wp,0.0_wp,1.0_wp]
        rotated_tensor = matmul(rotation,matmul(tensor,transpose(rotation)))
        gradient(:,1) = [2.0_wp,-1.0_wp,0.5_wp]
        rotated_gradient(:,1) = matmul(rotation,gradient(:,1))
        call material%init(darcy_anisotropic_t(10.0_wp,8.0_wp,6.0_wp, &
                           2.0_wp,-1.0_wp,1.5_wp,2.0_wp),1,status)
        call expect(status,cf_ok)
        call rotated%init(darcy_anisotropic_t(rotated_tensor(1,1), &
                          rotated_tensor(2,2),rotated_tensor(3,3), &
                          rotated_tensor(1,2),rotated_tensor(1,3), &
                          rotated_tensor(2,3),2.0_wp),1,status)
        call expect(status,cf_ok)
        call material%eval(gradient,velocity,status=status)
        call expect(status,cf_ok)
        call rotated%eval(rotated_gradient,rotated_velocity,status=status)
        call expect(status,cf_ok)
        do i = 1, ndim
            call close(rotated_velocity(i,1), &
                       dot_product(rotation(i,:),velocity(:,1)))
        end do

        call isotropic%init(darcy_isotropic_t(7.0_wp,2.0_wp),1,status)
        call expect(status,cf_ok)
        call diagonal%init(darcy_anisotropic_t(7.0_wp,7.0_wp,7.0_wp, &
                           0.0_wp,0.0_wp,0.0_wp,2.0_wp),1,status)
        call expect(status,cf_ok)
        call isotropic%eval(gradient,iso_velocity,status=status)
        call expect(status,cf_ok)
        call diagonal%eval(gradient,diag_velocity,status=status)
        call expect(status,cf_ok)
        do i = 1, ndim
            call close(iso_velocity(i,1),diag_velocity(i,1))
        end do
    end subroutine test_rotation_and_isotropic_limit

    subroutine test_sensitivities()
        type(material_t) :: tracked, upper, lower
        real(wp) :: values(nparameters), plus_values(nparameters), minus_values(nparameters)
        real(wp) :: gradient(ndim,2), velocity(ndim,2), dsdp(ndim,nparameters,2)
        real(wp) :: tangent(ndim,ndim,2), plus(ndim,2), minus(ndim,2)
        real(wp) :: bad(ndim,nparameters-1,2)
        real(wp) :: h, derivative, bound, relative_step
        integer, parameter :: ids(nparameters) = [7,6,5,4,3,2,1]
        integer :: status, iparam, ipoint, i, id

        values = [10.0_wp,8.0_wp,6.0_wp,2.0_wp,-1.0_wp,1.5_wp,2.0_wp]
        call tracked%init(model_from_values(values),2,status, &
            sensitivity_parameters=parameter_names)
        call expect(status,cf_ok)
        if (tracked%sensitivity_size() /= nparameters) error stop "Sensitivity count"
        if (tracked%state_storage_bytes() /= 0) error stop "Unnecessary derivative state"
        gradient(:,1) = [1.0_wp,-2.0_wp,0.5_wp]
        gradient(:,2) = [-1.2_wp,0.75_wp,0.4_wp]
        call tracked%eval(gradient,velocity,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(dsdp)) .or. .not. all(ieee_is_finite(tangent))) then
            error stop "Nonfinite anisotropic Darcy analytical derivatives"
        end if

        do iparam = 1, nparameters
            id = ids(iparam)
            if (id==nparameters) then
                h=values(id)*epsilon(1.0_wp)**(1.0_wp/3.0_wp)
            else
                ! Tensor coefficients enter affinely, and these finite
                ! perturbations retain positive definiteness.
                h=0.25_wp*abs(values(id))
            end if
            plus_values = values
            minus_values = values
            plus_values(id) = plus_values(id)+h
            minus_values(id) = minus_values(id)-h
            call upper%init(model_from_values(plus_values),2,status)
            call expect(status,cf_ok)
            call lower%init(model_from_values(minus_values),2,status)
            call expect(status,cf_ok)
            call upper%eval(gradient,plus,status=status)
            call expect(status,cf_ok)
            call lower%eval(gradient,minus,status=status)
            call expect(status,cf_ok)
            if (.not. all(ieee_is_finite(plus)) .or. .not. all(ieee_is_finite(minus))) then
                error stop "Nonfinite anisotropic Darcy parameter secant response"
            end if
            do ipoint = 1, 2
                do i = 1, ndim
                    derivative = (plus(i,ipoint)-minus(i,ipoint))/(2.0_wp*h)
                    bound=32.0_wp*epsilon(1.0_wp)* &
                        (abs(plus(i,ipoint))+abs(minus(i,ipoint)))/(2.0_wp*h)
                    if (id==nparameters) then
                        ! Exact central-difference bias for inverse viscosity.
                        relative_step=h/values(id)
                        bound=bound+abs(dsdp(i,iparam,ipoint))*relative_step**2/(1.0_wp-relative_step**2)
                    end if
                    if (.not. all(ieee_is_finite([derivative,bound]))) then
                        error stop "Nonfinite anisotropic Darcy parameter comparison"
                    end if
                    if (abs(derivative-dsdp(i,iparam,ipoint))>bound) then
                        error stop "Anisotropic Darcy parameter finite difference"
                    end if
                end do
            end do
        end do
        do iparam=1,nparameters
            do i=1,ndim
                call close(dsdp(i,iparam,1),reference_sensitivities(i,iparam))
            end do
        end do

        call tracked%commit(status)
        call expect(status,cf_ok)
        call tracked%eval(gradient,velocity,status=status)
        call expect(status,cf_invalid_argument)
        call tracked%eval(gradient,velocity,status=status,dsdp=bad)
        call expect(status,cf_invalid_shape)
        call tracked%eval(gradient,velocity,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
    end subroutine test_sensitivities

    subroutine test_common_parameter_scale(parameter_scale)
        !! Common scaling preserves the mobility, while all seven physical
        !! coefficient derivatives scale inversely with that common factor.
        real(wp), intent(in) :: parameter_scale
        type(material_t) :: material
        real(wp) :: tensor(ndim,ndim),gradient(ndim,1),velocity(ndim,1)
        real(wp) :: tangent(ndim,ndim,1),dsdp(ndim,nparameters,1)
        integer :: status,i,j

        call material%init(model_from_values(parameter_scale* &
            [10.0_wp,8.0_wp,6.0_wp,2.0_wp,-1.0_wp,1.5_wp,2.0_wp]),1,status, &
            sensitivity_parameters=parameter_names)
        call expect(status,cf_ok)
        call permeability_reference(tensor)
        gradient(:,1)=[1.0_wp,-2.0_wp,0.5_wp]
        call material%eval(gradient,velocity,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(velocity(1,1),-11.0_wp/4.0_wp)
        call close(velocity(2,1),53.0_wp/8.0_wp)
        call close(velocity(3,1),0.5_wp)
        do j=1,ndim
            do i=1,ndim
                call close(tangent(i,j,1),-0.5_wp*tensor(i,j))
            end do
        end do
        do j=1,nparameters
            do i=1,ndim
                call close(dsdp(i,j,1)*parameter_scale,reference_sensitivities(i,j))
            end do
        end do
    end subroutine test_common_parameter_scale

    subroutine test_cancelling_large_products()
        !! This SPD tensor has condition number five. Individual dot-product
        !! terms overflow although velocity, tangent and sensitivities do not.
        type(material_t) :: material
        real(wp), parameter :: tensor(ndim,ndim) = reshape( &
            [3.0_wp,2.0_wp,0.0_wp,2.0_wp,3.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp],[ndim,ndim])
        real(wp), parameter :: expected_derivative(ndim,nparameters) = reshape( &
            [2.0_wp,-2.0_wp,0.0_wp, 0.0_wp,0.0_wp,2.0_wp, 0.0_wp,0.0_wp,-2.0_wp, &
             2.0_wp,-2.0_wp,0.0_wp, 0.0_wp,0.0_wp,0.0_wp, &
             0.0_wp,2.0_wp,0.0_wp, -2.0_wp,0.0_wp,0.0_wp],[ndim,nparameters])
        real(wp) :: gradient(ndim,1),velocity(ndim,1),tangent(ndim,ndim,1)
        real(wp) :: dsdp(ndim,nparameters,1)
        real(wp) :: coefficient_scale
        integer :: status,i,j

        coefficient_scale=huge(1.0_wp)/4.0_wp
        call material%init(darcy_anisotropic_t(3.0_wp*coefficient_scale,3.0_wp*coefficient_scale, &
            coefficient_scale,2.0_wp*coefficient_scale,0.0_wp,0.0_wp,1.0_wp),1,status, &
            sensitivity_parameters=parameter_names)
        call expect(status,cf_ok)
        gradient(:,1)=[2.0_wp,-2.0_wp,0.0_wp]
        call material%eval(gradient,velocity,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(velocity(1,1)/coefficient_scale,-2.0_wp)
        call close(velocity(2,1)/coefficient_scale,2.0_wp)
        call close(velocity(3,1)/coefficient_scale,0.0_wp)
        do i=1,ndim
            call close(dsdp(i,1,1)/coefficient_scale,expected_derivative(i,1))
        end do
        do j=1,ndim
            do i=1,ndim
                call close(tangent(i,j,1)/coefficient_scale,-tensor(i,j))
            end do
        end do
        do j=2,nparameters
            do i=1,ndim
                call close(dsdp(i,j,1),expected_derivative(i,j))
            end do
        end do
    end subroutine test_cancelling_large_products

    subroutine test_cancellation_with_small_remainder()
        !! Overflowing opposite products cancel exactly, leaving a tiny but
        !! normal third contribution. Coordinate permutations must preserve it.
        integer, parameter :: norders = 6
        integer, parameter :: orders(ndim,norders) = reshape( &
            [1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1],[ndim,norders])
        type(material_t) :: material
        real(wp) :: base_tensor(ndim,ndim), tensor(ndim,ndim)
        real(wp) :: base_gradient(ndim), base_scale(ndim), base_reference(ndim)
        real(wp) :: gradient(ndim,1), velocity(ndim,1), tangent(ndim,ndim,1), dsdp(ndim,2,1)
        real(wp) :: large, small
        integer :: status, icase, i, j, coordinate

        large = scale(0.75_wp,maxexponent(1.0_wp))
        small = tiny(1.0_wp)
        base_tensor(:,1) = [large,0.5_wp*large,small]
        base_tensor(:,2) = [0.5_wp*large,0.375_wp*large,0.0_wp]
        base_tensor(:,3) = [small,0.0_wp,large]
        base_gradient = [2.0_wp,-4.0_wp,1.0_wp]
        base_scale = [small,large,large]
        base_reference = [-1.0_wp,0.5_wp,-1.0_wp]
        do icase = 1,norders
            tensor = base_tensor(orders(:,icase),orders(:,icase))
            gradient(:,1) = base_gradient(orders(:,icase))
            call material%init(darcy_anisotropic_t(tensor(1,1),tensor(2,2),tensor(3,3), &
                tensor(1,2),tensor(1,3),tensor(2,3),1.0_wp),1,status, &
                sensitivity_parameters=[character(len=18) :: "dynamic_viscosity","permeability_13"])
            call expect(status,cf_ok)
            call material%eval(gradient,velocity,tangent,status,dsdp=dsdp)
            call expect(status,cf_ok)
            if (.not. all(ieee_is_finite(velocity)) .or. .not. all(ieee_is_finite(tangent)) .or. &
                .not. all(ieee_is_finite(dsdp))) error stop "Nonfinite Darcy cancellation remainder"
            do i = 1,ndim
                coordinate = orders(i,icase)
                call close(velocity(i,1)/base_scale(coordinate),base_reference(coordinate))
                call close(dsdp(i,1,1)/base_scale(coordinate),-base_reference(coordinate))
                do j = 1,ndim
                    if (abs(tensor(i,j)) > 0.0_wp) then
                        call close(tangent(i,j,1)/tensor(i,j),-1.0_wp)
                    else
                        call close(tangent(i,j,1),0.0_wp)
                    end if
                end do
            end do
            call close(dsdp(1,2,1),-gradient(3,1))
            call close(dsdp(2,2,1),0.0_wp)
            call close(dsdp(3,2,1),-gradient(1,1))
        end do
    end subroutine test_cancellation_with_small_remainder

    subroutine test_sensitivity_range()
        type(material_t) :: tracked,plain
        real(wp) :: tensor(ndim,ndim),gradient(ndim,1),velocity(ndim,1)
        real(wp) :: tangent(ndim,ndim,1),dsdp(ndim,nparameters,1)
        real(wp) :: permeability_scale,viscosity
        integer :: status,i,j

        ! All velocities round to zero, but the coupled viscosity derivative
        ! must be computed from original factors and remains finite/nonzero.
        permeability_scale=tiny(1.0_wp)
        viscosity=sqrt(permeability_scale)
        call tracked%init(darcy_anisotropic_t(10.0_wp*permeability_scale,8.0_wp*permeability_scale, &
            6.0_wp*permeability_scale,2.0_wp*permeability_scale,-permeability_scale, &
            1.5_wp*permeability_scale,viscosity),1,status,sensitivity_parameters=parameter_names)
        call expect(status,cf_ok)
        call permeability_reference(tensor)
        gradient(:,1)=permeability_scale*[1.0_wp,-2.0_wp,0.5_wp]
        call tracked%eval(gradient,velocity,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        do i=1,ndim
            call close(velocity(i,1),0.0_wp)
            call close(dsdp(i,1,1)/permeability_scale,4.0_wp*reference_sensitivities(i,1))
        end do
        do j=1,ndim
            do i=1,ndim
                call close(tangent(i,j,1)/viscosity,-tensor(i,j))
            end do
        end do
        do j=2,nparameters
            do i=1,ndim
                call close(dsdp(i,j,1)/viscosity,2.0_wp*reference_sensitivities(i,j))
            end do
        end do

        viscosity=permeability_scale
        call plain%init(darcy_anisotropic_t(permeability_scale,permeability_scale,permeability_scale, &
            0.0_wp,0.0_wp,0.0_wp,viscosity),1,status)
        call expect(status,cf_ok)
        call tracked%init(darcy_anisotropic_t(permeability_scale,permeability_scale,permeability_scale, &
            0.0_wp,0.0_wp,0.0_wp,viscosity),1,status,sensitivity_parameters=parameter_names)
        call expect(status,cf_ok)
        gradient(:,1)=[8.0_wp,0.0_wp,0.0_wp]
        call plain%eval(gradient,velocity,status=status)
        call expect(status,cf_ok)
        call close(velocity(1,1),-8.0_wp)
        call tracked%eval(gradient,velocity,status=status,dsdp=dsdp)
        call expect(status,cf_evaluation_failed)
        if (tracked%last_failure_point()/=1) error stop "Anisotropic Darcy derivative overflow point"
        call tracked%commit(status)
        call expect(status,cf_no_trial)
        gradient(1,1)=1.0_wp
        call tracked%eval(gradient,velocity,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(dsdp(1,1,1)*viscosity,1.0_wp)
        call close(dsdp(1,nparameters,1)*viscosity,-1.0_wp)
    end subroutine test_sensitivity_range

    subroutine test_divided_mobility_underflow()
        !! A mobility entry may round to zero while its original permeability
        !! times the prescribed gradient, divided by viscosity, is nonzero.
        type(material_t) :: material
        real(wp) :: gradient(ndim,1),velocity(ndim,1),tangent(ndim,ndim,1)
        real(wp) :: dsdp(ndim,1,1),mixed_dsdp(ndim,2,1)
        real(wp) :: viscosity,large,small,gradient_over_viscosity,expected_small
        integer :: status

        large=huge(1.0_wp)
        viscosity=sqrt(large)
        small=tiny(1.0_wp)
        call material%init(darcy_anisotropic_t(small,1.0_wp,1.0_wp,0.0_wp,0.0_wp,0.0_wp,viscosity), &
            1,status,sensitivity_parameters=[character(len=15) :: "permeability_11"])
        call expect(status,cf_ok)
        gradient(:,1)=[viscosity,0.0_wp,0.0_wp]
        call material%eval(gradient,velocity,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(velocity(1,1)/small,-1.0_wp)
        call close(velocity(2,1),0.0_wp)
        call close(velocity(3,1),0.0_wp)
        call close(dsdp(1,1,1),-1.0_wp)
        call close(tangent(1,1,1),0.0_wp)
        call close(tangent(2,2,1)*viscosity,-1.0_wp)
        call close(tangent(3,3,1)*viscosity,-1.0_wp)

        ! The dominant entry in row one multiplies zero. Normalizing that
        ! entire row would erase the much smaller, physically active entry.
        small=nearest(0.0_wp,1.0_wp)
        call material%init(darcy_anisotropic_t(large/4.0_wp,1.0_wp,1.0_wp,small,0.0_wp,0.0_wp,viscosity), &
            1,status,sensitivity_parameters=[character(len=18) :: "dynamic_viscosity","permeability_12"])
        call expect(status,cf_ok)
        gradient(:,1)=[0.0_wp,large/4.0_wp,0.0_wp]
        gradient_over_viscosity=gradient(2,1)/viscosity
        expected_small=-small*gradient_over_viscosity
        call material%eval(gradient,velocity,tangent,status,dsdp=mixed_dsdp)
        call expect(status,cf_ok)
        call close(velocity(1,1)/expected_small,1.0_wp)
        call close(velocity(2,1)/gradient_over_viscosity,-1.0_wp)
        call close(velocity(3,1),0.0_wp)
        call close(mixed_dsdp(2,1,1),0.25_wp)
        call close(mixed_dsdp(1,2,1)/gradient_over_viscosity,-1.0_wp)
        call close(tangent(1,2,1),0.0_wp)
        call close(tangent(2,2,1)*viscosity,-1.0_wp)
    end subroutine test_divided_mobility_underflow

    subroutine test_subnormal_product_sum()
        !! The first row contributes two terms of 0.375d. Their sum rounds
        !! to d, although restoring the terms individually loses both.
        real(wp), parameter :: normal_tensor(ndim,ndim)=reshape( &
            [0.25_wp,0.25_wp,0.0_wp,0.25_wp,0.5_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp],[ndim,ndim])
        real(wp), parameter :: normal_velocity(ndim)=[-1.0_wp,-2.0_wp,0.0_wp]
        type(material_t) :: material
        real(wp) :: gradient(ndim,1),velocity(ndim,1),tangent(ndim,ndim,1)
        real(wp) :: dsdp(ndim,2,1)
        real(wp) :: smallest
        integer :: status,i,j

        smallest=nearest(0.0_wp,1.0_wp)
        call material%init(darcy_anisotropic_t(smallest,2.0_wp*smallest,4.0_wp*smallest, &
            smallest,0.0_wp,0.0_wp,4.0_wp),1,status, &
            sensitivity_parameters=[character(len=18) :: "dynamic_viscosity","permeability_12"])
        call expect(status,cf_ok)
        gradient(:,1)=[1.5_wp,1.5_wp,0.0_wp]
        call material%eval(gradient,velocity,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(velocity(1,1)/smallest,-1.0_wp)
        call close(velocity(2,1)/smallest,-1.0_wp)
        call close(velocity(3,1),0.0_wp)
        call close(tangent(1,1,1),0.0_wp)
        call close(tangent(2,2,1),0.0_wp)
        call close(tangent(3,3,1)/smallest,-1.0_wp)
        call close(dsdp(1,1,1),0.0_wp)
        call close(dsdp(2,1,1),0.0_wp)
        call close(dsdp(3,1,1),0.0_wp)
        call close(dsdp(1,2,1),-3.0_wp/8.0_wp)
        call close(dsdp(2,2,1),-3.0_wp/8.0_wp)
        call close(dsdp(3,2,1),0.0_wp)

        ! The supported binary kinds have gradual underflow. With normal
        ! mobility and subnormal gradients, the two first-row products are
        ! each m/2: the full row rounds to m, while separate products vanish.
        if (radix(1.0_wp) /= 2 .or. smallest >= tiny(1.0_wp)) return
        call material%init(darcy_anisotropic_t(0.25_wp,0.5_wp,1.0_wp,0.25_wp,0.0_wp,0.0_wp,1.0_wp), &
            1,status,sensitivity_parameters=[character(len=18) :: "dynamic_viscosity","permeability_12"])
        call expect(status,cf_ok)
        gradient(:,1)=[2.0_wp*smallest,2.0_wp*smallest,0.0_wp]
        call material%eval(gradient,velocity,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(velocity)) .or. .not. all(ieee_is_finite(tangent)) .or. &
            .not. all(ieee_is_finite(dsdp))) error stop "Nonfinite normal-mobility subnormal Darcy sum"
        do i=1,ndim
            call close(velocity(i,1)/smallest,normal_velocity(i))
            call close(dsdp(i,1,1)/smallest,-normal_velocity(i))
            call close(dsdp(i,2,1)/smallest,merge(-2.0_wp,0.0_wp,i<ndim))
            do j=1,ndim
                call close(tangent(i,j,1),-normal_tensor(i,j))
            end do
        end do
    end subroutine test_subnormal_product_sum

    subroutine test_failures()
        type(material_t) :: material
        real(wp) :: gradient(ndim,2), velocity(ndim,2), invalid(ndim-1,2)
        real(wp) :: extra(ndim,1,2)
        integer :: status

        call material%init(darcy_anisotropic_t(0.0_wp,8.0_wp,6.0_wp, &
                           0.0_wp,0.0_wp,0.0_wp,2.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(darcy_anisotropic_t(10.0_wp,8.0_wp,6.0_wp, &
                           10.0_wp,0.0_wp,0.0_wp,2.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(darcy_anisotropic_t(10.0_wp,8.0_wp,6.0_wp, &
                           0.0_wp,0.0_wp,0.0_wp,0.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(darcy_anisotropic_t(10.0_wp,8.0_wp,6.0_wp, &
                           ieee_value(0.0_wp,ieee_quiet_nan), &
                           0.0_wp,0.0_wp,2.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(darcy_anisotropic_t(huge(1.0_wp)/4.0_wp,huge(1.0_wp)/5.0_wp, &
                           huge(1.0_wp)/6.0_wp,0.0_wp,0.0_wp,0.0_wp,0.125_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(darcy_anisotropic_t(1.0e-12_wp,2.0e-12_wp, &
                           3.0e-12_wp,0.0_wp,0.0_wp,0.0_wp,1.0e-3_wp),2,status)
        call expect(status,cf_ok)
        gradient(:,1) = [1.0_wp,0.0_wp,0.0_wp]
        gradient(:,2) = gradient(:,1)
        call material%eval(gradient,velocity,status=status)
        call expect(status,cf_ok)
        call close(velocity(1,1)/1.0e-9_wp,-1.0_wp)

        call material%init(darcy_anisotropic_t(10.0_wp,8.0_wp,6.0_wp, &
                           2.0_wp,-1.0_wp,1.5_wp,2.0_wp),2,status)
        call expect(status,cf_ok)
        gradient(1,2) = ieee_value(0.0_wp,ieee_positive_inf)
        call material%eval(gradient,velocity,status=status)
        call expect(status,cf_evaluation_failed)
        if (material%last_failure_point() /= 2) error stop "Bad point index"
        call material%commit(status)
        call expect(status,cf_no_trial)
        gradient(1,2) = ieee_value(0.0_wp,ieee_quiet_nan)
        call material%eval(gradient,velocity,status=status)
        call expect(status,cf_evaluation_failed)
        gradient(1,2) = huge(1.0_wp)
        call material%eval(gradient,velocity,status=status)
        call expect(status,cf_evaluation_failed)
        gradient(1,2) = 1.0_wp
        call material%eval(gradient,velocity,status=status)
        call expect(status,cf_ok)
        call material%eval(invalid,velocity,status=status)
        call expect(status,cf_invalid_shape)
        call material%eval(gradient,velocity,status=status,dsdp=extra)
        call expect(status,cf_unsupported)
    end subroutine test_failures

end program test_darcy_anisotropic
