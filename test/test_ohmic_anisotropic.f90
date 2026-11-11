! SPDX-License-Identifier: Apache-2.0

program test_ohmic_anisotropic
    !! Test convention: electric field in V/m, current density in A/m2 and
    !! conductivity in S/m. Derived Joule power density is W/m3.
    !!
    !! Independent tensor, energy, covariance, Jacobian and sensitivity tests.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
    use constituf, only: wp, material_t, cf_ok, cf_invalid_argument, &
                        cf_invalid_shape, cf_evaluation_failed, cf_no_trial, cf_unsupported
    use constituf_ohmic_anisotropic, only: ohmic_anisotropic_t
    use constituf_ohmic_isotropic, only: ohmic_isotropic_t
    implicit none (type, external)

    integer, parameter :: ndim = 3, nparameters = 6

    call test_reference_and_schema()
    call test_tangent_and_batch()
    call test_rotation_and_isotropic_limit()
    call test_parameter_sensitivities()
    call test_scaled_conductivity(sqrt(tiny(1.0_wp)))
    call test_scaled_conductivity(sqrt(huge(1.0_wp)))
    call test_cancelling_large_products()
    call test_cancellation_with_small_remainder()
    call test_invalid_and_recovery()

contains

    function model_from_values(values) result(model)
        real(wp), intent(in) :: values(nparameters)
        type(ohmic_anisotropic_t) :: model

        model = ohmic_anisotropic_t(values(1),values(2),values(3), &
                                         values(4),values(5),values(6))
    end function model_from_values

    subroutine expect(actual, wanted)
        integer, intent(in) :: actual, wanted

        if (actual /= wanted) then
            print *, "actual status", actual, "wanted", wanted
            error stop "Anisotropic ohmic status mismatch"
        end if
    end subroutine expect

    subroutine close(actual, wanted)
        real(wp), intent(in) :: actual, wanted

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(wanted)) then
            error stop "Nonfinite anisotropic ohmic comparison"
        end if
        if (abs(actual-wanted) > 256.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(wanted))) then
            print *, "actual", actual, "wanted", wanted
            error stop "Anisotropic ohmic reference mismatch"
        end if
    end subroutine close

    subroutine tensor_reference(conductivity)
        real(wp), intent(out) :: conductivity(ndim,ndim)

        conductivity(:,1) = [10.0_wp,2.0_wp,-1.0_wp]
        conductivity(:,2) = [2.0_wp,8.0_wp,1.5_wp]
        conductivity(:,3) = [-1.0_wp,1.5_wp,6.0_wp]
    end subroutine tensor_reference

    subroutine test_reference_and_schema()
        type(material_t) :: material
        real(wp) :: conductivity(ndim,ndim)
        real(wp) :: electric(ndim,3), current_density(ndim,3), tangent(ndim,ndim,3)
        integer :: status, i, j, ipoint

        call material%init(ohmic_anisotropic_t(10.0_wp,8.0_wp,6.0_wp, &
                           2.0_wp,-1.0_wp,1.5_wp),3,status)
        call expect(status,cf_ok)
        if (material%input_size() /= ndim .or. material%output_size() /= ndim) then
            error stop "Anisotropic ohmic field sizes"
        end if
        if (material%state_size() /= 0 .or. material%state_storage_bytes() /= 0) then
            error stop "Anisotropic ohmic must be stateless"
        end if
        if (.not. material%field_matches("input", "electric.field", &
            "vector", "electric_field", "spatial", [ndim])) then
            error stop "Anisotropic electric field metadata"
        end if
        if (.not. material%field_matches("output", "electric.current_density", &
            "vector", "electric_current_density", "spatial", [ndim])) then
            error stop "Anisotropic electric current_density metadata"
        end if
        if (material%field_matches("input", "electric.field", &
            "vector", "electric_field", "reference", [ndim])) then
            error stop "Incompatible anisotropic frame accepted"
        end if
        if (material%field_matches("output", "electric.current_density", &
            "vector", "electric_displacement", "spatial", [ndim])) then
            error stop "Incompatible electric current measure accepted"
        end if

        call tensor_reference(conductivity)
        electric(:,1) = [1.0_wp,-2.0_wp,0.5_wp]
        electric(:,2) = [-3.0_wp,4.0_wp,1.25_wp]
        electric(:,3) = 0.0_wp
        call material%eval(electric,current_density,tangent,status)
        call expect(status,cf_ok)
        call close(current_density(1,1),5.5_wp)
        call close(current_density(2,1),-13.25_wp)
        call close(current_density(3,1),-1.0_wp)
        do ipoint = 1, 3
            do i = 1, ndim
                call close(current_density(i,ipoint), &
                           dot_product(conductivity(i,:),electric(:,ipoint)))
                do j = 1, ndim
                    call close(tangent(i,j,ipoint),conductivity(i,j))
                end do
            end do
            if (dot_product(electric(:,ipoint),current_density(:,ipoint)) < 0.0_wp) then
                error stop "Negative anisotropic Joule power density"
            end if
        end do

        call material%commit(status)
        call expect(status,cf_ok)
        call material%eval(electric,current_density,status=status)
        call expect(status,cf_ok)
        call material%revert()
        call material%commit(status)
        call expect(status,cf_no_trial)
    end subroutine test_reference_and_schema

    subroutine test_tangent_and_batch()
        type(material_t) :: batch, single
        real(wp) :: electric(ndim,2), current_density(ndim,2), tangent(ndim,ndim,2)
        real(wp) :: changed(ndim,2), positive(ndim,2), negative(ndim,2), alone(ndim,1)
        real(wp) :: derivative, h, error_bound
        integer :: status, ipoint, i, j

        call batch%init(ohmic_anisotropic_t(10.0_wp,8.0_wp,6.0_wp, &
                        2.0_wp,-1.0_wp,1.5_wp),2,status)
        call expect(status,cf_ok)
        call single%init(ohmic_anisotropic_t(10.0_wp,8.0_wp,6.0_wp, &
                         2.0_wp,-1.0_wp,1.5_wp),1,status)
        call expect(status,cf_ok)
        electric(:,1) = [1.2_wp,-2.3_wp,0.7_wp]
        electric(:,2) = [0.1_wp,3.4_wp,-1.1_wp]
        call batch%eval(electric,current_density,tangent,status)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(tangent))) error stop "Nonfinite anisotropic ohmic tangent"
        call single%eval(electric(:,1:1),alone,status=status)
        call expect(status,cf_ok)
        do i = 1, ndim
            call close(current_density(i,1),alone(i,1))
        end do

        do ipoint = 1, 2
            do j = 1, ndim
                ! Current is affine in each field component, so a finite
                ! central step avoids cancellation without truncation error.
                h = 0.25_wp*maxval(abs(electric(:,ipoint)))
                changed = electric
                changed(j,ipoint) = electric(j,ipoint)+h
                call batch%eval(changed,positive,status=status)
                call expect(status,cf_ok)
                changed(j,ipoint) = electric(j,ipoint)-h
                call batch%eval(changed,negative,status=status)
                call expect(status,cf_ok)
                if (.not. all(ieee_is_finite(positive)) .or. .not. all(ieee_is_finite(negative))) then
                    error stop "Nonfinite anisotropic ohmic tangent secant response"
                end if
                error_bound = 32.0_wp*epsilon(1.0_wp)* &
                    maxval(abs(positive(:,ipoint))+abs(negative(:,ipoint)))/(2.0_wp*h)
                if (.not. ieee_is_finite(error_bound)) error stop "Nonfinite anisotropic ohmic tangent bound"
                do i = 1, ndim
                    derivative = (positive(i,ipoint)-negative(i,ipoint))/(2.0_wp*h)
                    if (.not. ieee_is_finite(derivative) .or. &
                        abs(derivative-tangent(i,j,ipoint)) > error_bound) then
                        error stop "Anisotropic ohmic tangent finite difference"
                    end if
                end do
            end do
        end do
    end subroutine test_tangent_and_batch

    subroutine test_rotation_and_isotropic_limit()
        type(material_t) :: original, rotated, isotropic, diagonal
        real(wp) :: tensor(ndim,ndim), tensor_rotated(ndim,ndim), rotation(ndim,ndim)
        real(wp) :: electric(ndim,1), new_electric(ndim,1)
        real(wp) :: current_density(ndim,1), new_current_density(ndim,1)
        real(wp) :: isotropic_d(ndim,1), diagonal_d(ndim,1)
        integer :: status, i

        call tensor_reference(tensor)
        rotation(:,1) = [0.0_wp,1.0_wp,0.0_wp]
        rotation(:,2) = [-1.0_wp,0.0_wp,0.0_wp]
        rotation(:,3) = [0.0_wp,0.0_wp,1.0_wp]
        tensor_rotated = matmul(rotation,matmul(tensor,transpose(rotation)))
        electric(:,1) = [2.0_wp,-1.0_wp,0.5_wp]
        new_electric(:,1) = matmul(rotation,electric(:,1))

        call original%init(ohmic_anisotropic_t(10.0_wp,8.0_wp,6.0_wp, &
                           2.0_wp,-1.0_wp,1.5_wp),1,status)
        call expect(status,cf_ok)
        call rotated%init(ohmic_anisotropic_t(tensor_rotated(1,1), &
                          tensor_rotated(2,2),tensor_rotated(3,3), &
                          tensor_rotated(1,2),tensor_rotated(1,3), &
                          tensor_rotated(2,3)),1,status)
        call expect(status,cf_ok)
        call original%eval(electric,current_density,status=status)
        call expect(status,cf_ok)
        call rotated%eval(new_electric,new_current_density,status=status)
        call expect(status,cf_ok)
        do i = 1, ndim
            call close(new_current_density(i,1), &
                       dot_product(rotation(i,:),current_density(:,1)))
        end do

        call isotropic%init(ohmic_isotropic_t(7.0_wp),1,status)
        call expect(status,cf_ok)
        call diagonal%init(ohmic_anisotropic_t(7.0_wp,7.0_wp,7.0_wp, &
                           0.0_wp,0.0_wp,0.0_wp),1,status)
        call expect(status,cf_ok)
        call isotropic%eval(electric,isotropic_d,status=status)
        call expect(status,cf_ok)
        call diagonal%eval(electric,diagonal_d,status=status)
        call expect(status,cf_ok)
        do i = 1, ndim
            call close(isotropic_d(i,1),diagonal_d(i,1))
        end do
    end subroutine test_rotation_and_isotropic_limit

    subroutine test_parameter_sensitivities()
        type(material_t) :: tracked, upper, lower
        real(wp), parameter :: expected_derivative(ndim,nparameters) = reshape( &
            [0.0_wp,2.1_wp,-0.9_wp, 2.1_wp,0.0_wp,1.3_wp, &
             -0.9_wp,1.3_wp,0.0_wp, 0.0_wp,0.0_wp,2.1_wp, &
             0.0_wp,-0.9_wp,0.0_wp, 1.3_wp,0.0_wp,0.0_wp],[ndim,nparameters])
        real(wp) :: values(nparameters), plus_values(nparameters), minus_values(nparameters)
        real(wp) :: electric(ndim,2), current_density(ndim,2), dsdp(ndim,nparameters,2)
        real(wp) :: tangent(ndim,ndim,2), positive(ndim,2), negative(ndim,2)
        real(wp) :: wrong(ndim,nparameters-1,2)
        real(wp) :: derivative, h, error_bound
        integer, parameter :: parameter_ids(nparameters) = [6,5,4,3,2,1]
        integer :: status, iparam, ipoint, i, id

        values = [10.0_wp,8.0_wp,6.0_wp,2.0_wp,-1.0_wp,1.5_wp]
        call tracked%init(model_from_values(values),2,status, &
            sensitivity_parameters=[character(len=15) :: "conductivity_23", &
                                    "conductivity_13","conductivity_12", &
                                    "conductivity_33","conductivity_22", &
                                    "conductivity_11"])
        call expect(status,cf_ok)
        if (tracked%sensitivity_size() /= nparameters) then
            error stop "Missing anisotropic ohmic parameter sensitivities"
        end if
        if (tracked%state_storage_bytes() /= 0) error stop "Unexpected derivative history"
        electric(:,1) = [1.3_wp,-0.9_wp,2.1_wp]
        electric(:,2) = [-1.2_wp,0.75_wp,0.4_wp]
        call tracked%eval(electric,current_density,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(dsdp))) error stop "Nonfinite anisotropic ohmic sensitivity"
        if (.not. all(ieee_is_finite(tangent))) error stop "Nonfinite tracked anisotropic ohmic tangent"

        do iparam = 1, nparameters
            id = parameter_ids(iparam)
            ! Each independent tensor coefficient enters current affinely.
            ! These finite perturbations retain positive definiteness.
            h = 0.25_wp*abs(values(id))
            plus_values = values
            minus_values = values
            plus_values(id) = plus_values(id)+h
            minus_values(id) = minus_values(id)-h
            call upper%init(model_from_values(plus_values),2,status)
            call expect(status,cf_ok)
            call lower%init(model_from_values(minus_values),2,status)
            call expect(status,cf_ok)
            call upper%eval(electric,positive,status=status)
            call expect(status,cf_ok)
            call lower%eval(electric,negative,status=status)
            call expect(status,cf_ok)
            if (.not. all(ieee_is_finite(positive)) .or. .not. all(ieee_is_finite(negative))) then
                error stop "Nonfinite anisotropic ohmic parameter secant response"
            end if
            do ipoint = 1, 2
                error_bound = 32.0_wp*epsilon(1.0_wp)* &
                    maxval(abs(positive(:,ipoint))+abs(negative(:,ipoint)))/(2.0_wp*h)
                if (.not. ieee_is_finite(error_bound)) error stop "Nonfinite anisotropic ohmic parameter bound"
                do i = 1, ndim
                    derivative = (positive(i,ipoint)-negative(i,ipoint))/(2.0_wp*h)
                    if (.not. ieee_is_finite(derivative) .or. &
                        abs(derivative-dsdp(i,iparam,ipoint)) > error_bound) then
                        error stop "Anisotropic ohmic parameter derivative"
                    end if
                end do
            end do
        end do
        do iparam = 1,nparameters
            do i = 1,ndim
                call close(dsdp(i,iparam,1),expected_derivative(i,iparam))
            end do
        end do

        call tracked%commit(status)
        call expect(status,cf_ok)
        call tracked%eval(electric,current_density,status=status)
        call expect(status,cf_invalid_argument)
        call tracked%eval(electric,current_density,status=status,dsdp=wrong)
        call expect(status,cf_invalid_shape)
        call tracked%eval(electric,current_density,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
    end subroutine test_parameter_sensitivities

    subroutine test_scaled_conductivity(scale)
        !! Common coefficient scaling preserves definiteness. These scales
        !! would underflow/overflow raw cubic determinant products.
        real(wp), intent(in) :: scale
        type(material_t) :: material
        real(wp) :: conductivity(ndim,ndim), electric(ndim,1), current_density(ndim,1)
        real(wp) :: tangent(ndim,ndim,1), dsdp(ndim,1,1)
        integer :: status, i, j

        call material%init(model_from_values(scale*[10.0_wp,8.0_wp,6.0_wp, &
            2.0_wp,-1.0_wp,1.5_wp]),1,status, &
            sensitivity_parameters=[character(len=15) :: "conductivity_13"])
        call expect(status,cf_ok)
        call tensor_reference(conductivity)
        electric(:,1)=[1.0_wp,-2.0_wp,0.5_wp]
        call material%eval(electric,current_density,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(current_density(1,1)/scale,5.5_wp)
        call close(current_density(2,1)/scale,-13.25_wp)
        call close(current_density(3,1)/scale,-1.0_wp)
        call close(dot_product(current_density(:,1)/scale,electric(:,1)),31.5_wp)
        do j=1,ndim
            do i=1,ndim
                call close(tangent(i,j,1)/scale,conductivity(i,j))
            end do
        end do
        call close(dsdp(1,1,1),0.5_wp)
        call close(dsdp(2,1,1),0.0_wp)
        call close(dsdp(3,1,1),1.0_wp)
    end subroutine test_scaled_conductivity

    subroutine test_cancelling_large_products()
        !! This SPD tensor has condition number five. Individual products
        !! overflow, although their sums and every requested result are finite.
        type(material_t) :: material
        real(wp), parameter :: tensor(ndim,ndim) = reshape( &
            [3.0_wp,2.0_wp,0.0_wp,2.0_wp,3.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp],[ndim,ndim])
        real(wp), parameter :: expected_derivative(ndim,nparameters) = reshape( &
            [0.0_wp,0.0_wp,-2.0_wp, 0.0_wp,0.0_wp,2.0_wp, &
             -2.0_wp,2.0_wp,0.0_wp, 0.0_wp,0.0_wp,0.0_wp, &
             0.0_wp,-2.0_wp,0.0_wp, 2.0_wp,0.0_wp,0.0_wp],[ndim,nparameters])
        real(wp) :: electric(ndim,1), current_density(ndim,1), tangent(ndim,ndim,1)
        real(wp) :: dsdp(ndim,nparameters,1)
        real(wp) :: scale
        integer :: status, i, j

        scale=huge(1.0_wp)/4.0_wp
        call material%init(ohmic_anisotropic_t(3.0_wp*scale,3.0_wp*scale,scale, &
            2.0_wp*scale,0.0_wp,0.0_wp),1,status, &
            sensitivity_parameters=[character(len=15) :: "conductivity_23", &
                "conductivity_13","conductivity_12","conductivity_33", &
                "conductivity_22","conductivity_11"])
        call expect(status,cf_ok)
        electric(:,1)=[2.0_wp,-2.0_wp,0.0_wp]
        call material%eval(electric,current_density,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(current_density(1,1)/scale,2.0_wp)
        call close(current_density(2,1)/scale,-2.0_wp)
        call close(current_density(3,1)/scale,0.0_wp)
        do j=1,ndim
            do i=1,ndim
                call close(tangent(i,j,1)/scale,tensor(i,j))
            end do
        end do
        do j=1,nparameters
            do i=1,ndim
                call close(dsdp(i,j,1),expected_derivative(i,j))
            end do
        end do
    end subroutine test_cancelling_large_products

    subroutine test_cancellation_with_small_remainder()
        !! Equal overflowing terms cancel exactly before a tiny third
        !! term remains. Common scaling must not erase that finite current.
        type(material_t) :: material
        integer, parameter :: npermutations = 6
        integer, parameter :: permutations(ndim,npermutations) = reshape( &
            [1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1],[ndim,npermutations])
        real(wp) :: conductivity(ndim,ndim), permuted(ndim,ndim)
        real(wp) :: electric(ndim,1), current_density(ndim,1), tangent(ndim,ndim,1)
        real(wp) :: original_electric(ndim), expected(ndim), dsdp(ndim,1,1)
        real(wp) :: large, small
        integer :: status, ipermutation, i

        large = scale(0.75_wp,maxexponent(1.0_wp))
        small = tiny(1.0_wp)
        conductivity(:,1) = [large,0.5_wp*large,small]
        conductivity(:,2) = [0.5_wp*large,0.375_wp*large,0.0_wp]
        conductivity(:,3) = [small,0.0_wp,large]
        original_electric = [2.0_wp,-4.0_wp,1.0_wp]
        expected = [small,-0.5_wp*large,large]

        do ipermutation = 1, npermutations
            permuted = conductivity(permutations(:,ipermutation),permutations(:,ipermutation))
            electric(:,1) = original_electric(permutations(:,ipermutation))
            call material%init(ohmic_anisotropic_t(permuted(1,1),permuted(2,2), &
                permuted(3,3),permuted(1,2),permuted(1,3),permuted(2,3)),1,status, &
                sensitivity_parameters=[character(len=15) :: "conductivity_13"])
            call expect(status,cf_ok)
            call material%eval(electric,current_density,tangent,status,dsdp=dsdp)
            call expect(status,cf_ok)
            if (.not. all(ieee_is_finite(tangent))) error stop "Nonfinite cancellation tangent"
            do i = 1, ndim
                ! Normalize each component separately so the large currents
                ! cannot mask loss of the tiny exact remainder.
                call close(current_density(i,1)/expected(permutations(i,ipermutation)),1.0_wp)
            end do
            call close(dsdp(1,1,1),electric(3,1))
            call close(dsdp(2,1,1),0.0_wp)
            call close(dsdp(3,1,1),electric(1,1))
        end do
    end subroutine test_cancellation_with_small_remainder

    subroutine test_invalid_and_recovery()
        type(material_t) :: material
        real(wp) :: electric(ndim,2), current_density(ndim,2), wrong(ndim-1,2)
        real(wp) :: extra(ndim,1,2)
        integer :: status

        call material%init(ohmic_anisotropic_t(0.0_wp,8.0_wp,6.0_wp, &
                           0.0_wp,0.0_wp,0.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(ohmic_anisotropic_t(10.0_wp,8.0_wp,6.0_wp, &
                           10.0_wp,0.0_wp,0.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(ohmic_anisotropic_t(10.0_wp,8.0_wp,6.0_wp, &
                           ieee_value(0.0_wp,ieee_quiet_nan),0.0_wp,0.0_wp), &
                           2,status)
        call expect(status,cf_invalid_argument)
        call material%init(ohmic_anisotropic_t(1.0e-12_wp,1.0e-12_wp, &
                           1.0e-12_wp,0.0_wp,0.0_wp,0.0_wp),2,status)
        call expect(status,cf_ok)
        electric = 0.0_wp
        electric(1,:) = 3.0_wp
        call material%eval(electric,current_density,status=status)
        call expect(status,cf_ok)
        call close(current_density(1,1)/1.0e-12_wp,3.0_wp)
        call material%init(ohmic_anisotropic_t(10.0_wp,8.0_wp,6.0_wp, &
                           2.0_wp,-1.0_wp,1.5_wp),2,status)
        call expect(status,cf_ok)
        electric(1,2) = ieee_value(0.0_wp,ieee_positive_inf)
        call material%eval(electric,current_density,status=status)
        call expect(status,cf_evaluation_failed)
        if (material%last_failure_point() /= 2) error stop "Bad ohmic point"
        call material%commit(status)
        call expect(status,cf_no_trial)
        electric(1,2) = ieee_value(0.0_wp,ieee_quiet_nan)
        call material%eval(electric,current_density,status=status)
        call expect(status,cf_evaluation_failed)
        electric(1,2) = 2.0_wp
        call material%eval(electric,current_density,status=status)
        call expect(status,cf_ok)
        call material%eval(wrong,current_density,status=status)
        call expect(status,cf_invalid_shape)
        call material%eval(electric,current_density,status=status,dsdp=extra)
        call expect(status,cf_unsupported)

        electric(1,2) = huge(1.0_wp)/2.0_wp
        call material%eval(electric,current_density,status=status)
        call expect(status,cf_evaluation_failed)
        if (material%last_failure_point() /= 2) error stop "Bad overflowing ohmic point"
        call material%commit(status)
        call expect(status,cf_no_trial)
        electric(1,2) = 2.0_wp
        call material%eval(electric,current_density,status=status)
        call expect(status,cf_ok)
    end subroutine test_invalid_and_recovery

end program test_ohmic_anisotropic
