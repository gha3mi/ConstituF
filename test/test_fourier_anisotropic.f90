! SPDX-License-Identifier: Apache-2.0

program test_fourier_anisotropic
    !! Test convention: absolute temperature in K, length in m and time in s.
    !! Conductivity is W/(m*K), temperature gradient K/m and heat flux W/m2.
    !!
    !! Independent tensor, covariance, tangent, sensitivity and API tests.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
    use constituf, only: wp, material_t, cf_ok, cf_invalid_argument, &
                        cf_invalid_shape, cf_evaluation_failed, cf_no_trial, cf_unsupported
    use constituf_fourier_anisotropic, only: fourier_anisotropic_t
    use constituf_fourier_isotropic, only: fourier_isotropic_t
    implicit none (type, external)

    integer, parameter :: ndim = 3, npoints = 3, nparameters = 6

    call test_schema_and_flux()
    call test_tangent_and_batch()
    call test_covariance_and_isotropic_limit()
    call test_sensitivity()
    call test_scaled_conductivity(sqrt(tiny(1.0_wp)))
    call test_scaled_conductivity(sqrt(huge(1.0_wp)))
    call test_cancelling_large_products()
    call test_cancellation_with_small_remainder()
    call test_invalid_and_recovery()

contains

    function model_from_values(values) result(model)
        real(wp), intent(in) :: values(nparameters)
        type(fourier_anisotropic_t) :: model

        model = fourier_anisotropic_t(values(1),values(2),values(3), &
                                     values(4),values(5),values(6))
    end function model_from_values

    subroutine expect_status(value, expected)
        integer, intent(in) :: value, expected

        if (value /= expected) then
            print *, "Expected status", expected, "got", value
            error stop "Anisotropic Fourier status mismatch"
        end if
    end subroutine expect_status

    subroutine expect_close(value, expected)
        real(wp), intent(in) :: value, expected

        if (.not. ieee_is_finite(value) .or. .not. ieee_is_finite(expected)) then
            error stop "Nonfinite anisotropic Fourier comparison"
        end if
        if (abs(value-expected) > 256.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(expected))) then
            print *, "Expected", expected, "got", value
            error stop "Anisotropic Fourier numerical mismatch"
        end if
    end subroutine expect_close

    subroutine tensor_reference(conductivity)
        real(wp), intent(out) :: conductivity(ndim,ndim)

        conductivity(:,1) = [10.0_wp,2.0_wp,-1.0_wp]
        conductivity(:,2) = [2.0_wp,8.0_wp,1.5_wp]
        conductivity(:,3) = [-1.0_wp,1.5_wp,6.0_wp]
    end subroutine tensor_reference

    subroutine test_schema_and_flux()
        type(material_t) :: material
        real(wp) :: conductivity(ndim,ndim)
        real(wp) :: gradient(ndim,npoints), flux(ndim,npoints)
        real(wp) :: tangent(ndim,ndim,npoints)
        integer :: status, ipoint, i, j

        call material%init(model_from_values([10.0_wp,8.0_wp,6.0_wp, &
                                             2.0_wp,-1.0_wp,1.5_wp]), npoints, status)
        call expect_status(status, cf_ok)
        if (material%input_size() /= ndim .or. material%output_size() /= ndim) then
            error stop "Anisotropic thermal field dimension mismatch"
        end if
        if (material%state_size() /= 0 .or. material%state_storage_bytes() /= 0) then
            error stop "Anisotropic Fourier material must be stateless"
        end if
        if (.not. material%field_matches("input", "temperature.gradient", &
            "vector", "temperature_gradient", "spatial", [ndim])) then
            error stop "Anisotropic temperature-gradient schema mismatch"
        end if
        if (.not. material%field_matches("output", "heat.flux", &
            "vector", "heat_flux", "spatial", [ndim])) then
            error stop "Anisotropic heat-flux schema mismatch"
        end if
        if (material%field_matches("output", "heat.flux", &
            "vector", "heat_flux", "reference", [ndim])) then
            error stop "Incompatible thermal frame accepted"
        end if
        if (material%field_matches("input", "temperature.gradient", &
            "vector", "pressure_gradient", "spatial", [ndim])) then
            error stop "Incompatible temperature gradient measure accepted"
        end if

        call tensor_reference(conductivity)
        gradient(:,1) = [1.0_wp,-2.0_wp,0.5_wp]
        gradient(:,2) = [-3.0_wp,4.0_wp,1.25_wp]
        gradient(:,3) = 0.0_wp
        call material%eval(gradient, flux, tangent, status)
        call expect_status(status, cf_ok)
        call expect_close(flux(1,1), -5.5_wp)
        call expect_close(flux(2,1), 13.25_wp)
        call expect_close(flux(3,1), 1.0_wp)
        do ipoint = 1, npoints
            do i = 1, ndim
                call expect_close(flux(i,ipoint), &
                                  -dot_product(conductivity(i,:),gradient(:,ipoint)))
                do j = 1, ndim
                    call expect_close(tangent(i,j,ipoint), -conductivity(i,j))
                end do
            end do
            if (dot_product(flux(:,ipoint),gradient(:,ipoint)) > &
                128.0_wp*epsilon(1.0_wp)) then
                error stop "Positive heat flux along temperature gradient"
            end if
        end do
        call material%commit(status)
        call expect_status(status, cf_ok)
        call material%eval(gradient, flux, status=status)
        call expect_status(status, cf_ok)
        call material%revert()
        call material%commit(status)
        call expect_status(status, cf_no_trial)
    end subroutine test_schema_and_flux

    subroutine test_tangent_and_batch()
        type(material_t) :: batch, single
        real(wp) :: gradient(ndim,2), flux(ndim,2), tangent(ndim,ndim,2)
        real(wp) :: changed(ndim,2), positive(ndim,2), negative(ndim,2)
        real(wp) :: alone(ndim,1)
        real(wp) :: h, derivative, error_bound
        integer :: status, i, j, ipoint

        call batch%init(fourier_anisotropic_t(10.0_wp,8.0_wp,6.0_wp, &
                        2.0_wp,-1.0_wp,1.5_wp), 2, status)
        call expect_status(status, cf_ok)
        call single%init(fourier_anisotropic_t(10.0_wp,8.0_wp,6.0_wp, &
                         2.0_wp,-1.0_wp,1.5_wp), 1, status)
        call expect_status(status, cf_ok)
        gradient(:,1) = [1.2_wp,-2.3_wp,0.7_wp]
        gradient(:,2) = [0.1_wp,3.4_wp,-1.1_wp]
        call batch%eval(gradient, flux, tangent, status)
        call expect_status(status, cf_ok)
        if (.not. all(ieee_is_finite(tangent))) error stop "Nonfinite anisotropic Fourier analytical tangent"
        call single%eval(gradient(:,1:1), alone, status=status)
        call expect_status(status, cf_ok)
        do i = 1, ndim
            call expect_close(flux(i,1), alone(i,1))
        end do

        do ipoint = 1, 2
            do j = 1, ndim
                ! Flux is affine in each gradient component, so a finite
                ! central step avoids cancellation without truncation error.
                h = 0.25_wp*maxval(abs(gradient(:,ipoint)))
                changed = gradient
                changed(j,ipoint) = changed(j,ipoint)+h
                call batch%eval(changed, positive, status=status)
                call expect_status(status, cf_ok)
                changed(j,ipoint) = gradient(j,ipoint)-h
                call batch%eval(changed, negative, status=status)
                call expect_status(status, cf_ok)
                if (.not. all(ieee_is_finite(positive)) .or. .not. all(ieee_is_finite(negative))) then
                    error stop "Nonfinite anisotropic Fourier tangent response"
                end if
                error_bound = 32.0_wp*epsilon(1.0_wp)* &
                    maxval(abs(positive(:,ipoint))+abs(negative(:,ipoint)))/(2.0_wp*h)
                if (.not. ieee_is_finite(error_bound)) error stop "Nonfinite anisotropic Fourier tangent bound"
                do i = 1, ndim
                    derivative = (positive(i,ipoint)-negative(i,ipoint))/(2.0_wp*h)
                    if (.not. ieee_is_finite(derivative) .or. abs(derivative-tangent(i,j,ipoint)) > error_bound) then
                        error stop "Anisotropic Fourier tangent finite difference mismatch"
                    end if
                end do
            end do
        end do
    end subroutine test_tangent_and_batch

    subroutine test_covariance_and_isotropic_limit()
        type(material_t) :: original, rotated, isotropic, diagonal
        real(wp) :: conductivity(ndim,ndim), new_tensor(ndim,ndim)
        real(wp) :: rotation(ndim,ndim), gradient(ndim,1), new_gradient(ndim,1)
        real(wp) :: flux(ndim,1), new_flux(ndim,1), iso_flux(ndim,1), diag_flux(ndim,1)
        integer :: status, i

        call tensor_reference(conductivity)
        rotation(:,1) = [0.0_wp,1.0_wp,0.0_wp]
        rotation(:,2) = [-1.0_wp,0.0_wp,0.0_wp]
        rotation(:,3) = [0.0_wp,0.0_wp,1.0_wp]
        new_tensor = matmul(rotation,matmul(conductivity,transpose(rotation)))
        gradient(:,1) = [2.0_wp,-1.0_wp,0.5_wp]
        new_gradient(:,1) = matmul(rotation,gradient(:,1))

        call original%init(fourier_anisotropic_t(10.0_wp,8.0_wp,6.0_wp, &
                           2.0_wp,-1.0_wp,1.5_wp), 1, status)
        call expect_status(status, cf_ok)
        call rotated%init(fourier_anisotropic_t(new_tensor(1,1),new_tensor(2,2), &
                          new_tensor(3,3),new_tensor(1,2),new_tensor(1,3), &
                          new_tensor(2,3)), 1, status)
        call expect_status(status, cf_ok)
        call original%eval(gradient, flux, status=status)
        call expect_status(status, cf_ok)
        call rotated%eval(new_gradient, new_flux, status=status)
        call expect_status(status, cf_ok)
        do i = 1, ndim
            call expect_close(new_flux(i,1),dot_product(rotation(i,:),flux(:,1)))
        end do

        call isotropic%init(fourier_isotropic_t(7.0_wp), 1, status)
        call expect_status(status, cf_ok)
        call diagonal%init(fourier_anisotropic_t(7.0_wp,7.0_wp,7.0_wp, &
                           0.0_wp,0.0_wp,0.0_wp), 1, status)
        call expect_status(status, cf_ok)
        call isotropic%eval(gradient, iso_flux, status=status)
        call expect_status(status, cf_ok)
        call diagonal%eval(gradient, diag_flux, status=status)
        call expect_status(status, cf_ok)
        do i = 1, ndim
            call expect_close(iso_flux(i,1),diag_flux(i,1))
        end do
    end subroutine test_covariance_and_isotropic_limit

    subroutine test_sensitivity()
        type(material_t) :: tracked, positive, negative
        real(wp), parameter :: expected_derivative(ndim,nparameters) = reshape( &
            [0.0_wp,-2.1_wp,0.9_wp, -2.1_wp,0.0_wp,-1.3_wp, &
             0.9_wp,-1.3_wp,0.0_wp, 0.0_wp,0.0_wp,-2.1_wp, &
             0.0_wp,0.9_wp,0.0_wp, -1.3_wp,0.0_wp,0.0_wp],[ndim,nparameters])
        real(wp) :: values(nparameters), plus_values(nparameters), minus_values(nparameters)
        real(wp) :: gradient(ndim,2), flux(ndim,2), dsdp(ndim,nparameters,2)
        real(wp) :: tangent(ndim,ndim,2), plus_flux(ndim,2), minus_flux(ndim,2)
        real(wp) :: malformed(ndim,nparameters-1,2)
        real(wp) :: numerical, h, error_bound
        integer, parameter :: parameter_ids(nparameters) = [6,5,4,3,2,1]
        integer :: status, iparam, ipoint, i, id

        values = [10.0_wp,8.0_wp,6.0_wp,2.0_wp,-1.0_wp,1.5_wp]
        call tracked%init(model_from_values(values), 2, status, &
            sensitivity_parameters=[character(len=15) :: "conductivity_23", &
                                    "conductivity_13", "conductivity_12", &
                                    "conductivity_33", "conductivity_22", &
                                    "conductivity_11"])
        call expect_status(status, cf_ok)
        if (tracked%sensitivity_size() /= nparameters) then
            error stop "Missing anisotropic conductivity sensitivities"
        end if
        if (tracked%state_storage_bytes() /= 0) then
            error stop "Stateless anisotropic sensitivity allocated history"
        end if
        gradient(:,1) = [1.3_wp,-0.9_wp,2.1_wp]
        gradient(:,2) = [-1.2_wp,0.75_wp,0.4_wp]
        call tracked%eval(gradient, flux, tangent, status, dsdp=dsdp)
        call expect_status(status, cf_ok)
        if (.not. all(ieee_is_finite(dsdp))) error stop "Nonfinite anisotropic Fourier analytical sensitivity"

        do iparam = 1, nparameters
            id = parameter_ids(iparam)
            ! Flux is affine in every independent tensor coefficient; these
            ! finite perturbations preserve positive definiteness.
            h = 0.25_wp*abs(values(id))
            plus_values = values
            minus_values = values
            plus_values(id) = plus_values(id)+h
            minus_values(id) = minus_values(id)-h
            call positive%init(model_from_values(plus_values), 2, status)
            call expect_status(status, cf_ok)
            call negative%init(model_from_values(minus_values), 2, status)
            call expect_status(status, cf_ok)
            call positive%eval(gradient, plus_flux, status=status)
            call expect_status(status, cf_ok)
            call negative%eval(gradient, minus_flux, status=status)
            call expect_status(status, cf_ok)
            if (.not. all(ieee_is_finite(plus_flux)) .or. .not. all(ieee_is_finite(minus_flux))) then
                error stop "Nonfinite anisotropic Fourier parameter response"
            end if
            do ipoint = 1, 2
                error_bound = 32.0_wp*epsilon(1.0_wp)* &
                    maxval(abs(plus_flux(:,ipoint))+abs(minus_flux(:,ipoint)))/(2.0_wp*h)
                if (.not. ieee_is_finite(error_bound)) error stop "Nonfinite anisotropic Fourier parameter bound"
                do i = 1, ndim
                    numerical = (plus_flux(i,ipoint)-minus_flux(i,ipoint))/(2.0_wp*h)
                    if (.not. ieee_is_finite(numerical) .or. abs(numerical-dsdp(i,iparam,ipoint)) > error_bound) then
                        error stop "Anisotropic Fourier parameter derivative mismatch"
                    end if
                end do
            end do
        end do
        do iparam = 1,nparameters
            do i = 1,ndim
                call expect_close(dsdp(i,iparam,1),expected_derivative(i,iparam))
            end do
        end do

        call tracked%commit(status)
        call expect_status(status, cf_ok)
        call tracked%eval(gradient, flux, status=status)
        call expect_status(status, cf_invalid_argument)
        call tracked%eval(gradient, flux, status=status, dsdp=malformed)
        call expect_status(status, cf_invalid_shape)
        call tracked%eval(gradient, flux, status=status, dsdp=dsdp)
        call expect_status(status, cf_ok)
        call tracked%revert()
        call tracked%commit(status)
        call expect_status(status, cf_no_trial)
    end subroutine test_sensitivity

    subroutine test_scaled_conductivity(scale)
        !! Common coefficient scaling preserves definiteness. These scales
        !! would underflow/overflow raw cubic determinant products.
        real(wp), intent(in) :: scale
        type(material_t) :: material
        real(wp) :: conductivity(ndim,ndim), gradient(ndim,1), flux(ndim,1)
        real(wp) :: tangent(ndim,ndim,1), dsdp(ndim,1,1)
        integer :: status, i, j

        call material%init(model_from_values(scale*[10.0_wp,8.0_wp,6.0_wp, &
            2.0_wp,-1.0_wp,1.5_wp]),1,status, &
            sensitivity_parameters=[character(len=15) :: "conductivity_13"])
        call expect_status(status,cf_ok)
        call tensor_reference(conductivity)
        gradient(:,1)=[1.0_wp,-2.0_wp,0.5_wp]
        call material%eval(gradient,flux,tangent,status,dsdp=dsdp)
        call expect_status(status,cf_ok)
        call expect_close(flux(1,1)/scale,-5.5_wp)
        call expect_close(flux(2,1)/scale,13.25_wp)
        call expect_close(flux(3,1)/scale,1.0_wp)
        call expect_close(-dot_product(flux(:,1)/scale,gradient(:,1)),31.5_wp)
        do j=1,ndim
            do i=1,ndim
                call expect_close(tangent(i,j,1)/scale,-conductivity(i,j))
            end do
        end do
        call expect_close(dsdp(1,1,1),-0.5_wp)
        call expect_close(dsdp(2,1,1),0.0_wp)
        call expect_close(dsdp(3,1,1),-1.0_wp)
    end subroutine test_scaled_conductivity

    subroutine test_cancelling_large_products()
        !! This SPD tensor has condition number five. Individual products
        !! overflow, although their sums and every requested result are finite.
        type(material_t) :: material
        real(wp), parameter :: tensor(ndim,ndim) = reshape( &
            [3.0_wp,2.0_wp,0.0_wp,2.0_wp,3.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp],[ndim,ndim])
        real(wp), parameter :: expected_derivative(ndim,nparameters) = reshape( &
            [0.0_wp,0.0_wp,2.0_wp, 0.0_wp,0.0_wp,-2.0_wp, &
             2.0_wp,-2.0_wp,0.0_wp, 0.0_wp,0.0_wp,0.0_wp, &
             0.0_wp,2.0_wp,0.0_wp, -2.0_wp,0.0_wp,0.0_wp],[ndim,nparameters])
        real(wp) :: gradient(ndim,1), flux(ndim,1), tangent(ndim,ndim,1)
        real(wp) :: dsdp(ndim,nparameters,1)
        real(wp) :: scale
        integer :: status, i, j

        scale=huge(1.0_wp)/4.0_wp
        call material%init(fourier_anisotropic_t(3.0_wp*scale,3.0_wp*scale,scale, &
            2.0_wp*scale,0.0_wp,0.0_wp),1,status, &
            sensitivity_parameters=[character(len=15) :: "conductivity_23", &
                "conductivity_13","conductivity_12","conductivity_33", &
                "conductivity_22","conductivity_11"])
        call expect_status(status,cf_ok)
        gradient(:,1)=[2.0_wp,-2.0_wp,0.0_wp]
        call material%eval(gradient,flux,tangent,status,dsdp=dsdp)
        call expect_status(status,cf_ok)
        call expect_close(flux(1,1)/scale,-2.0_wp)
        call expect_close(flux(2,1)/scale,2.0_wp)
        call expect_close(flux(3,1)/scale,0.0_wp)
        do j=1,ndim
            do i=1,ndim
                call expect_close(tangent(i,j,1)/scale,-tensor(i,j))
            end do
        end do
        do j=1,nparameters
            do i=1,ndim
                call expect_close(dsdp(i,j,1),expected_derivative(i,j))
            end do
        end do
    end subroutine test_cancelling_large_products

    subroutine test_cancellation_with_small_remainder()
        !! Equal overflowing terms cancel exactly before a tiny third
        !! term remains. Common scaling must not erase that finite flux.
        type(material_t) :: material
        integer, parameter :: npermutations = 6
        integer, parameter :: permutations(ndim,npermutations) = reshape( &
            [1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1],[ndim,npermutations])
        real(wp) :: conductivity(ndim,ndim), permuted(ndim,ndim)
        real(wp) :: gradient(ndim,1), flux(ndim,1), tangent(ndim,ndim,1)
        real(wp) :: original_gradient(ndim), expected(ndim), dsdp(ndim,1,1)
        real(wp) :: large, small
        integer :: status, ipermutation, i

        large = scale(0.75_wp,maxexponent(1.0_wp))
        small = tiny(1.0_wp)
        conductivity(:,1) = [large,0.5_wp*large,small]
        conductivity(:,2) = [0.5_wp*large,0.375_wp*large,0.0_wp]
        conductivity(:,3) = [small,0.0_wp,large]
        original_gradient = [2.0_wp,-4.0_wp,1.0_wp]
        expected = [-small,0.5_wp*large,-large]

        do ipermutation = 1, size(permutations,2)
            permuted = conductivity(permutations(:,ipermutation),permutations(:,ipermutation))
            gradient(:,1) = original_gradient(permutations(:,ipermutation))
            call material%init(fourier_anisotropic_t(permuted(1,1),permuted(2,2), &
                permuted(3,3),permuted(1,2),permuted(1,3),permuted(2,3)),1,status, &
                sensitivity_parameters=[character(len=15) :: "conductivity_13"])
            call expect_status(status,cf_ok)
            call material%eval(gradient,flux,tangent,status,dsdp=dsdp)
            call expect_status(status,cf_ok)
            if (.not. all(ieee_is_finite(tangent))) error stop "Nonfinite cancellation tangent"
            do i = 1, ndim
                ! Normalize each component separately so the large fluxes
                ! cannot mask loss of the tiny exact remainder.
                call expect_close(flux(i,1)/expected(permutations(i,ipermutation)),1.0_wp)
            end do
            call expect_close(dsdp(1,1,1),-gradient(3,1))
            call expect_close(dsdp(2,1,1),0.0_wp)
            call expect_close(dsdp(3,1,1),-gradient(1,1))
        end do
    end subroutine test_cancellation_with_small_remainder

    subroutine test_invalid_and_recovery()
        type(material_t) :: material
        type(fourier_anisotropic_t) :: invalid
        real(wp) :: gradient(ndim,2), flux(ndim,2), malformed(ndim-1,2)
        real(wp) :: untracked_dsdp(ndim,1,2)
        integer :: status

        invalid = fourier_anisotropic_t()
        call material%init(invalid, 2, status)
        call expect_status(status, cf_invalid_argument)
        invalid = fourier_anisotropic_t(1.0_wp,1.0_wp,1.0_wp, &
                                       1.0_wp,0.0_wp,0.0_wp)
        call material%init(invalid, 2, status)
        call expect_status(status, cf_invalid_argument)
        invalid = fourier_anisotropic_t(10.0_wp,8.0_wp,6.0_wp, &
                                       10.0_wp,-1.0_wp,1.5_wp)
        call material%init(invalid, 2, status)
        call expect_status(status, cf_invalid_argument)
        invalid = fourier_anisotropic_t(10.0_wp,8.0_wp,6.0_wp, &
                                       2.0_wp,-1.0_wp,1.5_wp)
        invalid%conductivity_33 = -6.0_wp
        call material%init(invalid, 2, status)
        call expect_status(status, cf_invalid_argument)
        invalid%conductivity_33 = ieee_value(0.0_wp,ieee_quiet_nan)
        call material%init(invalid, 2, status)
        call expect_status(status, cf_invalid_argument)
        invalid%conductivity_33 = ieee_value(0.0_wp,ieee_positive_inf)
        call material%init(invalid, 2, status)
        call expect_status(status, cf_invalid_argument)

        call material%init(fourier_anisotropic_t(10.0_wp,8.0_wp,6.0_wp, &
                           2.0_wp,-1.0_wp,1.5_wp), 2, status)
        call expect_status(status, cf_ok)
        gradient(:,1) = [1.0_wp,-2.0_wp,0.5_wp]
        gradient(:,2) = [0.5_wp,1.0_wp,-0.5_wp]
        gradient(3,2) = ieee_value(0.0_wp,ieee_quiet_nan)
        call material%eval(gradient, flux, status=status)
        call expect_status(status, cf_evaluation_failed)
        if (material%last_failure_point() /= 2) error stop "Wrong failed thermal point"
        call material%commit(status)
        call expect_status(status, cf_no_trial)

        gradient(3,2) = -0.5_wp
        call material%eval(gradient, flux, status=status)
        call expect_status(status, cf_ok)
        call expect_close(flux(1,1), -5.5_wp)
        malformed = 0.0_wp
        call material%eval(malformed, flux, status=status)
        call expect_status(status, cf_invalid_shape)
        call material%init(fourier_anisotropic_t(), 2, status)
        call expect_status(status, cf_invalid_argument)
        call material%eval(gradient, flux, status=status)
        call expect_status(status, cf_ok)
        call expect_close(flux(1,1), -5.5_wp)

        gradient(1,2) = huge(1.0_wp)/2.0_wp
        call material%eval(gradient, flux, status=status)
        call expect_status(status, cf_evaluation_failed)
        call material%commit(status)
        call expect_status(status, cf_no_trial)
        gradient(1,2) = 0.5_wp
        call material%eval(gradient, flux, status=status)
        call expect_status(status, cf_ok)
        call material%eval(gradient, flux, status=status, dsdp=untracked_dsdp)
        call expect_status(status, cf_unsupported)
        call material%commit(status)
        call expect_status(status, cf_no_trial)
    end subroutine test_invalid_and_recovery

end program test_fourier_anisotropic
