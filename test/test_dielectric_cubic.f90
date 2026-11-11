! SPDX-License-Identifier: Apache-2.0

program test_dielectric_cubic
    !! Test convention: electric field in V/m, displacement in C/m2 and
    !! permittivity in F/m. Derived energy density is J/m3.
    !!
    !! Cubic-field energy, exact Jacobian, parameter and invalid-trial tests.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
    use constituf, only: wp, material_t, cf_ok, cf_invalid_argument, &
                        cf_invalid_shape, cf_evaluation_failed, cf_no_trial, cf_unsupported
    use constituf_dielectric_cubic, only: dielectric_cubic_t
    use constituf_dielectric_isotropic, only: dielectric_isotropic_t
    implicit none (type, external)

    integer, parameter :: ndim = 3

    call test_reference_and_schema()
    call test_tangent_and_rotation()
    call test_sensitivity()
    call test_large_linear_field()
    call test_large_field_norm()
    call test_large_cubic_coefficient()
    call test_small_field_norm()
    call test_subnormal_coefficients()
    call test_requested_output_range()
    call test_invalid_and_recovery()

contains

    subroutine expect(actual, wanted)
        integer, intent(in) :: actual, wanted

        if (actual /= wanted) then
            print *, "actual status", actual, "wanted", wanted
            error stop "Cubic dielectric status mismatch"
        end if
    end subroutine expect

    subroutine close(actual, wanted)
        real(wp), intent(in) :: actual, wanted

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(wanted)) then
            error stop "Nonfinite cubic dielectric comparison"
        end if
        if (abs(actual-wanted) > 1024.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(wanted))) then
            print *, "actual", actual, "wanted", wanted
            error stop "Cubic dielectric reference mismatch"
        end if
    end subroutine close

    subroutine test_reference_and_schema()
        type(material_t) :: material, limiting
        real(wp), parameter :: exact_displacement(ndim) = [101.0_wp/10.0_wp,-101.0_wp/20.0_wp,101.0_wp/40.0_wp]
        real(wp), parameter :: exact_tangent(ndim,ndim) = reshape( &
            [133.0_wp/20.0_wp,-4.0_wp/5.0_wp,2.0_wp/5.0_wp, &
             -4.0_wp/5.0_wp,109.0_wp/20.0_wp,-1.0_wp/5.0_wp, &
             2.0_wp/5.0_wp,-1.0_wp/5.0_wp,103.0_wp/20.0_wp],[ndim,ndim])
        real(wp) :: electric(ndim,2), displacement(ndim,2), tangent(ndim,ndim,2)
        real(wp) :: baseline(ndim,2), half_displacement(ndim,2)
        real(wp) :: norm_squared, effective, expected, coenergy, stored_energy
        integer :: status, ipoint, i, j

        call material%init(dielectric_cubic_t(4.0_wp,0.2_wp),2,status)
        call expect(status,cf_ok)
        if (material%input_size() /= ndim .or. material%output_size() /= ndim) then
            error stop "Cubic dielectric field sizes"
        end if
        if (material%state_storage_bytes() /= 0) error stop "Cubic dielectric history"
        if (.not. material%field_matches("input", "electric.field", &
            "vector", "electric_field", "spatial", [ndim])) then
            error stop "Cubic dielectric input schema"
        end if
        if (.not. material%field_matches("output", "electric.displacement", &
            "vector", "electric_displacement", "spatial", [ndim])) then
            error stop "Cubic dielectric output schema"
        end if
        if (material%field_matches("input", "electric.field", &
            "vector", "electric_displacement", "spatial", [ndim])) then
            error stop "Cubic dielectric wrong field measure accepted"
        end if

        electric(:,1) = [2.0_wp,-1.0_wp,0.5_wp]
        electric(:,2) = 0.0_wp
        call material%eval(electric,displacement,tangent,status)
        call expect(status,cf_ok)
        do ipoint = 1, 2
            norm_squared = dot_product(electric(:,ipoint),electric(:,ipoint))
            effective = 4.0_wp+0.2_wp*norm_squared
            coenergy = 2.0_wp*norm_squared+0.05_wp*norm_squared**2
            stored_energy = dot_product(electric(:,ipoint),displacement(:,ipoint)) - &
                            coenergy
            if (stored_energy < 0.0_wp) error stop "Cubic energy must be nonnegative"
            call close(stored_energy,2.0_wp*norm_squared+0.15_wp*norm_squared**2)
            do i = 1, ndim
                call close(displacement(i,ipoint),effective*electric(i,ipoint))
                do j = 1, ndim
                    expected = 0.4_wp*electric(i,ipoint)*electric(j,ipoint)
                    if (i == j) expected = expected+effective
                    call close(tangent(i,j,ipoint),expected)
                end do
            end do
        end do
        do i=1,ndim
            call close(displacement(i,1),exact_displacement(i))
            do j=1,ndim
                call close(tangent(i,j,1),exact_tangent(i,j))
            end do
        end do
        call close(tangent(1,1,2),4.0_wp)
        ! Simpson integration is exact along this proportional path because
        ! displacement is cubic in the field. This independently checks the
        ! coenergy integral and its distinct stored-energy Legendre transform.
        call material%eval(0.5_wp*electric,half_displacement,status=status)
        call expect(status,cf_ok)
        coenergy=dot_product(electric(:,1),displacement(:,1)+4.0_wp*half_displacement(:,1))/6.0_wp
        stored_energy=dot_product(electric(:,1),displacement(:,1))-coenergy
        call close(coenergy,3801.0_wp/320.0_wp)
        call close(stored_energy,4683.0_wp/320.0_wp)
        call material%commit(status)
        call expect(status,cf_ok)

        call material%init(dielectric_cubic_t(4.0_wp,0.0_wp),2,status)
        call expect(status,cf_ok)
        call limiting%init(dielectric_isotropic_t(4.0_wp),2,status)
        call expect(status,cf_ok)
        call material%eval(electric,displacement,tangent,status)
        call expect(status,cf_ok)
        call limiting%eval(electric,baseline,status=status)
        call expect(status,cf_ok)
        do ipoint = 1, 2
            do i = 1, ndim
                call close(displacement(i,ipoint),baseline(i,ipoint))
            end do
        end do
        call close(tangent(1,1,1),4.0_wp)
        call close(tangent(1,2,1),0.0_wp)
    end subroutine test_reference_and_schema

    subroutine test_tangent_and_rotation()
        type(material_t) :: batch, single
        real(wp) :: x(ndim,2), displacement(ndim,2), tangent(ndim,ndim,2)
        real(wp) :: changed(ndim,2), plus(ndim,2), minus(ndim,2), alone(ndim,1)
        real(wp) :: rotation(ndim,ndim), rotated_tangent(ndim,ndim)
        real(wp) :: h, derivative, error_bound
        integer :: status, i, j, ipoint

        call batch%init(dielectric_cubic_t(3.0_wp,0.25_wp),2,status)
        call expect(status,cf_ok)
        call single%init(dielectric_cubic_t(3.0_wp,0.25_wp),1,status)
        call expect(status,cf_ok)
        rotation(:,1) = [0.0_wp,1.0_wp,0.0_wp]
        rotation(:,2) = [-1.0_wp,0.0_wp,0.0_wp]
        rotation(:,3) = [0.0_wp,0.0_wp,1.0_wp]
        x(:,1) = [1.2_wp,-0.7_wp,2.0_wp]
        x(:,2) = matmul(rotation,x(:,1))
        call batch%eval(x,displacement,tangent,status)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(tangent))) error stop "Nonfinite cubic dielectric tangent"
        call single%eval(x(:,1:1),alone,status=status)
        call expect(status,cf_ok)
        do i = 1, ndim
            call close(displacement(i,1),alone(i,1))
            call close(displacement(i,2), &
                       dot_product(rotation(i,:),displacement(:,1)))
        end do
        rotated_tangent = matmul(rotation,matmul(tangent(:,:,1),transpose(rotation)))
        do j = 1, ndim
            do i = 1, ndim
                call close(tangent(i,j,2),rotated_tangent(i,j))
            end do
        end do

        do ipoint = 1, 2
            do j = 1, ndim
                h = epsilon(1.0_wp)**(1.0_wp/3.0_wp)*maxval(abs(x(:,ipoint)))
                changed = x
                changed(j,ipoint) = x(j,ipoint)+h
                call batch%eval(changed,plus,status=status)
                call expect(status,cf_ok)
                changed(j,ipoint) = x(j,ipoint)-h
                call batch%eval(changed,minus,status=status)
                call expect(status,cf_ok)
                if (.not. all(ieee_is_finite(plus)) .or. .not. all(ieee_is_finite(minus))) then
                    error stop "Nonfinite cubic dielectric tangent secant response"
                end if
                do i = 1, ndim
                    derivative = (plus(i,ipoint)-minus(i,ipoint))/(2.0_wp*h)
                    ! The exact central-difference truncation is alpha*h**2
                    ! on the perturbed diagonal and zero off the diagonal.
                    error_bound = 32.0_wp*epsilon(1.0_wp)* &
                        maxval(abs(plus(:,ipoint))+abs(minus(:,ipoint)))/(2.0_wp*h)+ &
                        merge(0.25_wp*h*h,0.0_wp,i==j)
                    if (.not. all(ieee_is_finite([derivative,error_bound]))) then
                        error stop "Nonfinite cubic dielectric tangent comparison"
                    end if
                    if (abs(derivative-tangent(i,j,ipoint)) > error_bound) then
                        error stop "Cubic dielectric tangent finite difference"
                    end if
                    call close(tangent(i,j,ipoint),tangent(j,i,ipoint))
                end do
            end do
        end do
        call batch%revert()
        call batch%commit(status)
        call expect(status,cf_no_trial)
        call batch%eval(x,displacement,status=status)
        call expect(status,cf_ok)
        call batch%commit(status)
        call expect(status,cf_ok)
    end subroutine test_tangent_and_rotation

    subroutine test_sensitivity()
        type(material_t) :: tracked, upper, lower
        real(wp), parameter :: cubic_derivative(ndim) = [21.0_wp/2.0_wp,-21.0_wp/4.0_wp,21.0_wp/8.0_wp]
        real(wp) :: electric(ndim,1), displacement(ndim,1), dsdp(ndim,2,1)
        real(wp) :: positive(ndim,1), negative(ndim,1), wrong(ndim,1,1)
        real(wp) :: derivative, h, error_bound
        integer :: status, iparam, i

        call tracked%init(dielectric_cubic_t(3.0_wp,0.3_wp),1,status, &
            sensitivity_parameters=[character(len=17) :: &
            "cubic_coefficient","permittivity"])
        call expect(status,cf_ok)
        if (tracked%sensitivity_size() /= 2) error stop "Cubic sensitivity selection"
        if (tracked%state_storage_bytes() /= 0) error stop "Cubic sensitivity history"
        electric(:,1) = [2.0_wp,-1.0_wp,0.5_wp]
        call tracked%eval(electric,displacement,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
        do i = 1, ndim
            call close(dsdp(i,1,1),cubic_derivative(i))
            call close(dsdp(i,2,1),electric(i,1))
        end do

        do iparam = 1, 2
            ! Both coefficients enter displacement affinely, so finite
            ! perturbations avoid cancellation without truncation error.
            select case (iparam)
            case (1)
                h=0.25_wp*0.3_wp
                call upper%init(dielectric_cubic_t(3.0_wp,0.3_wp+h),1,status)
                call expect(status,cf_ok)
                call lower%init(dielectric_cubic_t(3.0_wp,0.3_wp-h),1,status)
            case (2)
                h=0.25_wp*3.0_wp
                call upper%init(dielectric_cubic_t(3.0_wp+h,0.3_wp),1,status)
                call expect(status,cf_ok)
                call lower%init(dielectric_cubic_t(3.0_wp-h,0.3_wp),1,status)
            case default
                error stop "Invalid cubic sensitivity index"
            end select
            call expect(status,cf_ok)
            call upper%eval(electric,positive,status=status)
            call expect(status,cf_ok)
            call lower%eval(electric,negative,status=status)
            call expect(status,cf_ok)
            if (.not. all(ieee_is_finite(positive)) .or. .not. all(ieee_is_finite(negative))) then
                error stop "Nonfinite cubic dielectric parameter secant response"
            end if
            do i = 1, ndim
                derivative = (positive(i,1)-negative(i,1))/(2.0_wp*h)
                error_bound = 32.0_wp*epsilon(1.0_wp)* &
                    (abs(positive(i,1))+abs(negative(i,1)))/(2.0_wp*h)
                if (.not. all(ieee_is_finite([derivative,error_bound]))) then
                    error stop "Nonfinite cubic dielectric parameter comparison"
                end if
                if (abs(derivative-dsdp(i,iparam,1)) > error_bound) then
                    error stop "Cubic dielectric parameter finite difference"
                end if
            end do
        end do
        call tracked%commit(status)
        call expect(status,cf_ok)
        call tracked%eval(electric,displacement,status=status)
        call expect(status,cf_invalid_argument)
        call tracked%eval(electric,displacement,status=status,dsdp=wrong)
        call expect(status,cf_invalid_shape)
        call tracked%eval(electric,displacement,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
        ! The cubic derivative stays nonzero at the admissible alpha=0
        ! boundary, even though the response itself is exactly linear.
        call tracked%init(dielectric_cubic_t(3.0_wp,0.0_wp),1,status, &
            sensitivity_parameters=[character(len=17) :: "cubic_coefficient","permittivity"])
        call expect(status,cf_ok)
        call tracked%eval(electric,displacement,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
        do i = 1, ndim
            call close(displacement(i,1),3.0_wp*electric(i,1))
            call close(dsdp(i,1,1),cubic_derivative(i))
            call close(dsdp(i,2,1),electric(i,1))
        end do
    end subroutine test_sensitivity

    subroutine check_diagonal_tangent(tangent,diagonal)
        real(wp), intent(in) :: tangent(ndim,ndim), diagonal(ndim)
        integer :: i, j

        do j=1,ndim
            do i=1,ndim
                if (i==j) then
                    call close(tangent(i,j),diagonal(i))
                else
                    call close(tangent(i,j),0.0_wp)
                end if
            end do
        end do
    end subroutine check_diagonal_tangent

    subroutine test_large_linear_field()
        !! A zero cubic coefficient must preserve the linear limit even when
        !! the unused squared field norm is outside floating-point range.
        type(material_t) :: material
        real(wp) :: electric(ndim,1), displacement(ndim,1), tangent(ndim,ndim,1)
        real(wp) :: dsdp(ndim,1,1)
        real(wp) :: field_scale
        integer :: status, i

        field_scale=huge(1.0_wp)
        call material%init(dielectric_cubic_t(2.0_wp,0.0_wp),1,status, &
            sensitivity_parameters=[character(len=12) :: "permittivity"])
        call expect(status,cf_ok)
        electric(:,1)=field_scale*[0.25_wp,-0.125_wp,0.0_wp]
        call material%eval(electric,displacement,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        do i=1,ndim
            call close(displacement(i,1)/field_scale,2.0_wp*(electric(i,1)/field_scale))
            call close(dsdp(i,1,1)/field_scale,electric(i,1)/field_scale)
        end do
        call check_diagonal_tangent(tangent(:,:,1),[2.0_wp,2.0_wp,2.0_wp])
    end subroutine test_large_linear_field

    subroutine test_large_field_norm()
        !! The squared norm overflows, but the weighted nonlinear response
        !! and tangent remain finite. Only the cubic-parameter derivative
        !! is genuinely unrepresentable for this prescribed field.
        type(material_t) :: material
        real(wp) :: electric(ndim,1), displacement(ndim,1), tangent(ndim,ndim,1)
        real(wp) :: dsdp(ndim,1,1)
        real(wp) :: field_scale, alpha
        integer :: status

        field_scale=sqrt(huge(1.0_wp))
        alpha=1.0_wp/huge(1.0_wp)
        electric(:,1)=[2.0_wp*field_scale,0.0_wp,0.0_wp]
        call material%init(dielectric_cubic_t(1.0_wp,alpha),1,status, &
            sensitivity_parameters=[character(len=12) :: "permittivity"])
        call expect(status,cf_ok)
        call material%eval(electric,displacement,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(displacement(1,1)/field_scale,10.0_wp)
        call close(displacement(2,1),0.0_wp)
        call close(displacement(3,1),0.0_wp)
        call close(dsdp(1,1,1)/field_scale,2.0_wp)
        call check_diagonal_tangent(tangent(:,:,1),[13.0_wp,5.0_wp,5.0_wp])

        call material%init(dielectric_cubic_t(1.0_wp,alpha),1,status, &
            sensitivity_parameters=[character(len=17) :: "cubic_coefficient"])
        call expect(status,cf_ok)
        call material%eval(electric,displacement,status=status,dsdp=dsdp)
        call expect(status,cf_evaluation_failed)
        if (material%last_failure_point()/=1) error stop "Cubic derivative overflow point"
        call material%commit(status)
        call expect(status,cf_no_trial)
        electric(:,1)=[1.0_wp,0.0_wp,0.0_wp]
        call material%eval(electric,displacement,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(dsdp(1,1,1),1.0_wp)
    end subroutine test_large_field_norm

    subroutine test_large_cubic_coefficient()
        !! Forming twice the coefficient first overflows although all
        !! response, tangent and parameter-derivative entries are finite.
        type(material_t) :: material
        real(wp) :: electric(ndim,1), displacement(ndim,1), tangent(ndim,ndim,1)
        real(wp) :: dsdp(ndim,2,1)
        real(wp) :: alpha
        integer :: status, i

        alpha=huge(1.0_wp)
        call material%init(dielectric_cubic_t(1.0_wp,alpha),1,status, &
            sensitivity_parameters=[character(len=17) :: "cubic_coefficient","permittivity"])
        call expect(status,cf_ok)
        electric(:,1)=[0.25_wp,0.0_wp,0.0_wp]
        call material%eval(electric,displacement,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(displacement(1,1)/alpha,1.0_wp/64.0_wp)
        call close(dsdp(1,1,1),1.0_wp/64.0_wp)
        call close(dsdp(1,2,1),0.25_wp)
        do i=2,ndim
            call close(displacement(i,1),0.0_wp)
            call close(dsdp(i,1,1),0.0_wp)
            call close(dsdp(i,2,1),0.0_wp)
        end do
        call check_diagonal_tangent(tangent(:,:,1)/alpha, &
            [3.0_wp/16.0_wp,1.0_wp/16.0_wp,1.0_wp/16.0_wp])
    end subroutine test_large_cubic_coefficient

    subroutine test_small_field_norm()
        !! Here the raw squared field underflows to zero, while its weighted
        !! cubic contribution and every requested output remain normal.
        type(material_t) :: material
        real(wp) :: electric(ndim,1), displacement(ndim,1), tangent(ndim,ndim,1)
        real(wp) :: dsdp(ndim,1,1)
        real(wp) :: roundoff, field_scale
        integer :: status

        roundoff=epsilon(1.0_wp)
        field_scale=sqrt(roundoff)/(2.0_wp*sqrt(huge(1.0_wp)))
        electric(:,1)=[field_scale,0.0_wp,0.0_wp]
        call material%init(dielectric_cubic_t(roundoff/16.0_wp,huge(1.0_wp)/4.0_wp),1,status, &
            sensitivity_parameters=[character(len=12) :: "permittivity"])
        call expect(status,cf_ok)
        call material%eval(electric,displacement,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(displacement(1,1)/field_scale/roundoff,1.0_wp/8.0_wp)
        call close(displacement(2,1),0.0_wp)
        call close(displacement(3,1),0.0_wp)
        call close(dsdp(1,1,1)/field_scale,1.0_wp)
        call check_diagonal_tangent(tangent(:,:,1)/roundoff, &
            [1.0_wp/4.0_wp,1.0_wp/8.0_wp,1.0_wp/8.0_wp])
    end subroutine test_small_field_norm

    subroutine test_subnormal_coefficients()
        !! Round the complete displacement and diagonal tangent only once.
        !! At E1=9/8, separately rounded linear and cubic terms sum to 2d,
        !! although the exact displacement is 2.548828125d and rounds to 3d.
        type(material_t) :: material
        integer, parameter :: ncases = 3
        real(wp), parameter :: fields(ncases) = [1.5_wp,1.125_wp,0.0_wp]
        real(wp), parameter :: expected_displacement(ncases) = [5.0_wp,3.0_wp,0.0_wp]
        real(wp), parameter :: expected_longitudinal(ncases) = [8.0_wp,5.0_wp,1.0_wp]
        real(wp), parameter :: expected_transverse(ncases) = [3.0_wp,2.0_wp,1.0_wp]
        real(wp) :: electric(ndim,1), displacement(ndim,1), tangent(ndim,ndim,1)
        real(wp) :: dsdp(ndim,2,1)
        real(wp) :: smallest
        integer :: status, icase, i

        smallest = nearest(0.0_wp,1.0_wp)
        call material%init(dielectric_cubic_t(smallest,smallest),1,status, &
            sensitivity_parameters=[character(len=17) :: "cubic_coefficient","permittivity"])
        call expect(status,cf_ok)
        do icase = 1, ncases
            electric(:,1) = [fields(icase),0.0_wp,0.0_wp]
            call material%eval(electric,displacement,tangent,status,dsdp=dsdp)
            call expect(status,cf_ok)
            call close(displacement(1,1)/smallest,expected_displacement(icase))
            call check_diagonal_tangent(tangent(:,:,1)/smallest, &
                [expected_longitudinal(icase),expected_transverse(icase),expected_transverse(icase)])
            call close(dsdp(1,1,1),fields(icase)**3)
            call close(dsdp(1,2,1),fields(icase))
            do i = 2, ndim
                call close(displacement(i,1),0.0_wp)
                call close(dsdp(i,1,1),0.0_wp)
                call close(dsdp(i,2,1),0.0_wp)
            end do
        end do
    end subroutine test_subnormal_coefficients

    subroutine test_requested_output_range()
        !! The intermediate effective coefficient exceeds the real range,
        !! although displacement and both parameter derivatives are finite.
        !! The actual tangent is unrepresentable and must fail when requested.
        type(material_t) :: plain, tracked
        real(wp) :: electric(ndim,1), displacement(ndim,1), tangent(ndim,ndim,1)
        real(wp) :: dsdp(ndim,2,1)
        real(wp) :: coefficient_scale
        integer :: status, i

        coefficient_scale=huge(1.0_wp)
        call plain%init(dielectric_cubic_t(1.0_wp,0.75_wp*coefficient_scale),1,status)
        call expect(status,cf_ok)
        electric=0.75_wp
        call plain%eval(electric,displacement,status=status)
        call expect(status,cf_ok)
        do i=1,ndim
            call close(displacement(i,1)/coefficient_scale,243.0_wp/256.0_wp)
        end do

        call tracked%init(dielectric_cubic_t(1.0_wp,0.75_wp*coefficient_scale),1,status, &
            sensitivity_parameters=[character(len=17) :: "cubic_coefficient","permittivity"])
        call expect(status,cf_ok)
        call tracked%eval(electric,displacement,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
        do i=1,ndim
            call close(displacement(i,1)/coefficient_scale,243.0_wp/256.0_wp)
            call close(dsdp(i,1,1),81.0_wp/64.0_wp)
            call close(dsdp(i,2,1),0.75_wp)
        end do
        call tracked%eval(electric,displacement,tangent,status,dsdp=dsdp)
        call expect(status,cf_evaluation_failed)
        if (tracked%last_failure_point()/=1) error stop "Cubic tangent overflow point"
        call tracked%commit(status)
        call expect(status,cf_no_trial)
        electric=0.0_wp
        call tracked%eval(electric,displacement,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(displacement)) .or. .not. all(ieee_is_finite(dsdp))) then
            error stop "Nonfinite cubic dielectric zero-field recovery"
        end if
        call close(maxval(abs(displacement)),0.0_wp)
        call close(maxval(abs(dsdp)),0.0_wp)
        call check_diagonal_tangent(tangent(:,:,1),[1.0_wp,1.0_wp,1.0_wp])
    end subroutine test_requested_output_range

    subroutine test_invalid_and_recovery()
        type(material_t) :: material
        real(wp) :: electric(ndim,2), displacement(ndim,2), wrong(ndim-1,2)
        real(wp) :: extra(ndim,1,2)
        integer :: status

        call material%init(dielectric_cubic_t(0.0_wp,0.2_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(dielectric_cubic_t(1.0_wp,-0.1_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(dielectric_cubic_t(1.0_wp, &
                           ieee_value(0.0_wp,ieee_quiet_nan)),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(dielectric_cubic_t(3.0_wp,0.2_wp),2,status)
        call expect(status,cf_ok)
        electric = 0.0_wp
        electric(1,:) = 1.0_wp
        electric(2,2) = ieee_value(0.0_wp,ieee_positive_inf)
        call material%eval(electric,displacement,status=status)
        call expect(status,cf_evaluation_failed)
        if (material%last_failure_point() /= 2) error stop "Cubic failure point"
        call material%commit(status)
        call expect(status,cf_no_trial)
        electric(2,2) = ieee_value(0.0_wp,ieee_quiet_nan)
        call material%eval(electric,displacement,status=status)
        call expect(status,cf_evaluation_failed)
        electric(2,2) = huge(1.0_wp)/4.0_wp
        call material%eval(electric,displacement,status=status)
        call expect(status,cf_evaluation_failed)
        electric(2,2) = 2.0_wp
        call material%eval(electric,displacement,status=status)
        call expect(status,cf_ok)
        call close(displacement(1,2),4.0_wp)
        call material%eval(wrong,displacement,status=status)
        call expect(status,cf_invalid_shape)
        call material%eval(electric,displacement,status=status,dsdp=extra)
        call expect(status,cf_unsupported)
    end subroutine test_invalid_and_recovery

end program test_dielectric_cubic
