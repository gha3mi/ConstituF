! SPDX-License-Identifier: Apache-2.0

module constituf_verify
    !! Reusable material-point verification outside production integration.
    !!
    !! A central difference of the *discrete* trial response is compared
    !! against the advertised current-input algorithmic tangent, holding
    !! committed history, previous input and `dt` fixed. No trial is committed.
    !! The verification operation allocates workspace and makes repeated
    !! evaluations. It must not be called in a performance-critical FEM loop.
    !!
    !! Finite differences here are solely independent checks of analytical
    !! derivatives. They never implement a material tangent or sensitivity.
    !! Do not verify at yield/contact switches or other nondifferentiable
    !! points. Supply per-component perturbation step sizes when inputs have
    !! non-uniform physical units or scales. If a trial fails, its status is
    !! returned and no comparison is claimed.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_kinds, only: wp
    use constituf_material, only: material_t
    use constituf_status, only: cf_ok, cf_invalid_argument, cf_invalid_shape, cf_evaluation_failed, cf_allocation_failed
    implicit none (type, external)
    private

    public :: verify_tangent

contains

    subroutine verify_tangent(material, input, max_scaled_error, status, dt, step, relative_tolerance, absolute_tolerance)
        !! Compare the exact discrete tangent against centered perturbations.
        !!
        !! `input(nin,npoints)` must be finite and match the initialized model.
        !! The optional positive `step(nin)` gives input-component perturbations
        !! in each input's caller-chosen units. Otherwise, for input component
        !! \(x_j\), the numerical step \(h_j=\epsilon_{\mathrm{mach}}^{1/3}
        !! \max(1,\lvert x_j\rvert)\) is used, where \(\epsilon_{\mathrm{mach}}\)
        !! is the working-kind machine precision. The caller may override positive
        !! `relative_tolerance` and nonnegative `absolute_tolerance` (in output
        !! units per input units, suitable only for compatible tangent blocks).
        !! The default dimensionless relative tolerance is \(256\epsilon_{\mathrm{mach}}^{2/3}\).
        !! The reference magnitude 1 is in the caller's numerical representation.
        !! Select step and tolerances for that representation when changing units
        !! or checking fields with different scales. No unit metadata is consulted.
        !!
        !! For a finite-difference Jacobian entry \(J_{\mathrm{fd}}\) and its
        !! analytical counterpart \(J_{\mathrm{exact}}\), `max_scaled_error` is
        !! the maximum over points, outputs and inputs of
        !! \[\frac{\lvert J_{\mathrm{fd}}-J_{\mathrm{exact}}\rvert}
        !! {a_{\mathrm{tol}}+r_{\mathrm{tol}}\max(1,\lvert J_{\mathrm{fd}}\rvert,
        !! \lvert J_{\mathrm{exact}}\rvert)}.\]
        !! Here \(a_{\mathrm{tol}}\) and \(r_{\mathrm{tol}}\) denote the absolute
        !! and relative tolerances. Return `cf_ok` if this maximum is at most
        !! \(1\), otherwise `cf_evaluation_failed`. Other statuses propagate from the evaluator.
        !! A successful or failed call discards any pending trial, without
        !! altering the committed material history. Analytical parameter
        !! sensitivities, if enabled, are requested but not verified here.
        type(material_t), intent(inout) :: material
            !! Material handle whose configuration and history are queried or updated.
        real(wp), intent(in) :: input(:,:)
            !! Prescribed packed input values at which derivatives are checked.
        real(wp), intent(out) :: max_scaled_error
            !! Largest observed error normalized by the specified tolerances.
        integer, intent(out) :: status
            !! Operation status code; `cf_ok` indicates success.
        real(wp), intent(in), optional :: dt
            !! Prescribed time increment; strictly positive when required by the model.
        real(wp), intent(in), optional :: step(:)
            !! Positive absolute perturbations in each input component's physical units.
        real(wp), intent(in), optional :: relative_tolerance
            !! Relative tolerance for scaled derivative checks or calibration convergence.
        real(wp), intent(in), optional :: absolute_tolerance
            !! Nonnegative absolute term in the scaled-error denominator, in derivative units.
        real(wp), allocatable :: y(:,:)
        real(wp), allocatable :: positive(:,:)
        real(wp), allocatable :: negative(:,:)
        real(wp), allocatable :: jacobian(:,:,:)
        real(wp), allocatable :: xtrial(:,:)
        real(wp), allocatable :: dsdp(:,:,:)
        real(wp) :: rtol
        real(wp) :: atol
        real(wp) :: h
        real(wp) :: approximate
        real(wp) :: error
        real(wp) :: scale
        integer :: nin
        integer :: nout
        integer :: npoints
        integer :: i
        integer :: j
        integer :: k
        integer :: alloc_status
        character(256) :: allocation_message

        max_scaled_error = huge(1.0_wp)
        status = cf_invalid_shape
        nin = material%input_size()
        nout = material%output_size()
        npoints = size(input,2)
        if (nin <= 0 .or. nout <= 0 .or. npoints <= 0) return
        if (size(input,1) /= nin) return
        if (present(step)) then
            if (size(step) /= nin) return
        end if
        status = cf_invalid_argument
        if (.not. all(ieee_is_finite(input))) return
        if (present(step)) then
            if (.not. all(ieee_is_finite(step))) return
            if (any(step <= 0.0_wp)) return
        end if
        rtol = 256.0_wp*epsilon(1.0_wp)**(2.0_wp/3.0_wp)
        atol = 0.0_wp
        if (present(relative_tolerance)) rtol = relative_tolerance
        if (present(absolute_tolerance)) atol = absolute_tolerance
        if (.not. ieee_is_finite(rtol) .or. .not. ieee_is_finite(atol)) return
        if (rtol <= 0.0_wp .or. atol < 0.0_wp) return

        allocate(y(nout,npoints), stat=alloc_status, errmsg=allocation_message)
        if (alloc_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        allocate(positive(nout,npoints), stat=alloc_status, errmsg=allocation_message)
        if (alloc_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        allocate(negative(nout,npoints), stat=alloc_status, errmsg=allocation_message)
        if (alloc_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        allocate(jacobian(nout,nin,npoints), stat=alloc_status, errmsg=allocation_message)
        if (alloc_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        allocate(xtrial(nin,npoints), stat=alloc_status, errmsg=allocation_message)
        if (alloc_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        allocate(dsdp(nout,material%sensitivity_size(),npoints), stat=alloc_status, errmsg=allocation_message)
        if (alloc_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        if (material%sensitivity_size() > 0) then
            call material%eval(input, y, jacobian, status, dt=dt, dsdp=dsdp)
        else
            call material%eval(input, y, jacobian, status, dt=dt)
        end if
        if (status /= cf_ok) then
            call material%revert()
            return
        end if
        max_scaled_error = 0.0_wp
        input_loop: do j = 1, nin
            xtrial = input
            do i = 1, npoints
                h = epsilon(1.0_wp)**(1.0_wp/3.0_wp)*max(1.0_wp,abs(input(j,i)))
                if (present(step)) h = step(j)
                xtrial(j,i) = input(j,i) + h
            end do
            if (material%sensitivity_size() > 0) then
                call material%eval(xtrial, positive, status=status, dt=dt, dsdp=dsdp)
            else
                call material%eval(xtrial, positive, status=status, dt=dt)
            end if
            if (status /= cf_ok) exit input_loop
            do i = 1, npoints
                h = epsilon(1.0_wp)**(1.0_wp/3.0_wp)*max(1.0_wp,abs(input(j,i)))
                if (present(step)) h = step(j)
                xtrial(j,i) = input(j,i) - h
            end do
            if (material%sensitivity_size() > 0) then
                call material%eval(xtrial, negative, status=status, dt=dt, dsdp=dsdp)
            else
                call material%eval(xtrial, negative, status=status, dt=dt)
            end if
            if (status /= cf_ok) exit input_loop
            point_loop: do i = 1, npoints
                h = epsilon(1.0_wp)**(1.0_wp/3.0_wp)*max(1.0_wp,abs(input(j,i)))
                if (present(step)) h = step(j)
                output_loop: do k = 1, nout
                    approximate = (positive(k,i)-negative(k,i))/(2.0_wp*h)
                    scale = atol + rtol*max(1.0_wp,abs(approximate),abs(jacobian(k,j,i)))
                    error = abs(approximate-jacobian(k,j,i))/scale
                    if (.not. ieee_is_finite(error)) then
                        status = cf_evaluation_failed
                        exit output_loop
                    end if
                    max_scaled_error = max(max_scaled_error,error)
                end do output_loop
                if (status /= cf_ok) exit point_loop
            end do point_loop
            if (status /= cf_ok) exit input_loop
        end do input_loop
        call material%revert()
        if (status == cf_ok .and. max_scaled_error > 1.0_wp) status = cf_evaluation_failed
    end subroutine verify_tangent

end module constituf_verify
