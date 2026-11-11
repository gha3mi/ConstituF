! SPDX-License-Identifier: Apache-2.0

module constituf_potential_verify
    !! Numerical verification of constitutive energy/potential conjugacy.
    !!
    !! A stored energy, enthalpy or other smooth scalar potential \(\Phi\)
    !! may have mixed-sign physical conjugate outputs \(y_{c_j}\) for inputs \(x_j\):
    !! \[\frac{\partial\Phi}{\partial x_j}=s_j y_{c_j},
    !! \qquad s_j\in\{-1,+1\}.\]
    !! At fixed committed history and material parameters, this routine
    !! independently compares centered numerical differences of the
    !! declared scalar potential response against the supplied conjugate
    !! outputs. It does not use or generate a production material tangent
    !! or parameter derivative. For example: hyperelastic Helmholtz energy
    !! is work-conjugate to first Piola stress (all signs \(+1\)), while an
    !! electric enthalpy has strain-stress signs \(+1\) and electric-field/displacement
    !! \((\boldsymbol E,\boldsymbol D)\) signs \(-1\).
    !!
    !! Verification allocates workspace and performs repeated trial evals.
    !! It is intended for tests/scientific debugging, not FEM kernels.
    !! The routine never commits any trial and discards pending trials on
    !! completion. Committed history is unchanged. Avoid switches where
    !! the classical derivative is not unique. Input components with
    !! different units should use an explicitly dimensioned step vector.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_kinds, only: wp
    use constituf_material, only: material_t
    use constituf_status, only: cf_ok, cf_invalid_argument, cf_invalid_shape, cf_evaluation_failed, cf_allocation_failed
    implicit none (type, external)
    private

    public :: verify_potential

contains

    subroutine verify_potential(material,input,potential_index,conjugate_indices,conjugate_signs,max_scaled_error,status, dt,step,relative_tolerance,absolute_tolerance)
        !! Test an energy-conjugate constitutive gradient at a fixed history.
        !!
        !! `input` has shape `(nin,npoints)`. `potential_index` is one-based and
        !! selects the scalar output response for the potential.
        !! `conjugate_indices(nin)` selects output component \(c_j\) for each input \(x_j\),
        !! with `conjugate_signs(nin)` providing \(s_j\in\{-1,+1\}\). For example,
        !! the electric-field/displacement pair \((\boldsymbol E,\boldsymbol D)\) has sign \(-1\).
        !! The derivative is \(\partial\Phi/\partial x_j\) in the
        !! corresponding potential/input units. `step(nin)` sets positive
        !! physical perturbation sizes. The default is
        !! \(h_j=\epsilon_{\mathrm{mach}}^{1/3}\max(1,\lvert x_j\rvert)\),
        !! where \(\epsilon_{\mathrm{mach}}\) is the working-kind machine precision.
        !! This default is suitable only for uniformly scaled fields.
        !! The numerical reference magnitude 1 uses the caller's representation.
        !! Changing units may require different step and tolerance values.
        !! No unit metadata is consulted.
        !!
        !! Let \(g_j^{\mathrm{fd}}\) denote the centered finite-difference
        !! approximation to \(\partial\Phi/\partial x_j\). `max_scaled_error` is
        !! the maximum over \(j\) and all points of
        !! \[\frac{\lvert g_j^{\mathrm{fd}}-s_j y_{c_j}\rvert}
        !! {a_{\mathrm{tol}}+r_{\mathrm{tol}}\max(1,\lvert g_j^{\mathrm{fd}}\rvert,
        !! \lvert s_j y_{c_j}\rvert)}.\]
        !! The default relative tolerance is \(r_{\mathrm{tol}}=256\epsilon_{\mathrm{mach}}^{2/3}\)
        !! and the default absolute tolerance is \(a_{\mathrm{tol}}=0\).
        !! Success requires the maximum normalized error to be at most \(1\).
        !! This criterion requires the user to
        !! choose physically appropriate scales for mixed-unit systems.
        type(material_t), intent(inout) :: material
            !! Material handle whose configuration and history are queried or updated.
        real(wp), intent(in) :: input(:,:)
            !! Prescribed packed input values at which derivatives are checked.
        integer, intent(in) :: potential_index
            !! One-based response component containing the scalar potential.
        integer, intent(in) :: conjugate_indices(:)
            !! One-based output indices conjugate to the selected input components.
        real(wp), intent(in) :: conjugate_signs(:)
            !! Work-conjugacy signs used when comparing energy derivatives with outputs.
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
        real(wp), allocatable :: predicted(:,:)
        real(wp), allocatable :: positive(:,:)
        real(wp), allocatable :: negative(:,:)
        real(wp), allocatable :: test_input(:,:)
        real(wp), allocatable :: dsdp(:,:,:)
        real(wp) :: rtol
        real(wp) :: atol
        real(wp) :: h
        real(wp) :: finite_difference
        real(wp) :: exact
        real(wp) :: scale
        real(wp) :: discrepancy
        integer :: nin
        integer :: nout
        integer :: npoints
        integer :: j
        integer :: i
        integer :: allocation_status
        character(256) :: allocation_message

        max_scaled_error = huge(1.0_wp)
        status = cf_invalid_shape
        nin = material%input_size()
        nout = material%output_size()
        npoints = material%point_count()
        if (nin <= 0 .or. nout <= 0 .or. npoints <= 0) return
        if (size(input,1) /= nin .or. size(input,2) /= npoints) return
        if (size(conjugate_indices) /= nin .or. size(conjugate_signs) /= nin) return
        if (present(step)) then
            if (size(step) /= nin) return
        end if
        status = cf_invalid_argument
        if (potential_index < 1 .or. potential_index > nout) return
        if (any(conjugate_indices < 1) .or. any(conjugate_indices > nout)) return
        if (.not. all(ieee_is_finite(input))) return
        if (.not. all(ieee_is_finite(conjugate_signs))) return
        if (any(abs(abs(conjugate_signs)-1.0_wp) > 32.0_wp*epsilon(1.0_wp))) return
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

        allocate(predicted(nout,npoints),stat=allocation_status,errmsg=allocation_message)
        if (allocation_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        allocate(positive(nout,npoints),stat=allocation_status,errmsg=allocation_message)
        if (allocation_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        allocate(negative(nout,npoints),stat=allocation_status,errmsg=allocation_message)
        if (allocation_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        allocate(test_input(nin,npoints),stat=allocation_status,errmsg=allocation_message)
        if (allocation_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        allocate(dsdp(nout,material%sensitivity_size(),npoints),stat=allocation_status,errmsg=allocation_message)
        if (allocation_status /= 0) then
            status = cf_allocation_failed
            return
        end if

        if (material%sensitivity_size() > 0) then
            call material%eval(input,predicted,status=status,dt=dt,dsdp=dsdp)
        else
            call material%eval(input,predicted,status=status,dt=dt)
        end if
        if (status /= cf_ok) then
            call material%revert()
            return
        end if

        max_scaled_error = 0.0_wp
        input_loop: do j = 1,nin
            test_input = input
            do i = 1,npoints
                h = epsilon(1.0_wp)**(1.0_wp/3.0_wp)*max(1.0_wp,abs(input(j,i)))
                if (present(step)) h = step(j)
                test_input(j,i) = input(j,i)+h
            end do
            if (material%sensitivity_size() > 0) then
                call material%eval(test_input,positive,status=status,dt=dt,dsdp=dsdp)
            else
                call material%eval(test_input,positive,status=status,dt=dt)
            end if
            if (status /= cf_ok) exit input_loop
            do i = 1,npoints
                h = epsilon(1.0_wp)**(1.0_wp/3.0_wp)*max(1.0_wp,abs(input(j,i)))
                if (present(step)) h = step(j)
                test_input(j,i) = input(j,i)-h
            end do
            if (material%sensitivity_size() > 0) then
                call material%eval(test_input,negative,status=status,dt=dt,dsdp=dsdp)
            else
                call material%eval(test_input,negative,status=status,dt=dt)
            end if
            if (status /= cf_ok) exit input_loop

            do i = 1,npoints
                h = epsilon(1.0_wp)**(1.0_wp/3.0_wp)*max(1.0_wp,abs(input(j,i)))
                if (present(step)) h = step(j)
                finite_difference = (positive(potential_index,i)-negative(potential_index,i))/(2.0_wp*h)
                exact = conjugate_signs(j)*predicted(conjugate_indices(j),i)
                scale = atol+rtol*max(1.0_wp,abs(finite_difference),abs(exact))
                discrepancy = abs(finite_difference-exact)/scale
                if (.not. ieee_is_finite(discrepancy)) then
                    status = cf_evaluation_failed
                    exit input_loop
                end if
                max_scaled_error = max(max_scaled_error,discrepancy)
            end do
        end do input_loop
        call material%revert()
        if (status == cf_ok .and. max_scaled_error > 1.0_wp) status = cf_evaluation_failed
    end subroutine verify_potential

end module constituf_potential_verify
