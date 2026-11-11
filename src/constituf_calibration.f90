! SPDX-License-Identifier: Apache-2.0

module constituf_calibration
    !! Fortran-native damped least-squares calibration from material responses.
    !!
    !! A host-provided pure-in-effect callback constructs and evaluates its
    !! chosen constitutive model over the *entire prescribed experiment*,
    !! including all accepted history steps when the law is path dependent.
    !! It returns the predicted observations and the analytical material
    !! parameter Jacobian at the current candidate parameters. ConstituF's
    !! `material_t%eval(dsdp=...)` provides the required exact derivatives
    !! where the selected model supports them. This module never invents a
    !! finite-difference derivative, introspects arbitrary derived-type
    !! constructors, or assumes a particular mechanical stress measure.
    !!
    !! Given parameters \(p_j\), observations \(y_i^{\mathrm{obs}}\), known positive
    !! observation standard deviations \(s_i\) and callback predictions \(y_i(p)\), minimize
    !! \[\mathcal L(p)=\frac12\sum_i
    !! \left(\frac{y_i(p)-y_i^{\mathrm{obs}}}{s_i}\right)^2.\]
    !! Use dimensionless scaled variables \(z_j=p_j/a_j\) with user-supplied
    !! or default positive parameter scales \(a_j\). The weighted residual is
    !! \(r_i=(y_i(p)-y_i^{\mathrm{obs}})/s_i\), the weighted local
    !! Jacobian is \(J_{ij}=(a_j/s_i)\,\mathrm{d}y_i/\mathrm{d}p_j\) and each damped
    !! Gauss-Newton step solves
    !! \[(J^{\mathsf T}J+\lambda D)\delta z=-J^{\mathsf T}r,\qquad
    !! D_{jj}=\max(1,(J^{\mathsf T}J)_{jj}).\]
    !! A pivoted linear solve, finite-value checks, diagonal damping,
    !! bounded parameter projection and monotone objective acceptance
    !! prevent invalid trials from replacing the best estimate. Damping
    !! is adjusted after accepted/rejected steps. A rejected candidate
    !! never changes the best accepted parameters.
    !!
    !! This is a small dense, *local* Levenberg--Marquardt specialization,
    !! not a statistical inference library. It neither guarantees global
    !! convergence nor estimates confidence bounds/covariance. Parameter
    !! non-identifiability or severely ill-conditioned models may fail.
    !! Use the existing sensitivity/Fisher analysis tools to assess
    !! identifiability. The supplied observation scales must reflect
    !! physically meaningful measurement units/uncertainties.
    !!
    !! The callback is invoked serially, with candidate parameters and
    !! must be free of externally persistent side effects. Allocate and
    !! reinitialize material handles inside it as needed, commit only
    !! accepted *local loading-history* steps for that trial and never
    !! mutate the host FEM's current material state. This routine allocates
    !! optimizer work arrays and is intended only for offline calibration.
    !!
    !! Reference: D. W. Marquardt, "An Algorithm for Least-Squares
    !! Estimation of Nonlinear Parameters", Journal of the Society for
    !! Industrial and Applied Mathematics 11 (1963), 431--441,
    !! DOI 10.1137/0111030. The bounded projection, scaling, damping and
    !! stopping rules below are original engineering specializations.
    !! Code provenance: independently authored. No third-party code reused.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_kinds, only: wp
    use constituf_status, only: cf_ok, cf_invalid_argument, cf_invalid_shape, cf_evaluation_failed, cf_allocation_failed
    implicit none (type, external)
    private

    public :: calibrate_parameters

    abstract interface
        subroutine calibration_forward(parameters,predicted,jacobian,status)
            !! The callback returns `predicted(nobservations)` and the
            !! exact `jacobian(nobservations,nparameters)`. The caller must
            !! evaluate a fresh, reproducible complete experiment.
            import :: wp
            implicit none (type, external)
            real(wp), intent(in) :: parameters(:)
                !! Candidate material parameters in their declared order and physical units.
            real(wp), intent(out) :: predicted(:)
                !! Model-predicted observations in the prescribed experiment.
            real(wp), intent(out) :: jacobian(:,:)
                !! Exact observation-by-parameter Jacobian `(nobservations,nparameters)`.
            integer, intent(out) :: status
                !! Operation status code; `cf_ok` indicates success.
        end subroutine calibration_forward
    end interface

contains

    subroutine calibrate_parameters(initial,observed,noise_scale,lower_bound,upper_bound,evaluate,fitted,objective,status, parameter_scale,max_iterations,relative_tolerance,iterations)
        !! Fit selected scalar material parameters using exact model Jacobians.
        !!
        !! initial, lower_bound, upper_bound, fitted and optional
        !! parameter_scale have shape (nparameters). Bounds and all
        !! parameter scales share each parameter's physical unit.
        !! observed,noise_scale have shape (nobservations), possibly
        !! heterogeneous physical units. Each noise_scale entry is
        !! strictly positive in the corresponding observation's unit.
        !! evaluate(p,predicted,jacobian,status) supplies model predictions
        !! and exact derivatives, with jacobian(nobservations,nparameters)
        !! at fixed prescribed input histories. Observations must be
        !! sufficient in number to identify nparameters, but rank and
        !! identifiability require independent scientific examination.
        !!
        !! The optional max_iterations defaults to 60. The relative
        !! stopping tolerance defaults to \(10^{-8}\) in *scaled* variables.
        !! fitted receives the best accepted estimate even on failure.
        !! Objective is the corresponding dimensionless half chi-square.
        !! iterations reports accepted optimization steps. Return
        !! cf_ok after an interior stopping criterion is reached,
        !! cf_invalid_shape or cf_invalid_argument for bad inputs,
        !! cf_allocation_failed for allocation errors and
        !! cf_evaluation_failed for failed/undefined candidate models,
        !! failed linear solves, or exhausted iteration budget.
        real(wp), intent(in) :: initial(:)
            !! Initial candidate parameter values in their physical units.
        real(wp), intent(in) :: observed(:)
            !! Measured observations used as calibration targets.
        real(wp), intent(in) :: noise_scale(:)
            !! Positive observation standard deviations used for residual weighting.
        real(wp), intent(in) :: lower_bound(:)
            !! Inclusive lower bounds for each fitted parameter.
        real(wp), intent(in) :: upper_bound(:)
            !! Inclusive upper bounds for each fitted parameter.
        procedure(calibration_forward) :: evaluate
            !! Callback that computes predictions and their exact parameter Jacobian for each candidate.
        real(wp), intent(out) :: fitted(:)
            !! Best accepted parameter estimates in the original physical units.
        real(wp), intent(out) :: objective
            !! Accepted half sum of squared, uncertainty-weighted residuals.
        integer, intent(out) :: status
            !! Operation status code; `cf_ok` indicates success.
        real(wp), intent(in), optional :: parameter_scale(:)
            !! Positive parameter scales for dimensionless optimization coordinates.
        integer, intent(in), optional :: max_iterations
            !! Maximum number of calibration iterations; defaults to 60.
        real(wp), intent(in), optional :: relative_tolerance
            !! Relative tolerance for scaled derivative checks or calibration convergence.
        integer, intent(out), optional :: iterations
            !! Number of accepted calibration iterations.
        real(wp), allocatable :: predicted(:)
        real(wp), allocatable :: candidate_predicted(:)
        real(wp), allocatable :: jacobian(:,:)
        real(wp), allocatable :: candidate_jacobian(:,:)
        real(wp), allocatable :: weighted_jacobian(:,:)
        real(wp), allocatable :: residual(:)
        real(wp), allocatable :: scales(:)
        real(wp), allocatable :: normal(:,:)
        real(wp), allocatable :: system(:,:)
        real(wp), allocatable :: gradient(:)
        real(wp), allocatable :: candidate(:)
        real(wp), allocatable :: step(:)
        real(wp), allocatable :: diagonal(:)
        real(wp) :: damping
        real(wp) :: tol
        real(wp) :: cost
        real(wp) :: candidate_cost
        real(wp) :: improvement
        real(wp) :: relative_step
        integer :: nparameters
        integer :: nobservations
        integer :: nmax
        integer :: iteration
        integer :: attempt
        integer :: i
        integer :: allocation_status
        integer :: eval_status
        character(256) :: allocation_message
        logical :: solved
        logical :: accepted

        if (present(iterations)) iterations = 0
        objective = huge(1.0_wp)
        status = cf_invalid_shape
        nparameters = size(initial)
        nobservations = size(observed)
        if (nparameters < 1 .or. nobservations < nparameters) return
        if (size(fitted) /= nparameters .or. size(noise_scale) /= nobservations) return
        if (size(lower_bound) /= nparameters .or. size(upper_bound) /= nparameters) return
        if (present(parameter_scale)) then
            if (size(parameter_scale) /= nparameters) return
        end if
        fitted = initial

        status = cf_invalid_argument
        if (.not. all(ieee_is_finite(initial)) .or. .not. all(ieee_is_finite(observed))) return
        if (.not. all(ieee_is_finite(noise_scale)) .or. any(noise_scale <= 0.0_wp)) return
        if (.not. all(ieee_is_finite(lower_bound)) .or. .not. all(ieee_is_finite(upper_bound))) return
        if (any(upper_bound <= lower_bound)) return
        if (any(initial < lower_bound) .or. any(initial > upper_bound)) return
        nmax = 60
        if (present(max_iterations)) nmax = max_iterations
        if (nmax < 1) return
        tol = 1.0e-8_wp
        if (present(relative_tolerance)) tol = relative_tolerance
        if (.not. ieee_is_finite(tol) .or. tol <= 0.0_wp) return

        allocate(predicted(nobservations),stat=allocation_status,errmsg=allocation_message)
        if (allocation_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        allocate(candidate_predicted(nobservations),stat=allocation_status,errmsg=allocation_message)
        if (allocation_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        allocate(jacobian(nobservations,nparameters),stat=allocation_status,errmsg=allocation_message)
        if (allocation_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        allocate(candidate_jacobian(nobservations,nparameters),stat=allocation_status,errmsg=allocation_message)
        if (allocation_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        allocate(weighted_jacobian(nobservations,nparameters),stat=allocation_status,errmsg=allocation_message)
        if (allocation_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        allocate(residual(nobservations),stat=allocation_status,errmsg=allocation_message)
        if (allocation_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        allocate(scales(nparameters),stat=allocation_status,errmsg=allocation_message)
        if (allocation_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        allocate(normal(nparameters,nparameters),stat=allocation_status,errmsg=allocation_message)
        if (allocation_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        allocate(system(nparameters,nparameters),stat=allocation_status,errmsg=allocation_message)
        if (allocation_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        allocate(gradient(nparameters),stat=allocation_status,errmsg=allocation_message)
        if (allocation_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        allocate(candidate(nparameters),stat=allocation_status,errmsg=allocation_message)
        if (allocation_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        allocate(step(nparameters),stat=allocation_status,errmsg=allocation_message)
        if (allocation_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        allocate(diagonal(nparameters),stat=allocation_status,errmsg=allocation_message)
        if (allocation_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        if (present(parameter_scale)) then
            scales = parameter_scale
        else
            scales = abs(initial)
            do i = 1,nparameters
                if (scales(i) <= tiny(1.0_wp)) scales(i) = 0.5_wp*(upper_bound(i)-lower_bound(i))
            end do
        end if
        status = cf_invalid_argument
        if (.not. all(ieee_is_finite(scales))) return
        if (any(scales <= 0.0_wp)) return

        status = cf_evaluation_failed
        call evaluate(fitted,predicted,jacobian,eval_status)
        if (eval_status /= cf_ok) return
        if (.not. all(ieee_is_finite(predicted)) .or. .not. all(ieee_is_finite(jacobian))) return
        residual = (predicted-observed)/noise_scale
        cost = 0.5_wp*dot_product(residual,residual)
        if (.not. ieee_is_finite(cost)) return
        objective = cost
        damping = 1.0e-3_wp

        iteration_loop: do iteration = 1,nmax
            do i = 1,nparameters
                weighted_jacobian(:,i) = (jacobian(:,i)/noise_scale)*scales(i)
            end do
            normal = matmul(transpose(weighted_jacobian),weighted_jacobian)
            gradient = matmul(transpose(weighted_jacobian),residual)
            if (.not. all(ieee_is_finite(normal)) .or. .not. all(ieee_is_finite(gradient))) return
            if (maxval(abs(gradient)) <= tol*max(1.0_wp,sqrt(2.0_wp*cost))) then
                status = cf_ok
                return
            end if
            do i = 1,nparameters
                diagonal(i) = max(1.0_wp,normal(i,i))
            end do
            accepted = .false.
            damping_loop: do attempt = 1,16
                system = normal
                do i = 1,nparameters
                    system(i,i) = system(i,i)+damping*diagonal(i)
                end do
                call solve_pivoted(system,-gradient,step,solved)
                if (.not. solved) exit damping_loop
                candidate = min(upper_bound,max(lower_bound,fitted+scales*step))
                relative_step = maxval(abs((candidate-fitted)/scales))
                if (relative_step <= tol) then
                    status = cf_ok
                    return
                end if
                call evaluate(candidate,candidate_predicted,candidate_jacobian,eval_status)
                if (eval_status == cf_ok) then
                    if (all(ieee_is_finite(candidate_predicted)) .and. all(ieee_is_finite(candidate_jacobian))) then
                        residual = (candidate_predicted-observed)/noise_scale
                        candidate_cost = 0.5_wp*dot_product(residual,residual)
                        if (ieee_is_finite(candidate_cost)) then
                            if (candidate_cost < cost) then
                                accepted = .true.
                                exit damping_loop
                            end if
                        end if
                    end if
                end if
                if (damping > huge(damping)/10.0_wp) exit damping_loop
                damping = 10.0_wp*damping
            end do damping_loop
            if (.not. accepted) return

            improvement = cost-candidate_cost
            fitted = candidate
            predicted = candidate_predicted
            jacobian = candidate_jacobian
            cost = candidate_cost
            residual = (predicted-observed)/noise_scale
            objective = cost
            damping = max(1.0e-12_wp,damping/3.0_wp)
            if (present(iterations)) iterations = iteration
            if (improvement <= tol*max(1.0_wp,cost) .or. relative_step <= tol) then
                status = cf_ok
                return
            end if
        end do iteration_loop
    end subroutine calibrate_parameters

    subroutine solve_pivoted(matrix,rhs,solution,ok)
        !! Dense Gaussian elimination with partial pivoting outside kernels.
        real(wp), intent(in) :: matrix(:,:)
            !! Dense coefficient matrix for the pivoted linear solve.
        real(wp), intent(in) :: rhs(:)
            !! Right-hand side of the dense linear system.
        real(wp), intent(out) :: solution(:)
            !! Solution vector returned by partial-pivot Gaussian elimination.
        logical, intent(out) :: ok
            !! True if the requested numerical operation completed successfully.
        real(wp) :: a(size(rhs),size(rhs))
        real(wp) :: b(size(rhs))
        real(wp) :: swap_row(size(rhs))
        real(wp) :: factor
        real(wp) :: threshold
        real(wp) :: tmp
        integer :: n
        integer :: i
        integer :: j
        integer :: k
        integer :: pivot

        n = size(rhs)
        a = matrix
        b = rhs
        ok = .false.
        threshold = 64.0_wp*epsilon(1.0_wp)*maxval(abs(a))
        elimination: do k = 1,n
            pivot = k-1+maxloc(abs(a(k:n,k)),dim=1)
            if (abs(a(pivot,k)) <= threshold) return
            if (pivot /= k) then
                swap_row = a(k,:)
                a(k,:) = a(pivot,:)
                a(pivot,:) = swap_row
                tmp = b(k)
                b(k) = b(pivot)
                b(pivot) = tmp
            end if
            do i = k+1,n
                factor = a(i,k)/a(k,k)
                a(i,k:n) = a(i,k:n)-factor*a(k,k:n)
                b(i) = b(i)-factor*b(k)
            end do
        end do elimination
        back_substitution: do i = n,1,-1
            tmp = b(i)
            do j = i+1,n
                tmp = tmp-a(i,j)*solution(j)
            end do
            solution(i) = tmp/a(i,i)
        end do back_substitution
        ok = all(ieee_is_finite(solution))
    end subroutine solve_pivoted

end module constituf_calibration
