! SPDX-License-Identifier: Apache-2.0

module constituf_isotropic_elastic
    !! Three-dimensional, infinitesimal-strain, isotropic, linear elasticity.
    !!
    !! The isothermal stored energy per reference volume is
    !! \[\psi(\boldsymbol{\varepsilon}) =
    !! \frac{K}{2}[\mathrm{tr}(\boldsymbol{\varepsilon})]^2 +
    !! G\,\boldsymbol{\varepsilon}_{\mathrm{dev}}:
    !! \boldsymbol{\varepsilon}_{\mathrm{dev}},\]
    !! with the symmetric infinitesimal strain \(\boldsymbol{\varepsilon}\),
    !! \(\boldsymbol{\varepsilon}_{\mathrm{dev}} =
    !! \boldsymbol{\varepsilon} -
    !! \mathrm{tr}(\boldsymbol{\varepsilon})\boldsymbol{I}/3\),
    !! bulk modulus \(K>0\) and shear modulus \(G>0\), both in the chosen stress unit.
    !! The stress is
    !! \[\boldsymbol{\sigma} =
    !! K\,\mathrm{tr}(\boldsymbol{\varepsilon})\boldsymbol{I}
    !! + 2G\boldsymbol{\varepsilon}_{\mathrm{dev}}.\]
    !! The exact discrete response is the same algebraic map at each trial.
    !! There is no integration algorithm, rate dependence, or internal state.
    !! The exact material Jacobian at fixed parameters is
    !! \[\mathbb{C}_{\mathrm{alg}} =
    !! K\boldsymbol{I}\otimes\boldsymbol{I}
    !! + 2G\mathbb{I}_{\mathrm{dev}}.\]
    !!
    !! Strain and Cauchy stress use the same work-conjugate 3D Mandel map:
    !! \[\boldsymbol{A}_M=
    !! [A_{11},A_{22},A_{33},\sqrt{2}A_{23},\sqrt{2}A_{13},\sqrt{2}A_{12}]^T,\]
    !! where \(\boldsymbol{A}\) denotes either symmetric tensor. The fourth-order
    !! Jacobian is a \(6\times6\) Mandel matrix, with shear diagonals \(2G\).
    !! The sign convention is positive strain/stress in tension. All fields
    !! use a shared spatial Cartesian frame. This small-deformation model
    !! must not be used for finite rotations, large strains or thermomechanics.
    !!
    !! The optional analytical derivatives at fixed prescribed strain are
    !! \[\frac{\partial\boldsymbol{\sigma}}{\partial K}
    !! = \mathrm{tr}(\boldsymbol{\varepsilon})\boldsymbol{I},\qquad
    !! \frac{\partial\boldsymbol{\sigma}}{\partial G}
    !! = 2\boldsymbol{\varepsilon}_{\mathrm{dev}}.\]
    !! Because the model has no state, these are also its exact discrete-path
    !! sensitivities. They are not sensitivities of a parameter-dependent
    !! global equilibrium solution.
    !!
    !! References:
    !!
    !! J. C. Simo and T. J. R. Hughes, *Computational Inelasticity*,
    !! Interdisciplinary Applied Mathematics 7, Springer, 1998,
    !! DOI: 10.1007/b98904. ISBN 978-0-387-97520-7.
    !!
    !! T. Nagel, U.-J. Görke, K. M. Moerman and O. Kolditz,
    !! *On advantages of the Kelvin mapping in finite element implementations
    !! of deformation processes*, Environmental Earth Sciences 75, 937 (2016),
    !! DOI: 10.1007/s12665-016-5429-4. Equations (19), (20), (25) and (26)
    !! give the orthonormal tensor mapping, isotropic stiffness and stress.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: nmandel = 6

    public :: isotropic_elastic_t

    type, extends(material_model_t) :: isotropic_elastic_t
        !! Stateless 3D Hooke material, specified by positive bulk and shear moduli.
        !!
        !! Both moduli are in the chosen stress unit. Construct with
        !! `isotropic_elastic_t(bulk_modulus, shear_modulus)`, then pass the
        !! instance to `material_t%init`. The handle clones these parameters
        !! during setup. Changing the caller's instance cannot change it.
        !! Unsupported temperature, time and finite-strain effects are absent.
        real(wp) :: bulk_modulus = 0.0_wp
            !! Bulk modulus \(K\) in the chosen stress unit, finite and strictly positive.
        real(wp) :: shear_modulus = 0.0_wp
            !! Shear modulus \(G\) in the chosen stress unit, finite and strictly positive.
    contains
        procedure :: describe => describe_isotropic_elastic
            !! Declare the constitutive fields, parameter metadata and supported capabilities.
        procedure :: eval_batch => evaluate_isotropic_elastic
            !! Evaluate the batched constitutive response, optional tangent and trial history.
        procedure :: eval_sensitivity_batch => sensitivity_isotropic_elastic
            !! Evaluate exact parameter-path derivatives and corresponding trial history.
    end type isotropic_elastic_t

contains

    subroutine describe_isotropic_elastic(self, spec, status)
        !! Declare one 6-component Mandel strain and one Mandel Cauchy stress.
        !!
        !! Input: `strain.small`, dimensionless. Output: `stress.cauchy`, in the chosen stress unit.
        !! Both have physical shape `[3,3]`, spatial frame and 6 stored values.
        !! No history or time-step input is required. Invalid or numerically
        !! unrepresentable positive moduli return `cf_invalid_argument`.
        class(isotropic_elastic_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema receiving field layouts, parameters and supported operations.
        integer, intent(out) :: status
            !! Operation status. `cf_ok` indicates a valid model or successful initialization.

        status = cf_invalid_argument
        if (.not. ieee_is_finite(self%bulk_modulus)) return
        if (.not. ieee_is_finite(self%shear_modulus)) return
        if (self%bulk_modulus <= 0.0_wp .or. self%shear_modulus <= 0.0_wp) return
        if (.not. ieee_is_finite(2.0_wp*self%shear_modulus)) return
        if (.not. ieee_is_finite(self%bulk_modulus + (4.0_wp/3.0_wp)*self%shear_modulus)) return

        call spec%add_input("strain.small", 6, "mandel", "small_strain", "spatial", [3,3])
        call spec%add_output("stress.cauchy", 6, "mandel", "cauchy_stress", "spatial", [3,3])
        call spec%add_parameter("bulk_modulus")
        call spec%add_parameter("shear_modulus")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_isotropic_elastic

    subroutine evaluate_isotropic_elastic(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Evaluate all points independently. No allocations or persistent writes.
        !!
        !! `x_new(6,npoints)` is infinitesimal Mandel strain. `response(6,npoints)`
        !! is Cauchy stress in the chosen stress unit. `tangent(6,6,npoints)`, if requested, is the
        !! exact derivative of stress with respect to the Mandel strain.
        !! Invalid/nonfinite inputs or responses are flagged per point and
        !! the caller must discard all outputs after any failed point.
        class(isotropic_elastic_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed input `(ninput,npoints)` if required, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `strain.small` input `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued state `(nstate,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Trial real-valued state `(nstate,npoints)` returned by this evaluation.
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.cauchy` response `(noutput,npoints)` for the evaluated material points.
        real(wp), intent(in) :: dt
            !! Prescribed time increment. It is positive when the material requires time integration.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete response Jacobian `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(inout) :: point_status(:)
            !! Per-point status values; successful evaluations must set each entry to `cf_ok`.

        state_trial = state_old
        call evaluate_points(self, x_new, response, point_status, tangent)
    end subroutine evaluate_isotropic_elastic

    subroutine sensitivity_isotropic_elastic(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Fused stress, exact Mandel tangent and analytical parameter derivative.
        !!
        !! `selected_parameters` contains setup-resolved parameter indices:
        !! 1 for `bulk_modulus`, 2 for `shear_modulus`. `dsdp` has shape
        !! `(6,nselected,npoints)`, ordered as requested at initialization,
        !! dimensionless. The zero-sized stateless sensitivity history is
        !! unchanged. Input history is held fixed.
        class(isotropic_elastic_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed input `(ninput,npoints)` if required, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `strain.small` input `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued state `(nstate,npoints)`.
        integer, intent(in) :: selected_parameters(:)
            !! One-based IDs of selected model parameters, in sensitivity-axis order.
        real(wp), intent(in) :: sensitivity_old(:,:,:)
            !! Committed real-state derivatives `(nstate,nselected,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Trial real-valued state `(nstate,npoints)` returned by this evaluation.
        real(wp), intent(out) :: sensitivity_trial(:,:,:)
            !! Trial real-state derivatives `(nstate,nselected,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.cauchy` response `(noutput,npoints)` for the evaluated material points.
        real(wp), intent(out) :: dsdp(:,:,:)
            !! Exact total response derivatives `(noutput,nselected,npoints)` at fixed prescribed input history.
        real(wp), intent(in) :: dt
            !! Prescribed time increment. It is positive when the material requires time integration.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete response Jacobian `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(inout) :: point_status(:)
            !! Per-point status values; successful evaluations must set each entry to `cf_ok`.

        state_trial = state_old
        sensitivity_trial = sensitivity_old
        call evaluate_points(self, x_new, response, point_status, tangent, selected_parameters, dsdp)
    end subroutine sensitivity_isotropic_elastic

    subroutine evaluate_points(self, x_new, response, point_status, tangent, selected_parameters, dsdp)
        !! Shared stateless kernel: no parameter lookups or heap allocation.
        class(isotropic_elastic_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `strain.small` input `(ninput,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.cauchy` response `(noutput,npoints)` for the evaluated material points.
        integer, intent(inout) :: point_status(:)
            !! Per-point status values; successful evaluations must set each entry to `cf_ok`.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete response Jacobian `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(in), optional :: selected_parameters(:)
            !! One-based IDs of selected model parameters, in sensitivity-axis order.
        real(wp), intent(out), optional :: dsdp(:,:,:)
            !! Exact total response derivatives `(noutput,nselected,npoints)` at fixed prescribed input history.
        real(wp) :: tangent_matrix(nmandel,nmandel)
        real(wp) :: trace_strain
        integer :: i
        integer :: j

        if (present(tangent)) call elastic_tangent(self%bulk_modulus, self%shear_modulus, tangent_matrix)

        point_loop: do i = 1, size(x_new,2)
            point_status(i) = cf_invalid_argument
            if (.not. all(ieee_is_finite(x_new(:,i)))) cycle point_loop
            call elastic_point(self%bulk_modulus, self%shear_modulus, x_new(:,i), response(:,i), trace_strain)
            if (.not. all(ieee_is_finite(response(:,i)))) cycle point_loop
            if (present(tangent)) tangent(:,:,i) = tangent_matrix

            if (present(selected_parameters)) then
                parameter_loop: do j = 1, size(selected_parameters)
                    select case (selected_parameters(j))
                    case (1)
                        dsdp(1:3,j,i) = trace_strain
                        dsdp(4:6,j,i) = 0.0_wp
                    case (2)
                        dsdp(1:3,j,i) = 2.0_wp * (x_new(1:3,i) - trace_strain/3.0_wp)
                        dsdp(4:6,j,i) = 2.0_wp*x_new(4:6,i)
                    case default
                        exit parameter_loop
                    end select
                end do parameter_loop
                if (j <= size(selected_parameters)) cycle point_loop
                if (.not. all(ieee_is_finite(dsdp(:,:,i)))) cycle point_loop
            end if
            point_status(i) = cf_ok
        end do point_loop
    end subroutine evaluate_points

    pure subroutine elastic_point(bulk, shear, strain, stress, trace_strain)
        !! Hooke response and volumetric strain trace in 3D Mandel components.
        real(wp), intent(in) :: bulk
            !! Bulk modulus \(K\), in the chosen stress unit.
        real(wp), intent(in) :: shear
            !! Shear modulus \(G\), in the chosen stress unit.
        real(wp), intent(in) :: strain(nmandel)
            !! Prescribed small-strain components in the material's tensor convention.
        real(wp), intent(out) :: stress(nmandel)
            !! Cauchy stress components returned in the material's tensor convention.
        real(wp), intent(out) :: trace_strain
            !! Trace of the small-strain tensor.

        trace_strain = sum(strain(1:3))
        stress(1:3) = bulk*trace_strain + 2.0_wp*shear*(strain(1:3)-trace_strain/3.0_wp)
        stress(4:6) = 2.0_wp*shear*strain(4:6)
    end subroutine elastic_point

    pure subroutine elastic_tangent(bulk, shear, tangent)
        !! Exact \(6\times6\) Mandel stiffness. Shear diagonal entries are \(2G\).
        real(wp), intent(in) :: bulk
            !! Bulk modulus \(K\), in the chosen stress unit.
        real(wp), intent(in) :: shear
            !! Shear modulus \(G\), in the chosen stress unit.
        real(wp), intent(out) :: tangent(nmandel,nmandel)
            !! Optional exact discrete response Jacobian `(noutput,ninput,npoints)` at fixed committed history.
        real(wp) :: lambda
        integer :: i

        lambda = bulk - (2.0_wp/3.0_wp)*shear
        tangent = 0.0_wp
        tangent(1:3,1:3) = lambda
        do i = 1, 3
            tangent(i,i) = bulk + (4.0_wp/3.0_wp)*shear
        end do
        do i = 4, 6
            tangent(i,i) = 2.0_wp*shear
        end do
    end subroutine elastic_tangent

end module constituf_isotropic_elastic
