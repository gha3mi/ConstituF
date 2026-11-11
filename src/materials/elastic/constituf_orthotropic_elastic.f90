! SPDX-License-Identifier: Apache-2.0

module constituf_orthotropic_elastic
    !! Three-dimensional linear orthotropic elasticity in material principal axes.
    !!
    !! Nine independent engineering constants define the symmetric normal
    !! compliance \(\boldsymbol{S}\) and three shear moduli:
    !! \[\boldsymbol{S}=
    !! \begin{bmatrix}
    !! 1/E_1&-\nu_{12}/E_1&-\nu_{13}/E_1\\[2pt]
    !! -\nu_{12}/E_1&1/E_2&-\nu_{23}/E_2\\[2pt]
    !! -\nu_{13}/E_1&-\nu_{23}/E_2&1/E_3
    !! \end{bmatrix}.\]
    !! The reciprocal Poisson ratios satisfy
    !! \(\nu_{21}/E_2=\nu_{12}/E_1\),
    !! \(\nu_{31}/E_3=\nu_{13}/E_1\) and
    !! \(\nu_{32}/E_3=\nu_{23}/E_2\).
    !! The normal stiffness is \(\boldsymbol{C}_n=\boldsymbol{S}^{-1}\).
    !! In Mandel components, shear stresses satisfy
    !! \(\sigma^{\mathrm M}_{23}=2G_{23}\varepsilon^{\mathrm M}_{23}\)
    !! and similarly for 13 and 12. The energy per unit reference volume is
    !! \[\psi=\tfrac12\boldsymbol{\varepsilon}^{\mathrm M}\cdot
    !! \boldsymbol{C}^{\mathrm M}\boldsymbol{\varepsilon}^{\mathrm M},\]
    !! where \(\boldsymbol{C}^{\mathrm M}\) is the full \(6\times6\) Mandel stiffness.
    !!
    !! The material axes 1, 2, 3 are fixed and coincide with the spatial
    !! coordinate basis supplied by the solver. Neither this model nor the
    !! core rotates tensors into another material orientation. Both inputs
    !! and outputs use the symmetric 3D Mandel map
    !! \[\boldsymbol{A}^{\mathrm M}=
    !! [A_{11},A_{22},A_{33},\sqrt{2}A_{23},\sqrt{2}A_{13},\sqrt{2}A_{12}]^T,\]
    !! where \(\boldsymbol{A}\) denotes either strain or Cauchy stress.
    !! Input `strain.small` is dimensionless and output `stress.cauchy`
    !! is in the chosen stress unit. Positive values denote tension.
    !!
    !! \(E_1,E_2,E_3,G_{23},G_{13},G_{12}>0\) and all Poisson ratios are
    !! finite. Positive definiteness of \(\boldsymbol{S}\) is additionally
    !! required and verified by leading principal minors after diagonal
    !! scaling of the compliance. No arbitrary
    !! bound such as \(|\nu|<0.5\) is imposed on anisotropic materials.
    !! Any non-positive-definite or unrepresentable stiffness is rejected.
    !!
    !! The tangent \(\partial\boldsymbol{\sigma}^{\mathrm M}/
    !! \partial\boldsymbol{\varepsilon}^{\mathrm M}
    !! =\boldsymbol{C}^{\mathrm M}\) is
    !! constant and exact. Parameter sensitivities at fixed strain use
    !! \[\frac{\partial\boldsymbol{C}_n}{\partial p}
    !! =-\boldsymbol{C}_n\frac{\partial\boldsymbol{S}}{\partial p}
    !! \boldsymbol{C}_n,\qquad
    !! \frac{\partial\boldsymbol{\sigma}^{\mathrm M}}{\partial p}
    !! =\frac{\partial\boldsymbol{C}^{\mathrm M}}{\partial p}
    !! \boldsymbol{\varepsilon}^{\mathrm M}.\]
    !! The shear stiffness derivatives are \(2\) in the matching diagonal.
    !! These are analytical *local* derivatives, not global solution
    !! sensitivities. The model is stateless, infinitesimal-strain,
    !! isothermal and has no damage or finite-rotation effects.
    !!
    !! References:
    !!
    !! Dassault Systèmes, *Abaqus 2026 Material Reference*,
    !! "Linear Elastic Behavior", "Defining Orthotropic Elasticity by
    !! Specifying the Engineering Constants":
    !! https://docs.software.vt.edu/abaqusv2026/English/SIMACAEMATRefMap/simamat-c-linearelastic.htm.
    !!
    !! A. F. Bower, *Applied Mechanics of Solids*,
    !! CRC Press, 2010, ISBN 978-1-4398-0247-2, Section 3.2.13.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: ndim = 3
    integer, parameter :: nmandel = 6
    integer, parameter :: nparameters = 9

    public :: orthotropic_elastic_t

    type, extends(material_model_t) :: orthotropic_elastic_t
        !! Orthotropic solid defined in its aligned principal axes.
        !!
        !! Construct with nine engineering constants in the order below.
        !! All moduli are in the chosen stress unit. Poisson ratios are dimensionless.
        !! The material handle clones this definition during initialization.
        real(wp) :: young_1 = 0.0_wp
            !! Principal Young's modulus \(E_1>0\).
        real(wp) :: young_2 = 0.0_wp
            !! Principal Young's modulus \(E_2>0\).
        real(wp) :: young_3 = 0.0_wp
            !! Principal Young's modulus \(E_3>0\).
        real(wp) :: poisson_12 = 0.0_wp
            !! Negative transverse/axial strain ratio \(-\varepsilon_{22}/\varepsilon_{11}\) under uniaxial stress in axis 1.
        real(wp) :: poisson_13 = 0.0_wp
            !! Negative transverse/axial strain ratio \(-\varepsilon_{33}/\varepsilon_{11}\) under uniaxial stress in axis 1.
        real(wp) :: poisson_23 = 0.0_wp
            !! Negative transverse/axial strain ratio \(-\varepsilon_{33}/\varepsilon_{22}\) under uniaxial stress in axis 2.
        real(wp) :: shear_23 = 0.0_wp
            !! Principal engineering shear modulus \(G_{23}>0\).
        real(wp) :: shear_13 = 0.0_wp
            !! Principal engineering shear modulus \(G_{13}>0\).
        real(wp) :: shear_12 = 0.0_wp
            !! Principal engineering shear modulus \(G_{12}>0\).
    contains
        procedure :: describe => describe_orthotropic
            !! Declare the constitutive fields, parameter metadata and supported capabilities.
        procedure :: eval_batch => evaluate_orthotropic
            !! Evaluate the batched constitutive response, optional tangent and trial history.
        procedure :: eval_sensitivity_batch => sensitivity_orthotropic
            !! Evaluate exact parameter-path derivatives and corresponding trial history.
    end type orthotropic_elastic_t

contains

    subroutine describe_orthotropic(self, spec, status)
        !! Validate all engineering constants and declare material axes.
        class(orthotropic_elastic_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema receiving field layouts, parameters and supported operations.
        integer, intent(out) :: status
            !! Operation status. `cf_ok` indicates a valid model or successful initialization.
        real(wp) :: stiffness(nmandel,nmandel)

        call build_stiffness(self, stiffness, status)
        if (status /= cf_ok) return

        call spec%add_input("strain.small", nmandel, "mandel", "small_strain", "spatial", [ndim,ndim])
        call spec%add_output("stress.cauchy", nmandel, "mandel", "cauchy_stress", "spatial", [ndim,ndim])
        call spec%add_parameter("young_1")
        call spec%add_parameter("young_2")
        call spec%add_parameter("young_3")
        call spec%add_parameter("poisson_12")
        call spec%add_parameter("poisson_13")
        call spec%add_parameter("poisson_23")
        call spec%add_parameter("shear_23")
        call spec%add_parameter("shear_13")
        call spec%add_parameter("shear_12")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
    end subroutine describe_orthotropic

    subroutine evaluate_orthotropic(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Compute the exact Mandel stress and (optional) \(6\times6\) tangent.
        !!
        !! The state is empty. Failure invalidates all caller-owned outputs.
        class(orthotropic_elastic_t), intent(in) :: self
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
    end subroutine evaluate_orthotropic

    subroutine sensitivity_orthotropic(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Fused analytical stress and parameter derivatives.
        !!
        !! `dsdp(6,nselected,npoints)` follows the selected parameter names.
        !! The empty state and sensitivity histories remain empty.
        class(orthotropic_elastic_t), intent(in) :: self
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
    end subroutine sensitivity_orthotropic

    subroutine evaluate_points(self, x_new, response, point_status, tangent, selected_parameters, dsdp)
        !! Form stiffness once per batch. No allocations in point loops.
        class(orthotropic_elastic_t), intent(in) :: self
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
        real(wp) :: stiffness(nmandel,nmandel)
        real(wp) :: derivatives(nmandel,nmandel,nparameters)
        integer :: ipoint
        integer :: iparam
        integer :: status

        point_status = cf_invalid_argument
        call build_stiffness(self, stiffness, status)
        if (status /= cf_ok) return
        if (present(selected_parameters)) then
            do iparam = 1, size(selected_parameters)
                call stiffness_derivative(self, stiffness, selected_parameters(iparam), derivatives(:,:,iparam), status)
                if (status /= cf_ok) return
            end do
        end if

        point_loop: do ipoint = 1, size(x_new,2)
            if (.not. all(ieee_is_finite(x_new(:,ipoint)))) cycle point_loop
            response(:,ipoint) = matmul(stiffness,x_new(:,ipoint))
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop
            if (present(selected_parameters)) then
                do iparam = 1, size(selected_parameters)
                    dsdp(:,iparam,ipoint) = matmul(derivatives(:,:,iparam),x_new(:,ipoint))
                end do
                if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
            end if
            if (present(tangent)) tangent(:,:,ipoint) = stiffness
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

    subroutine build_stiffness(self, stiffness, status)
        !! Invert the diagonally scaled SPD normal compliance by cofactors.
        !!
        !! With \(\boldsymbol{D}=\operatorname{diag}(\sqrt{E_i})\), invert \(\boldsymbol{D}\boldsymbol{S}\boldsymbol{D}\) and recover \(\boldsymbol{C}=\boldsymbol{D}(\boldsymbol{D}\boldsymbol{S}\boldsymbol{D})^{-1}\boldsymbol{D}\).
        !! Reject nonpositive leading principal minors and nonfinite results.
        class(orthotropic_elastic_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp), intent(out) :: stiffness(nmandel,nmandel)
            !! Elastic stiffness matrix in the packed strain/stress convention.
        integer, intent(out) :: status
            !! Operation status. `cf_ok` indicates a valid model or successful initialization.
        real(wp) :: s(ndim,ndim)
        real(wp) :: young_roots(ndim)
        real(wp) :: determinant
        real(wp) :: minor

        status = cf_invalid_argument
        if (.not. all(ieee_is_finite([self%young_1,self%young_2,self%young_3, self%poisson_12,self%poisson_13,self%poisson_23, self%shear_23,self%shear_13,self%shear_12]))) return
        if (min(self%young_1,self%young_2,self%young_3, self%shear_23,self%shear_13,self%shear_12) <= 0.0_wp) return

        young_roots = sqrt([self%young_1,self%young_2,self%young_3])
        s = 0.0_wp
        s(1,1) = 1.0_wp
        s(2,2) = 1.0_wp
        s(3,3) = 1.0_wp
        s(1,2) = (-self%poisson_12/young_roots(1))*young_roots(2)
        s(1,3) = (-self%poisson_13/young_roots(1))*young_roots(3)
        s(2,3) = (-self%poisson_23/young_roots(2))*young_roots(3)
        s(2,1) = s(1,2)
        s(3,1) = s(1,3)
        s(3,2) = s(2,3)
        if (.not. all(ieee_is_finite(s))) return
        minor = s(1,1)*s(2,2)-s(1,2)*s(1,2)
        if (.not. ieee_is_finite(minor)) return
        if (minor <= 0.0_wp) return
        determinant = s(1,1)*(s(2,2)*s(3,3)-s(2,3)*s(2,3)) - s(1,2)*(s(1,2)*s(3,3)-s(1,3)*s(2,3)) + s(1,3)*(s(1,2)*s(2,3)-s(1,3)*s(2,2))
        if (.not. ieee_is_finite(determinant)) return
        if (determinant <= 0.0_wp) return

        stiffness = 0.0_wp
        stiffness(1,1) = ((s(2,2)*s(3,3)-s(2,3)*s(2,3))/determinant)*self%young_1
        stiffness(2,2) = ((s(1,1)*s(3,3)-s(1,3)*s(1,3))/determinant)*self%young_2
        stiffness(3,3) = (minor/determinant)*self%young_3
        stiffness(1,2) = (((s(1,3)*s(2,3)-s(1,2)*s(3,3))/determinant)*young_roots(1))*young_roots(2)
        stiffness(1,3) = (((s(1,2)*s(2,3)-s(1,3)*s(2,2))/determinant)*young_roots(1))*young_roots(3)
        stiffness(2,3) = (((s(1,2)*s(1,3)-s(1,1)*s(2,3))/determinant)*young_roots(2))*young_roots(3)
        stiffness(2,1) = stiffness(1,2)
        stiffness(3,1) = stiffness(1,3)
        stiffness(3,2) = stiffness(2,3)
        stiffness(4,4) = 2.0_wp*self%shear_23
        stiffness(5,5) = 2.0_wp*self%shear_13
        stiffness(6,6) = 2.0_wp*self%shear_12
        if (.not. all(ieee_is_finite(stiffness))) return
        status = cf_ok
    end subroutine build_stiffness

    subroutine stiffness_derivative(self, stiffness, parameter_id, derivative, status)
        !! Exact reciprocal-compliance chain rule for all nine parameters.
        !! Normalize the affected stiffness columns before forming outer products.
        class(orthotropic_elastic_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp), intent(in) :: stiffness(nmandel,nmandel)
            !! Elastic stiffness used to differentiate the compliance-based orthotropic law.
        integer, intent(in) :: parameter_id
            !! One-based scalar parameter identifier for differentiation.
        real(wp), intent(out) :: derivative(nmandel,nmandel)
            !! Exact derivative of the Mandel elastic stiffness matrix with respect to the requested parameter.
        integer, intent(out) :: status
            !! Operation status. `cf_ok` indicates a valid model or successful initialization.
        real(wp) :: column(ndim)
        real(wp) :: coupling(ndim)
        integer :: j

        status = cf_ok
        derivative = 0.0_wp
        select case (parameter_id)
        case (1)
            ! C(:,1)-nu12*C(:,2)-nu13*C(:,3)=E1*e1 follows from C*S=I.
            column = stiffness(1:3,1)/self%young_1
            do j = 1, ndim
                derivative(1:3,j) = -column*column(j)
            end do
            derivative(1:3,1) = derivative(1:3,1)+column
            derivative(1,1:3) = derivative(1,1:3)+column
        case (2)
            column = stiffness(1:3,2)/self%young_2
            coupling = (self%poisson_23*stiffness(1:3,3))/self%young_2
            do j = 1, ndim
                derivative(1:3,j) = column*(column(j)-coupling(j))-coupling*column(j)
            end do
        case (3)
            column = stiffness(1:3,3)/self%young_3
            do j = 1, ndim
                derivative(1:3,j) = column*column(j)
            end do
        case (4)
            column = stiffness(1:3,1)/self%young_1
            coupling = stiffness(1:3,2)
        case (5)
            column = stiffness(1:3,1)/self%young_1
            coupling = stiffness(1:3,3)
        case (6)
            column = stiffness(1:3,2)/self%young_2
            coupling = stiffness(1:3,3)
        case (7)
            derivative(4,4) = 2.0_wp
        case (8)
            derivative(5,5) = 2.0_wp
        case (9)
            derivative(6,6) = 2.0_wp
        case default
            status = cf_invalid_argument
            return
        end select
        if (parameter_id >= 4 .and. parameter_id <= 6) then
            do j = 1, ndim
                derivative(1:3,j) = column*coupling(j)+coupling*column(j)
            end do
        end if
        if (.not. all(ieee_is_finite(derivative))) status = cf_invalid_argument
    end subroutine stiffness_derivative

end module constituf_orthotropic_elastic
