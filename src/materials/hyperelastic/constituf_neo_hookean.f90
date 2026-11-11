! SPDX-License-Identifier: Apache-2.0

module constituf_neo_hookean
    !! Compressible three-dimensional neo-Hookean hyperelasticity.
    !!
    !! The deformation gradient \(\boldsymbol{F}\) maps reference to spatial
    !! coordinates, \(J=\det\boldsymbol{F}>0\) and the strain energy per
    !! reference volume is
    !! \[\Psi(\boldsymbol{F}) =
    !! \frac{\mu}{2}(\boldsymbol{F}:\boldsymbol{F}-3-2\ln J)
    !! +\frac{\lambda}{2}(\ln J)^2.\]
    !! The first Piola--Kirchhoff stress, work-conjugate to \(\boldsymbol{F}\)
    !! and the material tangent are
    !! \[\boldsymbol{P}=\mu(\boldsymbol{F}-\boldsymbol{F}^{-T})
    !! +\lambda\ln J\,\boldsymbol{F}^{-T},\]
    !! \[\frac{\partial P_{iJ}}{\partial F_{kL}}=
    !! \mu\delta_{ik}\delta_{JL}
    !! +\lambda(F^{-T})_{iJ}(F^{-T})_{kL}
    !! +(\mu-\lambda\ln J)(F^{-T})_{iL}(F^{-T})_{kJ}.\]
    !! Here \(\mu>0\) and \(\lambda\geq 0\) are Lamé moduli in the chosen stress unit.
    !! The restricted nonnegative \(\lambda\) domain is deliberate. Materials
    !! have infinitesimal Poisson ratio
    !! \(\nu=\lambda/[2(\lambda+\mu)]\), so \(0\leq\nu<1/2\).
    !! Auxetic calibrations are outside this implementation.
    !! The corresponding infinitesimal bulk modulus is
    !! \(K=\lambda+2\mu/3\).
    !!
    !! Inputs and outputs are full \(3\times3\) tensors in Fortran column-major component-index order
    !! \((11,21,31,12,22,32,13,23,33)\), not Mandel or Voigt.
    !! Both \(\boldsymbol{F}\) and \(\boldsymbol{P}\) are two-point tensors
    !! with spatial-reference component frames. Positive stress is tensile.
    !! The model is isothermal, rate-independent and stateless. No energy
    !! or Cauchy stress is returned by this interface. Invalid or inverted
    !! deformations, nonfinite inputs and nonfinite results fail by point
    !! status without any accepted history mutation.
    !!
    !! Optional local, exact parameter sensitivities at fixed deformation:
    !! \[\frac{\partial\boldsymbol{P}}{\partial\lambda}
    !! =\ln J\,\boldsymbol{F}^{-T},\qquad
    !! \frac{\partial\boldsymbol{P}}{\partial\mu}
    !! =\boldsymbol{F}-\boldsymbol{F}^{-T}.\]
    !! Neither these nor the tangent include global equilibrium derivatives.
    !!
    !! References: J. Bonet and R. D. Wood, *Nonlinear Continuum Mechanics
    !! for Finite Element Analysis*, 2nd ed., Cambridge University Press,
    !! 2008, DOI: 10.1017/CBO9780511755446.
    !! FEBio theory manual, "Neo-Hookean Hyperelasticity", sec. 5.2,
    !! energy and Cauchy stress equations (422)--(423):
    !! https://febiosoftware.github.io/febio-docs/theory/chapter5/5.2-compressible-materials/#neo-hookean-hyperelasticity.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: ndim = 3
    integer, parameter :: nfull = ndim*ndim

    public :: neo_hookean_t

    type, extends(material_model_t) :: neo_hookean_t
        !! Stateless 3D finite-strain constitutive law.
        !!
        !! Construct as `neo_hookean_t(lame_lambda, shear_modulus)`, both in the chosen stress unit,
        !! then pass to `material_t%init`. The material handle clones the model.
        !! The author does not own or modify integration-point state.
        real(wp) :: lame_lambda = 0.0_wp
            !! Nonnegative first Lamé modulus \(\lambda\), in the chosen stress unit.
        real(wp) :: shear_modulus = 0.0_wp
            !! Positive shear modulus \(\mu\), in the chosen stress unit.
    contains
        procedure :: describe => describe_neo_hookean
            !! Declare input and output fields, parameter names and constitutive capabilities.
        procedure :: eval_batch => evaluate_neo_hookean
            !! Evaluate the batched constitutive response and optional consistent tangent.
        procedure :: eval_sensitivity_batch => sensitivity_neo_hookean
            !! Evaluate exact scalar-parameter derivatives of the constitutive response.
    end type neo_hookean_t

contains

    subroutine describe_neo_hookean(self, spec, status)
        !! Declare full tensor \(\boldsymbol{F}\) input and first Piola stress output.
        !!
        !! Fields each have 9 stored components and physical shape `[3,3]`.
        !! The tangent is \(\partial\boldsymbol{P}/\partial\boldsymbol{F}\) in the same flattened order. Named scalar
        !! parameters `lame_lambda` and `shear_modulus` can be selected at init.
        class(neo_hookean_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        type(material_spec_t), intent(out) :: spec
            !! Packed field schema, material parameter metadata and declared capabilities.
        integer, intent(out) :: status
            !! Operation status; `cf_ok` indicates success.

        status = cf_invalid_argument
        if (.not. ieee_is_finite(self%lame_lambda)) return
        if (.not. ieee_is_finite(self%shear_modulus)) return
        if (self%lame_lambda < 0.0_wp .or. self%shear_modulus <= 0.0_wp) return
        if (.not. ieee_is_finite(2.0_wp*self%shear_modulus)) return

        call spec%add_input("deformation.f", nfull, "tensor_full", "deformation_gradient", "spatial_reference", [ndim,ndim])
        call spec%add_output("stress.piola", nfull, "tensor_full", "first_piola_kirchhoff", "spatial_reference", [ndim,ndim])
        call spec%add_parameter("lame_lambda")
        call spec%add_parameter("shear_modulus")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_neo_hookean

    subroutine evaluate_neo_hookean(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Evaluate \(\boldsymbol{P}(\boldsymbol{F})\) without modifying any committed state.
        !!
        !! `x_new(9,npoints)`, `response(9,npoints)`, `tangent(9,9,npoints)`.
        !! Failure invalidates all caller outputs, even successful points.
        class(neo_hookean_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed input `(ninput,npoints)` when needed, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `deformation.f` input `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued material history `(nstate,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Trial real-valued material history `(nstate,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.piola` output `(noutput,npoints)`.
        real(wp), intent(in) :: dt
            !! Prescribed time increment; required to be positive only for time-dependent models.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic response Jacobian `(noutput,ninput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Status of each point; a successful update must set its code to `cf_ok`.

        state_trial = state_old
        call evaluate_points(self, x_new, response, point_status, tangent)
    end subroutine evaluate_neo_hookean

    subroutine sensitivity_neo_hookean(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Exact local derivatives of first Piola stress with respect to
        !! selected Lamé parameters in the original selection order.
        !!
        !! `dsdp(9,nselected,npoints)` is dimensionless. Parameter IDs are
        !! 1 for `lame_lambda` and 2 for `shear_modulus`. Both sensitivities
        !! are local at fixed prescribed deformation history. Stateless
        !! sensitivity buffers have zero history components.
        class(neo_hookean_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed input `(ninput,npoints)` when needed, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `deformation.f` input `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued material history `(nstate,npoints)`.
        integer, intent(in) :: selected_parameters(:)
            !! One-based IDs of scalar parameters in the requested sensitivity-axis order.
        real(wp), intent(in) :: sensitivity_old(:,:,:)
            !! Committed derivatives of real material history `(nstate,nselected,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Trial real-valued material history `(nstate,npoints)`.
        real(wp), intent(out) :: sensitivity_trial(:,:,:)
            !! Trial derivatives of real material history `(nstate,nselected,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.piola` output `(noutput,npoints)`.
        real(wp), intent(out) :: dsdp(:,:,:)
            !! Exact discrete response derivatives `(noutput,nselected,npoints)` at fixed input histories.
        real(wp), intent(in) :: dt
            !! Prescribed time increment; required to be positive only for time-dependent models.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic response Jacobian `(noutput,ninput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Status of each point; a successful update must set its code to `cf_ok`.

        state_trial = state_old
        sensitivity_trial = sensitivity_old
        call evaluate_points(self, x_new, response, point_status, tangent, selected_parameters, dsdp)
    end subroutine sensitivity_neo_hookean

    subroutine evaluate_points(self, x_new, response, point_status, tangent, selected_parameters, dsdp)
        !! One model dispatch per batch, no heap allocation in the hot loop.
        class(neo_hookean_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `deformation.f` input `(ninput,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.piola` output `(noutput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Status of each point; a successful update must set its code to `cf_ok`.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic response Jacobian `(noutput,ninput,npoints)`.
        integer, intent(in), optional :: selected_parameters(:)
            !! One-based IDs of scalar parameters in the requested sensitivity-axis order.
        real(wp), intent(out), optional :: dsdp(:,:,:)
            !! Exact discrete response derivatives `(noutput,nselected,npoints)` at fixed input histories.
        real(wp) :: deformation(ndim,ndim)
        real(wp) :: cofactor(ndim,ndim)
        real(wp) :: inverse_t(ndim,ndim)
        real(wp) :: stress(ndim,ndim)
        real(wp) :: jacobian
        real(wp) :: log_jacobian
        real(wp) :: value
        integer :: ipoint
        integer :: i
        integer :: j
        integer :: k
        integer :: l
        integer :: row
        integer :: column

        point_loop: do ipoint = 1, size(x_new,2)
            point_status(ipoint) = cf_invalid_argument
            if (.not. all(ieee_is_finite(x_new(:,ipoint)))) cycle point_loop
            do j = 1, ndim
                deformation(:,j) = x_new(ndim*(j-1)+1:ndim*j,ipoint)
            end do
            call cofactor_matrix(deformation, cofactor)
            jacobian = sum(deformation(1,:)*cofactor(1,:))
            if (.not. ieee_is_finite(jacobian)) cycle point_loop
            if (jacobian <= 0.0_wp) cycle point_loop
            inverse_t = cofactor/jacobian
            if (.not. all(ieee_is_finite(inverse_t))) cycle point_loop
            log_jacobian = log(jacobian)
            stress = self%shear_modulus*(deformation-inverse_t) + self%lame_lambda*log_jacobian*inverse_t
            if (.not. all(ieee_is_finite(stress))) cycle point_loop

            if (present(tangent)) then
                do j = 1, ndim
                    do i = 1, ndim
                        row = i + ndim*(j-1)
                        do l = 1, ndim
                            do k = 1, ndim
                                column = k + ndim*(l-1)
                                value = self%lame_lambda*inverse_t(i,j)*inverse_t(k,l) + (self%shear_modulus - self%lame_lambda*log_jacobian)* inverse_t(i,l)*inverse_t(k,j)
                                if (i == k .and. j == l) value = value + self%shear_modulus
                                tangent(row,column,ipoint) = value
                            end do
                        end do
                    end do
                end do
                if (.not. all(ieee_is_finite(tangent(:,:,ipoint)))) cycle point_loop
            end if

            if (present(selected_parameters)) then
                parameter_loop: do k = 1, size(selected_parameters)
                    do j = 1, ndim
                        row = ndim*(j-1)
                        select case (selected_parameters(k))
                        case (1)
                            dsdp(row+1:row+ndim,k,ipoint) = log_jacobian*inverse_t(:,j)
                        case (2)
                            dsdp(row+1:row+ndim,k,ipoint) = deformation(:,j)-inverse_t(:,j)
                        case default
                            exit parameter_loop
                        end select
                    end do
                end do parameter_loop
                if (k <= size(selected_parameters)) cycle point_loop
                if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
            end if

            do j = 1, ndim
                response(ndim*(j-1)+1:ndim*j,ipoint) = stress(:,j)
            end do
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

    pure subroutine cofactor_matrix(deformation, cofactor)
        !! Cofactors are \((\det\boldsymbol{F})\boldsymbol{F}^{-T}\). Explicit unrolled 3D arithmetic.
        real(wp), intent(in) :: deformation(ndim,ndim)
            !! Deformation gradient \(\boldsymbol{F}\) in full three-dimensional tensor form.
        real(wp), intent(out) :: cofactor(ndim,ndim)
            !! Cofactor matrix \(\operatorname{cof}\boldsymbol{F}=J\boldsymbol{F}^{-T}\).

        cofactor(1,1) = deformation(2,2)*deformation(3,3) - deformation(2,3)*deformation(3,2)
        cofactor(1,2) = deformation(2,3)*deformation(3,1) - deformation(2,1)*deformation(3,3)
        cofactor(1,3) = deformation(2,1)*deformation(3,2) - deformation(2,2)*deformation(3,1)
        cofactor(2,1) = deformation(1,3)*deformation(3,2) - deformation(1,2)*deformation(3,3)
        cofactor(2,2) = deformation(1,1)*deformation(3,3) - deformation(1,3)*deformation(3,1)
        cofactor(2,3) = deformation(1,2)*deformation(3,1) - deformation(1,1)*deformation(3,2)
        cofactor(3,1) = deformation(1,2)*deformation(2,3) - deformation(1,3)*deformation(2,2)
        cofactor(3,2) = deformation(1,3)*deformation(2,1) - deformation(1,1)*deformation(2,3)
        cofactor(3,3) = deformation(1,1)*deformation(2,2) - deformation(1,2)*deformation(2,1)
    end subroutine cofactor_matrix

end module constituf_neo_hookean
