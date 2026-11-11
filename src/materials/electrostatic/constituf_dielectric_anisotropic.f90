! SPDX-License-Identifier: Apache-2.0

module constituf_dielectric_anisotropic
    !! Homogeneous linear symmetric anisotropic electrostatic dielectric in 3D.
    !!
    !! The prescribed spatial electric field \(\boldsymbol{E}\) (electric potential per length) and
    !! displacement \(\boldsymbol{D}\) (charge per area) obey
    !! \[\boldsymbol{D}=\boldsymbol{\varepsilon}\boldsymbol{E},\quad
    !! \boldsymbol{\varepsilon}=
    !! \begin{bmatrix}
    !! \varepsilon_{11}&\varepsilon_{12}&\varepsilon_{13}\\[2pt]
    !! \varepsilon_{12}&\varepsilon_{22}&\varepsilon_{23}\\[2pt]
    !! \varepsilon_{13}&\varepsilon_{23}&\varepsilon_{33}
    !! \end{bmatrix}.\]
    !! The six independent coefficients are *absolute* permittivities,
    !! with dimensions of electric displacement per electric field.
    !! Relative permittivities are not accepted. The tensor must be finite
    !! and strictly symmetric positive definite in the declared spatial axes.
    !! In the absolute-permittivity convention, \(\boldsymbol\varepsilon=\varepsilon_0\boldsymbol\varepsilon_r\)
    !! includes the vacuum contribution \(\varepsilon_0\boldsymbol I\).
    !! The supplied absolute tensor is used directly, with no additional
    !! vacuum-permittivity multiplier.
    !! Positive definiteness guarantees positive stored field energy
    !! \[w(\boldsymbol{E})=\tfrac12\boldsymbol{E}^{\mathsf T}
    !! \boldsymbol{\varepsilon}\boldsymbol{E}>0
    !! \quad \text{for }\boldsymbol{E}\ne\boldsymbol{0}.\]
    !! This reversible energy has units of energy per volume and is not an output field.
    !! Symmetry gives \(\boldsymbol{D}=\partial w/\partial\boldsymbol{E}\) and makes the linear
    !! polarization work independent of the loading path.
    !! The tangent is the exact local derivative
    !! \(\partial\boldsymbol{D}/\partial\boldsymbol{E}
    !! =\boldsymbol{\varepsilon}\). For each named coefficient
    !! \(p=\varepsilon_{ab}\) the optional local parameter sensitivity is
    !! \[\frac{\partial D_i}{\partial p}=
    !! \begin{cases}
    !! \delta_{ia}E_a,&a=b,\\[2pt]
    !! \delta_{ia}E_b+\delta_{ib}E_a,&a\ne b.
    !! \end{cases}\]
    !! Coefficients are specified in tensor-index order \((11,22,33,12,13,23)\).
    !!
    !! Input `electric.field` and output `electric.displacement` are
    !! spatial 3-vectors with component order \([1,2,3]\). Only a caller-
    !! prescribed \(\boldsymbol{E}\) is evaluated. No unit or frame transformation occurs.
    !! To rotate spatial axes, the caller must rotate both \(\boldsymbol{E}\) and the
    !! permittivity tensor. The law is stateless, reversible and isothermal,
    !! with no `dt` requirement, ferroelectric history, nonlinear effects,
    !! currents, Maxwell stress or field-equation solver. The trial fails
    !! on nonfinite inputs/results. Near-singular tensors may fail the
    !! scaled Cholesky check. Sensitivities hold prescribed \(\boldsymbol{E}\) fixed and
    !! exclude global electrostatic-equilibrium derivatives.
    !!
    !! Scientific foundation (not imported source code): J. D. Jackson,
    !! *Classical Electrodynamics*, 3rd ed., Wiley (1998),
    !! ISBN 978-0-471-30932-1, macroscopic dielectrics, Chapter 4.
    !! Direct tensor and energy-symmetry derivation: A. Authier,
    !! "Introduction to the properties of tensors", *International Tables
    !! for Crystallography* (2006), Vol. D, Ch. 1.1, pp. 10--31,
    !! section 1.1.4.5.2.1, example (2):
    !! https://onlinelibrary.wiley.com/iucr/itc/Da/ch1o1v0001/sec1o1o4/ .
    !! This implementation uses Cartesian coordinates and positive tensor
    !! eigenvalues. Code provenance: the displacement and its analytical
    !! derivatives directly implement the stated law. No third-party
    !! software code is reused in this revision.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: ndim = 3
    integer, parameter :: nparameters = 6

    public :: dielectric_anisotropic_t

    type, extends(material_model_t) :: dielectric_anisotropic_t
        !! Create with `dielectric_anisotropic_t(e11,e22,e33,e12,e13,e23)`.
        !! The coefficients are in electric displacement per electric field. `material_t%init` clones the model.
        real(wp) :: permittivity_11 = 0.0_wp
            !! Spatial tensor entry \(\varepsilon_{11}\), electric displacement per electric field.
        real(wp) :: permittivity_22 = 0.0_wp
            !! Spatial tensor entry \(\varepsilon_{22}\), electric displacement per electric field.
        real(wp) :: permittivity_33 = 0.0_wp
            !! Spatial tensor entry \(\varepsilon_{33}\), electric displacement per electric field.
        real(wp) :: permittivity_12 = 0.0_wp
            !! Symmetric entry \(\varepsilon_{12}=\varepsilon_{21}\), electric displacement per electric field.
        real(wp) :: permittivity_13 = 0.0_wp
            !! Symmetric entry \(\varepsilon_{13}=\varepsilon_{31}\), electric displacement per electric field.
        real(wp) :: permittivity_23 = 0.0_wp
            !! Symmetric entry \(\varepsilon_{23}=\varepsilon_{32}\), electric displacement per electric field.
    contains
        procedure :: describe => describe_dielectric_anisotropic
            !! Declare input and output fields, parameter names and constitutive capabilities.
        procedure :: eval_batch => evaluate_dielectric_anisotropic
            !! Evaluate the batched constitutive response and optional consistent tangent.
        procedure :: eval_sensitivity_batch => sensitivity_dielectric_anisotropic
            !! Evaluate exact scalar-parameter derivatives of the constitutive response.
    end type dielectric_anisotropic_t

contains

    subroutine describe_dielectric_anisotropic(self, spec, status)
        !! Define tensor material with six selectable scalar sensitivities.
        class(dielectric_anisotropic_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        type(material_spec_t), intent(out) :: spec
            !! Packed field schema, material parameter metadata and declared capabilities.
        integer, intent(out) :: status
            !! Operation status; `cf_ok` indicates success.
        real(wp) :: permittivity(ndim,ndim)

        call build_tensor(self,permittivity,status)
        if (status /= cf_ok) return

        call spec%add_input("electric.field", ndim, "vector", "electric_field", "spatial", [ndim])
        call spec%add_output("electric.displacement", ndim, "vector", "electric_displacement", "spatial", [ndim])
        call spec%add_parameter("permittivity_11")
        call spec%add_parameter("permittivity_22")
        call spec%add_parameter("permittivity_33")
        call spec%add_parameter("permittivity_12")
        call spec%add_parameter("permittivity_13")
        call spec%add_parameter("permittivity_23")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_dielectric_anisotropic

    subroutine evaluate_dielectric_anisotropic(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Evaluate \(\boldsymbol{D}=\boldsymbol{\varepsilon}\boldsymbol{E}\) over a caller-owned batch, without state.
        class(dielectric_anisotropic_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed input `(ninput,npoints)` when needed, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `electric.field` input `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued material history `(nstate,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Trial real-valued material history `(nstate,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `electric.displacement` output `(noutput,npoints)`.
        real(wp), intent(in) :: dt
            !! Prescribed time increment; required to be positive only for time-dependent models.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic response Jacobian `(noutput,ninput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Status of each point; a successful update must set its code to `cf_ok`.

        state_trial = state_old
        call evaluate_points(self,x_new,response,point_status,tangent)
    end subroutine evaluate_dielectric_anisotropic

    subroutine sensitivity_dielectric_anisotropic(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Exact local derivatives at fixed field in requested parameter order.
        !! `dsdp` has shape `(3,nselected,npoints)`. History buffers are empty.
        class(dielectric_anisotropic_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed input `(ninput,npoints)` when needed, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `electric.field` input `(ninput,npoints)`.
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
            !! Packed `electric.displacement` output `(noutput,npoints)`.
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
        call evaluate_points(self,x_new,response,point_status,tangent, selected_parameters,dsdp)
    end subroutine sensitivity_dielectric_anisotropic

    subroutine evaluate_points(self, electric, response, point_status, tangent, selected_parameters, dsdp)
        !! Validate tensor once per batch. No point-local allocations or I/O.
        class(dielectric_anisotropic_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        real(wp), intent(in) :: electric(:,:)
            !! Prescribed electric field `electric.field` at each point.
        real(wp), intent(out) :: response(:,:)
            !! Packed `electric.displacement` output `(noutput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Status of each point; a successful update must set its code to `cf_ok`.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic response Jacobian `(noutput,ninput,npoints)`.
        integer, intent(in), optional :: selected_parameters(:)
            !! One-based IDs of scalar parameters in the requested sensitivity-axis order.
        real(wp), intent(out), optional :: dsdp(:,:,:)
            !! Exact discrete response derivatives `(noutput,nselected,npoints)` at fixed input histories.
        real(wp) :: permittivity(ndim,ndim)
        integer :: ipoint
        integer :: iparam
        integer :: id
        integer :: i
        integer :: status

        point_status = cf_invalid_argument
        call build_tensor(self,permittivity,status)
        if (status /= cf_ok) return
        if (present(selected_parameters)) then
            if (any(selected_parameters < 1) .or. any(selected_parameters > nparameters)) return
        end if

        point_loop: do ipoint = 1, size(electric,2)
            if (.not. all(ieee_is_finite(electric(:,ipoint)))) cycle point_loop
            do i = 1, ndim
                response(i,ipoint) = dot_product(permittivity(i,:),electric(:,ipoint))
                if (.not. ieee_is_finite(response(i,ipoint))) then
                    response(i,ipoint) = scaled_displacement_component(permittivity(i,:),electric(:,ipoint))
                end if
            end do
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop

            if (present(selected_parameters)) then
                do iparam = 1, size(selected_parameters)
                    id = selected_parameters(iparam)
                    dsdp(:,iparam,ipoint) = 0.0_wp
                    select case (id)
                    case (1:3)
                        dsdp(id,iparam,ipoint) = electric(id,ipoint)
                    case (4)
                        dsdp(1,iparam,ipoint) = electric(2,ipoint)
                        dsdp(2,iparam,ipoint) = electric(1,ipoint)
                    case (5)
                        dsdp(1,iparam,ipoint) = electric(3,ipoint)
                        dsdp(3,iparam,ipoint) = electric(1,ipoint)
                    case (6)
                        dsdp(2,iparam,ipoint) = electric(3,ipoint)
                        dsdp(3,iparam,ipoint) = electric(2,ipoint)
                    case default
                        cycle point_loop
                    end select
                end do
                if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
            end if
            if (present(tangent)) tangent(:,:,ipoint) = permittivity
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

    pure function scaled_displacement_component(permittivity_row, electric) result(displacement)
        !! Retry an overflowing dot product using bounded product fractions.
        !! Sum larger terms first so cancellation preserves small finite terms.
        real(wp), intent(in) :: permittivity_row(:)
            !! Finite coefficients of one spatial absolute-permittivity row.
        real(wp), intent(in) :: electric(ndim)
            !! Finite prescribed spatial electric field.
        real(wp) :: displacement
        real(wp) :: term_fraction(ndim)
        real(wp) :: sum_fraction
        integer :: term_exponent(ndim)
        integer :: sum_exponent
        integer :: common_exponent
        integer :: largest
        integer :: i
        integer :: j

        term_fraction = 0.0_wp
        term_exponent = 0
        product_loop: do i = 1, ndim
            if (abs(permittivity_row(i)) <= 0.0_wp .or. abs(electric(i)) <= 0.0_wp) cycle product_loop
            term_fraction(i) = fraction(permittivity_row(i))*fraction(electric(i))
            term_exponent(i) = exponent(permittivity_row(i))+exponent(electric(i))+exponent(term_fraction(i))
            term_fraction(i) = fraction(term_fraction(i))
        end do product_loop

        sum_fraction = 0.0_wp
        sum_exponent = 0
        sum_loop: do i = 1, ndim
            largest = 0
            largest_loop: do j = 1, ndim
                if (abs(term_fraction(j)) <= 0.0_wp) cycle largest_loop
                if (largest == 0) then
                    largest = j
                    cycle largest_loop
                end if
                if (term_exponent(j) > term_exponent(largest)) then
                    largest = j
                else if (term_exponent(j) == term_exponent(largest)) then
                    if (abs(term_fraction(j)) > abs(term_fraction(largest))) largest = j
                end if
            end do largest_loop
            if (largest == 0) exit sum_loop

            if (abs(sum_fraction) <= 0.0_wp) then
                sum_fraction = term_fraction(largest)
                sum_exponent = term_exponent(largest)
            else
                common_exponent = max(sum_exponent,term_exponent(largest))
                sum_fraction = scale(sum_fraction,sum_exponent-common_exponent)+scale(term_fraction(largest),term_exponent(largest)-common_exponent)
                sum_exponent = common_exponent+exponent(sum_fraction)
                sum_fraction = fraction(sum_fraction)
            end if
            term_fraction(largest) = 0.0_wp
        end do sum_loop
        displacement = scale(sum_fraction,sum_exponent)
    end function scaled_displacement_component

    subroutine build_tensor(self, permittivity, status)
        !! Check SPD by Cholesky pivots after positive coefficient scaling.
        !! No numerical clipping or regularization is applied.
        class(dielectric_anisotropic_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        real(wp), intent(out) :: permittivity(ndim,ndim)
            !! Symmetric spatial dielectric permittivity tensor assembled from model coefficients.
        integer, intent(out) :: status
            !! Operation status; `cf_ok` indicates success.
        real(wp) :: scaled(ndim,ndim)
        real(wp) :: scale
        real(wp) :: l11
        real(wp) :: l21
        real(wp) :: l31
        real(wp) :: l22
        real(wp) :: l32
        real(wp) :: pivot

        status = cf_invalid_argument
        permittivity(1,:) = [self%permittivity_11,self%permittivity_12, self%permittivity_13]
        permittivity(2,:) = [self%permittivity_12,self%permittivity_22, self%permittivity_23]
        permittivity(3,:) = [self%permittivity_13,self%permittivity_23, self%permittivity_33]
        if (.not. all(ieee_is_finite(permittivity))) return
        scale = maxval(abs(permittivity))
        if (scale <= 0.0_wp) return
        scaled = permittivity/scale

        pivot = scaled(1,1)
        if (pivot <= 0.0_wp) return
        l11 = sqrt(pivot)
        l21 = scaled(2,1)/l11
        l31 = scaled(3,1)/l11
        pivot = scaled(2,2)-l21*l21
        if (.not. ieee_is_finite(pivot)) return
        if (pivot <= 0.0_wp) return
        l22 = sqrt(pivot)
        l32 = (scaled(3,2)-l31*l21)/l22
        pivot = scaled(3,3)-l31*l31-l32*l32
        if (.not. ieee_is_finite(pivot)) return
        if (pivot <= 0.0_wp) return
        status = cf_ok
    end subroutine build_tensor

end module constituf_dielectric_anisotropic
