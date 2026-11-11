! SPDX-License-Identifier: Apache-2.0

module constituf_piezoelectric_linear
    !! Isothermal, reversible, small-strain piezoelectric stress-charge law.
    !!
    !! The 3D spatial infinitesimal strain uses the six-component Mandel vector
    !! \[\boldsymbol{\varepsilon}_M=
    !! [\varepsilon_{11},\varepsilon_{22},\varepsilon_{33},\sqrt{2}\varepsilon_{23},
    !! \sqrt{2}\varepsilon_{13},\sqrt{2}\varepsilon_{12}]^T.\]
    !! Cauchy stress uses the same map. The electric field is
    !! \(\boldsymbol{E}=[E_1,E_2,E_3]^T\). Inputs are dimensionless `strain.small` and
    !! `electric.field` (electric potential per length), packed as nine
    !! components. Outputs are `stress.cauchy` (stress units),
    !! `electric.displacement` (charge per area) and
    !! `energy.electric_enthalpy` (energy per volume), in that order,
    !! ten components total. Tensile stresses are positive.
    !!
    !! A 6mm-style transversely isotropic *piezoelectric/dielectric*
    !! coupling matrix \(\boldsymbol{e}\) is combined with isotropic elastic moduli \(K,G\) at fixed
    !! electric field:
    !! \[\boldsymbol\sigma_M=\boldsymbol{C}^E\boldsymbol\varepsilon_M-\boldsymbol{e}^T\boldsymbol E,
    !! \qquad \boldsymbol D=\boldsymbol{e}\boldsymbol\varepsilon_M+
    !! \boldsymbol\kappa^S\boldsymbol E.\]
    !! The nonzero entries of the \(3\times6\) stress-charge matrix are
    !! \[e_{3,1}=e_{3,2}=e_{31},\quad e_{3,3}=e_{33},\quad
    !! e_{1,5}=e_{2,4}=\sqrt{2}e_{15}.\]
    !! The factors \(\sqrt{2}\) arise because the mechanical variables use
    !! Mandel stress and strain, *not* engineering Voigt shear.
    !! All components refer to fixed axes with the polarization/symmetry
    !! direction along the positive third axis. Reversing the polarization in these axes reverses
    !! all three signed \(\boldsymbol{e}\) coefficients and leaves the elastic/dielectric
    !! blocks unchanged. No material-axis rotation is performed here.
    !! The constant dielectric tensor is diagonal with positive
    !! transverse/axial entries \((\kappa_t,\kappa_l)\), with dimensions of
    !! electric displacement per electric field. These are absolute
    !! permittivities at fixed strain, including the vacuum contribution.
    !! This is an isotropic *elastic* specialization, not a full arbitrary
    !! 6mm anisotropic stiffness tensor.
    !!
    !! The coupled potential is the **electric enthalpy density**
    !! \[\mathcal H(\boldsymbol{\varepsilon}_M,\boldsymbol{E})=
    !! \tfrac12\boldsymbol{\varepsilon}_M^T\boldsymbol{C}^E\boldsymbol{\varepsilon}_M-
    !! \boldsymbol{E}^T\boldsymbol{e}\boldsymbol{\varepsilon}_M-\tfrac12\boldsymbol{E}^T\boldsymbol{\kappa}^S\boldsymbol{E},\]
    !! yielding \(\boldsymbol{\sigma}_M=\partial\mathcal H/\partial\boldsymbol{\varepsilon}_M\)
    !! and \(\boldsymbol{D}=-\partial\mathcal H/\partial\boldsymbol{E}\).
    !! The exact \(10\times9\) derivative matrix has cross blocks
    !! \(\partial\boldsymbol{\sigma}_M/\partial\boldsymbol{E}=-\boldsymbol{e}^T\) and
    !! \(\partial\boldsymbol{D}/\partial\boldsymbol{\varepsilon}_M=\boldsymbol{e}\), with a sign difference
    !! because \(\boldsymbol{D}\) is the *negative* enthalpy derivative.
    !! The enthalpy tangent row is \([\boldsymbol{\sigma}_M^T,-\boldsymbol{D}^T]\).
    !! This saddle-point energy is not a positive Helmholtz potential in
    !! independent \(\boldsymbol{E}\). Positivity of the elastic and dielectric blocks
    !! is the physically relevant restricted stability condition.
    !!
    !! Selected scalar-parameter sensitivities are exact at prescribed
    !! strain and field: `bulk_modulus` (stress units), `shear_modulus` (stress units),
    !! `e31`, `e33` and `e15` (charge per area), `transverse_permittivity` and
    !! `axial_permittivity` (electric displacement per electric field).
    !! The chosen mechanical and electrical units must be coherent, so
    !! charge per area also equals stress divided by electric field.
    !! \(K,G\) and both dielectric constants must be positive.
    !! The signed piezoelectric coefficients may be any finite
    !! values. No history or `dt` is needed. No hysteresis, ferroelectric
    !! switching, finite deformations, or field-dependent coefficients.
    !!
    !! Historical reference: ANSI/IEEE Std 176-1987,
    !! "IEEE Standard on Piezoelectricity", published 1988,
    !! DOI 10.1109/IEEESTD.1988.79638, Section 2.4, Eqs. (10)-(18)
    !! for electric enthalpy and the stress-charge law, Eqs. (23)-(25)
    !! and Table 2 for the compressed-component and engineering-shear
    !! conventions. https://standards.ieee.org/ieee/176/356/ .
    !! This standard has since been
    !! withdrawn. It is cited only for established stress-charge forms,
    !! not as a current compliance claim. The isotropic elastic block is
    !! the specialization used here. No third-party software code is reused
    !! in this revision.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: nstrain = 6
    integer, parameter :: nelectric = 3
    integer, parameter :: ninput = 9
    integer, parameter :: noutput = 10
    public :: piezoelectric_linear_t

    type, extends(material_model_t) :: piezoelectric_linear_t
        !! 3D linear stress-charge piezoelectric law with isotropic elasticity.
        !! Construct with \(K,G,e_{31},e_{33},e_{15},\kappa_t,\kappa_l\) in a coherent unit system.
        real(wp) :: bulk_modulus = 0.0_wp
            !! \(K>0\) in the chosen stress unit.
        real(wp) :: shear_modulus = 0.0_wp
            !! \(G>0\) in the chosen stress unit.
        real(wp) :: e31 = 0.0_wp
            !! Signed transverse extensional stress-charge coefficient, charge per area.
        real(wp) :: e33 = 0.0_wp
            !! Signed axial extensional stress-charge coefficient, charge per area.
        real(wp) :: e15 = 0.0_wp
            !! Signed transverse shear stress-charge coefficient, charge per area.
        real(wp) :: transverse_permittivity = 0.0_wp
            !! Positive clamped transverse permittivity, electric displacement per electric field.
        real(wp) :: axial_permittivity = 0.0_wp
            !! Positive clamped axial permittivity, electric displacement per electric field.
    contains
        procedure :: describe => describe_piezoelectric_linear
            !! Describe the material input, response, state and scalar parameter schema.
        procedure :: eval_batch => evaluate_piezoelectric_linear
            !! Evaluate the batched material response and trial history.
        procedure :: eval_sensitivity_batch => sensitivity_piezoelectric_linear
            !! Evaluate exact selected-parameter sensitivities through the discrete material update.
    end type piezoelectric_linear_t

contains

    subroutine describe_piezoelectric_linear(self, spec, status)
        !! Validate model parameters and register `strain.small` and `electric.field` input and `stress.cauchy` and `electric.displacement` and `energy.electric_enthalpy` response fields.
        class(piezoelectric_linear_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema receiving the packed field descriptions and supported capabilities.
        integer, intent(out) :: status
            !! Status code; `cf_ok` indicates success.

        status = cf_invalid_argument
        if (.not. all(ieee_is_finite([self%bulk_modulus,self%shear_modulus,self%e31,self%e33,self%e15, self%transverse_permittivity,self%axial_permittivity]))) return
        if (self%bulk_modulus <= 0.0_wp .or. self%shear_modulus <= 0.0_wp) return
        if (self%transverse_permittivity <= 0.0_wp .or. self%axial_permittivity <= 0.0_wp) return
        if (.not. ieee_is_finite(self%bulk_modulus+(4.0_wp/3.0_wp)*self%shear_modulus)) return

        call spec%add_input("strain.small", nstrain, "mandel", "small_strain", "spatial", [3,3])
        call spec%add_input("electric.field", nelectric, "vector", "electric_field", "spatial", [3])
        call spec%add_output("stress.cauchy", nstrain, "mandel", "cauchy_stress", "spatial", [3,3])
        call spec%add_output("electric.displacement", nelectric, "vector", "electric_displacement", "spatial", [3])
        call spec%add_output("energy.electric_enthalpy", 1, "scalar", "electric_enthalpy", "material")
        call spec%add_parameter("bulk_modulus")
        call spec%add_parameter("shear_modulus")
        call spec%add_parameter("e31")
        call spec%add_parameter("e33")
        call spec%add_parameter("e15")
        call spec%add_parameter("transverse_permittivity")
        call spec%add_parameter("axial_permittivity")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_piezoelectric_linear

    subroutine evaluate_piezoelectric_linear(self,x_old,x_new,state_old,state_trial,response,dt,tangent,point_status)
        !! Evaluate `stress.cauchy` and `electric.displacement` and `energy.electric_enthalpy` and update the trial material state for the prescribed inputs.
        class(piezoelectric_linear_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed inputs `(ninput,npoints)` when required, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `strain.small` and `electric.field` inputs `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued history `(nstate,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Computed uncommitted real-valued history `(nstate,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.cauchy` and `electric.displacement` and `energy.electric_enthalpy` responses `(noutput,npoints)`.
        real(wp), intent(in) :: dt
            !! Prescribed time increment; a positive value is required for time-dependent laws.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic tangent `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(inout) :: point_status(:)
            !! Per-point outcome; write `cf_ok` only for a successful constitutive update.

        state_trial = state_old
        call evaluate_points(self,x_new,response,point_status,tangent)
    end subroutine evaluate_piezoelectric_linear

    subroutine sensitivity_piezoelectric_linear(self,x_old,x_new,state_old,state_trial,selected_parameters, sensitivity_old,sensitivity_trial,response,dsdp,dt,tangent,point_status)
        !! Differentiate the discrete `stress.cauchy` and `electric.displacement` and `energy.electric_enthalpy` update through the selected material history.
        class(piezoelectric_linear_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed inputs `(ninput,npoints)` when required, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `strain.small` and `electric.field` inputs `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued history `(nstate,npoints)`.
        integer, intent(in) :: selected_parameters(:)
            !! One-based parameter IDs in the caller's selected derivative-axis order.
        real(wp), intent(in) :: sensitivity_old(:,:,:)
            !! Committed real-history parameter derivatives `(nstate,nselected,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Computed uncommitted real-valued history `(nstate,npoints)`.
        real(wp), intent(out) :: sensitivity_trial(:,:,:)
            !! Updated real-history parameter derivatives `(nstate,nselected,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.cauchy` and `electric.displacement` and `energy.electric_enthalpy` responses `(noutput,npoints)`.
        real(wp), intent(out) :: dsdp(:,:,:)
            !! Exact total response derivatives `(noutput,nselected,npoints)` at fixed prescribed driving histories.
        real(wp), intent(in) :: dt
            !! Prescribed time increment; a positive value is required for time-dependent laws.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic tangent `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(inout) :: point_status(:)
            !! Per-point outcome; write `cf_ok` only for a successful constitutive update.

        state_trial = state_old
        sensitivity_trial = sensitivity_old
        call evaluate_points(self,x_new,response,point_status,tangent,selected_parameters,dsdp)
    end subroutine sensitivity_piezoelectric_linear

    subroutine evaluate_points(self,x,response,point_status,tangent,selected_parameters,dsdp)
        !! Evaluate `stress.cauchy` and `electric.displacement` and `energy.electric_enthalpy` with its optional exact tangent and selected-parameter derivatives.
        class(piezoelectric_linear_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp), intent(in) :: x(:,:)
            !! Current prescribed `strain.small` and `electric.field` input values at each point.
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.cauchy` and `electric.displacement` and `energy.electric_enthalpy` responses `(noutput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Per-point outcome; write `cf_ok` only for a successful constitutive update.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic tangent `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(in), optional :: selected_parameters(:)
            !! One-based parameter IDs in the caller's selected derivative-axis order.
        real(wp), intent(out), optional :: dsdp(:,:,:)
            !! Exact total response derivatives `(noutput,nselected,npoints)` at fixed prescribed driving histories.
        real(wp), parameter :: sqrt_two = sqrt(2.0_wp)
        real(wp) :: elastic(nstrain,nstrain)
        real(wp) :: coupling(nelectric,nstrain)
        real(wp) :: permittivity(nelectric,nelectric)
        real(wp) :: strain(nstrain)
        real(wp) :: deviatoric_strain(nstrain)
        real(wp) :: electric_field(nelectric)
        real(wp) :: stress(nstrain)
        real(wp) :: piezoelectric_displacement(nelectric)
        real(wp) :: dielectric_displacement(nelectric)
        real(wp) :: electric_displacement(nelectric)
        real(wp) :: trace_strain
        real(wp) :: deviatoric_scale
        real(wp) :: scaled_norm_squared
        real(wp) :: deviatoric_energy
        real(wp) :: dielectric_energy
        real(wp) :: dielectric_energy_fraction(nelectric)
        real(wp) :: lambda
        real(wp) :: energy
        integer :: dielectric_energy_exponent(nelectric)
        integer :: common_exponent
        integer :: ipoint
        integer :: i
        integer :: k

        if (present(tangent)) then
            lambda = self%bulk_modulus-(2.0_wp/3.0_wp)*self%shear_modulus
            elastic = 0.0_wp
            elastic(1:3,1:3) = lambda
            do i = 1,3
                elastic(i,i) = self%bulk_modulus+(4.0_wp/3.0_wp)*self%shear_modulus
            end do
            do i = 4,nstrain
                elastic(i,i) = 2.0_wp*self%shear_modulus
            end do
            coupling = 0.0_wp
            coupling(1,5) = sqrt_two*self%e15
            coupling(2,4) = sqrt_two*self%e15
            coupling(3,1:2) = self%e31
            coupling(3,3) = self%e33
            permittivity = 0.0_wp
            permittivity(1,1) = self%transverse_permittivity
            permittivity(2,2) = self%transverse_permittivity
            permittivity(3,3) = self%axial_permittivity
        end if

        point_status = cf_invalid_argument
        point_loop: do ipoint = 1,size(x,2)
            if (.not. all(ieee_is_finite(x(:,ipoint)))) cycle point_loop
            strain = x(1:nstrain,ipoint)
            electric_field = x(nstrain+1:ninput,ipoint)
            trace_strain = sum(strain(1:3))
            ! Pairwise differences preserve a zero deviator for uniform dilation.
            deviatoric_strain = strain
            deviatoric_strain(1) = (strain(1)-strain(2))/3.0_wp+(strain(1)-strain(3))/3.0_wp
            deviatoric_strain(2) = (strain(2)-strain(1))/3.0_wp+(strain(2)-strain(3))/3.0_wp
            deviatoric_strain(3) = (strain(3)-strain(1))/3.0_wp+(strain(3)-strain(2))/3.0_wp
            stress = scaled_product(self%shear_modulus,deviatoric_strain,2.0_wp)
            stress(1:3) = stress(1:3)+self%bulk_modulus*trace_strain
            stress(1:2) = stress(1:2)-self%e31*electric_field(3)
            stress(3) = stress(3)-self%e33*electric_field(3)
            stress(4) = stress(4)-scaled_product(self%e15,electric_field(2),sqrt_two)
            stress(5) = stress(5)-scaled_product(self%e15,electric_field(1),sqrt_two)
            ! Retain the Mandel factor through the coefficient/input product,
            ! avoiding overflow or premature subnormal rounding.
            piezoelectric_displacement(1) = scaled_product(self%e15,strain(5),sqrt_two)
            piezoelectric_displacement(2) = scaled_product(self%e15,strain(4),sqrt_two)
            piezoelectric_displacement(3) = self%e31*strain(1)+self%e31*strain(2)+self%e33*strain(3)
            dielectric_displacement(1:2) = self%transverse_permittivity*electric_field(1:2)
            dielectric_displacement(3) = self%axial_permittivity*electric_field(3)
            electric_displacement = piezoelectric_displacement+dielectric_displacement
            deviatoric_energy = 0.0_wp
            deviatoric_scale = maxval(abs(deviatoric_strain))
            if (deviatoric_scale > 0.0_wp) then
                scaled_norm_squared = sum((deviatoric_strain/deviatoric_scale)**2)
                deviatoric_energy = scaled_product(self%shear_modulus,deviatoric_scale,deviatoric_scale,scaled_norm_squared)
            end if
            ! Sum the dielectric quadratic form before restoring its exponent;
            ! rounded displacement components can lose subnormal energy terms.
            dielectric_energy = 0.0_wp
            if (any(abs(electric_field) > 0.0_wp)) then
                dielectric_energy_fraction(1:2) = 0.5_wp*fraction(self%transverse_permittivity)*fraction(electric_field(1:2))**2
                dielectric_energy_fraction(3) = 0.5_wp*fraction(self%axial_permittivity)*fraction(electric_field(3))**2
                dielectric_energy_exponent(1:2) = exponent(self%transverse_permittivity)+2*exponent(electric_field(1:2))
                dielectric_energy_exponent(3) = exponent(self%axial_permittivity)+2*exponent(electric_field(3))
                common_exponent = maxval(dielectric_energy_exponent,mask=abs(electric_field) > 0.0_wp)
                dielectric_energy = scale(sum(scale(dielectric_energy_fraction,dielectric_energy_exponent-common_exponent)),common_exponent)
            end if
            energy = scaled_product(self%bulk_modulus,trace_strain,trace_strain,0.5_wp)+deviatoric_energy - dot_product(electric_field,piezoelectric_displacement) - dielectric_energy
            response(1:nstrain,ipoint) = stress
            response(nstrain+1:ninput,ipoint) = electric_displacement
            response(noutput,ipoint) = energy
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop
            if (present(tangent)) then
                tangent(:,:,ipoint) = 0.0_wp
                tangent(1:nstrain,1:nstrain,ipoint) = elastic
                tangent(1:nstrain,nstrain+1:ninput,ipoint) = -transpose(coupling)
                tangent(nstrain+1:ninput,1:nstrain,ipoint) = coupling
                tangent(nstrain+1:ninput,nstrain+1:ninput,ipoint) = permittivity
                tangent(noutput,1:nstrain,ipoint) = stress
                tangent(noutput,nstrain+1:ninput,ipoint) = -electric_displacement
                if (.not. all(ieee_is_finite(tangent(:,:,ipoint)))) cycle point_loop
            end if

            if (present(selected_parameters)) then
                parameter_loop: do k = 1,size(selected_parameters)
                    dsdp(:,k,ipoint) = 0.0_wp
                    select case (selected_parameters(k))
                    case (1)
                        dsdp(1:3,k,ipoint) = trace_strain
                        dsdp(noutput,k,ipoint) = (0.5_wp*trace_strain)*trace_strain
                    case (2)
                        dsdp(1:nstrain,k,ipoint) = 2.0_wp*deviatoric_strain
                        dsdp(noutput,k,ipoint) = dot_product(deviatoric_strain,deviatoric_strain)
                    case (3)
                        dsdp(1:2,k,ipoint) = -electric_field(3)
                        dsdp(9,k,ipoint) = sum(strain(1:2))
                        dsdp(noutput,k,ipoint) = -electric_field(3)*sum(strain(1:2))
                    case (4)
                        dsdp(3,k,ipoint) = -electric_field(3)
                        dsdp(9,k,ipoint) = strain(3)
                        dsdp(noutput,k,ipoint) = -electric_field(3)*strain(3)
                    case (5)
                        dsdp(4,k,ipoint) = -sqrt_two*electric_field(2)
                        dsdp(5,k,ipoint) = -sqrt_two*electric_field(1)
                        dsdp(7,k,ipoint) = sqrt_two*strain(5)
                        dsdp(8,k,ipoint) = sqrt_two*strain(4)
                        dsdp(noutput,k,ipoint) = -sqrt_two*(strain(5)*electric_field(1)+strain(4)*electric_field(2))
                    case (6)
                        dsdp(7:8,k,ipoint) = electric_field(1:2)
                        dsdp(noutput,k,ipoint) = -dot_product(0.5_wp*electric_field(1:2),electric_field(1:2))
                    case (7)
                        dsdp(9,k,ipoint) = electric_field(3)
                        dsdp(noutput,k,ipoint) = -(0.5_wp*electric_field(3))*electric_field(3)
                    case default
                        exit parameter_loop
                    end select
                end do parameter_loop
                if (k <= size(selected_parameters)) cycle point_loop
                if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
            end if
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

    pure elemental function scaled_product(coefficient, component, multiplier, weight) result(product)
        !! Form three-factor products with one final exponent restoration.
        !! A bounded optional weight preserves weighted norms before rounding.
        real(wp), intent(in) :: coefficient
            !! Finite elastic or piezoelectric coefficient.
        real(wp), intent(in) :: component
            !! Finite prescribed strain or electric-field component.
        real(wp), intent(in) :: multiplier
            !! Second input component or a Mandel or shear factor.
        real(wp), intent(in), optional :: weight
            !! Bounded normalized squared norm or the volumetric one-half factor.
        real(wp) :: product
        real(wp) :: product_fraction

        product_fraction = (fraction(coefficient)*fraction(component))*fraction(multiplier)
        if (present(weight)) product_fraction = product_fraction*weight
        product = scale(product_fraction,exponent(coefficient)+exponent(component)+exponent(multiplier))
    end function scaled_product

end module constituf_piezoelectric_linear
