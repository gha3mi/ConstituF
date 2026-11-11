! SPDX-License-Identifier: Apache-2.0

module constituf_biot_poroelastic
    !! Three-dimensional, small-strain, isotropic Biot poroelastic response.
    !!
    !! The prescribed infinitesimal strain \(\boldsymbol{\varepsilon}\) is
    !! dimensionless (six spatial Mandel components). The prescribed pore
    !! pressure \(p\) is the current signed excess pressure relative to a
    !! fixed zero reference (stress units), positive in compression. Negative excess
    !! pressures are admitted by this linear specialization. Total Cauchy stress
    !! \(\boldsymbol{\sigma}\) is positive in tension (stress units). The incremental
    !! fluid content \(\zeta\) is dimensionless, positive for fluid added
    !! per unit reference bulk volume. This local constitutive specialization is
    !! \[\boldsymbol{\sigma}
    !! =K\,\mathrm{tr}(\boldsymbol{\varepsilon})\boldsymbol{I}
    !! +2G\,\mathrm{dev}(\boldsymbol{\varepsilon})
    !! -b p\boldsymbol{I},\]
    !! \[\zeta=b\,\mathrm{tr}(\boldsymbol{\varepsilon})+\frac{p}{M}.\]
    !! Here \(K>0\) and \(G>0\) are drained bulk and shear moduli (stress units).
    !! \(0\le b\le1\) is the dimensionless Biot coefficient and \(M>0\)
    !! is the Biot storage modulus (stress units). The strain and stress use the Mandel map
    !! \[\boldsymbol{A}_M=
    !! [A_{11},A_{22},A_{33},\sqrt{2}A_{23},\sqrt{2}A_{13},\sqrt{2}A_{12}]^T,\]
    !! where \(\boldsymbol{A}\) denotes either symmetric tensor.
    !!
    !! Inputs are `strain.small` (six dimensionless components) and `pressure.pore`
    !! (one component, stress units), in that packed order. Outputs are
    !! `stress.total` (six components, stress units) and `fluid.content`
    !! (one dimensionless component). Thus the exact multi-field local Jacobian
    !! has total shape \(7\times7\):
    !! \[\frac{\partial\boldsymbol{\sigma}}
    !! {\partial\boldsymbol{\varepsilon}}
    !! =K\boldsymbol{I}\otimes\boldsymbol{I}
    !! +2G\mathbb{I}_{\mathrm{dev}},\qquad
    !! \frac{\partial\boldsymbol{\sigma}}{\partial p}
    !! =-b\boldsymbol{I},\]
    !! \[\frac{\partial\zeta}{\partial\boldsymbol{\varepsilon}}
    !! =b\boldsymbol{I},\qquad
    !! \frac{\partial\zeta}{\partial p}=\frac1M.\]
    !! Here \(\boldsymbol{I}\) is the second-order identity tensor and
    !! \(\mathbb{I}_{\mathrm{dev}}\) is the symmetric deviatoric
    !! fourth-order projector. The two cross-blocks have distinct signs
    !! because total stress and fluid content use the declared sign
    !! conventions. The tangent is a
    !! local constitutive derivative, not a global FEM assembly or
    !! time-discrete diffusion operator.
    !!
    !! With \(e=\mathrm{tr}(\boldsymbol{\varepsilon})\), the mixed strain-pressure
    !! potential density (energy per volume) is
    !! \[\mathcal G(\boldsymbol{\varepsilon},p)=\tfrac12 K e^2+
    !! G\,\mathrm{dev}(\boldsymbol{\varepsilon}):\mathrm{dev}(\boldsymbol{\varepsilon})
    !! -bpe-\frac{p^2}{2M}.\]
    !! It gives \(\boldsymbol{\sigma}=\partial_{\boldsymbol{\varepsilon}}\mathcal G\) and
    !! \(\zeta=-\partial_p\mathcal G\), explaining the cross-block signs.
    !! The corresponding stored energy with independent fluid content is
    !! \[\Psi(\boldsymbol{\varepsilon},\zeta)=\tfrac12 K e^2+
    !! G\,\mathrm{dev}(\boldsymbol{\varepsilon}):\mathrm{dev}(\boldsymbol{\varepsilon})
    !! +\tfrac12 M(\zeta-be)^2.\]
    !! Positive \(K,G,M\) make this energy positive definite. The interval
    !! \(b\in[0,1]\) is the conventional physical subclass selected here.
    !! Neither potential is a returned field. Drained loading has \(p=0\).
    !! The undrained condition \(\zeta=0\) gives \(p=-bMe\) and
    !! the bulk modulus \(K_u=K+b^2M\).
    !!
    !! The model is stateless: zero strain and zero pore pressure yield zero
    !! stress and fluid content. No separately prescribed prestress, fluid
    !! reference offset, Darcy flux, or mass-balance time integration is
    !! implemented. Its analytical derivatives with respect to the selected
    !! scalar parameters \(K\), \(G\), \(b\) and \(M\) hold both input fields fixed.
    !! `dsdp` follows the user's selection order at initialization. Nonfinite
    !! parameters, trial inputs, responses or requested derivatives fail.
    !! The host solver owns equilibrium
    !! and conservation equations.
    !!
    !! Reference: M. A. Biot, *General Theory of Three-Dimensional
    !! Consolidation*, Journal of Applied Physics **12** (1941), 155--164,
    !! DOI: 10.1063/1.1712886. Section 2, Eqs. (2.11)-(2.12), p. 158,
    !! give the coupled isotropic relations. Biot's \(\alpha\) and \(Q\) map to
    !! the present \(b\) and \(M\), respectively. His fluid-content variation
    !! \(\theta\) maps to \(\zeta\). The mixed potentials and undrained limit above
    !! follow from these equations and the stated sign conventions.
    !! No third-party software code is reused in this revision.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: nmandel = 6
    integer, parameter :: nfields = 7

    public :: biot_poroelastic_t

    type, extends(material_model_t) :: biot_poroelastic_t
        !! Construct as `biot_poroelastic_t(K,G,b,M)` with moduli in the chosen stress unit.
        real(wp) :: bulk_modulus = 0.0_wp
            !! Drained bulk modulus \(K>0\), in the chosen stress unit.
        real(wp) :: shear_modulus = 0.0_wp
            !! Drained shear modulus \(G>0\), in the chosen stress unit.
        real(wp) :: biot_coefficient = 0.0_wp
            !! Biot coefficient \(b\in[0,1]\).
        real(wp) :: biot_modulus = 0.0_wp
            !! Biot storage modulus \(M>0\), in the chosen stress unit.
    contains
        procedure :: describe => describe_biot
            !! Describe the material input, response, state and scalar parameter schema.
        procedure :: eval_batch => evaluate_biot
            !! Evaluate the batched material response and trial history.
        procedure :: eval_sensitivity_batch => sensitivity_biot
            !! Evaluate exact selected-parameter sensitivities through the discrete material update.
    end type biot_poroelastic_t

contains

    subroutine describe_biot(self, spec, status)
        !! Validate model parameters and register `strain.small` and `pressure.pore` input and `stress.total` and `fluid.content` response fields.
        class(biot_poroelastic_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema receiving the packed field descriptions and supported capabilities.
        integer, intent(out) :: status
            !! Status code; `cf_ok` indicates success.

        status = cf_invalid_argument
        if (.not. all(ieee_is_finite([self%bulk_modulus,self%shear_modulus, self%biot_coefficient,self%biot_modulus]))) return
        if (self%bulk_modulus <= 0.0_wp .or. self%shear_modulus <= 0.0_wp) return
        if (self%biot_modulus <= 0.0_wp) return
        if (self%biot_coefficient < 0.0_wp .or. self%biot_coefficient > 1.0_wp) return
        if (.not. ieee_is_finite(self%bulk_modulus+(4.0_wp/3.0_wp)*self%shear_modulus)) return

        call spec%add_input("strain.small", nmandel, "mandel", "small_strain", "spatial", [3,3])
        call spec%add_input("pressure.pore", 1, "scalar", "pore_pressure", "material")
        call spec%add_output("stress.total", nmandel, "mandel", "cauchy_stress", "spatial", [3,3])
        call spec%add_output("fluid.content", 1, "scalar", "fluid_content", "material")
        call spec%add_parameter("bulk_modulus")
        call spec%add_parameter("shear_modulus")
        call spec%add_parameter("biot_coefficient")
        call spec%add_parameter("biot_modulus")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_biot

    subroutine evaluate_biot(self,x_old,x_new,state_old,state_trial,response,dt,tangent,point_status)
        !! Batched stateless coupled evaluation. Output shape `(7,npoints)`.
        class(biot_poroelastic_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed inputs `(ninput,npoints)` when required, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `strain.small` and `pressure.pore` inputs `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued history `(nstate,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Computed uncommitted real-valued history `(nstate,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.total` and `fluid.content` responses `(noutput,npoints)`.
        real(wp), intent(in) :: dt
            !! Prescribed time increment; a positive value is required for time-dependent laws.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic tangent `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(inout) :: point_status(:)
            !! Per-point outcome; write `cf_ok` only for a successful constitutive update.

        state_trial = state_old
        call evaluate_points(self,x_new,response,point_status,tangent)
    end subroutine evaluate_biot

    subroutine sensitivity_biot(self,x_old,x_new,state_old,state_trial,selected_parameters,sensitivity_old,sensitivity_trial,response,dsdp, dt,tangent,point_status)
        !! Exact local scalar-parameter derivatives, with zero history.
        class(biot_poroelastic_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed inputs `(ninput,npoints)` when required, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `strain.small` and `pressure.pore` inputs `(ninput,npoints)`.
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
            !! Packed `stress.total` and `fluid.content` responses `(noutput,npoints)`.
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
    end subroutine sensitivity_biot

    subroutine evaluate_points(self,x,response,point_status,tangent,selected_parameters,dsdp)
        !! Evaluate `stress.total` and `fluid.content` with its optional exact tangent and selected-parameter derivatives.
        class(biot_poroelastic_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp), intent(in) :: x(:,:)
            !! Current prescribed `strain.small` and `pressure.pore` input values at each point.
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.total` and `fluid.content` responses `(noutput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Per-point outcome; write `cf_ok` only for a successful constitutive update.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic tangent `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(in), optional :: selected_parameters(:)
            !! One-based parameter IDs in the caller's selected derivative-axis order.
        real(wp), intent(out), optional :: dsdp(:,:,:)
            !! Exact total response derivatives `(noutput,nselected,npoints)` at fixed prescribed driving histories.
        real(wp) :: dev(nmandel)
        real(wp) :: tr
        real(wp) :: pressure
        real(wp) :: vol
        real(wp) :: matrix(nfields,nfields)
        real(wp) :: fluid_fraction(2)
        integer :: fluid_exponent(2)
        integer :: common_exponent
        integer :: ipoint
        integer :: i
        integer :: k

        if (present(tangent)) then
            matrix = 0.0_wp
            matrix(1:3,1:3) = self%bulk_modulus - (2.0_wp/3.0_wp)*self%shear_modulus
            do i = 1,3
                matrix(i,i) = self%bulk_modulus+(4.0_wp/3.0_wp)*self%shear_modulus
            end do
            do i = 4,nmandel
                matrix(i,i) = 2.0_wp*self%shear_modulus
            end do
            matrix(1:3,nfields) = -self%biot_coefficient
            matrix(nfields,1:3) = self%biot_coefficient
            matrix(nfields,nfields) = 1.0_wp/self%biot_modulus
            if (.not. all(ieee_is_finite(matrix))) then
                point_status = cf_invalid_argument
                return
            end if
        end if

        point_loop: do ipoint = 1,size(x,2)
            point_status(ipoint) = cf_invalid_argument
            if (.not. all(ieee_is_finite(x(:,ipoint)))) cycle point_loop
            tr = sum(x(1:3,ipoint))
            pressure = x(nfields,ipoint)
            dev = x(1:nmandel,ipoint)
            ! Pairwise differences preserve an exact zero deviator in uniform dilation.
            dev(1) = (x(1,ipoint)-x(2,ipoint))/3.0_wp+(x(1,ipoint)-x(3,ipoint))/3.0_wp
            dev(2) = (x(2,ipoint)-x(1,ipoint))/3.0_wp+(x(2,ipoint)-x(3,ipoint))/3.0_wp
            dev(3) = (x(3,ipoint)-x(1,ipoint))/3.0_wp+(x(3,ipoint)-x(2,ipoint))/3.0_wp
            vol = self%bulk_modulus*tr-self%biot_coefficient*pressure
            ! Include the factor two before rounding a subnormal shear response.
            response(1:nmandel,ipoint) = scale(fraction(self%shear_modulus)*fraction(dev),exponent(self%shear_modulus)+exponent(dev)+1)
            response(1:3,ipoint) = response(1:3,ipoint)+vol
            response(nfields,ipoint) = self%biot_coefficient*tr+pressure/self%biot_modulus
            if (abs(response(nfields,ipoint))<tiny(1.0_wp)) then
                ! Combine subnormal contributions before their final rounding.
                fluid_fraction = [fraction(self%biot_coefficient)*fraction(tr),fraction(pressure)/fraction(self%biot_modulus)]
                fluid_exponent = [exponent(self%biot_coefficient)+exponent(tr),exponent(pressure)-exponent(self%biot_modulus)]
                if (any(abs(fluid_fraction)>0.0_wp)) then
                    common_exponent = maxval(fluid_exponent,mask=abs(fluid_fraction)>0.0_wp)
                    response(nfields,ipoint) = scale(sum(scale(fluid_fraction,fluid_exponent-common_exponent)),common_exponent)
                end if
            end if
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop
            if (present(tangent)) tangent(:,:,ipoint) = matrix

            if (present(selected_parameters)) then
                do k = 1,size(selected_parameters)
                    dsdp(:,k,ipoint) = 0.0_wp
                    select case(selected_parameters(k))
                    case(1)
                        dsdp(1:3,k,ipoint) = tr
                    case(2)
                        dsdp(1:nmandel,k,ipoint) = 2.0_wp*dev
                    case(3)
                        dsdp(1:3,k,ipoint) = -pressure
                        dsdp(nfields,k,ipoint) = tr
                    case(4)
                        ! Keep the original pressure: a subnormal p/M can lose
                        ! precision before division by a small M restores scale.
                        dsdp(nfields,k,ipoint) = scale(-fraction(pressure)/fraction(self%biot_modulus)**2, exponent(pressure)-2*exponent(self%biot_modulus))
                    case default
                        cycle point_loop
                    end select
                end do
                if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
            end if
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

end module constituf_biot_poroelastic
