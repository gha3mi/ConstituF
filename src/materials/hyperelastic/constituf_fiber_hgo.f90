! SPDX-License-Identifier: Apache-2.0

module constituf_fiber_hgo
    !! Two-family fiber-reinforced compressible hyperelastic solid.
    !!
    !! For the deformation gradient \(\boldsymbol{F}\), \(J=\det\boldsymbol{F}>0\) and unit reference
    !! fiber directions \(\boldsymbol{a}_1,\boldsymbol{a}_2\), define \(I_{4a}=\|\boldsymbol{F}\boldsymbol{a}_a\|^2\) and
    !! \(z_a=\max(I_{4a}-1,0)\). The *adapted* HGO-type Helmholtz energy
    !! per reference volume is
    !! \[\Psi(\boldsymbol{F})=\frac{\mu}{2}(\boldsymbol{F}:\boldsymbol{F}-3-2\ln J)
    !! +\frac{\lambda}{2}(\ln J)^2+
    !! \sum_{a=1}^2\frac{k_1}{2k_2}
    !! [\exp(k_2z_a^2)-1].\]
    !! Matrix shear modulus \(\mu>0\), first Lamé modulus \(\lambda\geq0\) and fiber stiffness \(k_1\geq0\)
    !! share the chosen stress unit. Fiber curvature \(k_2>0\) is dimensionless.
    !! The underlying matrix is compressible neo-Hookean. The fibers
    !! resist stretch only, not compression.
    !! The fiber invariants use the full deformation gradient, so volume
    !! expansion can activate fibers. There is no fiber dispersion.
    !!
    !! Let \(\boldsymbol{B}=\boldsymbol{F}^{-T}\), \(l=\ln J\) and \(\boldsymbol{v}_a=\boldsymbol{F}\boldsymbol{a}_a\). The work-conjugate
    !! first Piola stress is the *analytic* derivative
    !! \[\boldsymbol{P}=\mu(\boldsymbol{F}-\boldsymbol{B})+\lambda l\boldsymbol{B}+
    !! \sum_a 2 k_1 z_a e^{k_2z_a^2}\boldsymbol{v}_a\otimes\boldsymbol{a}_a.\]
    !! The exact material tangent is
    !! \[\frac{\partial P_{iJ}}{\partial F_{kL}}=
    !! \mu\delta_{ik}\delta_{JL}+
    !! (\mu-\lambda l)B_{iL}B_{kJ}+\lambda B_{iJ}B_{kL}
    !! +\sum_{a:z_a>0}\left[
    !! 2k_1z_ae^{k_2z_a^2}\delta_{ik}a_{a,J}a_{a,L}
    !! +4k_1e^{k_2z_a^2}(1+2k_2z_a^2)
    !! v_{a,i}v_{a,k}a_{a,J}a_{a,L}\right].\]
    !! The stretch/compression switch at \(I_{4a}=1\) has a one-sided
    !! algorithmic tangent. Exactly at the switch the inactive tangent
    !! is used. Avoid tangent verification at this non-smooth surface.
    !!
    !! Input `deformation.f` and first output `stress.piola` are full
    !! \(3\times3\) spatial-reference tensors in Fortran column-major component-index order
    !! \((11,21,31,12,22,32,13,23,33)\) (not Mandel), dimensionless
    !! and stress units respectively. Second output `energy.helmholtz` is the
    !! scalar strain energy per reference volume in energy-per-volume units. The complete
    !! \(10\times9\) tangent's final row is \(\boldsymbol{P}=\partial\Psi/\partial\boldsymbol{F}\) in the same packed order.
    !!
    !! Analytical sensitivities are available for \(\mu\), \(\lambda\), \(k_1\) and \(k_2\)
    !! at fixed \(\boldsymbol{F}\) and fixed reference fiber directions, including the
    !! energy derivative. No `dt` or history is required. Rotating the
    !! current frame by \(\boldsymbol{R}\in\mathrm{SO}(3)\) gives \(\boldsymbol{P}(\boldsymbol{R}\boldsymbol{F})=\boldsymbol{R}\boldsymbol{P}(\boldsymbol{F})\) and leaves
    !! \(\Psi\) unchanged. Inverted/nonfinite deformations and overflow
    !! reject the material point.
    !!
    !! References: G. A. Holzapfel, T. C. Gasser and R. W. Ogden,
    !! "A New Constitutive Framework for Arterial Wall Mechanics and a
    !! Comparative Study of Material Models", Journal of Elasticity
    !! 61 (2000), 1--48, DOI 10.1023/A:1010835316564.
    !! The compressible matrix, two-fiber energy and no-dispersion
    !! assumptions are the present adaptation, not a complete layered
    !! arterial model or experimentally validated tissue calibration.
    !!
    !! D. R. Nolan, A. L. Gower, M. Destrade, R. W. Ogden and J. P. McGarry,
    !! "A robust anisotropic hyperelastic formulation for the modelling of
    !! soft tissue", Journal of the Mechanical Behavior of Biomedical
    !! Materials 39 (2014), 48--60, DOI 10.1016/j.jmbbm.2014.06.016.
    !! Section 2.4, Eq. (16), gives the full-invariant fiber contribution:
    !! https://maths.nuigalway.ie/~arturgower/Documents/Nolan_Gower_2014.pdf.
    !! Code provenance: independent original implementation.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: ndim = 3
    integer, parameter :: nfibers = 2
    integer, parameter :: nfull = 9
    integer, parameter :: noutput = 10
    public :: fiber_hgo_t

    type, extends(material_model_t) :: fiber_hgo_t
        !! Construct as `fiber_hgo_t(mu,lambda,k1,k2,directions)`.
        real(wp) :: shear_modulus = 0.0_wp
            !! \(\mu>0\) in the chosen stress unit.
        real(wp) :: lame_lambda = 0.0_wp
            !! \(\lambda\geq0\) in the chosen stress unit.
        real(wp) :: fiber_k1 = 0.0_wp
            !! \(k_1\geq0\) in the chosen stress unit.
        real(wp) :: fiber_k2 = 0.0_wp
            !! \(k_2>0\), dimensionless.
        real(wp), allocatable :: fiber_directions(:,:)
            !! Optional two unit reference-frame direction vectors with array shape `[3,2]`.
            !! If absent, use \([1,0,0]^T\) and \([0,1,0]^T\) for the two families.
            !! Supplied at model initialization, never changed in kernels.
    contains
        procedure :: describe => describe_fiber_hgo
            !! Declare input and output fields, parameter names and constitutive capabilities.
        procedure :: eval_batch => evaluate_fiber_hgo
            !! Evaluate the batched constitutive response and optional consistent tangent.
        procedure :: eval_sensitivity_batch => sensitivity_fiber_hgo
            !! Evaluate exact scalar-parameter derivatives of the constitutive response.
    end type fiber_hgo_t

contains

    subroutine describe_fiber_hgo(self,spec,status)
        !! Validate constitutive parameters and register `deformation.f` input and `stress.piola` and `energy.helmholtz` response fields.
        class(fiber_hgo_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        type(material_spec_t), intent(out) :: spec
            !! Packed field schema, material parameter metadata and declared capabilities.
        integer, intent(out) :: status
            !! Operation status; `cf_ok` indicates success.
        integer :: a
        real(wp) :: norm_squared

        status = cf_invalid_argument
        if (.not. all(ieee_is_finite([self%shear_modulus,self%lame_lambda,self%fiber_k1,self%fiber_k2]))) return
        if (self%shear_modulus <= 0.0_wp .or. self%lame_lambda < 0.0_wp) return
        if (self%fiber_k1 < 0.0_wp .or. self%fiber_k2 <= 0.0_wp) return
        if (allocated(self%fiber_directions)) then
            if (any(shape(self%fiber_directions) /= [ndim,nfibers])) return
            if (.not. all(ieee_is_finite(self%fiber_directions))) return
            do a = 1,nfibers
                norm_squared = dot_product(self%fiber_directions(:,a),self%fiber_directions(:,a))
                if (abs(norm_squared-1.0_wp) > 64.0_wp*epsilon(1.0_wp)) return
            end do
        end if

        call spec%add_input("deformation.f", nfull, "tensor_full", "deformation_gradient", "spatial_reference", [ndim,ndim])
        call spec%add_output("stress.piola", nfull, "tensor_full", "first_piola_kirchhoff", "spatial_reference", [ndim,ndim])
        call spec%add_output("energy.helmholtz", 1, "scalar", "helmholtz_energy", "material")
        call spec%add_parameter("shear_modulus")
        call spec%add_parameter("lame_lambda")
        call spec%add_parameter("fiber_k1")
        call spec%add_parameter("fiber_k2")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_fiber_hgo

    subroutine evaluate_fiber_hgo(self,x_old,x_new,state_old,state_trial,response,dt,tangent,point_status)
        !! Evaluate `stress.piola` and `energy.helmholtz` from `deformation.f` over the current point batch.
        class(fiber_hgo_t), intent(in) :: self
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
            !! Packed `stress.piola` and `energy.helmholtz` output `(noutput,npoints)`.
        real(wp), intent(in) :: dt
            !! Prescribed time increment; required to be positive only for time-dependent models.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic response Jacobian `(noutput,ninput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Status of each point; a successful update must set its code to `cf_ok`.

        state_trial = state_old
        call evaluate_points(self,x_new,response,point_status,tangent)
    end subroutine evaluate_fiber_hgo

    subroutine sensitivity_fiber_hgo(self,x_old,x_new,state_old,state_trial,selected_parameters,sensitivity_old,sensitivity_trial, response,dsdp,dt,tangent,point_status)
        !! Evaluate exact derivatives of `stress.piola` and `energy.helmholtz` with respect to the selected scalar parameters.
        class(fiber_hgo_t), intent(in) :: self
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
            !! Packed `stress.piola` and `energy.helmholtz` output `(noutput,npoints)`.
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
        call evaluate_points(self,x_new,response,point_status,tangent,selected_parameters,dsdp)
    end subroutine sensitivity_fiber_hgo

    subroutine evaluate_points(self,x,response,point_status,tangent,selected_parameters,dsdp)
        !! Compute `stress.piola` and `energy.helmholtz` and optional exact input and parameter derivatives for each point.
        class(fiber_hgo_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        real(wp), intent(in) :: x(:,:)
            !! Current prescribed `deformation.f` inputs for the point batch.
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.piola` and `energy.helmholtz` output `(noutput,npoints)`.
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
        real(wp) :: directions(ndim,nfibers)
        real(wp) :: fiber_vector(ndim,nfibers)
        real(wp) :: activation(nfibers)
        real(wp) :: exponential(nfibers)
        real(wp) :: energy_ratio(nfibers)
        real(wp) :: curvature_ratio(nfibers)
        real(wp) :: jacobian
        real(wp) :: log_jacobian
        real(wp) :: i1
        real(wp) :: value
        real(wp) :: energy
        real(wp) :: fiber_factor
        real(wp) :: z
        real(wp) :: weight
        real(wp) :: energy_derivative
        real(wp) :: series_term
        integer :: ipoint
        integer :: i
        integer :: j
        integer :: k
        integer :: l
        integer :: a
        integer :: row
        integer :: column
        integer :: iparam
        integer :: n
        logical :: valid
        logical :: needs_fiber_term

        directions = reshape([1.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp,0.0_wp],[ndim,nfibers])
        if (allocated(self%fiber_directions)) directions = self%fiber_directions
        needs_fiber_term = self%fiber_k1 > 0.0_wp
        if (present(selected_parameters)) needs_fiber_term = needs_fiber_term .or. any(selected_parameters == 3)

        point_status = cf_invalid_argument
        point_loop: do ipoint = 1,size(x,2)
            if (.not. all(ieee_is_finite(x(:,ipoint)))) cycle point_loop
            do j = 1,ndim
                deformation(:,j) = x(ndim*(j-1)+1:ndim*j,ipoint)
            end do
            call cofactor_matrix(deformation,cofactor)
            jacobian = dot_product(deformation(1,:),cofactor(1,:))
            if (.not. ieee_is_finite(jacobian)) cycle point_loop
            if (jacobian <= 0.0_wp) cycle point_loop
            inverse_t = cofactor/jacobian
            if (.not. all(ieee_is_finite(inverse_t))) cycle point_loop

            log_jacobian = log(jacobian)
            i1 = sum(deformation**2)
            stress = self%shear_modulus*(deformation-inverse_t)+ self%lame_lambda*log_jacobian*inverse_t
            energy = 0.5_wp*self%shear_modulus*(i1-3.0_wp-2.0_wp*log_jacobian)+ 0.5_wp*self%lame_lambda*log_jacobian**2
            valid = .true.
            activation = 0.0_wp
            fiber_stress_loop: do a = 1,nfibers
                if (.not. needs_fiber_term) exit fiber_stress_loop
                fiber_vector(:,a) = matmul(deformation,directions(:,a))
                activation(a) = max(0.0_wp,dot_product(fiber_vector(:,a),fiber_vector(:,a))-1.0_wp)
                if (activation(a) <= 0.0_wp) cycle fiber_stress_loop
                z = self%fiber_k2*activation(a)**2
                if (.not. ieee_is_finite(z)) then
                    valid = .false.
                    exit fiber_stress_loop
                end if
                exponential(a) = exp(z)
                if (.not. ieee_is_finite(exponential(a))) then
                    valid = .false.
                    exit fiber_stress_loop
                end if
                if (z < 0.5_wp) then
                    ! Positive Taylor sums for (exp(z)-1)/z and ((z-1)*exp(z)+1)/z**2.
                    energy_ratio(a) = 1.0_wp
                    curvature_ratio(a) = 0.5_wp
                    series_term = 1.0_wp
                    series_loop: do n = 1,digits(1.0_wp)
                        series_term = series_term*z/real(n,wp)
                        energy_ratio(a) = energy_ratio(a)+series_term/real(n+1,wp)
                        curvature_ratio(a) = curvature_ratio(a)+series_term/real(n+2,wp)
                        if (series_term <= epsilon(1.0_wp)) exit series_loop
                    end do series_loop
                else
                    energy_ratio(a) = (exponential(a)-1.0_wp)/z
                    curvature_ratio(a) = (exponential(a)-energy_ratio(a))/z
                end if
                energy = energy+0.5_wp*self%fiber_k1*activation(a)**2*energy_ratio(a)
                fiber_factor = 2.0_wp*self%fiber_k1*activation(a)*exponential(a)
                do j = 1,ndim
                    stress(:,j) = stress(:,j)+fiber_factor*fiber_vector(:,a)*directions(j,a)
                end do
            end do fiber_stress_loop
            if (.not. valid) cycle point_loop
            if (.not. all(ieee_is_finite(stress)) .or. .not. ieee_is_finite(energy)) cycle point_loop
            do j = 1,ndim
                response(ndim*(j-1)+1:ndim*j,ipoint) = stress(:,j)
            end do
            response(noutput,ipoint) = energy

            if (present(tangent)) then
                do j = 1,ndim
                    do i = 1,ndim
                        row = i+ndim*(j-1)
                        do l = 1,ndim
                            do k = 1,ndim
                                column = k+ndim*(l-1)
                                value = (self%shear_modulus-self%lame_lambda*log_jacobian)*inverse_t(i,l)*inverse_t(k,j) + self%lame_lambda*inverse_t(i,j)*inverse_t(k,l)
                                if (i == k .and. j == l) value = value+self%shear_modulus
                                fiber_tangent_loop: do a = 1,nfibers
                                    if (activation(a) <= 0.0_wp) cycle fiber_tangent_loop
                                    weight = directions(j,a)*directions(l,a)
                                    if (i == k) value = value+2.0_wp*self%fiber_k1* activation(a)*exponential(a)*weight
                                    value = value+4.0_wp*self%fiber_k1*exponential(a)* (1.0_wp+2.0_wp*self%fiber_k2*activation(a)**2)* fiber_vector(i,a)*fiber_vector(k,a)*weight
                                end do fiber_tangent_loop
                                tangent(row,column,ipoint) = value
                            end do
                        end do
                    end do
                end do
                tangent(noutput,:,ipoint) = response(1:nfull,ipoint)
                if (.not. all(ieee_is_finite(tangent(:,:,ipoint)))) cycle point_loop
            end if

            if (present(selected_parameters)) then
                parameter_loop: do iparam = 1,size(selected_parameters)
                    dsdp(:,iparam,ipoint) = 0.0_wp
                    select case (selected_parameters(iparam))
                    case (1)
                        do j = 1,ndim
                            row = ndim*(j-1)
                            dsdp(row+1:row+ndim,iparam,ipoint) = deformation(:,j)-inverse_t(:,j)
                        end do
                        dsdp(noutput,iparam,ipoint) = 0.5_wp*(i1-3.0_wp-2.0_wp*log_jacobian)
                    case (2)
                        do j = 1,ndim
                            row = ndim*(j-1)
                            dsdp(row+1:row+ndim,iparam,ipoint) = log_jacobian*inverse_t(:,j)
                        end do
                        dsdp(noutput,iparam,ipoint) = 0.5_wp*log_jacobian**2
                    case (3,4)
                        fiber_parameter_loop: do a = 1,nfibers
                            if (activation(a) <= 0.0_wp) cycle fiber_parameter_loop
                            if (selected_parameters(iparam) == 3) then
                                fiber_factor = 2.0_wp*activation(a)*exponential(a)
                                energy_derivative = 0.5_wp*activation(a)**2*energy_ratio(a)
                            else
                                fiber_factor = 2.0_wp*self%fiber_k1*activation(a)**3*exponential(a)
                                energy_derivative = 0.5_wp*self%fiber_k1*activation(a)**4*curvature_ratio(a)
                            end if
                            do j = 1,ndim
                                row = ndim*(j-1)
                                dsdp(row+1:row+ndim,iparam,ipoint) = dsdp(row+1:row+ndim,iparam,ipoint) + fiber_factor*fiber_vector(:,a)*directions(j,a)
                            end do
                            dsdp(noutput,iparam,ipoint) = dsdp(noutput,iparam,ipoint)+energy_derivative
                        end do fiber_parameter_loop
                    case default
                        exit parameter_loop
                    end select
                end do parameter_loop
                if (iparam <= size(selected_parameters)) cycle point_loop
                if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
            end if
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

    pure subroutine cofactor_matrix(deformation,cofactor)
        !! Cofactor matrix \(\operatorname{cof}\boldsymbol{F}=(\det\boldsymbol{F})\boldsymbol{F}^{-T}\).
        real(wp), intent(in) :: deformation(ndim,ndim)
            !! Deformation gradient \(\boldsymbol{F}\) in full three-dimensional tensor form.
        real(wp), intent(out) :: cofactor(ndim,ndim)
            !! Cofactor matrix \(\operatorname{cof}\boldsymbol{F}=J\boldsymbol{F}^{-T}\).

        cofactor(1,1) = deformation(2,2)*deformation(3,3)-deformation(2,3)*deformation(3,2)
        cofactor(1,2) = deformation(2,3)*deformation(3,1)-deformation(2,1)*deformation(3,3)
        cofactor(1,3) = deformation(2,1)*deformation(3,2)-deformation(2,2)*deformation(3,1)
        cofactor(2,1) = deformation(1,3)*deformation(3,2)-deformation(1,2)*deformation(3,3)
        cofactor(2,2) = deformation(1,1)*deformation(3,3)-deformation(1,3)*deformation(3,1)
        cofactor(2,3) = deformation(1,2)*deformation(3,1)-deformation(1,1)*deformation(3,2)
        cofactor(3,1) = deformation(1,2)*deformation(2,3)-deformation(1,3)*deformation(2,2)
        cofactor(3,2) = deformation(1,3)*deformation(2,1)-deformation(1,1)*deformation(2,3)
        cofactor(3,3) = deformation(1,1)*deformation(2,2)-deformation(1,2)*deformation(2,1)
    end subroutine cofactor_matrix

end module constituf_fiber_hgo
