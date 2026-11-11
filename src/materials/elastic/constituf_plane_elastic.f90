! SPDX-License-Identifier: Apache-2.0

module constituf_plane_elastic
    !! Isotropic small-strain 2D plane-strain or plane-stress elasticity.
    !!
    !! Input `strain.small` and output `stress.cauchy` are the in-plane \(2\times2\)
    !! symmetric tensors. Either tensor \(\boldsymbol{A}\) is stored as the Mandel
    !! vector \(\boldsymbol{A}_M=[A_{11},A_{22},\sqrt{2}A_{12}]^T\).
    !! They are work conjugate. Strain is dimensionless and stress is in the chosen stress unit.
    !! These distinct constraints are derived from 3D Hooke elasticity.
    !! Plane strain restricts the 3D strain. Plane stress condenses the
    !! out-of-plane strain to enforce zero out-of-plane stress.
    !!
    !! For plane strain, \(\varepsilon_{33}=\varepsilon_{13}=\varepsilon_{23}=0\), while the
    !! out-of-plane stress \(\sigma_{33}\) generally is nonzero. With 3D Lamé
    !! coefficient \(\lambda=K-2G/3\), the in-plane law is
    !! \[\sigma_{\alpha\beta}=
    !! \lambda(\varepsilon_{11}+\varepsilon_{22})\delta_{\alpha\beta}
    !! +2G\varepsilon_{\alpha\beta}.\]
    !! For plane stress, \(\sigma_{33}=\sigma_{13}=\sigma_{23}=0\) and \(\varepsilon_{33}\) is
    !! eliminated from three-dimensional Hooke elasticity. The condensed
    !! in-plane coefficient is
    !! \[\bar\lambda=\frac{2G(3K-2G)}{3K+4G},\qquad
    !! \boldsymbol C_\mathrm{2D}=
    !! \begin{bmatrix}\bar\lambda+2G&\bar\lambda&0\\
    !! \bar\lambda&\bar\lambda+2G&0\\0&0&2G\end{bmatrix}.\]
    !! Plane strain uses the same matrix with \(\bar\lambda\) replaced by \(\lambda\).
    !! The 3D stored-energy restriction/condensation yields these symmetric,
    !! positive-definite exact tangents for \(K>0\) and \(G>0\). Both states are
    !! stateless and have no time-step requirement.
    !!
    !! Analytical scalar parameter sensitivities are provided for \(K,G\).
    !! For either constraint let \(a\) be the appropriate 2D Lamé coefficient and
    !! \(t=\varepsilon_{11}+\varepsilon_{22}\). Then
    !! \[\partial_p\sigma_{\alpha\alpha}
    !! =(\partial_p a)t+2(\partial_pG)\varepsilon_{\alpha\alpha},
    !! \qquad \partial_p\sigma_{12}^{M}
    !! =2(\partial_pG)\varepsilon_{12}^{M}.\]
    !! In plane strain \(\partial a/\partial K=1\) and \(\partial a/\partial G=-2/3\). In plane stress,
    !! \[\frac{\partial a}{\partial K}
    !! =\frac{36G^2}{(3K+4G)^2},\quad
    !! \frac{\partial a}{\partial G}
    !! =\frac{18K^2-24KG-16G^2}{(3K+4G)^2}.\]
    !! The 2D response does not include \(\sigma_{33}\). Request a 3D model if
    !! out-of-plane stress is needed. This model is not a nonlinear
    !! plane-stress adapter for plasticity or finite deformations.
    !!
    !! References:
    !!
    !! J. C. Simo and T. J. R. Hughes, Computational
    !! Inelasticity, Springer (1998), DOI 10.1007/b98904. The algebraic
    !! plane-stress condensation is independently derived from 3D Hooke law.
    !!
    !! A. F. Bower, *Applied Mechanics of Solids*, online text,
    !! Sections 3.2.3 and 3.2.5 give the plane-deformation laws and elastic
    !! modulus conversions: https://solidmechanics.org/Text/Chapter3_2/Chapter3_2.php.
    !! Code provenance: original implementation from stated equations.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    public :: plane_elastic_t

    type, extends(material_model_t) :: plane_elastic_t
        !! 2D isotropic elasticity. `plane_stress=.false.` selects plane strain.
        !!
        !! Construct with `(bulk_modulus, shear_modulus, plane_stress)`.
        !! Initial input and `dt` are not required. The handle copies moduli
        !! and the constraint during initialization. Both moduli must be finite and positive.
        real(wp) :: bulk_modulus = 0.0_wp
            !! Three-dimensional bulk modulus \(K\) in the chosen stress unit.
        real(wp) :: shear_modulus = 0.0_wp
            !! Three-dimensional shear modulus \(G\) in the chosen stress unit.
        logical :: plane_stress = .false.
            !! If true, enforce \(\sigma_{33}=0\). Otherwise \(\varepsilon_{33}=0\).
    contains
        procedure :: describe => describe_plane_elastic
            !! Declare the constitutive fields, parameter metadata and supported capabilities.
        procedure :: eval_batch => evaluate_plane_elastic
            !! Evaluate the batched constitutive response, optional tangent and trial history.
        procedure :: eval_sensitivity_batch => sensitivity_plane_elastic
            !! Evaluate exact parameter-path derivatives and corresponding trial history.
    end type plane_elastic_t

contains

    subroutine describe_plane_elastic(self, spec, status)
        !! Validate the parameters and declare `strain.small` input and `stress.cauchy` output with their physical metadata.
        class(plane_elastic_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema receiving field layouts, parameters and supported operations.
        integer, intent(out) :: status
            !! Operation status. `cf_ok` indicates a valid model or successful initialization.

        status = cf_invalid_argument
        if (.not. ieee_is_finite(self%bulk_modulus)) return
        if (.not. ieee_is_finite(self%shear_modulus)) return
        if (self%bulk_modulus <= 0.0_wp .or. self%shear_modulus <= 0.0_wp) return
        if (.not. ieee_is_finite(3.0_wp*self%bulk_modulus + 4.0_wp*self%shear_modulus)) return

        call spec%add_input("strain.small", 3, "mandel", "small_strain", "spatial", [2,2])
        call spec%add_output("stress.cauchy", 3, "mandel", "cauchy_stress", "spatial", [2,2])
        call spec%add_parameter("bulk_modulus")
        call spec%add_parameter("shear_modulus")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_plane_elastic

    subroutine evaluate_plane_elastic(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Evaluate `stress.cauchy` from `strain.small` over the point batch and update trial history.
        class(plane_elastic_t), intent(in) :: self
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
    end subroutine evaluate_plane_elastic

    subroutine sensitivity_plane_elastic(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Evaluate exact selected-parameter derivatives of `stress.cauchy` through the discrete material update.
        class(plane_elastic_t), intent(in) :: self
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
    end subroutine sensitivity_plane_elastic

    subroutine evaluate_points(self, x_new, response, point_status, tangent, selected_parameters, dsdp)
        !! Compute `stress.cauchy` and the optional consistent tangent and parameter derivatives for each point.
        class(plane_elastic_t), intent(in) :: self
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
        real(wp) :: a
        real(wp) :: da_k
        real(wp) :: da_g
        real(wp) :: denominator
        real(wp) :: bulk_ratio
        real(wp) :: shear_ratio
        real(wp) :: trace_strain
        real(wp) :: shear
        integer :: ipoint
        integer :: j

        shear = self%shear_modulus
        if (self%plane_stress) then
            denominator = 3.0_wp*self%bulk_modulus + 4.0_wp*shear
            bulk_ratio = self%bulk_modulus/denominator
            shear_ratio = shear/denominator
            a = 2.0_wp*shear*((3.0_wp*self%bulk_modulus-2.0_wp*shear)/denominator)
            da_k = 36.0_wp*shear_ratio**2
            da_g = 18.0_wp*bulk_ratio**2 - 24.0_wp*bulk_ratio*shear_ratio - 16.0_wp*shear_ratio**2
        else
            a = self%bulk_modulus - 2.0_wp*shear/3.0_wp
            da_k = 1.0_wp
            da_g = -2.0_wp/3.0_wp
        end if
        point_status = cf_invalid_argument
        point_loop: do ipoint = 1, size(x_new,2)
            if (.not. all(ieee_is_finite(x_new(:,ipoint)))) cycle point_loop
            trace_strain = x_new(1,ipoint) + x_new(2,ipoint)
            response(1:2,ipoint) = a*trace_strain + 2.0_wp*shear*x_new(1:2,ipoint)
            response(3,ipoint) = 2.0_wp*shear*x_new(3,ipoint)
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop
            if (present(tangent)) then
                tangent(:,:,ipoint) = 0.0_wp
                tangent(1:2,1:2,ipoint) = a
                tangent(1,1,ipoint) = a + 2.0_wp*shear
                tangent(2,2,ipoint) = a + 2.0_wp*shear
                tangent(3,3,ipoint) = 2.0_wp*shear
            end if
            if (present(selected_parameters)) then
                parameter_loop: do j = 1, size(selected_parameters)
                    select case (selected_parameters(j))
                    case (1)
                        dsdp(1:2,j,ipoint) = da_k*trace_strain
                        dsdp(3,j,ipoint) = 0.0_wp
                    case (2)
                        dsdp(1:2,j,ipoint) = da_g*trace_strain + 2.0_wp*x_new(1:2,ipoint)
                        dsdp(3,j,ipoint) = 2.0_wp*x_new(3,ipoint)
                    case default
                        exit parameter_loop
                    end select
                end do parameter_loop
                if (j <= size(selected_parameters)) cycle point_loop
                if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
            end if
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

end module constituf_plane_elastic
