! SPDX-License-Identifier: Apache-2.0

module constituf_axisymmetric_elastic
    !! Torsion-free axisymmetric small-strain, three-dimensional isotropic elasticity.
    !!
    !! This specialization uses the four independent nonzero components of
    !! a symmetric three-dimensional spatial tensor under axisymmetry with
    !! zero circumferential displacement \(u_\theta=0\):
    !! \([\varepsilon_{rr},\varepsilon_{zz},
    !! \varepsilon_{\theta\theta},\sqrt{2}\varepsilon_{rz}]\).
    !! The output Cauchy stress uses exactly the same order and Mandel
    !! shear scaling, ensuring the work product
    !! \(\boldsymbol\sigma:\boldsymbol\varepsilon
    !! =\boldsymbol\sigma_M^T\boldsymbol\varepsilon_M\).
    !! The 3D physical tensor shape is `[3,3]`, but only four independent
    !! components are retained. The custom representation name
    !! `custom.axisymmetric_mandel` explicitly distinguishes these from
    !! plane strain and conventional 3/6-component Mandel layouts.
    !!
    !! For \(K>0\) and \(G>0\), with \(\lambda=K-2G/3\),
    !! \[\sigma_{ii}=\lambda(\varepsilon_{rr}+
    !! \varepsilon_{zz}+\varepsilon_{\theta\theta})
    !! +2G\varepsilon_{ii},\quad i=r,z,\theta,\qquad
    !! \sigma_{rz}^M=2G\varepsilon_{rz}^M.\]
    !! The exact consistent tangent is
    !! \[C=\begin{bmatrix}
    !! \lambda+2G&\lambda&\lambda&0\\
    !! \lambda&\lambda+2G&\lambda&0\\
    !! \lambda&\lambda&\lambda+2G&0\\
    !! 0&0&0&2G\end{bmatrix}.\]
    !! Selected analytical sensitivities to \(K\) and \(G\) use the same
    !! three-dimensional Hooke law. At prescribed strain,
    !! \(\partial_K\sigma_{ii}=\mathrm{tr}\varepsilon\),
    !! \(\partial_G\sigma_{ii}=2\varepsilon_{ii}-
    !! \tfrac23\mathrm{tr}\varepsilon\) and
    !! \(\partial_G\sigma_{rz}^M=2\varepsilon_{rz}^M\).
    !!
    !! The material is stateless and requires no time step. Strain uses
    !! dimensionless small strain, stresses/moduli in the chosen stress unit and positive stress
    !! means tension. The surrounding axisymmetric FEM solver must compute
    !! hoop strain from its displacement field and treat the \(r=0\)
    !! regularity condition itself. For \(r>0\),
    !! \(\varepsilon_{\theta\theta}=u_r/r\). At the axis a regular radial
    !! field gives the limit \(\partial u_r/\partial r\).
    !! This model does not compute geometry.
    !! Not a shell, generalized plane strain, or axisymmetric finite-strain
    !! law. No rotated or material-axis anisotropy is assumed.
    !!
    !! References:
    !!
    !! T. J. R. Hughes, The Finite Element Method: Linear
    !! Static and Dynamic Finite Element Analysis, Dover, 2000,
    !! ISBN 978-0-486-41181-1. The 4-component matrix is derived directly
    !! from 3D isotropic Hooke elasticity and axisymmetric kinematics.
    !!
    !! R. A. Regueiro, *Introduction to Finite Element Analysis*,
    !! University of Colorado Boulder, Fall 2024, Section 5.5, Eq. (5.40),
    !! gives the torsion-free axisymmetric kinematics and elastic matrix:
    !! https://ceae.colorado.edu/~regueiro/courses/CVEN4511_5511_notes.pdf.
    !! Code provenance: independently implemented. No third-party code.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: ncomponents = 4
    public :: axisymmetric_elastic_t

    type, extends(material_model_t) :: axisymmetric_elastic_t
        !! Axisymmetric 4-component Hooke material, construct with \(K,G\) in the chosen stress unit.
        real(wp) :: bulk_modulus = 0.0_wp
            !! Positive bulk modulus.
        real(wp) :: shear_modulus = 0.0_wp
            !! Positive shear modulus.
    contains
        procedure :: describe => describe_axisymmetric_elastic
            !! Declare the constitutive fields, parameter metadata and supported capabilities.
        procedure :: eval_batch => evaluate_axisymmetric_elastic
            !! Evaluate the batched constitutive response, optional tangent and trial history.
        procedure :: eval_sensitivity_batch => sensitivity_axisymmetric_elastic
            !! Evaluate exact parameter-path derivatives and corresponding trial history.
    end type axisymmetric_elastic_t

contains

    subroutine describe_axisymmetric_elastic(self, spec, status)
        !! Validate the parameters and declare `strain.small` input and `stress.cauchy` output with their physical metadata.
        class(axisymmetric_elastic_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema receiving field layouts, parameters and supported operations.
        integer, intent(out) :: status
            !! Operation status. `cf_ok` indicates a valid model or successful initialization.

        status = cf_invalid_argument
        if (.not. all(ieee_is_finite([self%bulk_modulus,self%shear_modulus]))) return
        if (self%bulk_modulus <= 0.0_wp .or. self%shear_modulus <= 0.0_wp) return
        if (.not. ieee_is_finite(2.0_wp*self%shear_modulus)) return
        if (.not. ieee_is_finite(self%bulk_modulus+(4.0_wp/3.0_wp)*self%shear_modulus)) return
        call spec%add_input("strain.small", ncomponents, "custom.axisymmetric_mandel", "small_strain", "spatial", [3,3])
        call spec%add_output("stress.cauchy", ncomponents, "custom.axisymmetric_mandel", "cauchy_stress", "spatial", [3,3])
        call spec%add_parameter("bulk_modulus")
        call spec%add_parameter("shear_modulus")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_axisymmetric_elastic

    subroutine evaluate_axisymmetric_elastic(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Evaluate `stress.cauchy` from `strain.small` over the point batch and update trial history.
        class(axisymmetric_elastic_t), intent(in) :: self
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
        call evaluate_points(self,x_new,response,point_status,tangent)
    end subroutine evaluate_axisymmetric_elastic

    subroutine sensitivity_axisymmetric_elastic(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Evaluate exact selected-parameter derivatives of `stress.cauchy` through the discrete material update.
        class(axisymmetric_elastic_t), intent(in) :: self
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
        call evaluate_points(self,x_new,response,point_status,tangent,selected_parameters,dsdp)
    end subroutine sensitivity_axisymmetric_elastic

    subroutine evaluate_points(self,x,response,point_status,tangent,selected_parameters,dsdp)
        !! Compute `stress.cauchy` and the optional consistent tangent and parameter derivatives for each point.
        class(axisymmetric_elastic_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp), intent(in) :: x(:,:)
            !! Prescribed `strain.small` input components for the point batch.
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
        real(wp) :: lambda
        real(wp) :: trace_strain
        real(wp) :: matrix(ncomponents,ncomponents)
        integer :: ipoint
        integer :: i
        integer :: j

        lambda = self%bulk_modulus-(2.0_wp/3.0_wp)*self%shear_modulus
        if (present(tangent)) then
            matrix = 0.0_wp
            matrix(1:3,1:3) = lambda
            do i = 1, 3
                matrix(i,i) = self%bulk_modulus+(4.0_wp/3.0_wp)*self%shear_modulus
            end do
            matrix(4,4) = 2.0_wp*self%shear_modulus
        end if
        point_status = cf_invalid_argument
        point_loop: do ipoint = 1, size(x,2)
            if (.not. all(ieee_is_finite(x(:,ipoint)))) cycle point_loop
            trace_strain = sum(x(1:3,ipoint))
            response(1:3,ipoint) = self%bulk_modulus*trace_strain+2.0_wp*self%shear_modulus*(x(1:3,ipoint)-trace_strain/3.0_wp)
            response(4,ipoint) = 2.0_wp*self%shear_modulus*x(4,ipoint)
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop
            if (present(tangent)) tangent(:,:,ipoint) = matrix
            if (present(selected_parameters)) then
                parameter_loop: do j = 1, size(selected_parameters)
                    select case (selected_parameters(j))
                    case (1)
                        dsdp(:,j,ipoint) = 0.0_wp
                        dsdp(1:3,j,ipoint) = trace_strain
                    case (2)
                        dsdp(1:3,j,ipoint) = 2.0_wp*(x(1:3,ipoint)-trace_strain/3.0_wp)
                        dsdp(4,j,ipoint) = 2.0_wp*x(4,ipoint)
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

end module constituf_axisymmetric_elastic
