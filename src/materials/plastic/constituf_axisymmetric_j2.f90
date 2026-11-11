! SPDX-License-Identifier: Apache-2.0

module constituf_axisymmetric_j2
    !! Four-component torsion-free axisymmetric small-strain \(J_2\) plasticity.
    !!
    !! With no circumferential displacement \(u_\theta=0\) and no
    !! \(\theta\)-dependence, the symmetric 3D spatial tensor in the local
    !! orthonormal \((r,z,\theta)\) basis has the independent Mandel components
    !! \[\boldsymbol\varepsilon^{\mathrm{axi}}_M=
    !! [\varepsilon_{rr},\varepsilon_{zz},
    !! \varepsilon_{\theta\theta},\sqrt2\varepsilon_{rz}]^T.\]
    !! The Cauchy stress is packed in the identical order, preserving
    !! the stress/strain work-conjugacy. The declared physical shape is
    !! \(3\times3\) and the *custom* representation is
    !! `custom.axisymmetric_mandel`, not a standard 2D Mandel tensor.
    !! A host FEM solver must compute the hoop strain
    !! \(\varepsilon_{\theta\theta}=u_r/r\) for \(r>0\), supply its regular
    !! limit \(\partial u_r/\partial r\) at \(r=0\) and assemble the
    !! axisymmetric geometrical weights.
    !!
    !! Constitutive dynamics are **exactly the 3D \(J_2\) radial return**
    !! described in [[constituf_j2_plasticity(module)]], with small
    !! strain, \(K>0\), \(G>0\), \(Y_0>0\), \(H\geq0\), isothermal, associative
    !! rate-independent flow and linear isotropic hardening:
    !! \[\Delta\gamma=
    !! \frac{\langle q_{\mathrm{tr}}-(Y_0+H\alpha_n)\rangle}{3G+H},
    !! \quad
    !! \boldsymbol\varepsilon^p_{n+1}=
    !! \boldsymbol\varepsilon^p_n+\Delta\gamma\boldsymbol n.\]
    !! 3D \(J_2\) internal state consists of six Mandel plastic strains
    !! and a scalar equivalent plastic strain \(\alpha\), all packed
    !! in the material-owned 7-component committed/trial history.
    !! There is no separate out-of-plane constitutive solve: the hoop
    !! strain is one of the four independent axisymmetric inputs.
    !!
    !! The packed axisymmetric input is embedded into standard 3D Mandel
    !! indices \((1,2,3,6)\), with 3D components \(4\) and \(5\) held at zero.
    !! The full 3D consistent *discrete* tangent is restricted to the
    !! same rows/columns:
    !! \[C^{\mathrm{axi}}_{ab}=C^{\mathrm{J2}}_{m_a m_b},
    !! \qquad m=[1,2,3,6].\]
    !! This exact submatrix is NOT a plane-stress condensation. All
    !! stress components and moduli are in the chosen stress unit, strains dimensionless.
    !!
    !! Optional analytical total parameter-path sensitivities with
    !! respect to `bulk_modulus`, `shear_modulus`, `yield_stress` and
    !! `hardening_modulus` are propagated by the same 3D \(J_2\) kernel through
    !! all committed loading steps. The 4-component response
    !! sensitivities retain the selected scalar parameter order.
    !! They hold the prescribed strain history fixed, not global
    !! FEM equilibrium displacements. A yield-surface branch switch
    !! may return `cf_nondifferentiable` when sensitivities are requested.
    !!
    !! Evaluation uses fixed-size point chunks and a single 3D \(J_2\)
    !! batch dispatch for each chunk. There are no per-point heap
    !! allocations, model registration lookups, I/O or global mutable
    !! caches. History is externally managed by `material_t` with normal
    !! trial/commit/revert and checkpoint semantics. No finite-strain,
    !! plane-stress reduction, anisotropy, or damage is included.
    !!
    !! Underlying 3D return and tangent: J. C. Simo and R. L. Taylor,
    !! "Consistent Tangent Operators for Rate-Independent Elastoplasticity",
    !! Computer Methods in Applied Mechanics and Engineering 48 (1985),
    !! 101--118, DOI: 10.1016/0045-7825(85)90070-2. Sections 3 and 4,
    !! Eqs. (3.8) and (4.12)--(4.13), specialized to zero back stress and
    !! linear isotropic hardening. J. C. Simo and T. J. R. Hughes,
    !! Computational Inelasticity, Springer (1998), DOI: 10.1007/b98904.
    !! Axisymmetric kinematics: R. A. Regueiro, Introduction to Finite
    !! Element Analysis, University of Colorado Boulder, Fall 2024,
    !! Section 5.5, pp. 133--134, discussion preceding Eq. (5.40):
    !! https://ceae.colorado.edu/~regueiro/courses/CVEN4511_5511_notes.pdf.
    !! The four-component embedding and Mandel scaling are independently
    !! derived from the 3D law and this torsion-free kinematic restriction.
    !! Code provenance: original Fortran adapter using ConstituF's
    !! own previously implemented 3D \(J_2\) return. No third-party code.
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok
    use constituf_j2_plasticity, only: j2_plasticity_t
    implicit none (type, external)
    private

    integer, parameter :: nmandel = 6
    integer, parameter :: naxis = 4
    integer, parameter :: max_selected = 4
    integer, parameter :: chunk_size = 64
    integer, parameter :: mapping(naxis) = [1,2,3,6]

    public :: axisymmetric_j2_t

    type, extends(material_model_t) :: axisymmetric_j2_t
        !! Construct with `axisymmetric_j2_t(bulk_modulus, shear_modulus, yield_stress, hardening_modulus)`,
        !! all in the chosen stress unit. \(H\) may be zero.
        real(wp) :: bulk_modulus = 0.0_wp
            !! \(K>0\).
        real(wp) :: shear_modulus = 0.0_wp
            !! \(G>0\).
        real(wp) :: yield_stress = 0.0_wp
            !! \(Y_0>0\).
        real(wp) :: hardening_modulus = 0.0_wp
            !! \(H\geq0\).
    contains
        procedure :: describe => describe_axisymmetric_j2
            !! Describe the material input, response, state and scalar parameter schema.
        procedure :: eval_batch => evaluate_axisymmetric_j2
            !! Evaluate the batched material response and trial history.
        procedure :: initialize_sensitivity => initialize_axisymmetric_j2_sensitivity
            !! Initialize analytical derivatives of the reference real history.
        procedure :: eval_sensitivity_batch => sensitivity_axisymmetric_j2
            !! Evaluate exact selected-parameter sensitivities through the discrete material update.
    end type axisymmetric_j2_t

contains

    subroutine describe_axisymmetric_j2(self,spec,status)
        !! Validate inherited \(J_2\) physical parameters and declare 4D fields.
        class(axisymmetric_j2_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema receiving the packed field descriptions and supported capabilities.
        integer, intent(out) :: status
            !! Status code; `cf_ok` indicates success.
        type(j2_plasticity_t) :: base
        type(material_spec_t) :: base_spec

        base = j2_plasticity_t(self%bulk_modulus,self%shear_modulus,self%yield_stress,self%hardening_modulus)
        call base%describe(base_spec,status)
        if (status /= cf_ok) return
        call spec%add_input("strain.small", naxis, "custom.axisymmetric_mandel", "small_strain", "spatial", [3,3])
        call spec%add_output("stress.cauchy", naxis, "custom.axisymmetric_mandel", "cauchy_stress", "spatial", [3,3])
        call spec%add_state("plastic.strain", nmandel, "mandel", "plastic_strain", "spatial", [3,3])
        call spec%add_state("plastic.accumulated", 1, "scalar", "equivalent_plastic_strain", "material")
        call spec%add_parameter("bulk_modulus")
        call spec%add_parameter("shear_modulus")
        call spec%add_parameter("yield_stress")
        call spec%add_parameter("hardening_modulus")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_axisymmetric_j2

    subroutine initialize_axisymmetric_j2_sensitivity(self,initial_input,initial_state,selected_parameters,sensitivity,status)
        !! Virgin \(J_2\) plastic strains and \(\alpha\) have zero parameter derivatives.
        class(axisymmetric_j2_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp), intent(in) :: initial_input(:,:)
            !! Prescribed input state used for material initialization.
        real(wp), intent(in) :: initial_state(:,:)
            !! Initial real-valued history `(nstate,npoints)`.
        integer, intent(in) :: selected_parameters(:)
            !! One-based parameter IDs in the caller's selected derivative-axis order.
        real(wp), intent(out) :: sensitivity(:,:,:)
            !! Initial exact history derivatives `(nstate,nselected,npoints)`.
        integer, intent(out) :: status
            !! Status code; `cf_ok` indicates success.

        sensitivity = 0.0_wp
        status = cf_ok
    end subroutine initialize_axisymmetric_j2_sensitivity

    subroutine evaluate_axisymmetric_j2(self,x_old,x_new,state_old,state_trial,response,dt,tangent,point_status)
        !! 3D \(J_2\) trial response restricted to four axisymmetric components.
        class(axisymmetric_j2_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed inputs `(ninput,npoints)` when required, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `strain.small` inputs `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued history `(nstate,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Computed uncommitted real-valued history `(nstate,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.cauchy` responses `(noutput,npoints)`.
        real(wp), intent(in) :: dt
            !! Prescribed time increment; a positive value is required for time-dependent laws.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic tangent `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(inout) :: point_status(:)
            !! Per-point outcome; write `cf_ok` only for a successful constitutive update.

        call evaluate_points(self,x_new,state_old,state_trial,response,point_status,tangent)
    end subroutine evaluate_axisymmetric_j2

    subroutine sensitivity_axisymmetric_j2(self,x_old,x_new,state_old,state_trial,selected_parameters,sensitivity_old,sensitivity_trial, response,dsdp,dt,tangent,point_status)
        !! Full-history analytical \(J_2\) derivatives with exact axisymmetric packing.
        class(axisymmetric_j2_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed inputs `(ninput,npoints)` when required, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `strain.small` inputs `(ninput,npoints)`.
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
            !! Packed `stress.cauchy` responses `(noutput,npoints)`.
        real(wp), intent(out) :: dsdp(:,:,:)
            !! Exact total response derivatives `(noutput,nselected,npoints)` at fixed prescribed driving histories.
        real(wp), intent(in) :: dt
            !! Prescribed time increment; a positive value is required for time-dependent laws.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic tangent `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(inout) :: point_status(:)
            !! Per-point outcome; write `cf_ok` only for a successful constitutive update.

        call evaluate_points(self,x_new,state_old,state_trial,response,point_status,tangent,selected_parameters, sensitivity_old,sensitivity_trial,dsdp)
    end subroutine sensitivity_axisymmetric_j2

    subroutine evaluate_points(self,x_new,state_old,state_trial,response,point_status,tangent, selected_parameters,sensitivity_old,sensitivity_trial,dsdp)
        !! Fixed-workspace chunks avoid dynamic allocation for arbitrarily
        !! large material batches, while dispatching only once per chunk.
        class(axisymmetric_j2_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `strain.small` inputs `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued history `(nstate,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Computed uncommitted real-valued history `(nstate,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.cauchy` responses `(noutput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Per-point outcome; write `cf_ok` only for a successful constitutive update.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic tangent `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(in), optional :: selected_parameters(:)
            !! One-based parameter IDs in the caller's selected derivative-axis order.
        real(wp), intent(in), optional :: sensitivity_old(:,:,:)
            !! Committed real-history parameter derivatives `(nstate,nselected,npoints)`.
        real(wp), intent(out), optional :: sensitivity_trial(:,:,:)
            !! Updated real-history parameter derivatives `(nstate,nselected,npoints)`.
        real(wp), intent(out), optional :: dsdp(:,:,:)
            !! Exact total response derivatives `(noutput,nselected,npoints)` at fixed prescribed driving histories.
        type(j2_plasticity_t) :: base
        real(wp) :: dummy_input(0,0)
        real(wp) :: full_input(nmandel,chunk_size)
        real(wp) :: full_stress(nmandel,chunk_size)
        real(wp) :: full_tangent(nmandel,nmandel,chunk_size)
        real(wp) :: full_dsdp(nmandel,max_selected,chunk_size)
        integer :: local_status(chunk_size)
        integer :: first
        integer :: last
        integer :: n
        integer :: p
        integer :: i
        integer :: j
        integer :: k
        integer :: nselected

        base = j2_plasticity_t(self%bulk_modulus,self%shear_modulus,self%yield_stress,self%hardening_modulus)
        nselected = 0
        if (present(selected_parameters)) nselected = size(selected_parameters)
        point_loop: do first = 1,size(x_new,2),chunk_size
            last = min(first+chunk_size-1,size(x_new,2))
            n = last-first+1
            full_input(:,1:n) = 0.0_wp
            full_input(mapping,1:n) = x_new(:,first:last)
            local_status(1:n) = cf_ok
            if (present(selected_parameters)) then
                if (present(tangent)) then
                    call base%eval_sensitivity_batch(dummy_input,full_input(:,1:n),state_old(:,first:last), state_trial(:,first:last),selected_parameters,sensitivity_old(:,:,first:last), sensitivity_trial(:,:,first:last),full_stress(:,1:n),full_dsdp(:,1:nselected,1:n), 0.0_wp,full_tangent(:,:,1:n),local_status(1:n))
                else
                    call base%eval_sensitivity_batch(dummy_input,full_input(:,1:n),state_old(:,first:last), state_trial(:,first:last),selected_parameters,sensitivity_old(:,:,first:last), sensitivity_trial(:,:,first:last),full_stress(:,1:n),full_dsdp(:,1:nselected,1:n), 0.0_wp,point_status=local_status(1:n))
                end if
            else
                if (present(tangent)) then
                    call base%eval_batch(dummy_input,full_input(:,1:n),state_old(:,first:last), state_trial(:,first:last),full_stress(:,1:n),0.0_wp,full_tangent(:,:,1:n),local_status(1:n))
                else
                    call base%eval_batch(dummy_input,full_input(:,1:n),state_old(:,first:last), state_trial(:,first:last),full_stress(:,1:n),0.0_wp,point_status=local_status(1:n))
                end if
            end if

            result_loop: do p = 1,n
                point_status(first+p-1) = local_status(p)
                if (local_status(p) /= cf_ok) cycle result_loop
                response(:,first+p-1) = full_stress(mapping,p)
                if (present(tangent)) then
                    do j = 1,naxis
                        do i = 1,naxis
                            tangent(i,j,first+p-1) = full_tangent(mapping(i),mapping(j),p)
                        end do
                    end do
                end if
                if (present(selected_parameters)) then
                    do k = 1,nselected
                        dsdp(:,k,first+p-1) = full_dsdp(mapping,k,p)
                    end do
                end if
            end do result_loop
        end do point_loop
    end subroutine evaluate_points

end module constituf_axisymmetric_j2
