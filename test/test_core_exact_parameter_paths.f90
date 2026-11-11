! SPDX-License-Identifier: Apache-2.0

program test_core_exact_parameter_paths
    !! Independent analytical identities for all six newly differentiated laws.
    !! No numerical differencing or perturbation of material parameters is used.
    !! Synthetic data use Pa, m and s. Cohesive stiffness uses Pa/m and viscosity Pa*s.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf, only: wp, material_t, cf_ok
    use constituf_chaboche, only: chaboche_t
    use constituf_cohesive_mixed_mode, only: cohesive_mixed_mode_t
    use constituf_drucker_prager, only: drucker_prager_t
    use constituf_ductile_damage, only: ductile_damage_t
    use constituf_j2_plane, only: j2_plane_t
    use constituf_j2_plasticity, only: j2_plasticity_t
    use constituf_perzyna, only: perzyna_t
    implicit none (type, external)

    integer, parameter :: nmandel=6, noutput_damage=8
    ! Independent formulas allow a small roundoff budget per accepted update.
    real(wp), parameter :: roundoff_factor=32.0_wp

    call verify_cohesive()
    call verify_drucker_prager()
    call verify_chaboche()
    call verify_perzyna()
    call verify_plane_j2()
    call verify_damage()

contains

    subroutine check(actual,expected,context,comparison_scale)
        real(wp),intent(in) :: actual,expected
        character(*),intent(in) :: context
        real(wp),intent(in),optional :: comparison_scale
        real(wp) :: scale,roundoff_scale,tolerance

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
            error stop "Nonfinite exact sensitivity reference"
        end if
        scale=max(1.0e-12_wp,abs(expected))
        roundoff_scale=abs(expected)
        if (present(comparison_scale)) then
            if (.not. ieee_is_finite(comparison_scale)) error stop "Nonfinite sensitivity comparison scale"
            roundoff_scale=max(roundoff_scale,abs(comparison_scale))
        end if
        ! Retain the original accuracy requirement when it exceeds working
        ! roundoff. Use uncancelled term scales where the identity subtracts.
        tolerance=max(1.0e-7_wp*scale,roundoff_factor*epsilon(1.0_wp)*roundoff_scale)
        if (abs(actual-expected)<=tolerance) return
        print "(a,2(1x,es17.9))",context,actual,expected
        error stop "Exact sensitivity reference mismatch"
    end subroutine check

    subroutine ok(status,context)
        integer,intent(in) :: status
        character(*),intent(in) :: context
        if (status==cf_ok) return
        print "(a,1x,i0)",context,status
        error stop "Analytical parameter sensitivity failed"
    end subroutine ok

    subroutine verify_cohesive()
        !! Exact loading and unloading r(kappa,beta), including history dκ/dβ.
        real(wp),parameter :: k=1000.0_wp,d0=0.01_wp,df=0.04_wp,beta=1.5_wp
        type(material_t) :: mat
        real(wp) :: x(3,1),y(3,1),d(3,2,1)
        real(wp) :: kappa,dkappa,r,dr
        integer :: status,step

        call mat%init(cohesive_mixed_mode_t(k,d0,df,beta),1,status, &
            sensitivity_parameters=[character(len=13) :: "stiffness","shear_weight"])
        call ok(status,"Mixed cohesive initialization")
        x(:,1)=[0.020_wp,0.004_wp,0.0_wp]
        kappa=sqrt(x(1,1)**2+(beta*x(2,1))**2)
        dkappa=beta*x(2,1)**2/kappa
        do step=1,2
            if (step==2) x(:,1)=[0.010_wp,0.001_wp,0.0_wp]
            r=d0*(df-kappa)/(kappa*(df-d0))
            dr=-d0*df*dkappa/((df-d0)*kappa**2)
            call mat%eval(x,y,status=status,dsdp=d)
            call ok(status,"Mixed cohesive parameter update")
            call check(d(1,1,1),y(1,1)/k,"Mixed cohesive stiffness derivative")
            call check(d(1,2,1),k*dr*x(1,1),"Mixed cohesive shear-weight normal derivative")
            call check(d(2,2,1),k*(dr*beta**2+2.0_wp*r*beta)*x(2,1),"Mixed cohesive shear-weight shear", &
                k*(abs(dr)*beta**2+2.0_wp*abs(r)*beta)*abs(x(2,1)))
            call mat%commit(status)
            call ok(status,"Mixed cohesive commit")
        end do
    end subroutine verify_cohesive

    subroutine verify_drucker_prager()
        !! Hydrostatic apex pressure-slope formula, then J2 limit with
        !! complete shear-modulus and yield-stress derivative histories.
        real(wp),parameter :: bulk=1000.0_wp,shear=400.0_wp,yield_stress=0.2_wp
        real(wp),parameter :: beta=0.3_wp,hard=50.0_wp
        type(material_t) :: mat,dp,j2
        real(wp) :: x(nmandel,1),y(nmandel,1),z(nmandel,1),d(nmandel,1,1),d_dp(nmandel,2,1),d_j2(nmandel,2,1)
        real(wp) :: mean_trial,den,delta_gamma,dgamma,expected,trial_q
        real(wp) :: trial_dev(nmandel)
        real(wp) :: stress_scale,strain_scale,roundoff,derivative_scale
        integer :: status,step,iparam

        call mat%init(drucker_prager_t(bulk,shear,yield_stress,beta,hard),1,status, &
            sensitivity_parameters=[character(len=14) :: "pressure_slope"])
        call ok(status,"Drucker-Prager apex setup")
        x(:,1)=[0.02_wp,0.02_wp,0.02_wp,0.0_wp,0.0_wp,0.0_wp]
        call mat%eval(x,y,status=status,dsdp=d)
        call ok(status,"Drucker-Prager apex derivative")
        mean_trial=bulk*sum(x(1:3,1))
        den=bulk*beta**2+hard
        delta_gamma=(beta*mean_trial-yield_stress)/den
        dgamma=(mean_trial-2.0_wp*bulk*beta*delta_gamma)/den
        expected=-bulk*(delta_gamma+beta*dgamma)
        call check(d(1,1,1),expected,"Drucker-Prager apex dmean/dslope", &
            bulk*(abs(delta_gamma)+abs(beta*dgamma)))
        call check(d(2,1,1),expected,"Drucker-Prager isotropic apex derivative", &
            bulk*(abs(delta_gamma)+abs(beta*dgamma)))

        ! Exact smooth-cone pressure-slope sensitivity. The yield
        ! direction does not depend on slope for a virgin state.
        call mat%init(drucker_prager_t(bulk,shear,yield_stress,beta,hard),1,status, &
            sensitivity_parameters=[character(len=14) :: "pressure_slope"])
        call ok(status,"Drucker-Prager smooth-cone setup")
        x(:,1)=[0.001_wp,-0.0005_wp,-0.0005_wp,0.02_wp,0.002_wp,0.0_wp]
        call mat%eval(x,y,status=status,dsdp=d)
        call ok(status,"Drucker-Prager smooth-cone derivative")
        trial_dev=2.0_wp*shear*x(:,1)
        trial_q=sqrt(1.5_wp*dot_product(trial_dev,trial_dev))
        mean_trial=bulk*sum(x(1:3,1))
        den=3.0_wp*shear+bulk*beta**2+hard
        delta_gamma=(trial_q+beta*mean_trial-yield_stress)/den
        if (trial_q<=3.0_wp*shear*delta_gamma) error stop "Smooth-cone oracle entered apex"
        dgamma=(mean_trial-2.0_wp*bulk*beta*delta_gamma)/den
        expected=-bulk*(delta_gamma+beta*dgamma)
        call check(sum(d(1:3,1,1))/3.0_wp,expected,"Drucker-Prager smooth-cone hydrostatic derivative", &
            bulk*(abs(delta_gamma)+abs(beta*dgamma)))

        call dp%init(drucker_prager_t(bulk,shear,yield_stress,0.0_wp,hard),1,status, &
            sensitivity_parameters=[character(len=13) :: "shear_modulus","yield_stress"])
        call ok(status,"Drucker-Prager J2-limit setup")
        call j2%init(j2_plasticity_t(bulk,shear,yield_stress,hard),1,status, &
            sensitivity_parameters=[character(len=13) :: "shear_modulus","yield_stress"])
        call ok(status,"J2 sensitivity reference setup")
        strain_scale=0.0_wp
        stress_scale=yield_stress
        do step=1,3
            x(:,1)=real(4-step,wp)*[0.005_wp,-0.0025_wp,-0.0025_wp,0.002_wp,0.0_wp,0.0_wp]
            if (step==2) x=-0.65_wp*x
            call dp%eval(x,y,status=status,dsdp=d_dp)
            call ok(status,"Drucker-Prager limit step")
            call j2%eval(x,z,status=status,dsdp=d_j2)
            call ok(status,"J2 analytic sensitivity limit step")
            if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(z)) .or. &
                .not. all(ieee_is_finite(d_dp)) .or. .not. all(ieee_is_finite(d_j2))) then
                error stop "Nonfinite Drucker-Prager J2 limit"
            end if
            ! Plastic corrections subtract elastic trial stresses. Retain their
            ! scale across the accepted path, rather than only the final stress.
            strain_scale=max(strain_scale,maxval(abs(x)))
            stress_scale=max(stress_scale,bulk*abs(sum(x(1:3,1))),2.0_wp*shear*strain_scale, &
                maxval(abs(y)),maxval(abs(z)))
            roundoff=roundoff_factor*real(step,wp)*epsilon(1.0_wp)
            if (maxval(abs(y-z))>max(1.0e-8_wp,roundoff*stress_scale)) error stop "Drucker-Prager J2 stress limit"
            do iparam=1,2
                derivative_scale=2.0_wp*strain_scale
                if (iparam==2) derivative_scale=1.0_wp
                derivative_scale=max(derivative_scale,maxval(abs(d_dp(:,iparam,:))),maxval(abs(d_j2(:,iparam,:))))
                if (maxval(abs(d_dp(:,iparam,:)-d_j2(:,iparam,:)))>max(1.0e-8_wp,roundoff*derivative_scale)) then
                    error stop "Drucker-Prager J2 derivative limit"
                end if
            end do
            call dp%commit(status)
            call ok(status,"Drucker-Prager commit")
            call j2%commit(status)
            call ok(status,"J2 commit")
        end do
    end subroutine verify_drucker_prager

    subroutine verify_chaboche()
        !! For C=gamma=0 the complete plastic-history derivatives must
        !! reproduce the independent full-history J2 analytical update.
        real(wp),parameter :: bulk=1000.0_wp,shear=400.0_wp,yield_stress=0.2_wp,hard=40.0_wp
        type(material_t) :: mat,j2
        integer, parameter :: nstate_chaboche=13
        real(wp), parameter :: kinematic=150.0_wp,recovery=5.0_wp
        real(wp) :: x(nmandel,1),y(nmandel,1),z(nmandel,1),d(nmandel,4,1),reference(nmandel,4,1)
        real(wp) :: d_active(nmandel,2,1),history(nstate_chaboche,1),history_sensitivity(nstate_chaboche,2,1)
        real(wp) :: trial(nmandel),normal(nmandel)
        real(wp) :: q,r,fac,den,dr_kinematic,dr_recovery
        real(wp) :: stress_scale,strain_scale,roundoff,derivative_scale
        integer :: step,status,iparam

        call mat%init(chaboche_t(bulk,shear,yield_stress,hard,0.0_wp,0.0_wp),1,status, &
            sensitivity_parameters=[character(len=17) :: "bulk_modulus","shear_modulus", &
                "yield_stress","hardening_modulus"])
        call ok(status,"Chaboche analytic setup")
        call j2%init(j2_plasticity_t(bulk,shear,yield_stress,hard),1,status, &
            sensitivity_parameters=[character(len=17) :: "bulk_modulus","shear_modulus", &
                "yield_stress","hardening_modulus"])
        call ok(status,"J2 analytic setup")
        strain_scale=0.0_wp
        stress_scale=yield_stress
        do step=1,4
            select case(step)
            case(1)
                x(:,1)=[0.01_wp,-0.005_wp,-0.005_wp,0.02_wp,0.005_wp,0.0_wp]
            case(2)
                x=-0.65_wp*x
            case(3)
                x=-1.5_wp*x
            case default
                x=0.35_wp*x
            end select
            call mat%eval(x,y,status=status,dsdp=d)
            call ok(status,"Chaboche exact sensitivity")
            call j2%eval(x,z,status=status,dsdp=reference)
            call ok(status,"J2 exact sensitivity reference")
            if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(z)) .or. &
                .not. all(ieee_is_finite(d)) .or. .not. all(ieee_is_finite(reference))) then
                error stop "Nonfinite Chaboche J2 limit"
            end if
            strain_scale=max(strain_scale,maxval(abs(x)))
            stress_scale=max(stress_scale,bulk*abs(sum(x(1:3,1))),2.0_wp*shear*strain_scale, &
                maxval(abs(y)),maxval(abs(z)))
            roundoff=roundoff_factor*real(step,wp)*epsilon(1.0_wp)
            if (maxval(abs(y-z))>max(1.0e-8_wp,roundoff*stress_scale)) error stop "Chaboche/J2 stress limit"
            do iparam=1,size(d,2)
                ! Modulus derivatives scale with strain; the yield-stress
                ! derivative has a dimensionless scale of order one.
                derivative_scale=2.0_wp*strain_scale
                if (iparam==1) derivative_scale=3.0_wp*strain_scale
                if (iparam==3) derivative_scale=1.0_wp
                derivative_scale=max(derivative_scale,maxval(abs(d(:,iparam,:))),maxval(abs(reference(:,iparam,:))))
                if (maxval(abs(d(:,iparam,:)-reference(:,iparam,:)))>max(2.0e-8_wp,roundoff*derivative_scale)) then
                    error stop "Chaboche/J2 analytic path limit"
                end if
            end do
            call mat%commit(status)
            call ok(status,"Chaboche history commit")
            call j2%commit(status)
            call ok(status,"J2 history commit")
        end do

        ! Exact nonzero kinematic/recovery parameter derivatives of the
        ! virgin implicit scalar residual R(r)=0. This exercises the
        ! nonlinear recovery denominator, not merely the J2 limit.
        call mat%init(chaboche_t(bulk,shear,yield_stress,hard,kinematic,recovery),1,status, &
            sensitivity_parameters=[character(len=17) :: "kinematic_modulus","recovery"])
        call ok(status,"Nonlinear Chaboche sensitivity setup")
        x(:,1)=[0.01_wp,-0.005_wp,-0.005_wp,0.02_wp,0.005_wp,0.0_wp]
        call mat%eval(x,y,status=status,dsdp=d_active)
        call ok(status,"Nonlinear Chaboche implicit derivatives")
        call mat%commit(status)
        call ok(status,"Nonlinear Chaboche commit")
        call mat%export_committed(history,status,sensitivity_state=history_sensitivity)
        call ok(status,"Nonlinear Chaboche history")
        r=history(nstate_chaboche,1)
        if (r<=0.0_wp) error stop "Chaboche analytical oracle needs plastic loading"
        fac=1.0_wp+recovery*r
        den=3.0_wp*shear+hard+kinematic/fac**2
        dr_kinematic=-r/(fac*den)
        dr_recovery=kinematic*r**2/(fac**2*den)
        trial=2.0_wp*shear*x(:,1)
        q=sqrt(1.5_wp*dot_product(trial,trial))
        normal=1.5_wp*trial/q
        call check(history_sensitivity(nstate_chaboche,1,1),dr_kinematic, &
            "Chaboche committed kinematic history derivative")
        call check(history_sensitivity(nstate_chaboche,2,1),dr_recovery, &
            "Chaboche committed recovery history derivative")
        call check(d_active(1,1,1),-2.0_wp*shear*dr_kinematic*normal(1), &
            "Chaboche nonzero kinematic derivative")
        call check(d_active(1,2,1),-2.0_wp*shear*dr_recovery*normal(1), &
            "Chaboche nonzero recovery derivative")
    end subroutine verify_chaboche

    subroutine verify_perzyna()
        !! Analytic derivative in the active overstress branch and exact
        !! equality of parameter-path derivatives across substep composition.
        real(wp),parameter :: bulk=160000.0_wp,shear=80000.0_wp,yield_stress=350.0_wp,hard=2000.0_wp
        real(wp),parameter :: viscosity=5000.0_wp,dt=0.2_wp
        integer,parameter :: nsub=4
        type(material_t) :: mat,automated,manual
        real(wp) :: initial(nmandel,1),x(nmandel,1),z(nmandel,1),y(nmandel,1),d(nmandel,1,1),d_auto(nmandel,1,1)
        real(wp) :: stress_trial(nmandel),normal(nmandel)
        real(wp) :: q,den,ddgamma,reference,derivative_scale,roundoff
        integer :: status,step

        initial=0.0_wp
        x(:,1)=[0.010_wp,-0.004_wp,-0.006_wp,0.001_wp,0.0_wp,0.001_wp]
        call mat%init(perzyna_t(bulk,shear,yield_stress,hard,viscosity),1,status,initial_input=initial, &
            sensitivity_parameters=[character(len=9) :: "viscosity"])
        call ok(status,"Perzyna analytic setup")
        call mat%eval(x,y,status=status,dt=dt,dsdp=d)
        call ok(status,"Perzyna analytical viscosity response")
        stress_trial=2.0_wp*shear*x(:,1)
        q=sqrt(1.5_wp*dot_product(stress_trial,stress_trial))
        normal=1.5_wp*stress_trial/q
        den=viscosity+dt*(3.0_wp*shear+hard)
        ddgamma=-dt*(q-yield_stress)/den**2
        reference=-2.0_wp*shear*ddgamma*normal(1)
        call check(d(1,1,1),reference,"Perzyna exact backward-Euler viscosity derivative")

        call automated%init(perzyna_t(bulk,shear,yield_stress,hard,viscosity,0.003_wp,nsub), &
            1,status,initial_input=initial,sensitivity_parameters=[character(len=9) :: "viscosity"])
        call ok(status,"Perzyna subdivision setup")
        call manual%init(perzyna_t(bulk,shear,yield_stress,hard,viscosity),1,status, &
            initial_input=initial,sensitivity_parameters=[character(len=9) :: "viscosity"])
        call ok(status,"Perzyna manual subdivision setup")
        call automated%eval(x,y,status=status,dt=dt,dsdp=d_auto)
        call ok(status,"Perzyna complete-increment derivative")
        do step=1,nsub
            z=real(step,wp)/real(nsub,wp)*x
            call manual%eval(z,y,status=status,dt=dt/real(nsub,wp),dsdp=d)
            call ok(status,"Perzyna substep reference")
            call manual%commit(status)
            call ok(status,"Perzyna manual substep commit")
        end do
        if (.not. all(ieee_is_finite(d_auto)) .or. .not. all(ieee_is_finite(d))) then
            error stop "Nonfinite Perzyna substep derivatives"
        end if
        derivative_scale=max(abs(reference),maxval(abs(d_auto)),maxval(abs(d)))
        roundoff=roundoff_factor*real(nsub,wp)*epsilon(1.0_wp)
        if (maxval(abs(d_auto-d))>max(2.0e-10_wp,roundoff*derivative_scale)) then
            error stop "Perzyna substep chain-rule mismatch"
        end if
    end subroutine verify_perzyna

    subroutine verify_plane_j2()
        !! An independently evaluated 3D analytic J2 parameter gradient
        !! is condensed with the same exact implicit stress constraint.
        real(wp),parameter :: bulk=150000.0_wp,shear=75000.0_wp,sy=350.0_wp,hard=2000.0_wp
        integer,parameter :: npar=4
        integer,parameter :: inplane(3)=[1,2,6]
        type(material_t) :: plane,three_d
        real(wp) :: x(3,1),y(3,1),d(3,npar,1),state(noutput_damage,1),x3(nmandel,1),z(nmandel,1)
        real(wp) :: history_sensitivity(noutput_damage,npar,1)
        real(wp) :: tangent(nmandel,nmandel,1),d3(nmandel,npar,1)
        real(wp) :: expected
        integer :: status,i,j

        call plane%init(j2_plane_t(bulk,shear,sy,hard,.true.),1,status, &
            sensitivity_parameters=[character(len=17) :: "bulk_modulus","shear_modulus", &
                "yield_stress","hardening_modulus"])
        call ok(status,"Plane J2 sensitivity setup")
        x(:,1)=[0.008_wp,-0.003_wp,0.0015_wp]
        call plane%eval(x,y,status=status,dsdp=d)
        call ok(status,"Plane-stress parameter derivatives")
        call plane%commit(status)
        call ok(status,"Plane-stress state commit")
        call plane%export_committed(state,status,sensitivity_state=history_sensitivity)
        call ok(status,"Plane-stress out-of-plane state")
        x3=0.0_wp
        x3(inplane,1)=x(:,1)
        x3(3,1)=state(8,1)
        call three_d%init(j2_plasticity_t(bulk,shear,sy,hard),1,status, &
            sensitivity_parameters=[character(len=17) :: "bulk_modulus","shear_modulus", &
                "yield_stress","hardening_modulus"])
        call ok(status,"3D sensitivity comparison setup")
        call three_d%eval(x3,z,tangent,status,dsdp=d3)
        call ok(status,"3D analytical reference")
        do j=1,npar
            do i=1,3
                expected=d3(inplane(i),j,1)-tangent(inplane(i),3,1)*d3(3,j,1)/tangent(3,3,1)
                call check(d(i,j,1),expected,"Plane-stress implicit-constraint derivative", &
                    abs(d3(inplane(i),j,1))+abs(tangent(inplane(i),3,1)*d3(3,j,1)/tangent(3,3,1)))
            end do
        end do
    end subroutine verify_plane_j2

    subroutine verify_damage()
        !! Damage kinetics parameters cannot affect effective J2 plasticity,
        !! so stress/energy/dissipation follow three exact closed forms.
        real(wp),parameter :: bulk=150000.0_wp,shear=75000.0_wp,sy=350.0_wp,hard=2000.0_wp
        real(wp),parameter :: alpha0=0.001_wp,rate=150.0_wp,cap=0.85_wp
        type(material_t) :: mat
        real(wp) :: x(nmandel,1),y(noutput_damage,1),d(noutput_damage,3,1),state(noutput_damage,1)
        real(wp) :: history_sensitivity(noutput_damage,3,1)
        real(wp) :: alpha,damage,offset,exponential,effective_stress,energy
        real(wp) :: derivative(3)
        integer :: status,j,i

        call mat%init(ductile_damage_t(bulk,shear,sy,hard,alpha0,rate,cap),1,status, &
            sensitivity_parameters=[character(len=16) :: "damage_threshold","damage_rate","maximum_damage"])
        call ok(status,"Ductile damage sensitivity setup")
        x(:,1)=[0.012_wp,-0.006_wp,-0.006_wp,0.002_wp,0.0_wp,0.0_wp]
        call mat%eval(x,y,status=status,dsdp=d)
        call ok(status,"Damage exact energy and dissipation sensitivities")
        call mat%commit(status)
        call ok(status,"Damage state commit")
        call mat%export_committed(state,status,sensitivity_state=history_sensitivity)
        call ok(status,"Damage state export")
        alpha=state(7,1)
        damage=state(8,1)
        if (alpha<=alpha0) error stop "Damage reference requires active damage"
        offset=alpha-alpha0
        exponential=exp(-rate*offset)
        derivative=[-cap*rate*exponential,cap*offset*exponential,1.0_wp-exponential]
        do j=1,3
            call check(history_sensitivity(noutput_damage,j,1),derivative(j), &
                "Damage committed history derivative")
        end do
        energy=y(7,1)/(1.0_wp-damage)
        do j=1,3
            do i=1,6
                effective_stress=y(i,1)/(1.0_wp-damage)
                call check(d(i,j,1),-effective_stress*derivative(j),"Damage stress parameter derivative")
            end do
            call check(d(7,j,1),-energy*derivative(j),"Damage stored-energy parameter derivative")
            call check(d(8,j,1),(energy-sy*alpha)*derivative(j),"Damage dissipation parameter derivative", &
                (abs(energy)+abs(sy*alpha))*abs(derivative(j)))
        end do
    end subroutine verify_damage
end program test_core_exact_parameter_paths
