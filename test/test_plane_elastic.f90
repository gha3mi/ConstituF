! SPDX-License-Identifier: Apache-2.0

program test_plane_elastic
    !! Test convention: stress and stress-like material parameters in Pa.
    !! Strain and deformation gradients are dimensionless. Explicit rescaling
    !! cases state their factors.
    !!
    !! Independent 3D constraint reduction and derivative verification.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf, only: wp, material_t, verify_tangent, cf_ok, cf_invalid_argument
    use constituf_isotropic_elastic, only: isotropic_elastic_t
    use constituf_plane_elastic, only: plane_elastic_t
    implicit none (type, external)

    integer, parameter :: nmandel = 6
    real(wp), parameter :: bulk = 120.0_wp, shear = 45.0_wp
    type(material_t) :: plane, reference
    real(wp) :: x(3,2), y(3,2), jac(3,3,2), dsdp(3,2,2)
    real(wp) :: x3(nmandel,2), y3(nmandel,2)
    real(wp) :: lam, error
    real(wp) :: yp(3,2), ym(3,2)
    real(wp) :: scaled_jac(3,3,2), scaled_dsdp(3,2,2)
    real(wp) :: h, tolerance, modulus_scale
    integer :: constraint, status, iscale

    x(:,1) = [0.012_wp,-0.004_wp,0.007_wp]
    x(:,2) = [-0.003_wp,0.009_wp,-0.005_wp]
    lam = bulk-2.0_wp*shear/3.0_wp
    do constraint = 1, 2
        call plane%init(plane_elastic_t(bulk,shear,constraint==2),2,status, &
                        sensitivity_parameters=[character(len=15) :: "bulk_modulus", "shear_modulus"])
        if (status /= cf_ok) error stop "2D elastic initialization"
        if (.not. plane%field_matches("input","strain.small","mandel","small_strain","spatial",[2,2])) then
            error stop "2D Mandel field contract"
        end if
        call reference%init(isotropic_elastic_t(bulk,shear),2,status)
        if (status /= cf_ok) error stop "3D elastic initialization"
        x3 = 0.0_wp
        x3(1:2,:) = x(1:2,:)
        x3(6,:) = x(3,:)
        if (constraint == 2) then
            x3(3,:) = -lam*(x(1,:)+x(2,:))/(lam+2.0_wp*shear)
        end if
        call reference%eval(x3,y3,status=status)
        if (status /= cf_ok) error stop "3D reference evaluation"
        call plane%eval(x,y,jac,status,dsdp=dsdp)
        if (status /= cf_ok) error stop "2D elastic evaluation"
        call compare(y(1,:),y3(1,:))
        call compare(y(2,:),y3(2,:))
        call compare(y(3,:),y3(6,:))
        ! E = 120 Pa and nu = 1/3 give independent in-plane reference values.
        if (constraint == 1) then
            call compare(y(:,1),[1.8_wp,0.36_wp,0.63_wp])
            call compare(y(:,2),[0.27_wp,1.35_wp,-0.45_wp])
        else
            call compare(y(:,1),[1.44_wp,0.0_wp,0.63_wp])
            call compare(y(:,2),[0.0_wp,1.08_wp,-0.45_wp])
        end if
        if (constraint == 2) call compare(y3(3,:),[0.0_wp,0.0_wp])
        if (.not. all(ieee_is_finite(jac))) error stop "Nonfinite 2D tangent"
        if (.not. all(ieee_is_finite(dsdp))) error stop "Nonfinite 2D sensitivity"
        if (maxval(abs(jac(1,2,:)-jac(2,1,:))) > 1024.0_wp*epsilon(1.0_wp)) then
            error stop "Asymmetric 2D tangent"
        end if
        call plane%revert()

        call reference%init(plane_elastic_t(bulk,shear,constraint==2),2,status)
        if (status /= cf_ok) error stop "Tangent verification setup"
        call verify_tangent(reference,x,error,status)
        if (status /= cf_ok .or. .not. ieee_is_finite(error) .or. error > 1.0_wp) then
            error stop "2D tangent verification"
        end if

        h = epsilon(1.0_wp)**(1.0_wp/3.0_wp)*bulk
        call reference%init(plane_elastic_t(bulk+h,shear,constraint==2),2,status)
        if (status /= cf_ok) error stop "Bulk plus initialization"
        call reference%eval(x,yp,status=status)
        if (status /= cf_ok) error stop "Parameter plus evaluation"
        call reference%init(plane_elastic_t(bulk-h,shear,constraint==2),2,status)
        if (status /= cf_ok) error stop "Bulk minus initialization"
        call reference%eval(x,ym,status=status)
        if (status /= cf_ok) error stop "Parameter minus evaluation"
        ! Balance central-difference truncation and roundoff for the selected precision.
        tolerance = 8.0_wp*(h/bulk)**2*maxval(abs(x)) + &
                    16.0_wp*epsilon(1.0_wp)*max(1.0_wp,maxval(abs(yp)),maxval(abs(ym)))/h
        if (.not. all(ieee_is_finite((yp-ym)/(2.0_wp*h)))) error stop "Nonfinite bulk difference"
        if (maxval(abs((yp-ym)/(2.0_wp*h)-dsdp(:,1,:))) > tolerance) error stop "d stress / d bulk modulus"

        h = epsilon(1.0_wp)**(1.0_wp/3.0_wp)*shear
        call reference%init(plane_elastic_t(bulk,shear+h,constraint==2),2,status)
        if (status /= cf_ok) error stop "Shear plus initialization"
        call reference%eval(x,yp,status=status)
        if (status /= cf_ok) error stop "Shear plus evaluation"
        call reference%init(plane_elastic_t(bulk,shear-h,constraint==2),2,status)
        if (status /= cf_ok) error stop "Shear minus initialization"
        call reference%eval(x,ym,status=status)
        if (status /= cf_ok) error stop "Shear minus evaluation"
        tolerance = 8.0_wp*(h/shear)**2*maxval(abs(x)) + &
                    16.0_wp*epsilon(1.0_wp)*max(1.0_wp,maxval(abs(yp)),maxval(abs(ym)))/h
        if (.not. all(ieee_is_finite((yp-ym)/(2.0_wp*h)))) error stop "Nonfinite shear difference"
        if (maxval(abs((yp-ym)/(2.0_wp*h)-dsdp(:,2,:))) > tolerance) error stop "d stress / d shear modulus"

        if (constraint == 2) then
            do iscale = 1, 2
                if (iscale == 1) then
                    modulus_scale = sqrt(tiny(1.0_wp))*epsilon(1.0_wp)
                else
                    modulus_scale = sqrt(huge(1.0_wp))
                end if
                call reference%init(plane_elastic_t(bulk*modulus_scale,shear*modulus_scale,.true.),2,status, &
                    sensitivity_parameters=[character(len=15) :: "bulk_modulus", "shear_modulus"])
                if (status /= cf_ok) error stop "Scaled plane-stress initialization"
                call reference%eval(x,yp,scaled_jac,status,dsdp=scaled_dsdp)
                if (status /= cf_ok) error stop "Scaled plane-stress evaluation"
                call compare(reshape(yp/modulus_scale,[6]),reshape(y,[6]))
                call compare(reshape(scaled_jac/modulus_scale,[18]),reshape(jac,[18]))
                call compare(reshape(scaled_dsdp,[12]),reshape(dsdp,[12]))
            end do
        end if
    end do

    call test_auxetic_response()
    call plane%init(plane_elastic_t(bulk,0.0_wp),2,status)
    if (status /= cf_invalid_argument) error stop "Invalid shear modulus accepted"

contains

    subroutine compare(actual, expected)
        real(wp), intent(in) :: actual(:), expected(:)

        if (.not. all(ieee_is_finite(actual)) .or. .not. all(ieee_is_finite(expected))) then
            error stop "Nonfinite plane-elastic comparison"
        end if
        if (any(abs(actual-expected) > 2000.0_wp*epsilon(1.0_wp)* &
            max(1.0_wp,maxval(abs(expected))))) error stop "3D and 2D elastic mismatch"
    end subroutine compare

    subroutine test_auxetic_response()
        type(material_t) :: material
        real(wp), parameter :: young = 45.0_wp, poisson = -0.25_wp
        real(wp) :: strain(3,2), stress(3,2)
        integer :: mode, result

        strain(:,1) = [0.001_wp,-poisson*0.001_wp,0.0_wp]
        ! Physical tensor shear epsilon_12 = 0.001 gives sigma_12 = 0.06 Pa.
        strain(:,2) = [0.0_wp,0.0_wp,sqrt(2.0_wp)*0.001_wp]
        do mode = 1, 2
            call material%init(plane_elastic_t(young/(3.0_wp*(1.0_wp-2.0_wp*poisson)), &
                young/(2.0_wp*(1.0_wp+poisson)),mode==2),2,result)
            if (result /= cf_ok) error stop "Auxetic plane-elastic initialization"
            call material%eval(strain,stress,status=result)
            if (result /= cf_ok) error stop "Auxetic plane-elastic evaluation"
            if (mode == 1) then
                call compare(stress(:,1),[0.0475_wp,0.0025_wp,0.0_wp])
            else
                call compare(stress(:,1),[0.045_wp,0.0_wp,0.0_wp])
            end if
            call compare(stress(:,2),[0.0_wp,0.0_wp,sqrt(2.0_wp)*0.06_wp])
        end do
    end subroutine test_auxetic_response
end program test_plane_elastic
