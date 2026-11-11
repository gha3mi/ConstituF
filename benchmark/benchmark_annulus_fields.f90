! SPDX-License-Identifier: Apache-2.0

module benchmark_annulus_fields
    !! Prescribed annular loading and radial response projection.
    !! Material response and derivatives are always evaluated by ConstituF.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use benchmark_annulus_mesh, only: r_inner, r_outer, interface_radius, nangle, cell_radius, cell_angle
    use benchmark_history, only: npath, path
    use benchmark_models, only: cases, history_case, interface_case, &
                              load_mandel, load_plane, load_axisym, load_F, load_cohesive_i, load_cohesive_mixed, &
                              load_thermoelastic, load_heat, load_heat_temperature, load_concentration, &
                              load_electric, load_piezo, load_pressure, load_pressure_pressure, load_biot
    use constituf, only: wp
    implicit none (type, external)
    private

    real(wp), parameter :: displacement_amplitude=2.5e-4_wp
    real(wp), parameter :: temperature_inner=360.0_wp, temperature_outer=300.0_wp
    real(wp), parameter :: concentration_inner=1000.0_wp, concentration_outer=500.0_wp
    real(wp), parameter :: electric_potential_inner=500.0_wp, electric_potential_outer=0.0_wp
    real(wp), parameter :: fluid_pressure_inner=2.0e6_wp, fluid_pressure_outer=1.0e6_wp
    real(wp), parameter :: root_two=sqrt(2.0_wp)

    public :: fill_annulus_inputs, project_annulus_response, annulus_displacement

contains

    subroutine fill_annulus_inputs(id,x,amplitude)
        !! Prescribe the same annular kinematic profile for all bulk
        !! mechanical materials. Scalar transport and coupled fields use
        !! prescribed logarithmic profiles, not solved boundary-value problems.
        integer, intent(in) :: id
        real(wp), intent(out) :: x(:,:)
        real(wp), intent(in) :: amplitude
        integer :: p, i, j
        real(wp) :: r, angle, c, s, rr, tt, exx, eyy, exy, radial_gradient
        real(wp) :: temperature, pressure, electric_field, grad_temperature
        real(wp) :: grad_concentration, grad_pressure, opening, logarithm

        logarithm=log(r_outer/r_inner)
        x=0.0_wp

        do p=1,size(x,2)
            j=mod(p-1,nangle)+1
            i=(p-1)/nangle+1
            r=cell_radius(i)
            if (interface_case(id)) r=interface_radius
            angle=cell_angle(j)
            c=cos(angle)
            s=sin(angle)
            radial_gradient=amplitude*displacement_amplitude*r_inner/r**2
            rr=-radial_gradient
            tt=radial_gradient
            exx=rr*c*c+tt*s*s
            eyy=rr*s*s+tt*c*c
            exy=(rr-tt)*c*s
            temperature=temperature_inner+(temperature_outer-temperature_inner)* &
                log(r/r_inner)/logarithm
            pressure=fluid_pressure_inner+(fluid_pressure_outer-fluid_pressure_inner)* &
                log(r/r_inner)/logarithm
            electric_field=(electric_potential_inner-electric_potential_outer)/(r*logarithm)
            grad_temperature=(temperature_outer-temperature_inner)/(r*logarithm)
            grad_concentration=(concentration_outer-concentration_inner)/(r*logarithm)
            grad_pressure=(fluid_pressure_outer-fluid_pressure_inner)/(r*logarithm)

            select case(cases(id)%input_kind)
            case(load_mandel)
                x(:,p)=[exx,eyy,0.0_wp,0.0_wp,0.0_wp,root_two*exy]
            case(load_plane)
                x(:,p)=[exx,eyy,root_two*exy]
            case(load_axisym)
                ! Custom axisymmetric Mandel order: rr, zz, theta-theta, sqrt(2)*rz.
                x(:,p)=[rr,0.0_wp,tt,0.0_wp]
            case(load_F)
                ! Fortran column-major, spatial/Cartesian F; no Mandel scaling.
                x(:,p)=[1.0_wp+exx,exy,0.0_wp,exy,1.0_wp+eyy, &
                         0.0_wp,0.0_wp,0.0_wp,1.0_wp]
            case(load_cohesive_i)
                opening=amplitude*3.0e-5_wp*(1.0_wp+0.15_wp*cos(2.0_wp*angle))
                x(1,p)=opening
            case(load_cohesive_mixed)
                opening=amplitude*3.0e-5_wp*(1.0_wp+0.15_wp*cos(2.0_wp*angle))
                x(:,p)=[opening,0.2_wp*opening,0.0_wp]
            case(load_thermoelastic)
                x(:,p)=[exx,eyy,0.0_wp,0.0_wp,0.0_wp,root_two*exy,temperature]
            case(load_heat)
                x(:,p)=[grad_temperature*c,grad_temperature*s,0.0_wp]
            case(load_heat_temperature)
                x(:,p)=[grad_temperature*c,grad_temperature*s,0.0_wp,temperature]
            case(load_concentration)
                x(:,p)=[grad_concentration*c,grad_concentration*s,0.0_wp]
            case(load_electric)
                x(:,p)=[electric_field*c,electric_field*s,0.0_wp]
            case(load_piezo)
                x(:,p)=[exx,eyy,0.0_wp,0.0_wp,0.0_wp,root_two*exy, &
                    electric_field*c,electric_field*s,0.0_wp]
            case(load_pressure)
                x(:,p)=[grad_pressure*c,grad_pressure*s,0.0_wp]
            case(load_pressure_pressure)
                x(:,p)=[grad_pressure*c,grad_pressure*s,0.0_wp,pressure]
            case(load_biot)
                x(:,p)=[exx,eyy,0.0_wp,0.0_wp,0.0_wp,root_two*exy,pressure]
            case default
                error stop "Annulus input case missing"
            end select
        end do
    end subroutine fill_annulus_inputs

    subroutine project_annulus_response(id,kind,response,field)
        !! Return a *physical* radial projection, not an arbitrary packed
        !! component. Piola and Cauchy projections remain distinctly labeled.
        integer, intent(in) :: id,kind
        real(wp), intent(in) :: response(:,:)
        real(wp), intent(out) :: field(:)
        real(wp) :: c,s,theta
        integer :: p,j

        if (size(field)/=size(response,2)) error stop "Annulus projection output shape mismatch"
        do p=1,size(field)
            j=mod(p-1,nangle)+1
            theta=cell_angle(j)
            c=cos(theta)
            s=sin(theta)
            select case(kind)
            case(1)
                field(p)=c*c*response(1,p)+s*s*response(2,p)+ &
                    root_two*c*s*response(6,p)
            case(2)
                field(p)=c*c*response(1,p)+s*s*response(2,p)+ &
                    root_two*c*s*response(3,p)
            case(3)
                field(p)=response(1,p)
            case(4)
                field(p)=c*c*response(1,p)+c*s*(response(2,p)+response(4,p))+ &
                    s*s*response(5,p)
            case(5)
                field(p)=c*response(1,p)+s*response(2,p)
            case(6)
                field(p)=response(1,p)
            case default
                error stop "Annulus radial projection missing"
            end select
        end do
        if (.not. all(ieee_is_finite(field))) error stop "Invalid annulus radial response"
    end subroutine project_annulus_response

    pure real(wp) function annulus_displacement(id,r) result(u)
        !! VTK displacement corresponds to the final accepted loading step.
        integer, intent(in) :: id
        real(wp), intent(in) :: r
        real(wp) :: amplitude

        u=0.0_wp
        if (.not.any(cases(id)%input_kind==[load_mandel,load_plane,load_axisym,load_F, &
                     load_thermoelastic,load_piezo,load_biot])) return
        amplitude=1.0_wp
        if (history_case(id)) amplitude=path(npath)
        u=amplitude*displacement_amplitude*r_inner/r
    end function annulus_displacement
end module benchmark_annulus_fields
