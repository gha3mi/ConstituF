! SPDX-License-Identifier: Apache-2.0

module benchmark_cavity3d_fields
    !! Analytically differentiated 3D prescribed kinematics and scalar fields.
    !! No global equilibrium, transport equation, or material law is solved here.
    use benchmark_cavity3d_mesh, only: cavity_radius, half_width, unit_normal, radial_coordinate
    use benchmark_models, only: cases, ncases, load_mandel, load_plane, load_axisym, load_F, &
                              load_cohesive_i, load_cohesive_mixed, load_thermoelastic, &
                              load_heat, load_heat_temperature, load_concentration, &
                              load_electric, load_piezo, load_pressure, load_pressure_pressure, load_biot
    use constituf, only: wp
    implicit none (type, external)
    private

    integer, parameter :: n_mandel=6, n_tensor=9
    real(wp), parameter :: sqrt_two=sqrt(2.0_wp)
    real(wp), parameter :: cavity_opening=30.0e-6_wp
    real(wp), parameter :: macro_strain(3,3)=reshape([ &
        1.2e-3_wp,0.35e-3_wp,-0.2e-3_wp, &
        0.35e-3_wp,-0.5e-3_wp,0.25e-3_wp, &
        -0.2e-3_wp,0.25e-3_wp,0.2e-3_wp],[3,3])

    public :: supports_cavity3d, fill_cavity_inputs, project_cavity_response
    public :: cavity_displacement

contains

    pure logical function supports_cavity3d(id) result(supported)
        integer, intent(in) :: id

        supported=.false.
        if (id<1.or.id>ncases) return
        supported=cases(id)%input_kind/=load_plane.and.cases(id)%input_kind/=load_axisym
    end function supports_cavity3d

    pure function cavity_displacement(position,amplitude) result(displacement)
        real(wp), intent(in) :: amplitude
        real(wp), intent(in) :: position(3)
        real(wp) :: radius
        real(wp) :: displacement(3)

        radius=radial_coordinate(position)
        displacement=amplitude*(matmul(macro_strain,position)+ &
            cavity_opening*cavity_radius**2*position/radius**3)
    end function cavity_displacement

    pure subroutine kinematics(position,amplitude,strain,deformation)
        real(wp), intent(in) :: amplitude
        real(wp), intent(in) :: position(3)
        real(wp), intent(out) :: strain(n_mandel),deformation(n_tensor)
        real(wp) :: factor, radius
        real(wp) :: gradient(3,3), identity(3,3), symmetric(3,3)
        integer :: i,j

        radius=radial_coordinate(position)
        factor=cavity_opening*cavity_radius**2/radius**3
        identity=0.0_wp
        do i=1,3
            identity(i,i)=1.0_wp
        end do
        gradient=amplitude*(macro_strain+factor*identity)
        do j=1,3
            do i=1,3
                gradient(i,j)=gradient(i,j)- &
                    3.0_wp*amplitude*factor*position(i)*position(j)/radius**2
            end do
        end do
        symmetric=0.5_wp*(gradient+transpose(gradient))
        strain=[symmetric(1,1),symmetric(2,2),symmetric(3,3), &
                sqrt_two*symmetric(2,3),sqrt_two*symmetric(1,3),sqrt_two*symmetric(1,2)]
        deformation=reshape(identity+gradient,[9])
    end subroutine kinematics

    subroutine fill_cavity_inputs(id,position,amplitude,input)
        integer, intent(in) :: id
        real(wp), intent(in) :: amplitude
        real(wp), intent(in) :: position(:,:)
        real(wp), intent(out) :: input(:,:)
        real(wp) :: radius
        real(wp) :: x(3), normal(3), strain(n_mandel), deformation(n_tensor)
        real(wp) :: temperature, pressure
        real(wp) :: electric(3)
        real(wp) :: temperature_gradient(3),concentration_gradient(3),pressure_gradient(3)
        real(wp) :: opening
        integer :: p,kind

        if (id<1.or.id>ncases.or.size(position,1)/=3.or.size(position,2)/=size(input,2)) then
            error stop "Invalid cavity input dimensions"
        end if
        kind=cases(id)%input_kind
        if (.not.supports_cavity3d(id)) error stop "2D material passed to three-dimensional benchmark"
        input=0.0_wp
        do p=1,size(input,2)
            x=position(:,p)
            radius=radial_coordinate(x)
            normal=unit_normal(x)
            call kinematics(x,amplitude,strain,deformation)

            temperature=300.0_wp+45.0_wp*cavity_radius/radius+15.0_wp*x(1)/half_width
            temperature_gradient=-45.0_wp*cavity_radius*x/radius**3
            temperature_gradient(1)=temperature_gradient(1)+15.0_wp/half_width
            concentration_gradient=-500.0_wp*cavity_radius*x/radius**3
            concentration_gradient(3)=concentration_gradient(3)+80.0_wp/half_width
            pressure=1.0e6_wp+1.0e6_wp*cavity_radius/radius+0.2e6_wp*x(2)/half_width
            pressure_gradient=-1.0e6_wp*cavity_radius*x/radius**3
            pressure_gradient(2)=pressure_gradient(2)+0.2e6_wp/half_width
            electric=400.0_wp*cavity_radius*x/radius**3
            electric(1)=electric(1)-70.0_wp/half_width

            select case(kind)
            case(load_mandel)
                input(:,p)=strain
            case(load_F)
                input(:,p)=deformation
            case(load_cohesive_i,load_cohesive_mixed)
                opening=amplitude*30.0e-6_wp*(1.0_wp+ &
                    0.15_wp*(3.0_wp*normal(3)**2-1.0_wp)/2.0_wp)
                input(1,p)=opening
                if (kind==load_cohesive_mixed) input(2,p)=0.2_wp*opening
            case(load_thermoelastic)
                input(:,p)=[strain,temperature]
            case(load_heat)
                input(:,p)=temperature_gradient
            case(load_heat_temperature)
                input(:,p)=[temperature_gradient,temperature]
            case(load_concentration)
                input(:,p)=concentration_gradient
            case(load_electric)
                input(:,p)=electric
            case(load_piezo)
                input(:,p)=[strain,electric]
            case(load_pressure)
                input(:,p)=pressure_gradient
            case(load_pressure_pressure)
                input(:,p)=[pressure_gradient,pressure]
            case(load_biot)
                input(:,p)=[strain,pressure]
            case default
                error stop "Unsupported three-dimensional material input"
            end select
        end do
    end subroutine fill_cavity_inputs

    subroutine project_cavity_response(id,position,response,field)
        !! Radial Cauchy, first-Piola, vector-flux or cohesive-normal component.
        !! A symmetric tensor \(\boldsymbol{A}\) uses the Mandel order
        !! \([A_{11},A_{22},A_{33},\sqrt{2}A_{23},\sqrt{2}A_{13},\sqrt{2}A_{12}]^T\).
        integer, intent(in) :: id
        real(wp), intent(in) :: position(:,:),response(:,:)
        real(wp), intent(out) :: field(:)
        integer :: p,kind
        real(wp) :: n(3),tensor(3,3)

        if (size(position,1)/=3.or.size(position,2)/=size(response,2).or. &
            size(field)/=size(response,2)) error stop "Cavity response projection shape mismatch"
        kind=cases(id)%response_component
        do p=1,size(field)
            n=unit_normal(position(:,p))
            select case(kind)
            case(1)
                field(p)=n(1)**2*response(1,p)+n(2)**2*response(2,p)+n(3)**2*response(3,p)+ &
                    sqrt_two*(n(2)*n(3)*response(4,p)+n(1)*n(3)*response(5,p)+ &
                              n(1)*n(2)*response(6,p))
            case(4)
                tensor=reshape(response(:,p),[3,3])
                field(p)=dot_product(n,matmul(tensor,n))
            case(5)
                field(p)=dot_product(n,response(1:3,p))
            case(6)
                field(p)=response(1,p)
            case default
                error stop "Unsupported three-dimensional response representation"
            end select
        end do
    end subroutine project_cavity_response

end module benchmark_cavity3d_fields
