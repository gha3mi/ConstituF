! SPDX-License-Identifier: Apache-2.0

module benchmark_pages
    !! Generate each study's material-results gallery from the shared registry.
    use benchmark_models, only: cases, names, ncases, load_mandel, load_plane, load_axisym, load_F, &
                                load_cohesive_i, load_cohesive_mixed, load_thermoelastic, load_heat, &
                                load_heat_temperature, load_concentration, load_electric, load_piezo, &
                                load_pressure, load_pressure_pressure, load_biot
    implicit none (type, external)
    private

    public :: write_benchmark_index

contains

    subroutine write_benchmark_index(study,title,included)
        !! An optional registry mask limits the gallery to compatible materials.
        !! Call only after a complete, unfiltered sensitivity run succeeds.
        character(*), intent(in) :: study,title
        logical, intent(in), optional :: included(:)
        integer :: u,id,family,i,letter
        character(:), allocatable :: slug,label
        character(44), parameter :: family_titles(4) = [character(len=44) :: &
            "Solid and interface mechanics", "Heat and species transport", &
            "Electrical and electromechanical response", "Porous flow and poromechanics"]

        if (present(included)) then
            if (size(included)/=ncases) error stop "Invalid benchmark gallery material mask"
        end if

        open(newunit=u,file="docs/pages/benchmark-"//study//"-results.md",status="replace",action="write")
        write(u,"(a)") "---"
        write(u,"(a)") "title: "//title//" material results"
        write(u,"(a)") "---"
        write(u,"(a)") ""
        write(u,"(a)") "Each compatible material has a constitutive response preview and"
        write(u,"(a)") "an analytical parameter-sensitivity heatmap. Select a material to open"
        write(u,"(a)") "its response, derivative metrics and numerical exports."
        write(u,"(a)") "See the [benchmark formulation](|page|/benchmark-"//study//".html)"
        write(u,"(a)") "and [all benchmarks](|page|/benchmarks.html)."
        write(u,"(a)") ""
        write(u,"(a)") "## Material results"
        write(u,"(a)") ""

        do family=1,size(family_titles)
            write(u,"(a)") "### "//trim(family_titles(family))
            write(u,"(a)") ""
            write(u,"(a)") '<div class="cf-benchmark-atlas">'
            material_index: do id=1,ncases
                if (present(included)) then
                    if (.not.included(id)) cycle material_index
                end if
                if (material_family(cases(id)%input_kind)/=family) cycle material_index
                slug=trim(names(id))
                label=slug
                do i=1,len(label)
                    if (label(i:i)=="_") label(i:i)=" "
                end do
                letter=iachar(label(1:1))
                if (letter>=iachar("a").and.letter<=iachar("z")) label(1:1)=achar(letter-32)
                write(u,"(a)") '  <a class="cf-benchmark-tile" href="benchmark-'//study//"-"//slug//'.html">'
                write(u,"(a)") '    <img src="../media/benchmarks/'//study//"/"//slug// &
                    '_response.svg" loading="lazy" alt="Constitutive response for '//label//'">'
                write(u,"(a)") "    <span>"//label//"</span>"
                write(u,"(a)") "  </a>"
            end do material_index
            write(u,"(a)") "</div>"
            write(u,"(a)") ""
        end do
        close(u)
    end subroutine write_benchmark_index

    integer function material_family(kind) result(family)
        integer, intent(in) :: kind

        select case (kind)
        case (load_mandel, load_plane, load_axisym, load_F, load_cohesive_i, load_cohesive_mixed, load_thermoelastic)
            family=1
        case (load_heat, load_heat_temperature, load_concentration)
            family=2
        case (load_electric, load_piezo)
            family=3
        case (load_pressure, load_pressure_pressure, load_biot)
            family=4
        case default
            error stop "Unsupported benchmark input kind"
        end select
    end function material_family

end module benchmark_pages
