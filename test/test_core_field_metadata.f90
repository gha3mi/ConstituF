! SPDX-License-Identifier: Apache-2.0

program test_core_field_metadata
    !! Exercise physical shape, Mandel packing, exact metadata matching and
    !! rejection of ambiguous or contradictory author specifications.
    !! Illustrative dimensional fields use SI units, including Pa, K, m and J.
    !! Their units are a caller convention and are absent from the schema.
    use constituf_author, only: material_spec_t
    implicit none (type, external)

    type(material_spec_t) :: spec
    character(:), allocatable :: convention, measure, frame
    integer, allocatable :: shape(:)
    logical :: found
    integer :: empty_shape(0)

    call spec%add_input("deformation.f", 9, "tensor_full", &
                        "deformation_gradient", "spatial_reference", [3,3])
    call spec%add_output("stress.cauchy", 6, "mandel", &
                         "cauchy_stress", "spatial", [3,3])
    call spec%add_state("damage", 1, "scalar", "damage", "material")
    if (.not. spec%valid()) error stop "Valid tensor schema rejected"
    if (spec%input_size() /= 9 .or. spec%output_size() /= 6) then
        error stop "Unexpected packed field sizes"
    end if
    if (.not. spec%matches_field("output", "stress.cauchy", "mandel", &
                                 "cauchy_stress", "spatial", [3,3])) then
        error stop "Compatible Mandel field rejected"
    end if
    if (spec%matches_field("output", "stress.missing", "mandel", &
                           "cauchy_stress", "spatial", [3,3])) then
        error stop "Unknown field name accepted"
    end if
    if (spec%matches_field("output", "stress.cauchy", "mandel", &
                           "piola_stress", "spatial", [3,3])) then
        error stop "Stress measure mismatch accepted"
    end if
    if (spec%matches_field("output", "stress.cauchy", "mandel", &
                           "cauchy_stress", "reference", [3,3])) then
        error stop "Reference frame mismatch accepted"
    end if
    if (spec%matches_field("output", "stress.cauchy", "tensor_full", &
                           "cauchy_stress", "spatial", [3,3])) then
        error stop "Storage convention mismatch accepted"
    end if
    if (spec%matches_field("output", "stress.cauchy", "mandel", &
                           "cauchy_stress", "spatial", [2,2])) then
        error stop "Physical shape mismatch accepted"
    end if
    call spec%field_info("output", "stress.cauchy", convention, measure, &
                         frame, shape, found)
    if (.not. found) error stop "Field lookup failed"
    if (convention /= "mandel" .or. &
        measure /= "cauchy_stress" .or. frame /= "spatial") then
        error stop "Field metadata corrupted"
    end if
    if (size(shape) /= 2) error stop "Incorrect physical tensor rank"
    if (any(shape /= [3,3])) error stop "Incorrect tensor shape"
    call spec%field_info("output", "missing", convention, measure, frame, shape, found)
    if (found .or. allocated(shape)) error stop "Missing field reported as present"
    if (allocated(convention) .or. allocated(measure) .or. allocated(frame)) then
        error stop "Missing field returned metadata"
    end if
    if (spec%matches_field("unexpected", "damage", "scalar", &
                           "damage", "material", empty_shape)) then
        error stop "Unknown field role accepted"
    end if

    ! Representation validation is scientific metadata, not a component count guess.
    call invalid_mandel()
    call invalid_vector()
    call invalid_tensor()
    call invalid_scalar()
    call invalid_custom()
    call valid_array()

contains

    subroutine invalid_mandel()
        type(material_spec_t) :: bad
        call bad%add_input("strain", 9, "mandel", &
                           "small_strain", "material", [3,3])
        call bad%add_output("stress", 6, "mandel", &
                            "cauchy_stress", "spatial", [3,3])
        if (bad%valid()) error stop "Mandel storage mismatch accepted"
    end subroutine invalid_mandel

    subroutine invalid_vector()
        type(material_spec_t) :: bad
        call bad%add_input("gradient", 3, "vector", &
                           "temperature_gradient", "reference", [2])
        call bad%add_output("flux", 3, "vector", &
                            "heat_flux", "reference", [3])
        if (bad%valid()) error stop "Vector shape mismatch accepted"
    end subroutine invalid_vector

    subroutine invalid_tensor()
        type(material_spec_t) :: bad
        call bad%add_input("f", 8, "tensor_full", &
                           "deformation_gradient", "spatial_reference", [3,3])
        call bad%add_output("p", 9, "tensor_full", &
                            "piola_stress", "reference", [3,3])
        if (bad%valid()) error stop "Full tensor component mismatch accepted"
    end subroutine invalid_tensor

    subroutine invalid_scalar()
        type(material_spec_t) :: bad
        call bad%add_input("temperature", 2, "scalar", &
                           "temperature", "none")
        call bad%add_output("energy", 1, "scalar", &
                            "energy", "none")
        if (bad%valid()) error stop "Invalid scalar accepted"
    end subroutine invalid_scalar

    subroutine invalid_custom()
        type(material_spec_t) :: bad
        call bad%add_input("tensor", 6, "unknown", &
                           "strain", "material", [3,3])
        call bad%add_output("stress", 1, "scalar", &
                            "axial_stress", "material")
        if (bad%valid()) error stop "Unqualified custom convention accepted"
    end subroutine invalid_custom

    subroutine valid_array()
        type(material_spec_t) :: good
        call good%add_input("slips", 12, "array", &
                            "slip_state", "crystal", [4,3])
        call good%add_output("resistance", 12, "array", &
                             "slip_resistance", "crystal", [4,3])
        if (.not. good%valid()) error stop "Generic array incorrectly rejected"
        if (.not. good%matches_field("input", "slips", "array", &
                                     "slip_state", "crystal", [4,3])) then
            error stop "Generic array compatibility mismatch"
        end if
    end subroutine valid_array

end program test_core_field_metadata
