! SPDX-License-Identifier: Apache-2.0

module benchmark_metrics
    !! Weighted analytical constitutive parameter sensitivity metrics.
    !! Geometry-specific callers supply positive cell areas or volumes.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf, only: wp
    implicit none (type, external)
    private

    type, public :: sensitivity_metric_t
        !! A model/parameter pair, with exact derivative and reference value.
        logical :: present=.false., relative_defined=.false., below_resolution=.false.
        character(48) :: parameter=""
        character(24) :: units=""
        character(12) :: scale_source=""
        real(wp) :: reference_value=0.0_wp, parameter_scale=0.0_wp
        real(wp) :: response_rms=0.0_wp, sensitivity_rms=0.0_wp
        real(wp) :: sensitivity_linf=0.0_wp, signed_mean=0.0_wp
        real(wp) :: relative_rms=0.0_wp, sign_coherence=0.0_wp
    end type sensitivity_metric_t

    public :: calculate_metric, parameter_scale, verify_metric_contract

contains

    pure real(wp) function weighted_rms(a,weights) result(value)
        !! Stable RMS: scale before squaring to avoid intermediate overflow
        !! and subnormal underflow in real32/real128 source compatibility.
        real(wp), intent(in) :: a(:),weights(:)
        real(wp) :: scale

        value=0.0_wp
        scale=maxval(abs(a))
        if (scale<=0.0_wp) return
        value=scale*sqrt(sum(weights*(a/scale)**2)/sum(weights))
    end function weighted_rms

    pure real(wp) function parameter_scale(index,values,units) result(scale)
        !! A zero or negative parameter has no logarithmic derivative.
        !! Use \(|p|\) for nonzero \(p\), or the largest \(|p|\) with the same SI unit
        !! for a zero coefficient (e.g., unused anisotropic off-diagonal).
        !! The chosen scale is explicitly recorded in the CSV.
        integer, intent(in) :: index
        real(wp), intent(in) :: values(:)
        character(*), intent(in) :: units(:)
        integer :: j

        scale=abs(values(index))
        if (scale>0.0_wp) return
        do j=1,size(values)
            if (trim(units(j))==trim(units(index))) scale=max(scale,abs(values(j)))
        end do
    end function parameter_scale

    subroutine calculate_metric(metric,field,derivative,weights,scale)
        !! Response and scaled derivative share physical units, so their
        !! weighted-RMS ratio is dimensionless across different laws.
        type(sensitivity_metric_t), intent(inout) :: metric
        real(wp), intent(in) :: field(:),derivative(:),weights(:)
        real(wp), intent(in) :: scale
        real(wp) :: max_abs,weighted_average
        real(wp) :: scaled(size(field))

        if (size(field)/=size(derivative).or.size(field)/=size(weights)) then
            error stop "Sensitivity metric field/weight shape mismatch"
        end if
        if (any(weights<=0.0_wp).or..not.all(ieee_is_finite(weights))) error stop "Invalid integration weights"
        if (.not.all(ieee_is_finite(field)).or..not.all(ieee_is_finite(derivative))) error stop "Nonfinite field"
        if (.not.ieee_is_finite(scale).or.scale<=0.0_wp) error stop "Nonpositive parameter reference scale"
        metric%present=.true.
        metric%parameter_scale=scale
        scaled=scale*derivative
        if (.not.all(ieee_is_finite(scaled))) error stop "Nonfinite scaled derivative"
        metric%response_rms=weighted_rms(field,weights)
        metric%sensitivity_rms=weighted_rms(scaled,weights)
        max_abs=maxval(abs(scaled))
        metric%sensitivity_linf=max_abs
        weighted_average=sum(weights*(scaled/max(max_abs,tiny(1.0_wp))))/sum(weights)*max_abs
        metric%signed_mean=weighted_average
        if (metric%sensitivity_rms>0.0_wp) then
            metric%sign_coherence=max(-1.0_wp,min(1.0_wp,metric%signed_mean/metric%sensitivity_rms))
        else
            metric%sign_coherence=0.0_wp
        end if
        metric%relative_defined=metric%response_rms>0.0_wp
        if (metric%relative_defined) metric%relative_rms=metric%sensitivity_rms/metric%response_rms
        metric%below_resolution=metric%relative_defined.and.metric%relative_rms>0.0_wp.and. &
            metric%relative_rms<=256.0_wp*epsilon(1.0_wp)
        if (.not.all(ieee_is_finite([metric%response_rms,metric%sensitivity_rms,metric%signed_mean, &
                                    metric%relative_rms,metric%sign_coherence]))) error stop "Nonfinite metric"
    end subroutine calculate_metric

    subroutine verify_metric_contract()
        !! Independent synthetic weighted identities; no parameter
        !! perturbations or finite-difference derivative calculations.
        type(sensitivity_metric_t) :: metric
        integer, parameter :: nsamples=2,nparameter=3
        real(wp), parameter :: samples(nsamples)=[2.0_wp,4.0_wp]
        real(wp), parameter :: derivatives(nsamples)=[1.0_wp,-1.0_wp]
        real(wp), parameter :: weights(nsamples)=[1.0_wp,3.0_wp]
        real(wp), parameter :: sample_parameters(nparameter)=[0.0_wp,4.0_wp,6.0_wp]
        character(2), parameter :: units(nparameter)=[character(len=2) :: "Pa","Pa","m"]
        real(wp) :: tolerance

        tolerance=512.0_wp*epsilon(1.0_wp)
        if (abs(parameter_scale(1,sample_parameters,units)-4.0_wp)>tolerance) then
            error stop "Same-unit zero-parameter scaling identity"
        end if
        if (abs(parameter_scale(3,sample_parameters,units)-6.0_wp)>tolerance) then
            error stop "Absolute-value parameter scaling identity"
        end if
        call calculate_metric(metric,samples,derivatives,weights,2.0_wp)
        if (abs(metric%response_rms-sqrt(13.0_wp))>tolerance) error stop "Weighted response RMS identity"
        if (abs(metric%sensitivity_rms-2.0_wp)>tolerance) error stop "Weighted sensitivity RMS identity"
        if (abs(metric%relative_rms-2.0_wp/sqrt(13.0_wp))>tolerance) then
            error stop "Dimensionless weighted sensitivity ratio identity"
        end if
        if (abs(metric%signed_mean+1.0_wp)>tolerance) error stop "Weighted signed mean identity"
        if (abs(metric%sign_coherence+0.5_wp)>tolerance) error stop "Signed coherence identity"

        call calculate_metric(metric,samples,0.0_wp*derivatives,weights,2.0_wp)
        if (.not.metric%relative_defined.or.metric%relative_rms>0.0_wp) error stop "Exact zero sensitivity identity"
        call calculate_metric(metric,0.0_wp*samples,derivatives,weights,2.0_wp)
        if (metric%relative_defined) error stop "Undefined normalization identity"
    end subroutine verify_metric_contract

end module benchmark_metrics
