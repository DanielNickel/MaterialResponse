using Tensors, ForwardDiff, Waveforms, Plots
using MaterialModelsBase
using Newton
import CairoMakie as CM

# Include the supplied files
include(joinpath(@__DIR__, "utils.jl"))
include(joinpath(@__DIR__, "perfect_plasticity.jl"))
include(joinpath(@__DIR__, "chaboche.jl"))
include(joinpath(@__DIR__, "uniaxial_stress.jl"))
include(joinpath(@__DIR__, "chabocheNorton.jl"))

amplitudes_in_percent=[0.4, 0.8]
amplitude_monotonic_percent=1.0
load_rates=[1.0, 2.0 ]
load_rate_monotonic=1.0
Triangular_Wave=true
Plot_over_Strain=true #if true-> plot over strain, else: plot over time
Which_Model=3 #1 for perfect_plasticity, 2 for chaboche, 3 for chabocheNorton
#because perfect_plasticity and chaboche are time independent, 
# plotting over time is only permitted for chabocheNorton

if Which_Model==1 || Which_Model==2
    Plot_over_Strain=true
end

# Helper function to get a triangle wave strain time history
function get_cycles(;amplitude=1.0, num_cycles=1, num_steps_per_cycle=100)
    return amplitude*trianglewave1.(range(0.0, num_cycles, num_steps_per_cycle*num_cycles))
end

function createstrainramp(m::AbstractMaterial)

    num_steps_rise = 100
    rise_values = range(0.0, stop = twice_yield_limit_in_strain(m), length = num_steps_rise)#
    constant_value = twice_yield_limit_in_strain(m)
    num_steps_constant = 100
    constant_values = repeat([constant_value], inner = num_steps_constant)

    ϵ11 = (vcat(rise_values, constant_values))

    fig2 = CM.Figure()
    ax2 = CM.Axis(fig2[1, 1])  # Fix here: Use fig2 instead of fig

    CM.lines!(ax2, ϵ11, markersize = 5, color = :blue)
    CM.save("strain_ramp.pdf",fig2)
    
    return ϵ11, twice_yield_limit_in_strain(m)
end

# Define the material
if Which_Model==1
    m_perfect = PerfectPlasticity(G=80.e3, K=160.e3, Y=250.0)
    title_model="Perfect Plasticity"
elseif Which_Model==2
    m_chaboche = Chaboche(G=80.e3, K=160.e3, Y=250.0, Hiso=100.e3, κ∞=400.0, Hkin=40.e3, β∞=50.0)
    title_model="Chaboche"
else
    m_chabocheNorton=ChabocheNorton(G=80.e3, K=160.e3, Y=250.0, Hiso=100e3, κ∞=400.0, Hkin=40e3, β∞=50.0, t_star=1.0, n=1.0)
    title_model="Chaboche Norton"
end

# Define function to plot results
function plot_results(m, ϵ11,ΔL_monotonic, loading_rate; kwargs...)
    fig = CM.Figure();
    if Triangular_Wave
        if Plot_over_Strain
            ax = CM.Axis(fig[1,1]; xlabel = "ϵ₁₁ [%]", ylabel = "σ₁₁ [MPa]", title="$title_model: σ₁₁ in relation to ϵ₁₁, Triangular Wave")
            plot_results!(ax, m, ϵ11,ΔL_monotonic, loading_rate; kwargs...)
        else
            ax = CM.Axis(fig[1,1]; xlabel = "t [s]", ylabel = "σ₁₁ [MPa]", title="$title_model: σ₁₁ in relation to time, Triangular Wave")
            plot_results!(ax, m, ϵ11,ΔL_monotonic, loading_rate; kwargs...)
        end
    else
        if Plot_over_Strain
            ax = CM.Axis(fig[1,1]; xlabel = "ϵ₁₁ [%]", ylabel = "σ₁₁ [MPa]", title="$title_model: σ₁₁ in relation to ϵ₁₁, Strain ramp")
            plot_results!(ax, m, ϵ11,ΔL_monotonic, loading_rate; kwargs...)
        else
            ax = CM.Axis(fig[1,1]; xlabel = "t [s]", ylabel = "σ₁₁ [MPa]", title="$title_model: σ₁₁ in relation to time, Strain ramp")
            plot_results!(ax, m, ϵ11,ΔL_monotonic, loading_rate; kwargs...)
        end
    end
    
    return fig, ax
end
function plot_results!(ax, m::AbstractMaterial, ϵ11, displacement, loading_rate; label="") #time history
    if Triangular_Wave
        t_end=displacement/loading_rate
    else
        t_end=displacement/(loading_rate) #Reason: if there is a strain ramp, the loading phase is only half of the time
        #the other half of the time there is no loading rate anymore, because the constant value is reached after half of the time
    end
    t = collect(range(0,t_end,length(ϵ11)))
    σ11 = uniaxial_stress(m, ϵ11, t, displacement)
    if Plot_over_Strain
        CM.lines!(ax, 100*ϵ11, σ11; label=label)
    else
        CM.lines!(ax, t, σ11; label=label)
    end
    return nothing
end

# Example how to plot 1 cycle, for different strain amplitudes and monotonic loading
function create_example_plot(m)
    if Triangular_Wave
        fig, ax = plot_results(m,collect(range(0, amplitude_monotonic_percent/100, 200)),amplitude_monotonic_percent/100,load_rate_monotonic; label="monotonic, ϵ₁₁ = $amplitude_monotonic_percent%, ϵ̇  = $load_rate_monotonic%s⁻¹")#createstrainramp()
        for loading_rate in load_rates
            for amplitude_percent in amplitudes_in_percent
                ϵ11 = get_cycles(;amplitude=amplitude_percent/100, num_steps_per_cycle=Int(round(500*amplitude_percent)))
                ΔL=4*(amplitude_percent/100)
                plot_results!(ax, m, ϵ11,ΔL,loading_rate; label="ϵ₁₁ = $amplitude_percent%, ϵ̇  = $loading_rate%s⁻¹")
            end
        end
    else
        ϵ11, ΔL=createstrainramp(m)
        ΔL_rounded=round(ΔL*100, digits=3)
        fig, ax = plot_results(m,collect(range(0,ΔL,200)),ΔL,load_rate_monotonic; label="monotonic, ϵ11=$ΔL_rounded%, ϵ̇ = $load_rate_monotonic%s⁻¹")
        for loading_rate in load_rates
            plot_results!(ax, m, ϵ11,2*ΔL,loading_rate; label="ϵ₁₁ = $ΔL_rounded%, ϵ̇  = $loading_rate%s⁻¹")
        end
    end
    
        CM.axislegend(ax; position=:rb)
    CM.save("Chaboche_Norton_triangular_over_time_2.pdf",fig)
    return fig, ax
end

# # for Task A4
# function create_example_plot(m)
#     # fig, ax = plot_results(m, collect(range(0, 0.0/100, 200)))# label="monotonic")
#     amplitude_percent=0.15
#     ϵ11 = get_cycles(;amplitude=amplitude_percent/100, num_steps_per_cycle=Int(round(500*amplitude_percent)))
#     fig, ax=plot_results(m, ϵ11; label="ϵ₁₁ = ±$amplitude_percent %")
#     amplitude_percent=0.9
#     ϵ11 = get_cycles(;amplitude=amplitude_percent/100, num_steps_per_cycle=Int(round(500*amplitude_percent)))
#     plot_results!(ax, m, ϵ11; label="ϵ₁₁ = ±$amplitude_percent %")
#     CM.axislegend(ax; position=:rb)
#     CM.save("CA2_TaskA5.pdf",fig)
#     return fig, ax
# end

# # for Task A5
# function create_example_plot(m)
#     # fig, ax = plot_results(m, collect(range(0, 0.0/100, 200)))# label="monotonic")
#     amplitude_percent=0.5
#     ϵ11 = get_cycles10(;amplitude=amplitude_percent/100, num_steps_per_cycle=10)
#     fig, ax=plot_results(m, ϵ11; label="10 time steps")

#     ϵ11 = get_cycles(;amplitude=amplitude_percent/100, num_steps_per_cycle=100)
#     plot_results!(ax, m, ϵ11; label="100 time steps")

#     ϵ11 = get_cycles1000(;amplitude=amplitude_percent/100, num_steps_per_cycle=1000)
#     plot_results!(ax, m, ϵ11; label="1000 time steps")

#     CM.axislegend(ax; position=:rb)
#     CM.save("CA2_TaskA5.pdf",fig)
#     return fig, ax
# end


fig, ax = create_example_plot(m_chabocheNorton)
fig

if Which_Model==1
    fig, ax = create_example_plot(m_perfect)
    fig
elseif Which_Model==2
    fig, ax = create_example_plot(m_chaboche)
    fig
else
    fig, ax = create_example_plot(m_chabocheNorton)
    fig
end
# specify strain rate as input