"""
The barrel vault (Scordelis-Lo) roof is one of the benchmarks for linear elastic
analysis of shells. 

The candidate element's usefulness in irregular geometries (and most practical
cases involve a high degree of geometric irregularity) is tested. As would be
expected,the irregular mesh results are not as good as those provided by a
regular meshwith the same number of variables. 

Problem description

The physical basis of the problem is a deeply arched roof supported only
bydiaphragms at its curved edges (an aircraft hanger), deforming under its own
weight. It is interesting to observe that the geometry is such that the
centerpoint of the roof moves upward under the self-weight(downwardly directed)
load. Perhaps this is one reason why the problem is not straightforward
numerically. 
"""
module scordelis_lo_benchmark

using LinearAlgebra
using FinEtools
using FinEtools.FTypesModule: FInt, FFlt, FFltMat, FFltVec
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellQ4RSModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

# analytical solution for the vertical deflection and the midpoint of the
# free edge 
const analyt_sol = -0.30202
# Parameters:
const E = 4.32e8
const nu = 0.0
const thickness = 0.25 # geometrical dimensions are in feet
const R = 25.0
const L = 50.0

cylindrical!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, feid::FInt, qpid::FInt) = begin
    r = vec(XYZ)
    r[2] = 0.0
    r[3] += R
    csmatout[:, 3] .= vec(r) / norm(vec(r))
    csmatout[:, 2] .= (0.0, 1.0, 0.0) #  this is along the axis
    cross3!(view(csmatout, :, 1), view(csmatout, :, 2), view(csmatout, :, 3))
    return csmatout
end

function _execute(formul, n = 8, visualize = true)
    tolerance = R / n / 10
    # fens, fes = Q4blockrand(40/360*2*pi,L/2,n,n);
    fens, fes = Q4block(40 / 360 * 2 * pi, L / 2, n, n)
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a = fens.xyz[i, 1]
        y = fens.xyz[i, 2]
        fens.xyz[i, :] .= (R * sin(a), y, R * (cos(a) - 1))
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    ocsys = CSys(3, 3, cylindrical!)

    sfes = FESetShellQ4()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, GaussRule(2, 2), thickness), mater)
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))

    # Apply EBC's
    # rigid diaphragm
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    for i in [1, 3, 5]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Y
    l1 = selectnode(fens; box = Float64[-Inf Inf L / 2 L / 2 -Inf Inf], inflate = tolerance)
    for i in [2, 4, 6]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to X
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1, 5, 6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi)

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi)

    # Midpoint of the free edge
    nl = selectnode(fens; box = Float64[sin(40 / 360 * 2 * pi) * 25 sin(40 / 360 * 2 * pi) * 25 L / 2 L / 2 -Inf Inf], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, GaussRule(2, 2)))
    fi = ForceIntensity(FFlt[0, 0, -90, 0, 0, 0])
    F = distribloads(lfemm, geom0, dchi, fi, 2)

    # Solve
    solve_blocked!(dchi, K, F)
    U = gathersysvec(dchi, DOF_KIND_ALL)
    result = dchi.values[nl, 3][1]
    @info "n=$(n): $(round(result, digits = 6)), error $(round(100 - result/analyt_sol*100, digits = 4))%"

    # Visualization
    if visualize

        # Generate a graphical display of resultants
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
            push!(scalars, ("m$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
            push!(scalars, ("em$nc", fld.values))
        end
        vtkwrite("scordelis_lo_benchmark-$(n)-m.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
            push!(scalars, ("n$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
            push!(scalars, ("en$nc", fld.values))
        end
        vtkwrite("scordelis_lo_benchmark-$(n)-n.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
        scalars = []
        for nc in 1:2
            fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            push!(scalars, ("q$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            push!(scalars, ("eq$nc", fld.values))
        end
        vtkwrite("scordelis_lo_benchmark-$(n)-q.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])

        vtkwrite("scordelis_lo_benchmark-$(n)-uur.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3]), ("ur", dchi.values[:, 4:6])])

        scattersysvec!(dchi, (L / 8) / maximum(abs.(U)) .* U, DOF_KIND_ALL)
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[0 0 -L / 2]; [L / 2 L / 2 L / 2]]),
            #plot_nodes(fens),
            plot_midsurface(fens, fes; x = geom0.values, facecolor = "rgb(12, 12, 123)"),
            plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
            dims = 1)
        pl = render(plots)
    end

    result
end

function test_convergence(ns = [4, 6, 8, 16, 32], visualize = false)
    formul = FEMMShellQ4RSModule
    @info "Scordelis-Lo shell benchmark"
    results = []
    for n in ns
        v = _execute(formul, n, visualize)
        push!(results, v / analyt_sol * 100)
    end
    return ns, results
end

end # module

using .scordelis_lo_benchmark

# Visualized internal resultants
ns, results = scordelis_lo_benchmark.test_convergence([2 * 16,], false)


# These results come from Table 9 of An efficient three‑node triangular
# Mindlin–Reissner flat shell element, Hosein Sangtarash1 · Hamed Ghohani
# Arab1 · Mohammad R. Sohrabi1 · Mohammad R. Ghasemi1
# they are all normalized relative to 0.3024# .   # 
# Mesh Subdivision
# 4, 8, 16, 32

Allman = [
    1.004
    0.987
    0.987
    0.988] .* 100
Cook = [
    0.907
    0.929
    0.950
    0.981] .* 100
Providas_Kattis = [
    0.734
    0.815
    0.873
    0.967] .* 100
Shin_Lee = [
    1.379
    1.023
    1.004
    missing] .* 100
MITC3plus = [
    0.669
    missing
    0.857
    0.955] .* 100
TMRFS = [
    0.924
    0.963
    0.974
    0.998
] .* 100

using PGFPlotsX

objects = []

ns, results = scordelis_lo_benchmark.test_convergence()

compensate(r) = r
all_results = [("Present", compensate(results), "*"), ("Allman", compensate(Allman), "x"), ("ProvKat", compensate(Providas_Kattis), "triangle"), ("Cook", compensate(Cook), "square"), ("SL", compensate(Shin_Lee), "o"), ("MITC3+", compensate(MITC3plus), "diamond")]

for r in all_results
    @pgf p = PGFPlotsX.Plot(
        {
            color = "black",
            line_width = 0.7,
            style = "solid",
            mark = "$(r[3])"
        },
        Coordinates([v for v in zip(ns, r[2]) if v[2] !== missing])
    )
    push!(objects, p)
    push!(objects, LegendEntry("$(r[1])"))
end


@pgf ax = Axis(
    {
        xlabel = "Number of Elements / side [ND]",
        ylabel = "Normalized Displacement [ND]",
        # xmin = range[1],
        # xmax = range[2],
        xmode = "linear",
        ymode = "linear",
        yminorgrids = "true",
        grid = "both",
        legend_pos = "north east"
    },
    objects...
)

display(ax)
pgfsave("scordelis_lo_benchmark-convergence.pdf", ax)

# ns, results = scordelis_lo_benchmark.test_convergence([4, 8, 16])
# q1, q2, q3 = results
# @show qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)

# ns, results = scordelis_lo_benchmark.test_convergence(2 .* [4, 8, 16])
# q1, q2, q3 = results
# @show qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)

ns, results = scordelis_lo_benchmark.test_convergence(2 .* [16, 32, 64])
q1, q2, q3 = results
@show qtrue = (q2^2 - q1 * q3) / (2 * q2 - q1 - q3) * scordelis_lo_benchmark.analyt_sol / 100
