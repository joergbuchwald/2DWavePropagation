

# basic definitions
@enum Sourcetype line=1 plane=2

@enum Boundarytype rigid=1 free=2 absorbing=3

struct Displacementfield
    v::Array{Float64, 3}
    zeta::Array{Float64, 2}
end

struct Boundarydisplacements
    vt::Array{Float64, 2}
    vb::Array{Float64, 2}
    vl::Array{Float64, 2}
    vr::Array{Float64, 2}
end

struct Boundaries
    top::Boundarytype
    bottom::Boundarytype
    left::Boundarytype
    right::Boundarytype
end

struct Params
    #LX,LZ: number of grid points in x and z direction (note: boundaries are at points with indices 2, lx-1,lz-1)
    lx::Int32
    lz::Int32
    #vp, eps, kappa: background velocity, std.dev. of random fluctuations, parameter for density fluctuations (3F10.0)
    vp::Float64 #//[vp]=km/s
    eps::Float64 # 0.05
    kappa::Float64 # 0.0;
    #source signal (KUEPPERS,1958):
    nqx::Int32 # NQX: x index of source, last config: 2000-pt lattice 1000
    nqz::Int32 # NQZ: z index of source, last config: 2000-pt lattice 100
    iquell::Sourcetype #IQUELL:=1: line source,  =2: plane wave starting at z=0  (F10.0,4I10)
    dx::Float64 # 0.12;//[dx]=km
    nsig::Int32 #2; //NSIG: number of exrema
    nmax::Int32 #16000; ?
end

# Einganspuls (Kuepperwavelet)
function signal(t::Float64, delta::Float64=6.0, nsig::Int=2)
    ampmax = zeros(10)
    ampmax[1]=1.333
    ampmax[2]=1.299
    ampmax[3]=1.6
    ampmax[4]=1.585
    ampmax[5]=1.714
    ampmax[6]=1.706
    ampmax[7]=1.778
    ampmax[8]=1.773
    ampmax[9]=1.818
    ampmax[10]=1.815
    fm = (1.0 * nsig + 2.0) / (1.0 * nsig)
    del = nsig * pi / delta
    signal = (t < 0 || t > delta) ? 0.0 : (sin(del * t)-((sin(fm* del *t))/fm)) / ampmax[nsig]
end

function timeloop(p::Params, bd::Boundarydisplacements, d::Displacementfield)
    for nt in 1:p.nmax
        saveboundarydisplacements(p, bd, d)
    end
end

function setinitialcondtions(field::Displacementfield, p::Params)
    
end

function saveboundarydisplacements!(p::Params, bd::Boundarydisplacements, d::Displacementfield)
    for i in 2:p.lx-1
        bd.vt[i][1] = d.v[i][2][1]
        bd.vt[i][2] = d.v[i][3][1]
        bd.vb[i][1] = d.v[i][p.lz-1][1]
        bd.vb[i][2] = d.v[i][p.lz-2][1]
    end
    for j in 2:p.lz-1
        bd.vl[j][1] = d.v[2][j][1]
        bd.vl[j][2] = d.v[3][j][1]
        bd.vr[j][1] = d.v[p.lx-1][j][1]
        bd.vr[j][2] = d.v[p.lx-2][j][1]
    end
end

function calculatedisplacement!(p::Params, d::Displacementfield)
    for i in 2:p.lx-1
        for j in 2:p.lz-1
            r5 = 1. / (1. + p.kappa * fij)
            f1 = 1. + r5 * (1. + p.kappa * d.zeta[i + 1][j])
            f2 = 1. + r5 * (1. + p.kappa * d.zeta[i - 1][j])
            f3 = 1. + r5 * (1. + p.kappa * d.zeta[i][j+1])
            f4 = 1. + r5 * (1. + p.kappa * d.zeta[i][j-1])
            v[i][j][1] = -d.v[i][j][1] + fact1 * (1. + d.zeta[i][j]) * (1. + d.zeta[i][j]) *
                ((d.v[i+1][j][2] - d.v[i][j][2]) / f1 - (d.v[i][j][2] - d.v[i-1][j][2]) / f2 +
                 (d.v[i][j+1][2] - d.v[i][j][2]) / f3 - (d.v[i][j][2] - d.v[i][j-1][2]) / f4) + 2. * d.v[i][j][2]
        end
    end
end

function applyboundaryconditions!()
end
