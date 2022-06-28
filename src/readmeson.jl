# 6-19
# kangjiayin
workdir = "/Users/kjy/Desktop/program/julia/DSE_BSE"
cd(workdir)
using TOML
using LinearAlgebra
using Dierckx
using JLD2
using Gaussquad # View in github: https://github.com/kangjiayin/Gaussquad.jl
using FastGaussQuadrature
using ChebyshevFun 

dataset = TOML.parsefile("src/config.toml")
# quark system
quarkrepoint = dataset["readsetting"]["repoint"]
quarkintstep = dataset["readsetting"]["quarkintstep"]
quarkm = dataset["readsetting"]["quarkmass"]

# data
logofcutoff = dataset["readsetting"]["logofcutoff"]
cutup = 10. ^logofcutoff
cutdown = 10. ^(-logofcutoff)

kstep = dataset["readsetting"]["kstep"]
zstep = dataset["readsetting"]["zstep"]
Pstep = dataset["readsetting"]["Pstep"]
# 取点方式
meshk,weightk= gausslegendremesh(cutdown,cutup,kstep,2);
meshz,weightz= gausschebyshev(zstep,2);

# using SPMinterpolation
Ppoints=Pstep
# kstep=128
# k=gausslegendremesh(10. ^-4,10. ^4,32,2);
P2=Array{Float64}(undef,Ppoints,1)
F1k=Array{Float64}(undef,Ppoints,kstep,zstep)
F2k=Array{Float64}(undef,Ppoints,kstep,zstep)
F3k=Array{Float64}(undef,Ppoints,kstep,zstep)
F4k=Array{Float64}(undef,Ppoints,kstep,zstep)
mP2=Array{Float64}(undef,Ppoints,1)
mF1k=Array{Float64}(undef,Ppoints,kstep,zstep)
mF2k=Array{Float64}(undef,Ppoints,kstep,zstep)
mF3k=Array{Float64}(undef,Ppoints,kstep,zstep)
mF4k=Array{Float64}(undef,Ppoints,kstep,zstep)
# F1k=Array{Float64}(undef,Ppoints,kstep)
# F2k=Array{Float64}(undef,Ppoints,kstep)
# F3k=Array{Float64}(undef,Ppoints,kstep)
# F4k=Array{Float64}(undef,Ppoints,kstep)
kname=Array{String}(undef,Ppoints,1)
if dataset["readsetting"]["readmode"] == 1
    if ispath("data/pseudo_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint") == true
        print("导入--","data/pseudo_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint\n")
        for i=1:Ppoints
            # local a, b
            global P2, F1k, F2k, F3k, F4k
            P2[i], F1k[i,:,:], F2k[i,:,:], F3k[i,:,:], F4k[i,:,:] =load("data/pseudo_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint/P&F1-4_$i-$Pstep.jld2","P2", "F1", "F2", "F3", "F4")
        end #for i
    else #ispath
        print("不存在该文件","data/pseudo_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint","\n")
    end
elseif dataset["readsetting"]["readmode"] == 2
    if ispath("data/scalar_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint") == true
        print("导入","data/scalar_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint\n")
        for i=1:Ppoints
            # local a, b
            global P2, F1k, F2k, F3k, F4k
            P2[i], F1k[i,:,:], F2k[i,:,:], F3k[i,:,:], F4k[i,:,:] =load("data/scalar_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint/P&F1-4_$i-$Pstep.jld2","P2", "F1", "F2", "F3", "F4")
        end #for i
    else #ispath
        print("不存在该文件,","data/scalar_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint","\n")
    end
end # if

F1pq=zeros(Pstep,kstep);
F2pq=zeros(Pstep,kstep);
F3pq=zeros(Pstep,kstep);
F4pq=zeros(Pstep,kstep);
for i=1:Pstep
    for j=1:kstep
        # i for P2[i], j for q2[j]
        for m=1:zstep
            F1pq[i,j] += weightz[m]*F1k[i,j,m]*2/pi
            F2pq[i,j] += weightz[m]*F2k[i,j,m]*2/pi
            F3pq[i,j] += weightz[m]*F3k[i,j,m]*2/pi
            F4pq[i,j] += weightz[m]*F4k[i,j,m]*2/pi
            # f_gammaPq[i,j] += weightz[m]*(P2[i]*A2[i,j,m]*((F1k[i,j,m] + 4*F4k[i,j,m]*q2[j])*A1[i,j,m] - 2*F2k[i,j,m]*B1[i,j,m]) - 2*F2k[i,j,m]*P2[i]*A1[i,j,m]*B2[i,j,m] + 4*sqrt(P2[i]*q2[j])*(F2k[i,j,m] + F3k[i,j,m]*q2[j])*meshz[m]*(A2[i,j,m]*B1[i,j,m] - A1[i,j,m]*B2[i,j,m]) - 2*P2[i]*q2[j]*meshz[m]^2*(2*F4k[i,j,m]*A1[i,j,m]*A2[i,j,m] + F3k[i,j,m]*A2[i,j,m]*B1[i,j,m] + F3k[i,j,m]*A1[i,j,m]*B2[i,j,m]) - 4*F1k[i,j,m]*(q2[j]*A1[i,j,m]*A2[i,j,m] + B1[i,j,m]*B2[i,j,m]))
        end
    end
end

function Inport()
    global z2, z4
    local A, B, k
    A, B, k, z2, z4=load("data/quark_gap_equation/quark-ABkz2z4-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint.jld2","A","B","k", "z2", "z4");
    print("确认导入的数据为--", "quark-ABkz2z4-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint.jld2\n")
    global AA
    global BB
    AA=Spline1D(k,A)
    BB=Spline1D(k,B)
    return true
end
Inport()

include("decayconst.jl")