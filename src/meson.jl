# workdir = "/Users/kjy/Desktop/program/julia/DSE_BSE"
# # 放在一个模块儿里运行，与主题函数隔开
# cd(workdir)

module mesonbse
using TOML
using LinearAlgebra
using Dierckx
using JLD2
using Gaussquad # View in github: https://github.com/kangjiayin/Gaussquad.jl
using FastGaussQuadrature
using ChebyshevFun 

dataset = TOML.parsefile("src/config.toml")
# quark system
quarkrepoint = dataset["quarkDSE"]["repoint"]
quarkintstep = dataset["quarkDSE"]["quarkintstep"]
quarkm = dataset["quarkDSE"]["quarkmass"]

# data
logofcutoff = dataset["mesonBSE"]["logofcutoff"]
mt = dataset["data"]["mt"]
τ = dataset["data"]["tau"]
Λ = dataset["data"]["lambda"]
ω = dataset["data"]["omega"]
dd = dataset["data"]["dd"]
Nf = dataset["data"]["Nf"]
rm = dataset["data"]["rm"]
cutup = 10. ^logofcutoff
cutdown = 10. ^(-logofcutoff)

kstep = dataset["mesonBSE"]["kstep"]
zstep = dataset["mesonBSE"]["zstep"]
Pstep = dataset["mesonBSE"]["Pstep"]
dim = kstep * zstep
# 取点方式
plist=[0.001+1/Pstep*(i-1) for i=1:Pstep]
meshk,weightk= gausslegendremesh(cutdown,cutup,kstep,2);
meshz,weightz= gausschebyshev(zstep,2);

function Inport()
    global z2, z4
    local A, B, k
    A, B, k, z2, z4=load("data/quark_gap_equation/quark-ABkz2z4-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint.jld2","A","B","k", "z2", "z4");
    global AA
    global BB
    AA=Spline1D(k,A)
    BB=Spline1D(k,B)
    return true
end
Inport()

D(t::Float64)=8*pi^2*(dd*exp(-t/(ω^2))/ω^4+rm*((-expm1(-t/(4*mt)^2))/t)/log(τ+(1+t/Λ^2)^2))
branchfunction(x::Float64)=(x*AA(x)^2+BB(x)^2)
#切比雪夫展开
include(joinpath(pwd(),"src/mesonfile/chebyshevD.jl"))

print("参数导入完毕,开始计算\n")
if dataset["mesonBSE"]["mesonmode"] == 1
    if ispath("data/pseudo_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint") == false
        mkdir("data/pseudo_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint")
        cp(joinpath(pwd(),"src/config.toml"),joinpath(pwd(),"data/pseudo_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint/log.toml"))
        for indexforp2=1:Pstep
            timetest1=time()
            global P2
            P2=plist[indexforp2]
            # 分配点与权重
            include(joinpath(pwd(),"src/mesonfile/setupkernel.jl"))
            # 计算kernel
            include(joinpath(pwd(),"src/mesonfile/mode1.jl"))
            # 求解函数
            include(joinpath(pwd(),"src/mesonfile/solve_kernel.jl"))
            # 保存文件
            jldsave("data/pseudo_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint/P&F1-4_$indexforp2-$Pstep.jld2";P2, F1, F2, F3, F4)
            print("$P2 for $indexforp2/$Pstep done, takes",round((time()-timetest1)*100)/100,"s \n")
        end
    else # 
        print("已存在文件")
    end
elseif dataset["mesonBSE"]["mesonmode"] == 2
    if ispath("data/scalar_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint") == false
        mkdir("data/scalar_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint")
        cp(joinpath(pwd(),"src/config.toml"),joinpath(pwd(),"data/scalar_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint/log.toml"))
        for indexforp2=1:Pstep
            timetest1=time()
            global P2
            P2=plist[indexforp2]
            # 分配点与权重
            include(joinpath(pwd(),"src/mesonfile/setupkernel.jl"))
            # 计算kernel
            include(joinpath(pwd(),"src/mesonfile/mode2.jl"))
            # 求解函数
            include(joinpath(pwd(),"src/mesonfile/solve_kernel.jl"))
            # 保存文件
            jldsave("data/scalar_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint/P&F1-4_$indexforp2-$Pstep.jld2";P2, F1, F2, F3, F4)
            print("$P2 for $indexforp2/$Pstep done, takes",round((time()-timetest1)*100)/100,"s \n")
        end # for
    else # 
        print("已存在文件")
    end # ispath
end # if for mode
end # module