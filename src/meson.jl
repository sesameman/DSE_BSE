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
include(joinpath(pwd(),"src/moudule_dse/inter.jl"))

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
cutdown = 10. ^(-4)

kstep = dataset["mesonBSE"]["kstep"]
zstep = dataset["mesonBSE"]["zstep"]
Pstep = dataset["mesonBSE"]["Pstep"]
dim = kstep * zstep
# 取点方式
plist=gausslegendremesh(cutdown,cutup,Pstep,2)[1]
if Pstep == 1
    plist = [0.001]
end 
meshk,weightk= gausslegendremesh(cutdown,cutup,kstep,2);
meshz,weightz= gausschebyshev(zstep,2);

D(t::Float64)=8*pi^2*(dd*exp(-t/(ω^2))/ω^4+rm*((-expm1(-t/(4*mt)^2))/t)/log(τ+(1+t/Λ^2)^2))
#切比雪夫展开
include(joinpath(pwd(),"src/mesonfile/chebyshevD.jl"))

consta1=0.
consta2=0.
constb1=0.
constb2=0.
function Inport()
    global z2, z4, AA1, BB1, consta1, consta2, constb1, constb2
    local A, B, k, mk
    A, B, k, z2, z4=load("data/quark_gap_equation/quark-ABkz2z4-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint.jld2","A","B","k", "z2", "z4");
    consta1, constb1 = Main.interDSE.ExtraA(k,A)
    consta2, constb2 = Main.interDSE.ExtraB(k,B)
    AA1=Spline1D(k,A)
    BB1=Spline1D(k,B)
    mk=k[length(k)]
    return mk
end
maxofk=Inport()
function AA(x)
    if x<maxofk
        return AA1(x)
    else
        return consta1/log(x/0.234^2)^constb1
    end
end
function BB(x)
    if x<maxofk
        return BB1(x)
    else
        return consta2/log(x/0.234^2)^constb2
    end
end
branchfunction(x::Float64)=(x*AA(x)^2+BB(x)^2)

print("参数导入完毕,开始计算mesonBSA\n")
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
            print("mesonBSA--$P2 for $indexforp2/$Pstep done, takes",round((time()-timetest1)*100)/100,"s \n")
        end
    else # 
        print("已存在文件--",joinpath(pwd(),"data/pseudo_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint"),"\n")
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
            print("mesonBSA--$P2 for $indexforp2/$Pstep done, takes",round((time()-timetest1)*100)/100,"s \n")
        end # for
    else # 
        print("已存在文件--",joinpath(pwd(),"data/scalar_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint"),"\n")
    end # ispath
end # if for mode


# if dataset["subterm"]["calculate"] == 1
#     bigm = dataset["subterm"]["bigm"]
#     function Inport()
#         global z2, z4
#         local A, B, k
#         A, B, k, z2, z4=load("data/quark_gap_equation/holdz2z4-$bigm-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint.jld2","A","B","k", "z2", "z4");
#         global AA
#         global BB
#         AA=Spline1D(k,A)
#         BB=Spline1D(k,B)
#         return true
#     end
#     Inport()
#     branchfunction(x::Float64)=(x*AA(x)^2+BB(x)^2)
    
#     print("开始计算减除项\n")
#     if dataset["mesonBSE"]["mesonmode"] == 1
#         if ispath("data/pseudo_BSE/submeson-$kstep-$zstep-$Pstep-$bigm-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint") == false
#             mkdir("data/pseudo_BSE/submeson-$kstep-$zstep-$Pstep-$bigm-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint")
#             cp(joinpath(pwd(),"src/config.toml"),joinpath(pwd(),"data/pseudo_BSE/submeson-$kstep-$zstep-$Pstep-$bigm-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint/log.toml"))
#             for indexforp2=1:Pstep
#                 timetest1=time()
#                 global P2
#                 P2=plist[indexforp2]
#                 # 分配点与权重
#                 include(joinpath(pwd(),"src/mesonfile/setupkernel.jl"))
#                 # 计算kernel
#                 include(joinpath(pwd(),"src/mesonfile/mode1.jl"))
#                 # 求解函数
#                 include(joinpath(pwd(),"src/mesonfile/solve_kernel.jl"))
#                 # 保存文件
#                 jldsave("data/pseudo_BSE/submeson-$kstep-$zstep-$Pstep-$bigm-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint/P&F1-4_$indexforp2-$Pstep.jld2";P2, F1, F2, F3, F4)
#                 print("subterm--$P2 for $indexforp2/$Pstep done, takes",round((time()-timetest1)*100)/100,"s \n")
#             end
#         else # 
#             print("已存在剪除项")
#         end
#     elseif dataset["mesonBSE"]["mesonmode"] == 2
#         if ispath("data/scalar_BSE/submeson-$kstep-$zstep-$Pstep-$bigm-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint") == false
#             mkdir("data/scalar_BSE/submeson-$kstep-$zstep-$Pstep-$bigm-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint")
#             cp(joinpath(pwd(),"src/config.toml"),joinpath(pwd(),"data/scalar_BSE/submeson-$kstep-$zstep-$Pstep-$bigm-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint/log.toml"))
#             for indexforp2=1:Pstep
#                 timetest1=time()
#                 global P2
#                 P2=plist[indexforp2]
#                 # 分配点与权重
#                 include(joinpath(pwd(),"src/mesonfile/setupkernel.jl"))
#                 # 计算kernel
#                 include(joinpath(pwd(),"src/mesonfile/mode2.jl"))
#                 # 求解函数
#                 include(joinpath(pwd(),"src/mesonfile/solve_kernel.jl"))
#                 # 保存文件
#                 jldsave("data/scalar_BSE/submeson-$kstep-$zstep-$Pstep-$bigm-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint/P&F1-4_$indexforp2-$Pstep.jld2";P2, F1, F2, F3, F4)
#                 print("subterm--$P2 for $indexforp2/$Pstep done, takes",round((time()-timetest1)*100)/100,"s \n")
#             end # for
#         else # 
#             print("已剪除项")
#         end # ispath
#     end # if for mode
# end # for subterm
end # module