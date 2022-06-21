workdir = "/Users/kjy/Desktop/program/julia/DSE_BSE"
# 放在一个模块儿里运行，与主题函数隔开
cd(workdir)
using TOML
using LinearAlgebra
using Dierckx
using JLD2
using Gaussquad # View in github: https://github.com/kangjiayin/Gaussquad.jl
using FastGaussQuadrature
using ChebyshevFun 


dataset = TOML.parsefile("src/config.toml")

# data
logofcutoff = dataset["subterm"]["masscutoff"]
cutup = 10. ^logofcutoff
cutdown = 10. ^(-logofcutoff)

# readsetting
kstep = dataset["readsetting"]["kstep"]
zstep = dataset["readsetting"]["zstep"]

Pstep = dataset["subterm"]["Pstep"]

# quark system
quarkrepoint = dataset["readsetting"]["repoint"]
quarkintstep = dataset["readsetting"]["quarkintstep"]
quarkm = dataset["readsetting"]["quarkmass"]

plist=[0.001+1/Pstep*(i-1) for i=1:Pstep]
meshk,weightk= gausslegendremesh(cutdown,cutup,kstep,2);
meshz,weightz= gausschebyshev(zstep,2);

function Inport()
    global z2, z4
    local A, B, k
    A, B, k, z2, z4=load("data/quark_gap_equation/quark-ABkz2z4-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint.jld2","A","B","k", "z2", "z4");
    print("确认导入的数据为--", "quark-ABkz2z4-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint.jld2\n")
    # global AA
    # global BB
    AA=Spline1D(k,A)
    BB=Spline1D(k,B)
    return true
end
Inport()

constant1 = 1
constant2 = -2
bigm1 = dataset["quarkDSE"]["quarkmass"]^2
bigm2 = dataset["quarkDSE"]["quarkmass"]^2

psub_form1_pqz=zeros(Pstep,kstep,zstep)
psub_form1_pq=zeros(Pstep, kstep)
psub_form1_p=zeros(Pstep)

psub_form2_pqz=zeros(Pstep,kstep,zstep)
psub_form2_pq=zeros(Pstep, kstep)
psub_form2_p=zeros(Pstep)

psub=zeros(Pstep)

# psub_form2_pqz=zeros(Pstep,kstep,zstep)
# psub_form2_pq=zeros(Pstep, kstep)
# psub_form2_p=zeros(Pstep)

Threads.@threads for i = 1:Pstep
    for j = 1:kstep
        # i for P2[i], j for q2[j]
        for m = 1:zstep
            P2i = plist[i]
            q2i = meshk[j]
            zq = meshz[m]
            psub_form1_pqz[i,j,m] = 4*( q2i + bigm1 - P2i/4)/(q2i+P2i+sqrt(q2i*P2i)*zq+bigm1)/(q2i+P2i-sqrt(q2i*P2i)*zq+bigm1)
            psub_form2_pqz[i,j,m] = 4*( q2i + bigm2 - P2i/4)/(q2i+P2i+sqrt(q2i*P2i)*zq+bigm2)/(q2i+P2i-sqrt(q2i*P2i)*zq+bigm2)
            # psub_form2_pqz[i,j,m] = 2/q2i+(-2*bigm1-3/2*P2i+sqrt(P2i)*zq)/q2i/q2i
            psub_form1_pq[i,j] += weightz[m]*psub_form1_pqz[i,j,m]
            psub_form2_pq[i,j] += weightz[m]*psub_form2_pqz[i,j,m]
            # psub_form2_pq[i,j] += weightz[m]*psub_form2_pqz[i,j,m]
            # f_gammaPq[i,j] += weightz[m]*(P2[i]*A2[i,j,m]*((F1k[i,j,m] + 4*F4k[i,j,m]*q2[j])*A1[i,j,m] - 2*F2k[i,j,m]*B1[i,j,m]) - 2*F2k[i,j,m]*P2[i]*A1[i,j,m]*B2[i,j,m] + 4*sqrt(P2[i]*q2[j])*(F2k[i,j,m] + F3k[i,j,m]*q2[j])*meshz[m]*(A2[i,j,m]*B1[i,j,m] - A1[i,j,m]*B2[i,j,m]) - 2*P2[i]*q2[j]*meshz[m]^2*(2*F4k[i,j,m]*A1[i,j,m]*A2[i,j,m] + F3k[i,j,m]*A2[i,j,m]*B1[i,j,m] + F3k[i,j,m]*A1[i,j,m]*B2[i,j,m]) - 4*F1k[i,j,m]*(q2[j]*A1[i,j,m]*A2[i,j,m] + B1[i,j,m]*B2[i,j,m]))
        end
    end
end


for i=1:Pstep
    for j=1:kstep
        # i for P2[i], j for q2[j]
        psub_form1_p[i] += weightk[j]*psub_form1_pq[i,j]*meshk[j]
        psub_form2_p[i] += weightk[j]*psub_form2_pq[i,j]*meshk[j]
        # psub_form2_p[i]+=weightk[j]*psub_form2_pq[i,j]*meshk[j]
    end
end

psub_form1_p = psub_form1_p / 8 / pi^3  * z4^2 / z2^2 * constant1 * 3
psub_form2_p = psub_form2_p / 8 / pi^3  * z4^2 / z2^2 * constant2 * 3

psub = psub_form1_p + psub_form2_p

print(psub_form1_p,"\n")
print(psub_form2_p,"\n")
print(psub)
# psub_form2_p = psub_form2_p / 32 / pi^2 * (-1) * z4 / z2

# print(psub_form2_p)