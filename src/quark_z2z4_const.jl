module quarksub
export k, A, B, z2, z4, M

# using TOML
# workdir = "/Users/kjy/Desktop/program/julia/DSE_BSE"
# # 放在一个模块儿里运行，与主题函数隔开
# cd(workdir)


using TOML
using Gaussquad
using Dierckx
using JLD2
using Main.quarkdse: z2, z4
#=
现在定义命名规范
quark-ABkz2z4-$quarkmass-$logofcutoff-$quarkintstep-$repoint.jld2
=#
print("----------------quarksub-------------------\n")
dataset=TOML.parsefile("src/config.toml")

repoint = dataset["quarkDSE"]["repoint"]
intstep = dataset["quarkDSE"]["quarkintstep"]
quarkm = dataset["quarkDSE"]["quarkmass"]
logofcutoff = dataset["mesonBSE"]["logofcutoff"]

bigm = dataset["subterm"]["bigm"]

# 检测是否有重复的文件
function testfile()
    isfile("data/quark_gap_equation/holdz2z4-$bigm-$quarkm-$logofcutoff-$intstep-$repoint.jld2")
end

if testfile() == false
print("文件不存在,正在计算相应于cutoff=10^$logofcutoff--z2=$z2--z4=$z4--的剪除项","\n")
mt = dataset["data"]["mt"]
τ = dataset["data"]["tau"]
Λ = dataset["data"]["lambda"]
ω = dataset["data"]["omega"]
dd = dataset["data"]["dd"]
Nf = dataset["data"]["Nf"]
rm = dataset["data"]["rm"]
cutup = 10. ^logofcutoff
cutdown = 10. ^(-logofcutoff)


A = Array{Float64}(undef,intstep,1)
B = Array{Float64}(undef,intstep,1)
# 点与权重
k,w = gausslegendremesh(cutdown,cutup,intstep,2)
delta(x,y) = ==(x,y)
F(x, mt) = (1-exp(-x/(4*mt)^2))/x
D(t, dd, ω, rm, mt, τ, Λ) = 8*pi^2*(dd*exp(-t/(ω^2))/ω^4+rm*F(t, mt)/log(τ+(1+t/Λ^2)^2))
# 先做角度积分
# 采用第二类切比雪夫
# i代表外动量，j代表内动量
IntA=Array{Float64}(undef,intstep+1,intstep)
IntB=Array{Float64}(undef,intstep+1,intstep)
IntBz4=Array{Float64}(undef,intstep)
Threads.@threads for i= 1:intstep+1
    for j = 1:intstep
        q2=k[j]
        if i==intstep+1
            k2=repoint ^2
        else
            k2=k[i]
        end
        kdotq(z)=sqrt(k2*q2)*z
        IntA[i,j]=(4/3)*gausschebyshevint512(z->D(k2+q2-2*kdotq(z), dd, ω, rm, mt, τ, Λ)*(3*k2*kdotq(z)+3*q2*kdotq(z)-2*k2*q2-4*(kdotq(z))^2)/k2/(k2+q2-2*kdotq(z)))/(2*pi)^3
        IntB[i,j]=(4/3)*gausschebyshevint512(z->D(k2+q2-2*kdotq(z), dd, ω, rm, mt, τ, Λ))/(2*pi)^3
    end
end
function FAf(i)
    sum=0.::Float64
    for j=1:intstep
        sum+=k[j]*A[j]/(k[j]*A[j]^2+B[j]^2)*IntA[i,j]*w[j]
    end
    result=A[i]-z2-z2^2*sum
    return result
end

function FBf(i)
    sum=0.::Float64
    for j=1:intstep
        sum+=3*k[j]*B[j]/(k[j]*A[j]^2+B[j]^2)*IntB[i,j]*w[j]
    end
    result=B[i]-z4*bigm-z2^2*sum
    return result
end

function renormpoint()
    sumA=0.::Float64
    sumB=0.::Float64
    for j=1:intstep
        sumA+=k[j]*A[j]/(k[j]*A[j]^2+B[j]^2)*IntA[intstep+1,j]*w[j]
        sumB+=3*k[j]*B[j]/(k[j]*A[j]^2+B[j]^2)*IntB[intstep+1,j]*w[j]
    end
    return sumA,sumB
end

jacobifAA(i,j)=delta(i,j)-z2^2*w[j]*k[j]*(B[j]^2-k[j]*A[j]^2)/(k[j]*A[j]^2+B[j]^2)^2*IntA[i,j]
jacobifAB(i,j)=-z2^2*w[j]*k[j]*(-2*A[j]*B[j])/(k[j]*A[j]^2+B[j]^2)^2*IntA[i,j]
jacobifBA(i,j)=-3*z2^2*w[j]*k[j]*(-2*A[j]*B[j]*k[j])/(k[j]*A[j]^2+B[j]^2)^2*IntB[i,j]
jacobifBB(i,j)=delta(i,j)-3*z2^2*w[j]*k[j]*(k[j]*A[j]^2-B[j]^2)/(k[j]*A[j]^2+B[j]^2)^2*IntB[i,j]

A=fill(1.7,intstep);
B=fill(1.2,intstep);
Δ=fill(1.,2*intstep)
st=0::Int
# z4old=0.

# k19sumA, k19sumB=repointalpoint()
while st<500 && maximum(abs.(Δ))>10^-8 
# while st<1
    global A, B, Δ, st, z2, z4, z2old, FA, FB, k19sumB
    st+=1
    k19sumA, k19sumB=renormpoint()
    FA=[FAf(i) for i=1:(intstep)]
    FB=[FBf(i) for i=1:(intstep)]
    jacobiAA=[jacobifAA(i,j) for i=1:(intstep),j=1:(intstep)]
    jacobiAB=[jacobifAB(i,j) for i=1:(intstep),j=1:(intstep)]
    jacobiBA=[jacobifBA(i,j) for i=1:(intstep),j=1:(intstep)]
    jacobiBB=[jacobifBB(i,j) for i=1:(intstep),j=1:(intstep)]
    jacobi=[jacobiAA jacobiAB; jacobiBA jacobiBB]
    Δ=jacobi\[FA; FB]
    A-=Δ[1:(intstep)]
    B-=Δ[(intstep+1):(2intstep)]
    if st == 499
        print("未达到要求，已循环500次，并未收敛")
    end
end
print("循环次数为 $st 次\n")
AA=Spline1D(k,A)
BB=Spline1D(k,B)
print("无穷远点的应取值A=", z2, ",","B=",z4*bigm,"\n")
print("无穷远点真实值为A=", AA(10^logofcutoff), ",","B=",BB(10^logofcutoff),"\n")
print("z2=$z2\nz4=$z4\n")

jldsave("data/quark_gap_equation/holdz2z4-$bigm-$quarkm-$logofcutoff-$intstep-$repoint.jld2";A, B, k, z2, z4)
# print("已保存到",joinpath(pwd(),"data/quark_gap_equation\n"))
else # if for testfile
print("有现有文件,已读取数据\n")
A, B, k=load("data/quark_gap_equation/holdz2z4-$bigm-$quarkm-$logofcutoff-$intstep-$repoint.jld2","A","B","k")
AA=Spline1D(k,A)
BB=Spline1D(k,B)
print("无穷远点点的理想值A=", z2, ",","B=",z4*bigm,"\n")
print("无穷远点真实值为A=", AA(10^logofcutoff), ",","B=",BB(10^logofcutoff),"\n")
print("z2=$z2\nz4=$z4\n")
end #if


# M=[B[i]/A[i] for i=1:intstep]


end #module
