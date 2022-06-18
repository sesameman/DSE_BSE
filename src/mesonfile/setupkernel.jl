# 初始化变量
kernel=Array{Float64}(undef,kstep , zstep, kstep, zstep, 4, 4)
A1=Array{Float64}(undef,kstep,zstep);
B1=Array{Float64}(undef,kstep,zstep);
A2=Array{Float64}(undef,kstep,zstep);
B2=Array{Float64}(undef,kstep,zstep);
for j=1:kstep
for m=1:zstep
    A1[j,m]=AA((P2/4+meshk[j]+sqrt(P2*meshk[j])*meshz[m]))
    B1[j,m]=BB((P2/4+meshk[j]+sqrt(P2*meshk[j])*meshz[m]))
    A2[j,m]=AA((P2/4+meshk[j]-sqrt(P2*meshk[j])*meshz[m]))
    B2[j,m]=BB((P2/4+meshk[j]-sqrt(P2*meshk[j])*meshz[m]))
end
end


# print("完成kernel准备,用时",round((time()-timetest)*100)/100,"s \n")
