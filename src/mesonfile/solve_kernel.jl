timetest=time()
# 压平二维数组
function flat(x::Array)
    output=Array{Float64}(undef,dim)
    for i=1:kstep
        for j=1:zstep
            output[j+(i-1)*zstep]=x[i,j]
        end
    end
    return output
end
function getkernel(n1,n2)
    kernelsetup=Array{Float64}(undef,dim,dim)
    for j=1:kstep
        for m=1:zstep
            kernelsetup[:,m+(j-1)*zstep]=flat(kernel[:,:,j,m,n1,n2])
        end
    end
    return kernelsetup
end

kernelsolve=[getkernel(1,1) getkernel(1,2) getkernel(1,3) getkernel(1,4);
             getkernel(2,1) getkernel(2,2) getkernel(2,3) getkernel(2,4);
             getkernel(3,1) getkernel(3,2) getkernel(3,3) getkernel(3,4);
             getkernel(4,1) getkernel(4,2) getkernel(4,3) getkernel(4,4)]

kernel=nothing
right=[z4*ones(dim) ;zeros(3*dim)]
kernelsolve=I-kernelsolve
solution=kernelsolve\right
right=nothing

# # 切比雪夫
# F1=zeros(kstep)
# F2=zeros(kstep)
# F3=zeros(kstep)
# F4=zeros(kstep)

# for u=1:4*dim
#     if u<=dim
#         u1=u
#         nk=((u1-1)÷zstep+1)
#         nz=((u1-1)%zstep+1)
#         F1[nk]+=solution[u]*weightz[nz]*2/pi
#     elseif u<=2*dim
#         u1=u-dim
#         nk=((u1-1)÷zstep+1)
#         nz=((u1-1)%zstep+1)
#         F2[nk]+=solution[u]*weightz[nz]*2/pi
#     elseif u<=3*dim
#         u1=u-2*dim
#         nk=((u1-1)÷zstep+1)
#         nz=((u1-1)%zstep+1)
#         F3[nk]+=solution[u]*weightz[nz]*2/pi
#     elseif u<=4*dim
#         u1=u-3*dim
#         nk=((u1-1)÷zstep+1)
#         nz=((u1-1)%zstep+1)
#         F4[nk]+=solution[u]*weightz[nz]*2/pi
#     end
# end

# 切比雪夫
F1=zeros(kstep,zstep)
F2=zeros(kstep,zstep)
F3=zeros(kstep,zstep)
F4=zeros(kstep,zstep)

for u=1:4*dim
    if u<=dim
        u1=u
        nk=((u1-1)÷zstep+1)
        nz=((u1-1)%zstep+1)
        F1[nk,nz]=solution[u]
    elseif u<=2*dim
        u1=u-dim
        nk=((u1-1)÷zstep+1)
        nz=((u1-1)%zstep+1)
        F2[nk,nz]=solution[u]
    elseif u<=3*dim
        u1=u-2*dim
        nk=((u1-1)÷zstep+1)
        nz=((u1-1)%zstep+1)
        F3[nk,nz]=solution[u]
    elseif u<=4*dim
        u1=u-3*dim
        nk=((u1-1)÷zstep+1)
        nz=((u1-1)%zstep+1)
        F4[nk,nz]=solution[u]
    end
end

# print("完成kernel求解，用时",round((time()-timetest)*100)/100,"s \n")
# print("预计误差", round(abs(B1[1,1]-F1[1]*0.003)/B1[1,1]*100000)/1000,"%\n")