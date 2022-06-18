function Inport()
    global z2, z4
    local A, B, k
    A, B, k, z2, z4=load("/Users/kjy/Desktop/program/julia/Gamma_23/data/quark_equation/ABkz2z4334.jld2","A","B","k", "z2", "z4");
    global AA
    global BB
    AA=Spline1D(k,A)
    BB=Spline1D(k,B)
    return true
end
# 导入数据
Inport();
A1=Array{Float64}(undef,Pstep,kstep,zstep);
B1=Array{Float64}(undef,Pstep,kstep,zstep);
A2=Array{Float64}(undef,Pstep,kstep,zstep);
B2=Array{Float64}(undef,Pstep,kstep,zstep);
for i=1:Pstep
for j=1:kstep
for m=1:zstep
    A1[i,j,m]=AA((P2[i]/4+meshk[j]+sqrt(P2[i]*meshk[j])*meshz[m]))
    B1[i,j,m]=BB((P2[i]/4+meshk[j]+sqrt(P2[i]*meshk[j])*meshz[m]))
    A2[i,j,m]=AA((P2[i]/4+meshk[j]-sqrt(P2[i]*meshk[j])*meshz[m]))
    B2[i,j,m]=BB((P2[i]/4+meshk[j]-sqrt(P2[i]*meshk[j])*meshz[m]))
end
end
end

# q2=meshk
f_gamma5Pqzq=zeros(Pstep,kstep,zstep)
f_gamma5Pq=zeros(Pstep, kstep)
f_gamma5P=zeros(Pstep)

f_gamma5muPqzq=zeros(Pstep,kstep,zstep)
f_gamma5muPq=zeros(Pstep, kstep)
f_gamma5muP=zeros(Pstep)

f_mo0rePqzq_why=zeros(Pstep,kstep,zstep)
f_mo0rePq_why=zeros(Pstep, kstep)
f_mo0reP_why=zeros(Pstep)

f_mo0rePqzq=zeros(Pstep,kstep,zstep)
f_mo0rePq=zeros(Pstep, kstep)
f_mo0reP=zeros(Pstep)

f_mo0imPqzq=zeros(Pstep,kstep,zstep)
f_mo0imPq=zeros(Pstep, kstep)
f_mo0imP=zeros(Pstep)

f_mo2rePqzq=zeros(Pstep,kstep,zstep)
f_mo2rePq=zeros(Pstep, kstep)
f_mo2reP=zeros(Pstep)

f_mo2imPqzq=zeros(Pstep,kstep,zstep)
f_mo2imPq=zeros(Pstep, kstep)
f_mo2imP=zeros(Pstep)

s_Pi_Pqzq=zeros(Pstep,kstep,zstep)
s_Pi_Pq=zeros(Pstep, kstep)
s_Pi_P=zeros(Pstep)

mz4 = 0.
Threads.@threads for i=1:Pstep
    for j=1:kstep
        # i for P2[i], j for q2[j]
        for m=1:zstep
            P2i=P2[i]
            q2i=meshk[j]
            zq=meshz[m]
            A1i=A1[i,j,m]
            A2i=A2[i,j,m]
            B1i=B1[i,j,m]
            B2i=B2[i,j,m]
            F1=F1k[i,j,m]
            F2=F2k[i,j,m]
            F3=F3k[i,j,m]
            F4=F4k[i,j,m]
            mF1=mF1k[i,j,m]
            mF2=mF2k[i,j,m]
            mF3=mF3k[i,j,m]
            mF4=mF4k[i,j,m]
            # F1=0.
            # F2=0.
            # F3=0.
            # F4=0.
            forget_div=1/(A1i^2*(P2i/4+q2i+sqrt(P2i*q2i)*zq)+B1i^2)/(A2i^2*(P2i/4+q2i-sqrt(P2i*q2i)*zq)+B2i^2)
            f_gamma5Pqzq[i,j,m]=(
                4*(F1 - mz4)*(B1i*B2i + A1i*A2i*q2i) + P2i*(2*A2i*B1i*(F2 - mz4) + 2*A1i*B2i*(F2 - mz4) - A1i*A2i*(F1 - mz4 + 4*(F4 - mz4)*q2i)) - 4*(A2i*B1i - A1i*B2i)*sqrt(P2i*q2i)*(F2 - mz4 + (F3 - mz4)*q2i)*zq + 2*(A2i*B1i*(F3 - mz4) + A1i*B2i*(F3 - mz4) + 2*A1i*A2i*(F4 - mz4))*P2i*q2i*zq^2
                )*forget_div
            f_gamma5muPqzq[i,j,m]=(
                2*A2i*B1i*F1 + 2*A1i*B2i*F1 - 4*B1i*B2i*F2 + A1i*A2i*F2*P2i + 4*(A1i*A2i*F2 + A2i*B1i*F4 + A1i*B2i*F4)*q2i + (4*(-(A2i*B1i) + A1i*B2i)*F1*q2i*zq)/sqrt(P2i*q2i) + A1i*A2i*F3*P2i*q2i*zq^2 - 4*q2i*(2*A1i*A2i*F2 + B1i*B2i*F3 + A2i*B1i*F4 + A1i*B2i*F4 + A1i*A2i*F3*q2i)*zq^2
                )*forget_div
            f_mo0rePqzq_why[i,j,m]=(
                A1i*A2i*F2*P2i + 2*(-2*B1i*B2i*F2 + A2i*B1i*(F1 + 2*F4*q2i - F4*sqrt(P2i*q2i)*zq) + A1i*B2i*(F1 + 2*F4*q2i + F4*sqrt(P2i*q2i)*zq) + A1i*A2i*q2i*(2*F2 + F3*P2i*zq^2))
                )*forget_div
            f_mo0rePqzq[i,j,m]=(
                A1i*A2i*F2*(P2i + 4*q2i) + A1i*A2i*(-8*F2 + F3*(P2i - 4*q2i))*q2i*zq^2 - 4*B1i*B2i*(F2 + F3*q2i*zq^2) + 2*A2i*B1i*(F1 - 2*F1*sqrt(q2i/P2i)*zq - 2*F4*q2i*(-1 + zq^2)) + 2*A1i*B2i*(F1 + 2*F1*sqrt(q2i/P2i)*zq - 2*F4*q2i*(-1 + zq^2))
                )*forget_div
            f_mo0imPqzq[i,j,m]=(
                sqrt(q2i/P2i)*sqrt(1 - zq^2)*(4*A2i*B1i*F1 - 4*A1i*B2i*F1 + 4*sqrt(P2i*q2i)*(2*A1i*A2i*F2 + B1i*B2i*F3 + A2i*B1i*F4 + A1i*B2i*F4 + A1i*A2i*F3*q2i)*zq + P2i*(-2*A2i*B1i*F4 + 2*A1i*B2i*F4 + A1i*A2i*F3*sqrt(P2i*q2i)*zq))
                )*forget_div
            f_mo2rePqzq[i,j,m]=(
                (q2i/P2i)^1.5*zq*(-3 + 4*zq^2)*(-4*A2i*B1i*F1 + 4*A1i*B2i*F1 - 4*sqrt(P2i*q2i)*(2*A1i*A2i*F2 + B1i*B2i*F3 + A2i*B1i*F4 + A1i*B2i*F4 + A1i*A2i*F3*q2i)*zq + P2i*(2*A2i*B1i*F4 - 2*A1i*B2i*F4 - A1i*A2i*F3*sqrt(P2i*q2i)*zq)) + (q2i*(-1 + 2*zq^2)*(-4*B1i*B2i*F2 + 2*A2i*B1i*(F1 + 2*F4*q2i - F4*sqrt(P2i*q2i)*zq) + 2*A1i*B2i*(F1 + 2*F4*q2i + F4*sqrt(P2i*q2i)*zq) + A1i*A2i*(F2*(P2i + 4*q2i) + 2*F3*P2i*q2i*zq^2)))/P2i
                )*forget_div
            f_mo2imPqzq[i,j,m]=(
                sqrt(1 - zq^2)*((q2i/P2i)^1.5*(1 - 4*zq^2)*(-4*A2i*B1i*F1 + 4*A1i*B2i*F1 - 4*sqrt(P2i*q2i)*(2*A1i*A2i*F2 + B1i*B2i*F3 + A2i*B1i*F4 + A1i*B2i*F4 + A1i*A2i*F3*q2i)*zq + P2i*(2*A2i*B1i*F4 - 2*A1i*B2i*F4 - A1i*A2i*F3*sqrt(P2i*q2i)*zq)) - (2*q2i*zq*(-4*B1i*B2i*F2 + 2*A2i*B1i*(F1 + 2*F4*q2i - F4*sqrt(P2i*q2i)*zq) + 2*A1i*B2i*(F1 + 2*F4*q2i + F4*sqrt(P2i*q2i)*zq) + A1i*A2i*(F2*(P2i + 4*q2i) + 2*F3*P2i*q2i*zq^2)))/P2i)
                )*forget_div
            s_Pi_Pqzq[i,j,m]=(
                -4*B1i*B2i*F1 + 4*(A1i*A2i*F1 + A2i*B1i*F3 + A1i*B2i*F3)*q2i - P2i*(A1i*A2i*(F1 + 4*F4*q2i) + 2*(A2i*B1i - A1i*B2i)*F2*sqrt(P2i*q2i)*zq) + 2*sqrt(P2i*q2i)*zq*(-(A2i*B1i*F3) + A1i*B2i*F3 + 2*(A2i*B1i*F2 + A1i*B2i*F2 + A1i*A2i*F4)*sqrt(P2i*q2i)*zq)
                )*forget_div
            f_gamma5Pq[i,j] += weightz[m]*f_gamma5Pqzq[i,j,m]
            f_gamma5muPq[i,j] += weightz[m]* f_gamma5muPqzq[i,j,m]
            f_mo0rePq_why[i,j] += weightz[m]* f_mo0rePqzq_why[i,j,m]
            f_mo0rePq[i,j] += weightz[m]* f_mo0rePqzq[i,j,m]
            f_mo0imPq[i,j] += weightz[m]* f_mo0imPqzq[i,j,m]
            f_mo2rePq[i,j] += weightz[m]* f_mo2rePqzq[i,j,m]
            f_mo2imPq[i,j] += weightz[m]* f_mo2imPqzq[i,j,m]
            s_Pi_Pq[i,j] += weightz[m]* s_Pi_Pqzq[i,j,m]
            # f_gammaPq[i,j] += weightz[m]*(P2[i]*A2[i,j,m]*((F1k[i,j,m] + 4*F4k[i,j,m]*q2[j])*A1[i,j,m] - 2*F2k[i,j,m]*B1[i,j,m]) - 2*F2k[i,j,m]*P2[i]*A1[i,j,m]*B2[i,j,m] + 4*sqrt(P2[i]*q2[j])*(F2k[i,j,m] + F3k[i,j,m]*q2[j])*meshz[m]*(A2[i,j,m]*B1[i,j,m] - A1[i,j,m]*B2[i,j,m]) - 2*P2[i]*q2[j]*meshz[m]^2*(2*F4k[i,j,m]*A1[i,j,m]*A2[i,j,m] + F3k[i,j,m]*A2[i,j,m]*B1[i,j,m] + F3k[i,j,m]*A1[i,j,m]*B2[i,j,m]) - 4*F1k[i,j,m]*(q2[j]*A1[i,j,m]*A2[i,j,m] + B1[i,j,m]*B2[i,j,m]))
        end
    end
end

for i=1:Pstep
    for j=1:kstep
        # i for P2[i], j for q2[j]
        f_gamma5P[i]+=weightk[j]*f_gamma5Pq[i,j]*meshk[j]
        # if  f_gamma5muPq[i,j]<0
        #     f_gamma5muPq[i,j]=0
        # end
        f_gamma5muP[i]+=weightk[j]*f_gamma5muPq[i,j]*meshk[j]
        f_mo0reP_why[i]+=weightk[j]*f_mo0rePq_why[i,j]*meshk[j]
        f_mo0reP[i]+=weightk[j]*f_mo0rePq[i,j]*meshk[j]
        f_mo0imP[i]+=weightk[j]*f_mo0imPq[i,j]*meshk[j]
        f_mo2reP[i]+=weightk[j]*f_mo2rePq[i,j]*meshk[j]
        f_mo2imP[i]+=weightk[j]*f_mo2imPq[i,j]*meshk[j]
        s_Pi_P[i]+=weightk[j]*s_Pi_Pq[i,j]*meshk[j]
        
    end
end
f_gamma5P=f_gamma5P *(2*pi)^-3*z4/2*3
f_gamma5muP=f_gamma5muP *(2*pi)^-3*z2/2*3
f_mo0reP_why=f_mo0reP_why *(2*pi)^-3*z2/2*3
f_mo0reP=f_mo0reP *(2*pi)^-3*z2/2*3
f_mo0imP=f_mo0imP *(2*pi)^-3*z2/2*3
f_mo2reP=f_mo2reP *(2*pi)^-3*z2/2*3*2^2
f_mo2imP=f_mo2imP *(2*pi)^-3*z2/2*3*2^2
s_Pi_P=s_Pi_P *(2*pi)^-3*z4/2*3
