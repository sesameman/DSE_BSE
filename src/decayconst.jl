

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
pseudo_vp_pqzq=zeros(Pstep,kstep,zstep)
pseudo_vp_pq=zeros(Pstep, kstep)
pseudo_vp_p=zeros(Pstep)

f_gamma5muPqzq=zeros(Pstep,kstep,zstep)
f_gamma5muPq=zeros(Pstep, kstep)
f_gamma5muP=zeros(Pstep)

f_gamma5Pqzq=zeros(Pstep,kstep,zstep)
f_gamma5Pq=zeros(Pstep, kstep)
f_gamma5P=zeros(Pstep)

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
            pseudo_vp_pqzq[i,j,m]=(
                4*F1*(B1i*B2i + A1i*A2i*q2i) + P2i*(2*A2i*B1i*F2 + 2*A1i*B2i*F2 - A1i*A2i*(F1 + 4*F4*q2i)) - 4*(A2i*B1i - A1i*B2i)*sqrt(P2i*q2i)*(F2 + F3*q2i)*zq + 2*(A2i*B1i*F3 + A1i*B2i*F3 + 2*A1i*A2i*F4)*P2i*q2i*zq^2
                )*forget_div
            # pseudo_vp_pqzq[i,j,m]=(
            #     4*F1*(B1i*B2i + A1i*A2i*q2i) + P2i*(A1i*A2i*F1) 
            #     )*forget_div
            f_gamma5Pqzq[i,j,m]=(
                4*F1*(B1i*B2i + A1i*A2i*q2i) + P2i*(2*A2i*B1i*F2 + 2*A1i*B2i*F2 - A1i*A2i*(F1 + 4*F4*q2i)) - 4*(A2i*B1i - A1i*B2i)*sqrt(P2i*q2i)*(F2 + F3*q2i)*zq + 2*(A2i*B1i*F3 + A1i*B2i*F3 + 2*A1i*A2i*F4)*P2i*q2i*zq^2
                )*forget_div
            f_gamma5muPqzq[i,j,m]=(
                2*A2i*B1i*F1 + 2*A1i*B2i*F1 - 4*B1i*B2i*F2 + A1i*A2i*F2*P2i + 4*(A1i*A2i*F2 + A2i*B1i*F4 + A1i*B2i*F4)*q2i + (4*(-(A2i*B1i) + A1i*B2i)*F1*q2i*zq)/sqrt(P2i*q2i) + A1i*A2i*F3*P2i*q2i*zq^2 - 4*q2i*(2*A1i*A2i*F2 + B1i*B2i*F3 + A2i*B1i*F4 + A1i*B2i*F4 + A1i*A2i*F3*q2i)*zq^2
                )*forget_div

            pseudo_vp_pq[i,j] += weightz[m]*pseudo_vp_pqzq[i,j,m]
            f_gamma5Pq[i,j] += weightz[m]*f_gamma5Pqzq[i,j,m]
            f_gamma5muPq[i,j] += weightz[m]* f_gamma5muPqzq[i,j,m]

            # f_gammaPq[i,j] += weightz[m]*(P2[i]*A2[i,j,m]*((F1k[i,j,m] + 4*F4k[i,j,m]*q2[j])*A1[i,j,m] - 2*F2k[i,j,m]*B1[i,j,m]) - 2*F2k[i,j,m]*P2[i]*A1[i,j,m]*B2[i,j,m] + 4*sqrt(P2[i]*q2[j])*(F2k[i,j,m] + F3k[i,j,m]*q2[j])*meshz[m]*(A2[i,j,m]*B1[i,j,m] - A1[i,j,m]*B2[i,j,m]) - 2*P2[i]*q2[j]*meshz[m]^2*(2*F4k[i,j,m]*A1[i,j,m]*A2[i,j,m] + F3k[i,j,m]*A2[i,j,m]*B1[i,j,m] + F3k[i,j,m]*A1[i,j,m]*B2[i,j,m]) - 4*F1k[i,j,m]*(q2[j]*A1[i,j,m]*A2[i,j,m] + B1[i,j,m]*B2[i,j,m]))
        end
    end
end

for i=1:Pstep
    for j=1:kstep
        pseudo_vp_p[i]+=weightk[j]*pseudo_vp_pq[i,j]*meshk[j]
        f_gamma5muP[i]+=weightk[j]*f_gamma5muPq[i,j]*meshk[j]
        f_gamma5muP[i]+=weightk[j]*f_gamma5muPq[i,j]*meshk[j]
    end
end
pseudo_vp_p = pseudo_vp_p *(2*pi)^-3*z4*3
f_gamma5P=f_gamma5P *(2*pi)^-3*z4/2*3
f_gamma5muP=f_gamma5muP *(2*pi)^-3*z2/2*3

print(pseudo_vp_p)