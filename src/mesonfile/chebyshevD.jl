timetest=time()
#  切比雪夫展开
chebyshev_D_step=3
chebyshevD=Array{Float64}(undef,kstep,kstep,chebyshev_D_step+1);
Threads.@threads for i=1:chebyshev_D_step
    for kn=1:kstep
        for  qn=1:kstep
            k2=meshk[kn]
            q2=meshk[qn]
            # @expand_D()
            chebyshevD[kn,qn,i]=gausschebyshevint64(z->D(k2+q2-2*sqrt(k2*q2)*z)chebyshevU(i-1,z))*2/pi
        end
    end
end

function D_chebyshev(kn,qn,x)
    result=0.::Float64
    for i=1:chebyshev_D_step
        result+=chebyshevD[kn,qn,i]*chebyshevU(i-1,x)
    end
    return result
end

# print("完成切比雪夫展开，用时",round((time()-timetest)*100)/100,"s \n")
