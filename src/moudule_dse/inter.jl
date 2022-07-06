module interDSE
export ExtraA, ExtraB
using Interpolations
using MathLink
function ExtraA(k,A)
    # IntA=Splines1D(k,A)
    lengthofk=length(k)
    a=weval(W`FindFit[Adata, {a /Log[x/0.234^2]^b}, {a, b}, x][[1]][[2]]`; Adata=[k A][Int(lengthofk-floor(1/10*lengthofk)):(lengthofk-1),:])
    b=weval(W`FindFit[Adata, {a /Log[x/0.234^2]^b}, {a, b}, x][[2]][[2]]`; Adata=[k A][Int(lengthofk-floor(1/10*lengthofk)):(lengthofk-1),:])
    return a, b
end

function ExtraB(k,B)
    # IntA=Splines1D(k,A)
    lengthofk=length(k)
    a=weval(W`FindFit[Adata, {a /Log[x/0.234^2]^b}, {a, b}, x][[1]][[2]]`; Adata=[k B][Int(lengthofk-floor(1/10*lengthofk)):(lengthofk-1),:])
    b=weval(W`FindFit[Adata, {a /Log[x/0.234^2]^b}, {a, b}, x][[2]][[2]]`; Adata=[k B][Int(lengthofk-floor(1/10*lengthofk)):(lengthofk-1),:])
    return a, b
end

end # module