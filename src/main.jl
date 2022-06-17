# 康嘉胤 2022-6-17
using TOML
workdir = "/Users/kjy/Desktop/program/julia/DSE_BSE"
# 放在一个函数里运行，与主题函数隔开
function innerplace(dir) 
    include(dir)
end
innerplace(joinpath(workdir, "src", "quark_equation_dse.jl"))