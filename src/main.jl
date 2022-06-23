# 康嘉胤 2022-6-17
using TOML
workdir = "/Users/kjy/Desktop/program/julia/DSE_BSE"
# 放在一个模块儿里运行，与主题函数隔开
cd(workdir)
include(joinpath(pwd(),"src/quark_equation_dse.jl"))
# include(joinpath(pwd(),"src/meson.jl"))