import pyomo.environ as pyo
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core.surrogate.pysmo_surrogate import PysmoSurrogate

# establish new surrogates
min_st_surrogate = PysmoSurrogate.load_from_file(
    "watertap/examples/flowsheets/pfas_treatment/min_st_pysmo_surr_spline.json",
)
throughput_surrogate = PysmoSurrogate.load_from_file(
    "watertap/examples/flowsheets/pfas_treatment/throughput_pysmo_surr_linear.json",
)

# set up solver
solver = get_solver()

# build models
m = pyo.ConcreteModel()

# min st
m.fs_min_st = FlowsheetBlock(dynamic=False)

m.fs_min_st.n_inv = pyo.Var(initialize=0.8)
m.fs_min_st.n_inv.fix()
m.fs_min_st.bi = pyo.Var(initialize=15)
m.fs_min_st.bi.fix()
m.fs_min_st.min_st = pyo.Var(initialize=1)

inputs = [m.fs_min_st.n_inv, m.fs_min_st.bi]
outputs = [m.fs_min_st.min_st]

m.fs_min_st.min_st_surrogate = SurrogateBlock(concrete=True)
m.fs_min_st.min_st_surrogate.build_model(
    min_st_surrogate,
    input_vars=inputs,
    output_vars=outputs,
)

# solve simulation
results = solver.solve(m, tee=False)

m.fs_min_st.display()

# throughput
m.fs_throughput = FlowsheetBlock(dynamic=False)

m.fs_throughput.n_inv = pyo.Var(initialize=0.8)
m.fs_throughput.n_inv.fix()
m.fs_throughput.bi = pyo.Var(initialize=15)
m.fs_throughput.bi.fix()
m.fs_throughput.conc_ratio = pyo.Var(initialize=0.5)
m.fs_throughput.conc_ratio.fix()
m.fs_throughput.throughput = pyo.Var(initialize=1)

inputs = [m.fs_throughput.n_inv, m.fs_throughput.bi, m.fs_throughput.conc_ratio]
outputs = [m.fs_throughput.throughput]

m.fs_throughput.throughput_surrogate = SurrogateBlock(concrete=True)
m.fs_throughput.throughput_surrogate.build_model(
    throughput_surrogate,
    input_vars=inputs,
    output_vars=outputs,
)

# solve simulation
results = solver.solve(m, tee=False)

m.fs_throughput.display()
