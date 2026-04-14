import pytest
from pyomo.environ import value, units as pyunits
import watertap.flowsheets.multistage_RO.multistage_RO as ro
import watertap.flowsheets.multistage_RO.utils as utils

lcow_results = {
    15: {
        1: {True: 0.35918616895308614, False: 0.39953201596014143},
        2: {True: 0.35585296050593535, False: 0.39545334229694845},
        3: {True: 0.3546436870590851, False: 0.39395906713287787},
    },
    20: {
        1: {True: 0.4137118661354158, False: 0.4672037381844089},
        2: {True: 0.40893472183682167, False: 0.46199031817129044},
        3: {True: 0.40704123711819895, False: 0.4599080338046219},
    },
    25: {
        1: {True: 0.46988686058001833, False: 0.5370253418725763},
        2: {True: 0.4627053541320692, False: 0.5302779753924234},
        3: {True: 0.4604311571208316, False: 0.5271988155949776},
    },
    30: {
        1: {True: 0.5277864238057798, False: 0.6090382027770745},
        2: {True: 0.5174849364328002, False: 0.5995756322894906},
        3: {True: 0.5149038293058242, False: 0.5960157066081466},
    },
    35: {
        1: {True: 0.5874835107244917, False: 0.683295808493933},
        2: {True: 0.5734098961102648, False: 0.6703883154448503},
        3: {True: 0.570554590714878, False: 0.6664482260120036},
    },
    40: {
        1: {True: 0.6490570480113436, False: 0.759870231339564},
        2: {True: 0.6305995763467285, False: 0.7428796932465613},
        3: {True: 0.6274863296198796, False: 0.7385992178391493},
    },
    45: {
        1: {True: 0.712593109685068, False: 0.8396191379911527},
        2: {True: 0.6891688665784267, False: 0.817835461490402},
        3: {True: 0.6858052859998074, False: 0.8129475491932157},
    },
    50: {
        1: {True: 0.7794856002583126, False: 0.925523129397549},
        2: {True: 0.7504576745371792, False: 0.8985526237344573},
        3: {True: 0.7465011163644455, False: 0.8926861130091104},
    },
}
sec_results = {
    15: {
        1: {True: 1.4180446083742657, False: 2.1732370966927297},
        2: {True: 1.401804923044571, False: 2.148033398732634},
        3: {True: 1.388910339403487, False: 2.1398330891356485},
    },
    20: {
        1: {True: 1.73161278158706, False: 2.676143580324469},
        2: {True: 1.6857474269427224, False: 2.6449911208729957},
        3: {True: 1.6773992646687006, False: 2.62314002865006},
    },
    25: {
        1: {True: 2.0558200783515703, False: 3.196418010828393},
        2: {True: 1.9779002762066002, False: 3.135613649171163},
        3: {True: 1.9707111642411612, False: 3.115800767800333},
    },
    30: {
        1: {True: 2.3902189597167856, False: 3.7327662200415332},
        2: {True: 2.276125025698808, False: 3.6389688560649662},
        3: {True: 2.268864992816272, False: 3.619353742548616},
    },
    35: {
        1: {True: 2.7345777099457305, False: 4.284484313277294},
        2: {True: 2.580207822473266, False: 4.154326654038163},
        3: {True: 2.5724652342584937, False: 4.1337234916542505},
    },
    40: {
        1: {True: 3.0888784824645197, False: 4.851361937131379},
        2: {True: 2.8905756830753813, False: 4.6815279519122805},
        3: {True: 2.882141344771226, False: 4.6595002378896995},
    },
    45: {
        1: {True: 3.4532722457086473, False: 5.508564371552928},
        2: {True: 3.2077579387893516, False: 5.29196391366613},
        3: {True: 3.198495367807448, False: 5.249774232212644},
    },
    50: {
        1: {True: 3.9031959770031985, False: 6.242193996873167},
        2: {True: 3.6081146831710393, False: 5.989304551989373},
        3: {True: 3.586079056107577, False: 5.940146385155745},
    },
}
salinity = [15, 20, 25, 30, 35, 40, 45, 50]
n_stages = [1, 2, 3]
add_erd = [True, False]

default_ro_op_dict = {
    "A_comp": 1.5 * pyunits.liter / pyunits.m**2 / pyunits.hour / pyunits.bar,
    "B_comp": 0.1 * pyunits.liter / pyunits.m**2 / pyunits.hour,
}


@pytest.mark.parametrize("n_stages", n_stages)
@pytest.mark.parametrize("add_erd", add_erd)
@pytest.mark.parametrize("salinity", salinity)
@pytest.mark.integration
def test_multistage_ro(salinity, n_stages, add_erd):
    m = ro.run_n_stage_system(
        n_stages=n_stages,
        salinity=salinity,
        water_recovery=0.5,
        pump_dict={1: True, 2: True, 3: False},
        ro_op_dict=default_ro_op_dict,
        add_erd=add_erd,
    )

    m = ro.set_system_recovery(m, 0.5)

    _ = utils.solve(model=m, tee=False)

    assert pytest.approx(lcow_results[salinity][n_stages][add_erd], rel=1e-3) == value(
        m.fs.costing.LCOW
    )
    assert pytest.approx(sec_results[salinity][n_stages][add_erd], rel=1e-3) == value(
        m.fs.costing.SEC
    )
