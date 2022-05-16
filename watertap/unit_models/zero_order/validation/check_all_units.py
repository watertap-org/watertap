import watertap.unit_models.zero_order as zo

from unit_validator import ZeroOrderUnitChecker

def run_and_print_result(name, file):
    try:
        checker = ZeroOrderUnitChecker(zero_order_model=name)
        error = checker.check_unit() 
        file.write(f"{name},{error}\n")
    except Exception as error:
        file.write(f"{name},{error}\n")

with open("test_results.csv", "w") as f:
    f.write("unit_name,result\n")
    for n in dir(zo):
        if isinstance(getattr(zo, n), type):
            run_and_print_result(n, f)
