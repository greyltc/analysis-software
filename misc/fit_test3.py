import lmfit
import sympy
from scipy import constants


I_0, I_02, I_ph, R_s, R_sh, n, n_2, I, V, V_th = sympy.symbols(
    'I_0 I_02 I_ph R_s R_sh n n_2 I V V_th'
    )

lhs = I
rhs = I_ph-((V+I*R_s)/R_sh)-I_0*(sympy.exp((V+I*R_s)/(n*V_th))-1)

# one diode model
one_diode_model = sympy.Eq(lhs, rhs)

print('Solving solar cell equation...')
current = sympy.solveset(one_diode_model, I)
voltage = sympy.solveset(one_diode_model, V)
dark_current = sympy.solveset(one_diode_model, I_0)
photocurrent = sympy.solveset(one_diode_model, I_ph)
series_resistance = sympy.solveset(one_diode_model, R_s)
shunt_resistance = sympy.solveset(one_diode_model, R_sh)
ideality_factor = sympy.solveset(one_diode_model, n)
print('Solving completed!')

solutions = [current, voltage, dark_current, photocurrent, series_resistance,
             shunt_resistance, ideality_factor]

for solution in solutions:
    sympy.pprint(solution)


def thermal_voltage(temp_in_c):
    """returns a thermal voltage given a temperature in celsius"""
    return constants.k*(temp_in_c+constants.zero_Celsius)/constants.e


assumed_cell_temperature = 29.0  # degrees celsius
vee_tee_aitch = thermal_voltage(assumed_cell_temperature)

pass
