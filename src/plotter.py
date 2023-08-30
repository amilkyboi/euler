# module plotter
'''
Plots the residual vector vs. the number of iterations.
'''

import matplotlib.pyplot as plt
import scienceplots # pylint: disable=unused-import
import pandas as pd

res_data = pd.read_csv('../data/res.csv')

plt.style.use(['science', 'grid'])

plt.plot(res_data['iter'], res_data['rho'],   label='$\\rho$')
plt.plot(res_data['iter'], res_data['rho_u'], label='$\\rho{}u$')
plt.plot(res_data['iter'], res_data['rho_v'], label='$\\rho{}v$')
plt.plot(res_data['iter'], res_data['rho_e'], label='$\\rho{}E$')

plt.plot(res_data['iter'], res_data['err_rho'],   label='$\\rho$ Error',    linestyle='dashdot')
plt.plot(res_data['iter'], res_data['err_rho_u'], label='$\\rho{}u$ Error', linestyle='dashdot')
plt.plot(res_data['iter'], res_data['err_rho_v'], label='$\\rho{}v$ Error', linestyle='dashdot')
plt.plot(res_data['iter'], res_data['err_rho_e'], label='$\\rho{}E$ Error', linestyle='dashdot')

plt.yscale('log')
plt.title('Residual vs. Iteration Number')
plt.xlabel('Iteration Number')
plt.ylabel('Residual')
plt.legend()

plt.show()
