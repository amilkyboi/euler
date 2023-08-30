# module plotter
'''
Plots the residual vector vs. the number of iterations.
'''

import matplotlib.pyplot as plt
import scienceplots # pylint: disable=unused-import
import pandas as pd

res_data = pd.read_csv('../data/res.csv')
frc_data = pd.read_csv('../data/frc.csv')

plt.style.use(['science', 'grid'])
plt.figure(figsize=(1920/96, 1080/96), dpi=96)

plt.plot(res_data['iter'], res_data['rho'],   label='$\\rho$')
plt.plot(res_data['iter'], res_data['rho_u'], label='$\\rho{}u$')
plt.plot(res_data['iter'], res_data['rho_v'], label='$\\rho{}v$')
plt.plot(res_data['iter'], res_data['rho_e'], label='$\\rho{}E$')

plt.plot(res_data['iter'], res_data['err_rho'],   label='$\\rho$ Error',    linestyle='dashdot')
plt.plot(res_data['iter'], res_data['err_rho_u'], label='$\\rho{}u$ Error', linestyle='dashdot')
plt.plot(res_data['iter'], res_data['err_rho_v'], label='$\\rho{}v$ Error', linestyle='dashdot')
plt.plot(res_data['iter'], res_data['err_rho_e'], label='$\\rho{}E$ Error', linestyle='dashdot')

plt.yscale('log')
plt.title('$R_{ij} - D_{ij}$ \\& Respective Error vs. Iteration Number')
plt.xlabel('Iteration Number')
plt.ylabel('Residual')
plt.legend()

# plt.savefig('../img/res.webp', dpi=96 * 2)
plt.show()

plt.plot(frc_data['dist'], frc_data['fp_x'], label='$F_{px}$')
plt.plot(frc_data['dist'], frc_data['fp_y'], label='$F_{py}$')
plt.title('Pressure Force on the Bottom Wall')
plt.xlabel('$x$-distance')
plt.ylabel('Pressure Force')
plt.legend()

# plt.savefig('../img/frc.webp', dpi=96 * 2)
plt.show()
