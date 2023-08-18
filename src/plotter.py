# module plotter
'''Plots residual.'''

import csv
import matplotlib.pyplot as plt
import scienceplots

x = []
y1 = []
y2 = []
y3 = []
y4 = []

with open('../data/res.csv', 'r', encoding='UTF-8') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    for row in lines:
        x.append(int(row[0]))
        y1.append(float(row[1]))
        y2.append(float(row[2]))
        y3.append(float(row[3]))
        y4.append(float(row[4]))

plt.style.use(['science', 'grid'])
plt.plot(x, y1, label='$\\rho$')
plt.plot(x, y2, label='$\\rho{}u$')
plt.plot(x, y3, label='$\\rho{}v$')
plt.plot(x, y4, label='$\\rho{}E$')
plt.yscale("log")
plt.title('Residual vs. Iteration Number')
plt.xlabel('Iteration Number')
plt.ylabel('Residual')
plt.legend()
plt.show()

# x = []
# y1 = []
# y2 = []

# with open('../data/force.csv', 'r', encoding='UTF-8') as csvfile:
#     lines = csv.reader(csvfile, delimiter=',')
#     for row in lines:
#         x.append(float(row[0]))
#         y1.append(float(row[1]))
#         y2.append(float(row[2]))

# plt.plot(x, y1, label='$F_x$')
# plt.plot(x, y2, label='$F_y$')
# plt.title('Force on Lower Bump')
# plt.xlabel('$x$-Direction Location')
# plt.ylabel('Force')
# plt.legend()
# plt.show()
