# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# read data
df_intrinsic = pd.read_csv('./intrinsic_velocities.csv', '\t')
df_initial = pd.read_csv('./velocities_at_x.csv', '\t')


# do plotting
fig, axes = plt.subplots(2, 2, figsize = (12, 7))


# sorted by coefficients of friction \mu
for key in df_intrinsic.columns[1:]:
    data = df_intrinsic[key].as_matrix().reshape(19, 3)
    lambdas = data[::,0]
    angular = data[::,1]
    transversal = data[::,2]
    axes[0][0].plot(lambdas, angular, 'o', label = r'$\mu = %s$' % key)
    axes[0][1].plot(lambdas, transversal, 'o', label = r'$\mu = %s$' % key)


# sorted by initial angular velocities \dot{\phi}
for key in (r'\frac{2}{3}\pi', r'\pi', r'\frac{3}{2}\pi', r'2\pi',
            r'\frac{5}{2}\pi'):
    data = df_initial[key].as_matrix().reshape(21, 3)
    pos = data[::,0]
    angular = data[::,1]
    transversal = data[::,2]
    label_template = r'$\dot{\phi}_0 = %s \frac{rad}{sec}$'
    axes[1][0].plot(pos, angular, 'o', label = label_template % key)
    axes[1][1].plot(pos, transversal, 'o', label = label_template % key)


axes[0][0].set(title = r'Spezifische Winkelgeschwindigkeiten',
               xlabel = r'$\lambda$/m',
               ylabel = r'$\dot{\phi} \frac{rad}{sec}$')
axes[0][1].set(title = r'Spezifische Transversalgeschwindigkeiten',
               xlabel = r'$\lambda$/m',
               ylabel = r'$V/\frac{m}{sec}$')

axes[1][0].set(title = r'Entwicklung der Winkelgeschwindigkeit - $\mu=0.2, \lambda=6cm$',
               xlabel = r'x/m',
               ylabel = r'$\dot{\phi} / \frac{rad}{sec}$')
axes[1][1].set(title = r'Entwicklung der Transversalgeschwindigkeit - $\mu=0.2, \lambda=6cm$',
               xlabel = r'x/m',
               ylabel = r'$V/\frac{m}{sec}$')

# show legend
for ax in axes.ravel():
    ax.legend()

# save figure
plt.tight_layout()
plt.savefig('data.png', dpi=300)

