import numpy as np
from matplotlib import rcParams, cycler
import matplotlib.pyplot as plt
from read_likwid import LikwidReader


plt.clf()# clear , only need if you have plotted earlier in the code
fig, ax = plt.subplots()


x_var, y_var = 'n', 'Runtime [s]'


data = LikwidReader('ICACHE/prof_sm/')
readings = data.make()

x_var, y_var = 'n', 'Cache misses'
xss, yss, l = data.select(readings, x_var, y_var)

data = LikwidReader('basic/ICACHE/cache_basic_sm/')
readings = data.make()
xss2, yss2, l = data.select(readings, x_var, y_var)


xss = [ [xss[i][0], xss2[i][0]] for i in range(len(xss))]
yss = [ [yss[i][0], yss2[i][0]] for i in range(len(xss))]

#xss=[[]]
#yss=[[]]
#l=[[]]


all_ = [(xss[i],yss[i],l[i]) for i in range(len(l))]
all_.sort(key=lambda x:x[2])
#all_ = all_[8:]+all_[:8]
xss, yss, l = [ a[0] for a in all_ ], [ a[1] for a in all_ ], [ a[2] for a in all_ ]

#pf = np.polyfit([x[0] for x in xss], [y[0] for y in yss], 2)
#xp = np.linspace(min([x[0] for x in xss]),max([x[0] for x in xss]), 50)
#out = np.poly1d(pf)

cmap = plt.cm.viridis # - quite like this one too: viridis
rcParams['axes.prop_cycle'] = cycler(color=cmap(np.linspace(0, 1, len(xss))))
for i in range(len(xss)):
    ax.scatter(xss[i], yss[i], color=cmap(i*(1/len(xss))), marker='x', label=l[i])#, markersize = 12)
ax.set_yscale("log")
#ax.plot(xp, out(xp), 'r--')#, label='poly reg. with: '+str(num1)+'$x^{2}$ '+str(num2)+'x')


plt.annotate(
    'Basic',
    xy=(4000, 100000), xytext=(10, 5),
    textcoords='offset points', ha='left', va='bottom',
    bbox=dict(boxstyle='round,pad=0.2', fc='grey', alpha=0.5),
    arrowprops=dict(arrowstyle = '-', connectionstyle='arc3,rad=0'))
plt.annotate(
    'Optimised',
    xy=(6000, 6000), xytext=(10, 5),
    textcoords='offset points', ha='left', va='bottom',
    bbox=dict(boxstyle='round,pad=0.2', fc='grey', alpha=0.5),
    arrowprops=dict(arrowstyle = '-', connectionstyle='arc3,rad=0'))




box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

#ax.legend(loc='right upper', fancybox=True, shadow=True)
plt.xlabel(x_var)
plt.ylabel("log "+y_var)
plt.title('n vs cache misses (small matrices)')
plt.plot()
#plt.savefig('filename.png',dpi=500, bbox_inches = 'tight')# dpi is quality dots per inch