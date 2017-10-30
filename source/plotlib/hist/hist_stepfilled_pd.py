import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
matplotlib.style.use('grayscale')

# ts = pd.Series(np.random.randn(1000), index=pd.date_range('1/1/2000', periods=1000))
# ts = ts.cumsum()
# ts.plot()

# df = pd.DataFrame(np.random.randn(1000, 4), index=ts.index, columns=list('ABCD'))
# df = df.cumsum()
# plt.figure(); df.plot();

# df3 = pd.DataFrame(np.random.randn(1000, 2), columns=['B', 'C']).cumsum()
# df3['A'] = pd.Series(list(range(len(df))))
# df3.plot(x='A', y='B')

# plt.figure()
# df.ix[5].plot(kind='bar', color='red'); plt.axhline(0, color='k')

# plt.figure()
# df2 = pd.DataFrame(np.random.rand(10, 4), columns=['a', 'b', 'c', 'd'])
# df2.plot(kind='bar')

# df4 = pd.DataFrame({'a': np.random.randn(1000) + 1, 'b': np.random.randn(1000), 'c': np.random.randn(1000) - 1}, columns=['a', 'b', 'c'])
# df4.plot(kind='hist', stacked=True, bins=20, alpha=0.8)

# # df = pd.DataFrame(np.random.rand(10, 4), columns=['a', 'b', 'c', 'd'])
# # df.plot(kind='bar', hatch='/', color='k');

# # df = pd.DataFrame(np.random.rand(10, 4), columns=['a', 'b', 'c', 'd'])
# ax = plt.figure(figsize=(10, 6)).add_subplot(111)
# # df.plot(ax=ax, kind='bar', legend=False)

# bars = ax.patches
# hatches = ''.join(h*len(df4) for h in 'x/O.')

# print hatches

# for bar, hatch in zip(bars, hatches):
#     bar.set_hatch(hatch)

# ax.legend(loc='center right', bbox_to_anchor=(1, 1), ncol=3)

df4 = pd.DataFrame({'a': np.random.randn(1000) + 1, 'b': np.random.randn(1000), 'c': np.random.randn(1000) - 1}, columns=['a', 'b', 'c'])
ax = plt.figure(figsize=(10, 6)).add_subplot(111)
df4.plot(ax=ax, kind='hist', stacked=True, bins=20, alpha=0.8)

bars = ax.patches
print len(bars)
print len(df4)
# hatches = ''.join(h*len(df4) for h in 'x/O.')
hatches = ''.join(h*len(df4) for h in 'o')
print hatches

for bar, hatch in zip(bars, hatches):
    bar.set_hatch(hatch)

# hatches = ('-', '+', 'x', '*', 'o', 'O', '.')
# for bar, pattern in zip(bars, hatches):
#      bar.set_hatch(pattern)

plt.title('Hist', fontsize=30)
plt.ylabel('$Ratio f = N_{HI}$/$N^*_{HI}$', fontsize=35)
plt.xlabel('$log_{10}(N^*_{HI}/10^{20} cm^{-2}$)', fontsize=35)
# plt.xlim(0.0, 2.0)
# plt.ylim(-1.0, 6.0)
plt.grid(True)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)
# plt.text(0.0, 3.2, 'adsaadas', color='blue', fontsize=20)
plt.legend(loc='upper right', fontsize=18)
# plt.savefig("test.png",bbox_inches='tight')
# for i in range(len(sc)):
# 	plt.annotate('('+str(sc[i])+')', xy=(xdata[i], ydata[i]), xycoords='data',
#             xytext=(-50.,30.), textcoords='offset points',
#             arrowprops=dict(arrowstyle="->"),fontsize=18,
#             )
plt.show()