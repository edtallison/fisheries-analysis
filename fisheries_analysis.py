from functions import *
from matplotlib import pyplot as plt

# carrying capacity
k = 1000000
# initial population
n0 = 750000
# nominal quota
f0 = 130000
# birth rate
r = 0.5

# step size 
h = 1
# initial time (years)
start = 0
# final time (years)
end = 50

# set up plot
f, ax1 = plt.subplots(1,1)

# get 50 year time-span values for original quota management system, for plotting
t1, y1 = solve_explicit_rk(dndt_quota, start, end, n0, h, 'rk4', [r, k, f0])
line1, = plt.plot(t1,y1)

# fraction that cannot be exceeded
fr = 0.20
# get 50 year time-span values for kaitiakitanga management system, for plotting
t2, y2 = solve_explicit_rk(dndt_kaitiakitanga, start, end, n0, h, 'rk4', [r, k, f0, fr])
line2, = plt.plot(t2,y2)

# rahui length (years)
x = 3
# get 50 year time-span values for rahui management system, for plotting
t3, y3 = solve_explicit_rk(dndt_rahui, start, end, n0, h, 'rk4', [r, k, f0, x])
line3, = plt.plot(t3,y3)

# plot settings
hline = ax1.axhline(y=k/2, color = 'k', linestyle='-')
ax1.set_title('Fish Population under Different Management Scenarios over a 50 Year Time Period')
ax1.set_xlabel('Time (years)')
ax1.set_ylabel('Fish population (number of fish)')
line1.set_label('Original Quota System')
line2.set_label('Kaitiakitanga System (max fraction = 20%)')
line3.set_label('Rahui System (3 year period)')
hline.set_label('Healthy Fishery Threshold (half of capacity)')
ax1.legend()
plt.show()
