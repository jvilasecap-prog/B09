import numpy as np

filename = 'ms313.dat'
with open(filename) as f:
    title = f.readline().strip()
    n_upper, n_lower = (int(float(v)) for v in f.readline().split())

# loadtxt skips the blank lines automatically
data = np.loadtxt(filename, skiprows=2)   # shape (90, 2)

upper = data[:n_upper]                    # x, y of the upper surface
lower = data[n_upper:n_upper + n_lower]   # x, y of the lower surface

x_upper, y_upper = upper[:, 0], upper[:, 1]
x_lower, y_lower = lower[:, 0], lower[:, 1]

print(title)
print(upper.shape, lower.shape)

c_LE = 0.2
c_TE = 0.8 

yLE= np.array([np.interp(c_LE, x_lower, y_lower), np.interp(c_LE, x_upper, y_upper)])
yTE= np.array([np.interp(c_TE, x_lower, y_lower), np.interp(c_TE, x_upper, y_upper)])

print(x_lower)
print(y_lower)

print(yLE, yTE)