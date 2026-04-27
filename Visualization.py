import numpy as np
import matplotlib.pyplot as plt

with open("data.txt") as file:
    N, M, M_vis = map(int, file.readline().split())
    x, y = [], []
    x_seg, y_seg = [], []
    k = 1
    for i in range(M):
        xi, yi = map(float, file.readline().split())
        x.append(xi)
        y.append(yi)
        if i == (N-1)*k:
            x_seg.append(xi)
            y_seg.append(yi)
            k += 1
    dot, func, lagr = [], [], []
    for i in range(M_vis):
        di, fi, li = map(float, file.readline().split())
        dot.append(di)
        func.append(fi)
        lagr.append(li)

plt.figure(figsize=(7, 7))
plt.plot(dot, func, 'b-', label='f(x)')
plt.plot(dot, lagr, 'r--', label='НП')
plt.scatter(x, y, c = 'g', marker = 'o')
plt.scatter(x_seg, y_seg, c = 'b', marker = 'o')
ax = plt.gca()  
ax.spines['left'].set_position('zero')   
ax.spines['bottom'].set_position('zero') 
ax.spines['right'].set_color('none')     
ax.spines['top'].set_color('none')

plt.grid(True)
plt.legend()
plt.show()
