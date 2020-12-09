import lsqlin
import numpy as np

C = np.array([[0.0372, 0.2869], [0.6861, 0.7071],[0.6233, 0.6245], [0.6344, 0.6170]])
d = np.array([0.8587, 0.1781, 0.0747, 0.8405])
ret = lsqlin.lsqnonneg(C, d, {'show_progress': False})
print ret['x'].T

C = np.array([[0.9501,0.7620,0.6153,0.4057],[0.2311,0.4564,0.7919,0.9354],[0.6068,0.0185,0.9218,0.9169],[0.4859,0.8214,0.7382,0.4102],[0.8912,0.4447,0.1762,0.8936]])

#A = np.array([[0.2027,0.2721,0.7467,0.4659],[0.1987,0.1988,0.4450,0.4186],[0.6037,0.0152,0.9318,0.8462]])
A = np.array([[0.1987,0.1988,0.4450,0.4186]])


d = np.array([0.0578, 0.3528, 0.8131, 0.0098, 0.1388])
b=np.array([0])
lb = np.array([-0.1] * 4)
ub = np.array([2] * 4)

#ret = lsqlin.lsqlin(C, d, 0, A, b, None, None,lb, ub, None, {'show_progress': False})
ret = lsqlin.lsqlin(C, d, 0, None, None, A, b, ub, None, {'show_progress': False})

print ret['x'].T

