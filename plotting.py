import matplotlib.pyplot as plt
import numpy as np

f = open("example.txt", "r")
xyz = f.readline()

values = {
    'x' : [],
    'y' : [],
    'z' : []
}

for line in f:
    curr_values = line.split(' ')

    values['x'].append(float(curr_values[0]))
    values['y'].append(float(curr_values[1]))
    values['z'].append(float(curr_values[2]))

# print(values)

# displaying my function
# plt.plot(values['x'], values['y'], color = "green", label = "my func")
plt.scatter(values['x'], values['y'], color = "green", marker = "*", s=30, label = "my func")
 
# displaying 'correct' function
x = np.arange(0, 84, 0.1)
# ans: y = e^(2 * x) * (-0.00341797*x^4 + 0.0625*x^3 - 0.25*x^2 + 0.75*x - 9.75)
# y = np.exp(2*x) * (-0.00341797*(x**4) + 0.0625*(x**3) - 0.25*(x**2) + 0.75*(x) - 0.75)
y = 643.788-277.162*np.log(np.power(np.tan(0.0297*x+1.89),2)+1)

print(x)
print(y)

plt.plot(x,y,color = "red", label = "wolf func")

plt.xlabel('x - axis')
plt.ylabel('y - axis')
plt.title('My func')

plt.legend()

plt.show()
