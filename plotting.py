import matplotlib.pyplot as plt
# import numpy as np

f = open("output.txt", "r")
xyz = f.readline() # first line is for visuals only

values = {
    'x' : [],
    'y' : [],
    'z' : []
}

# declare positions
xi_pos = 0 
yi_pos = 1
zi_pos = 2

for line in f:
    curr_values = line.split(' ')

    # stop writing values if the function is too steep or y becomes negative
    if float(curr_values[zi_pos]) < -100 or \
       float(curr_values[yi_pos]) < 0:
        break

    values['x'].append(float(curr_values[xi_pos]))
    values['y'].append(float(curr_values[yi_pos]))
    values['z'].append(float(curr_values[zi_pos]))

print(values)

# displaying function
plt.plot(values['x'], values['y'], color = "green")

plt.xlabel('x - axis')
plt.ylabel('y - axis')

plt.show()
