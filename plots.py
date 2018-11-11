import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

nr_objects = [1740, 2377, 3127, 4277, 6830, 9386]
rfc_performance = [94, 93, 92, 91, 89, 87]
plt.plot(nr_objects, rfc_performance, 'g')
plt.show()
