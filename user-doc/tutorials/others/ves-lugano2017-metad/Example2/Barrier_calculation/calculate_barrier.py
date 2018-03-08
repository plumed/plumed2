import numpy as np

# Total number of fes files in folder
total_files=101
# Min and max initial guesses
min_min=50
min_max=90
max_min=90
max_max=130

for i in range(total_files):
        file_name="fes_" + str(i) + ".dat"
        matrix=np.genfromtxt(file_name)
        minimum=np.amin(matrix[min_min:min_max,1])
        maximum=np.amax(matrix[max_min:max_max,1])
        print(str(i) + " " + str(minimum) + " " + str(maximum) + " " + str(maximum-minimum))
