import numpy as np


param_values = []
x = np.linspace(10, 20, 2)
param_values.append(x)
x = np.linspace(30, 60, 2)
param_values.append(x)
dims = 40
for k in range(dims):
    param_values.append(np.array([k + 3]))

print(param_values)
mix_idx = []
mix_mesh_arrs = []
single_arrs = []
single_idx = []
for i, arr in enumerate(param_values):
    if arr.shape[0] == 1:
        single_arrs.append(arr)
        single_idx.append(i)
    else:
        mix_mesh_arrs.append(arr)
        mix_idx.append(i)
temp_global_combo_array = np.array(np.meshgrid(*mix_mesh_arrs, indexing="ij"))
print("temp_global_combo_array", temp_global_combo_array)
temp_global_combo_array = temp_global_combo_array.reshape(len(mix_mesh_arrs), -1).T
print("temp_global_combo_array.T", temp_global_combo_array)
global_combo_array = np.zeros((temp_global_combo_array.shape[0], len(param_values)))
print("global_combo_array", global_combo_array)
for i, g_i in enumerate(single_idx):
    global_combo_array[:, g_i] = single_arrs[i][0]
print("global_combo_array", global_combo_array)
for i, g_i in enumerate(mix_idx):
    global_combo_array[:, g_i] = temp_global_combo_array[:, i]
print("global_combo_array", global_combo_array)
for i in range(temp_global_combo_array.shape[0]):
    print(i, global_combo_array[i:,])
