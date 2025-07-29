import os
import torch
import numpy as np
from src.util import beijing_time, load_norm_params, load_model, init_weights
from src.model import DNN
from warnings import simplefilter
simplefilter(action='ignore', category=FutureWarning)

print(f"Process starts at: [{beijing_time()}]")

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f" Using device: [{device}]")


inputs_npy = os.path.join('../nn_jacobian', 'inputs.npy')
params_npy = os.path.join('../nn_jacobian', 'params.npy')
info_npy = os.path.join('../nn_jacobian', 'info.npy')
conc = np.load(inputs_npy)
p = np.load(params_npy)   
i = np.load(info_npy)

remain_list = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 49, 51, 52, 54, 55, 66, 67, 70]
tensor_C = torch.tensor(conc[:, remain_list], dtype=torch.float32, device=device)
tensor_P = torch.tensor(p, dtype=torch.float32, device=device)
tensor_I = torch.tensor(i, dtype=torch.float32, device=device)

print(f"{tensor_C.shape = }") 
print(f"{tensor_P.shape = }")  
print(f"{tensor_I.shape = }")  

norm_path = os.path.join('../model', 'data/norm_values.pkl')
norm_args = load_norm_params(norm_fp=norm_path)

model = DNN(
    input_size=56,  
    param_size=26, 
    output_size=66, 
    norm_args=norm_args  
).to(device)

model_path = '../path/20250617-231942.pt'
if load_model(model, model_path, device):
    print("  Found model path.")
else:
    print("  No model found.")

J_all = [] 
remain_list = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 49, 51, 52, 54, 55, 60, 61, 64]
for i in range(tensor_C.shape[0]):   
    tensor_C_i = tensor_C[i].unsqueeze(0).requires_grad_(True) 
    tensor_P_i = tensor_P[i].unsqueeze(0).detach()  
    tensor_I_i = tensor_I[i].unsqueeze(0).detach()  



    def func_reduced(x):
        full_output = model(x, tensor_P_i, tensor_I_i[:,0].unsqueeze(0), tensor_I_i[:,1].unsqueeze(0), return_delta=True)[1]
        selected_output = full_output[:, remain_list] 
        return selected_output

    J_i = torch.autograd.functional.jacobian(func_reduced, tensor_C_i) 
    J_i = J_i.squeeze().detach().cpu().numpy() 
    J_all.append(J_i) 

    print(f"Processed sample {i+1}/{tensor_C.shape[0]}, J_i shape: {J_i.shape}")

J_all = np.stack(J_all)  
nn_jacobian_path = os.path.join('../nn_jacobian', 'nn_jacobian.npy')
np.save(nn_jacobian_path, J_all) 


print(f"Process ends at: [{beijing_time()}]")