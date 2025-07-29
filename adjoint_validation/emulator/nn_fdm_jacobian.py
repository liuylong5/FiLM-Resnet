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

model.eval()

sample_idx = 0 
tensor_C_i = tensor_C[sample_idx].unsqueeze(0)  
tensor_P_i = tensor_P[sample_idx].unsqueeze(0) 
tensor_I_i = tensor_I[sample_idx].unsqueeze(0)  
print(f"{tensor_C_i.shape = }") 
print(f"{tensor_P_i.shape = }")  
print(f"{tensor_I_i.shape = }")  
print(f"{tensor_I_i = }") 
print(tensor_C_i)

with torch.no_grad():
    base_output = model(tensor_C_i, tensor_P_i, tensor_I_i[:,0].unsqueeze(0), tensor_I_i[:,1].unsqueeze(0) , return_delta=True)
    base_output = base_output[1].squeeze(0).cpu().numpy() 
print(base_output)

num_features = 56
perturbation = 0.1

J_all = [] 

for sample_idx in range(tensor_C.shape[0]):  
    print(f'processing {sample_idx = }')
    tensor_C_i = tensor_C[sample_idx].unsqueeze(0)  
    tensor_P_i = tensor_P[sample_idx].unsqueeze(0) 
    tensor_I_i = tensor_I[sample_idx].unsqueeze(0) 
    with torch.no_grad():
        base_output = model(tensor_C_i, tensor_P_i, tensor_I_i[:,0].unsqueeze(0), tensor_I_i[:,1].unsqueeze(0), return_delta=True)
        base_output = base_output[1].squeeze(0).cpu().numpy() 

    jacobian_matrix = np.zeros((len(base_output), num_features))  

    for i in range(num_features):
        conc_perturbed = tensor_C_i.clone() 
        conc_perturbed[0, i] = tensor_C_i[0, i] * (1 + perturbation) 
        with torch.no_grad():
            perturbed_output = model(conc_perturbed, tensor_P_i, tensor_I_i[:,0].unsqueeze(0), tensor_I_i[:,1].unsqueeze(0), return_delta=True)
            perturbed_output = perturbed_output[1].squeeze(0).cpu().numpy()

        jacobian_matrix[:, i] = (perturbed_output - base_output) / (tensor_C_i[0, i].cpu().numpy() * perturbation)  

    J_all.append(jacobian_matrix) 

J_all = np.stack(J_all)  
remain_list = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 49, 51, 52, 54, 55, 60, 61, 64]
J_all = J_all[:, remain_list, :]
nn_fdm_jacobian_path = os.path.join('../nn_jacobian', 'nn_fwd_jacobian.npy')
np.save(nn_fdm_jacobian_path, J_all) 
print(J_all.shape)