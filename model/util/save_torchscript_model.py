import torch
import torch.nn as nn
from glob import glob
import os
import sys
import pickle
from src.util import beijing_time
from src.model import DNN
from warnings import simplefilter
simplefilter(action='ignore', category=FutureWarning)


def convert_dict_to_tensors(norm_dict, device):
    tensor_dict = {}
    for key, value in norm_dict.items():
        if isinstance(value, torch.Tensor):
            tensor_dict[key] = value.clone().detach().to(dtype=torch.float32, device=device)
        else:
            tensor_dict[key] = torch.tensor(value, dtype=torch.float32, device=device)
    return tensor_dict

def load_norm(norm_fp='data/norm_values.pkl'):
    if not os.path.exists(norm_fp):
        raise FileNotFoundError(f"  Normalization parameters file not found: {norm_fp}")
    
    with open(norm_fp, 'rb') as f:
        norm_dict = pickle.load(f)  
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    norm_dict = convert_dict_to_tensors(norm_dict, device)  # 转换为 PyTorch 张量
    
    print(" Loaded normalization parameters successfully.")
    return norm_dict

def load_norm_params(norm_fp='data/norm_values.pkl'):
    norm_dict = load_norm(norm_fp)
    return tuple(norm_dict[key] for key in ["Xmu", "Xstd", "Dmu", "Dstd", "Pmin", "Pmax"])



device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: [{device}]")

norm_args = load_norm_params('../data/norm_values.pkl')

model = DNN(
    input_size=56,  
    param_size=26, 
    output_size=66, 
    norm_args=norm_args  
).to("cpu")

checkpoint = torch.load("../path/20250617-231942.pt", map_location="cpu")
model.load_state_dict(checkpoint["model_state_dict"], strict=False)
model.eval()

class InferenceWrapper(nn.Module):
    def __init__(self, model):
        super().__init__()
        self.model = model

    def forward(self, init_conc, params, if_light, min_step):
        delta_conc, final_conc = self.model(
            init_conc, params, if_light, min_step, return_delta=True
        )
        return delta_conc, final_conc

wrapped_model = InferenceWrapper(model)
wrapped_model.eval()

example_inputs = (
    torch.rand(8, 56),                        # init_conc
    torch.rand(8, 26),                        # params
    torch.ones(8, 1),                         # if_light
    torch.full((8, 1), 5.0)                   # min_step
)

traced = torch.jit.trace(wrapped_model, example_inputs)
traced.save("CONC_MODEL.pt")
print("Traced model saved.")
