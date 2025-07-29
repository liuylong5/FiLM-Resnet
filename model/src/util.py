from datetime import datetime, timezone
import pytz
import torch
import os
import pickle
import torch.nn as nn

def beijing_time():
    utc_now = datetime.utcnow().replace(tzinfo=timezone.utc)
    beijing_tz = pytz.timezone('Asia/Shanghai')
    beijing_time = utc_now.astimezone(beijing_tz)
    formatted_time = beijing_time.strftime('%Y-%m-%d %H:%M:%S')
    
    return formatted_time

def init_weights(m):
    if isinstance(m, nn.Linear): 
        nn.init.xavier_uniform_(m.weight)
        if m.bias is not None:
            nn.init.zeros_(m.bias)


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


def load_model(model, model_path, device):
    if os.path.exists(model_path): 
        checkpoint = torch.load(model_path, map_location=device)
        
        if "model_state_dict" in checkpoint:
            model.load_state_dict(checkpoint["model_state_dict"]) 
            print(f"  Loaded  model from: {model_path}")
            return True  
        else:
            print("  Checkpoint does not contain 'model_state_dict'.")
            return False
    else:
        print(f"  Model not found: {model_path}")
        return False 

