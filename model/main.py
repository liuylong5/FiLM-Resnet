import torch
import torch.nn as nn
from glob import glob
import os
import pickle
from src.util import beijing_time, load_norm_params, load_model, init_weights
from src.model import DNN
from src.train import train_model
from warnings import simplefilter
simplefilter(action='ignore', category=FutureWarning)

def main():
    print(f"Training process starts at: [{beijing_time()}]")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f" Using device: [{device}]")
    norm_args = load_norm_params()
    
    model = DNN(
        input_size=56,  
        param_size=26, 
        output_size=66, 
        norm_args=norm_args  
    ).to(device)

    # for block in model.res_blocks[:4]:
    #     for param in block.parameters():
    #         param.requires_grad = False

    model_path = ""  
    if load_model(model, model_path, device):
        print("  Using existing model weights.")
    else:
        print("  Initializing model from scratch.")
        model.apply(init_weights)
    
    train_model(model, device, num_epochs=50)
    print(f"Training process ends at: [{beijing_time()}]")

if __name__ == "__main__":
    main()