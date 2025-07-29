import torch
import random
import time
from datetime import datetime
import os
from .dataset import HDF5FileIterator
from glob import glob
import json
import gc
import psutil
import torch.multiprocessing
torch.multiprocessing.set_sharing_strategy('file_system')

time_now = datetime.now().strftime("%Y%m%d-%H%M%S")
model_fp = os.path.join('path', f'{time_now}' + ".pt")
trace_model_fp = os.path.join('path', f'trace_{time_now}' + ".pt")

train_hdf5_files = sorted(glob('../merged_datasets/training_dataset_*.h5'))
val_hdf5_files = sorted(glob('../merged_datasets/validation_dataset_*.h5'))
test_hdf5_files = sorted(glob('../merged_datasets/test_dataset_*.h5'))


def print_memory_usage(epoch, message):
    process = psutil.Process(os.getpid())
    memory_info = process.memory_info()
    print(f"Epoch {epoch + 1} - Memory Usage (CPU): {memory_info.rss / 1024**2:.2f} MB, {message}")

def train_model(model, device, num_epochs=5, use_cosine_schedule=True):
    print(" Training Process Starting...")

    if use_cosine_schedule:
        optimizer = torch.optim.Adam(model.parameters(), lr=1e-4, weight_decay=3e-4)
        lr_scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
            optimizer, T_max=num_epochs, eta_min=1e-6
        )
    else:
        optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=3e-4)
        lr_scheduler = None 
    best_val_loss, patience_counter = float('inf'), 0


    for epoch in range(num_epochs):
        model.train()
        total_loss = 0.0
        total_samples = 0

        random.shuffle(train_hdf5_files)
        train_loader = HDF5FileIterator(train_hdf5_files, dataset_type='train', batch_size=256*16)
        for file_idx, (train_dataloader, train_dataset) in enumerate(train_loader, start=1):
            print(f"  Training: Processing file 【{file_idx}/{len(train_hdf5_files)}】")
            
            hdf5_loss = 0.0
            hdf5_samples = 0
            start_time = time.time()
            print_memory_usage(epoch, 'before batch')

            for batch_idx, (all_init_conc, final_conc, init_conc, params, if_light, min_step) in enumerate(train_dataloader):
                all_init_conc = all_init_conc.to(device)
                init_conc = init_conc.to(device)
                params = params.to(device)
                final_conc = final_conc.to(device)
                if_light = if_light.to(device)
                min_step = min_step.to(device)

                optimizer.zero_grad()

                d_pred_norm, d_norm = model(init_conc, params, final_conc, all_init_conc, if_light, min_step)

                loss = model.loss(d_pred_norm, d_norm)  

                loss.backward()
                optimizer.step()

                batch_loss = loss.item()
                batch_samples = init_conc.size(0)
                total_loss += batch_loss * batch_samples
                total_samples += batch_samples

                hdf5_loss += batch_loss * batch_samples
                hdf5_samples += batch_samples

                if batch_idx % 500 == 0:
                    print(f"    Epoch [{epoch + 1: >3}/{num_epochs: >3}], "
                        f"Step [{batch_idx: >5}/{len(train_dataloader): >5}], "
                        f"Loss: {batch_loss: >7.4f}")
            print_memory_usage(epoch, 'after batch')

            end_time = time.time()
            elapsed_time = end_time - start_time

            if hdf5_samples > 0:
                hdf5_loss /= hdf5_samples
            else:
                hdf5_loss = float('inf')

            print(f"  Average Loss for HDF5 file '{train_dataset.hdf5_filename}': {hdf5_loss:.6f}")
            print(f"  Time taken for HDF5 file '{train_dataset.hdf5_filename}' training: 【{elapsed_time:.2f} seconds】")

            del train_dataloader, train_dataset
            gc.collect()
            torch.cuda.empty_cache()

        lr_scheduler.step()
        total_loss /= total_samples
        print(f"  Epoch [{epoch + 1}/{num_epochs}] Total Loss: {total_loss:.6f}")
        
        current_lr = optimizer.param_groups[0]['lr']
        print(f"  Epoch [{epoch + 1}/{num_epochs}] Learning Rate: {current_lr:.5e}")

        print(" Validating...")
        val_loss = validate(model, device)

        best_val_loss, patience_counter, early_stopping = early_stopping_helper(
            val_loss, best_val_loss, patience_counter, epoch, model, optimizer, lr_scheduler)

        if early_stopping:
            break
        del train_loader
        print_memory_usage(epoch, 'end epoch')
    # Test the model after training
    print(" Testing the model...")
    test_model(model, device)

def validate(model, device):
    model.eval()
    val_loss = 0.0
    total_samples = 0

    with torch.no_grad():
        val_loader = HDF5FileIterator(val_hdf5_files, dataset_type='val', batch_size=256*8)
        for val_dataloader, val_dataset in val_loader:
            for (all_init_conc, final_conc, init_conc, params, if_light, min_step) in val_dataloader:
                all_init_conc = all_init_conc.to(device)
                init_conc = init_conc.to(device)
                params = params.to(device)
                final_conc = final_conc.to(device)
                if_light = if_light.to(device)
                min_step = min_step.to(device)
                d_pred_norm, d_norm = model(init_conc, params, final_conc, all_init_conc, if_light, min_step)
                loss = model.loss(d_pred_norm, d_norm)  
                batch_loss = loss.item()
                batch_samples = init_conc.size(0)
                val_loss += batch_loss * batch_samples
                total_samples += batch_samples

    val_loss = val_loss / total_samples if total_samples > 0 else float('inf')
    print(f"  Validation Loss (per sample): {val_loss:.6f}")
    return val_loss


def test_model(model, device):
    model.eval()
    test_loss = 0.0
    total_samples = 0

    with torch.no_grad():
        test_loader = HDF5FileIterator(test_hdf5_files, dataset_type='test', batch_size=256*8)
        for test_dataloader, test_dataset in test_loader:
            for (all_init_conc, final_conc, init_conc, params, if_light, min_step) in test_dataloader:
                all_init_conc = all_init_conc.to(device)
                init_conc = init_conc.to(device)
                params = params.to(device)
                final_conc = final_conc.to(device)
                if_light = if_light.to(device)
                min_step =  min_step.to(device)
                d_pred_norm, d_norm = model(init_conc, params, final_conc, all_init_conc, if_light, min_step)
                loss = model.loss(d_pred_norm, d_norm)  
                batch_loss = loss.item()
                batch_samples = init_conc.size(0)
                test_loss += batch_loss * batch_samples
                total_samples += batch_samples

    test_loss = test_loss / total_samples if total_samples > 0 else float('inf')
    print(f"  Test Loss (per sample): {test_loss:.6f}")


def early_stopping_helper(val_loss, best_val_loss, patience_counter, epoch, model, optimizer, lr_scheduler):
    early_stopping_patience = 25
    if val_loss < best_val_loss:
        best_val_loss = val_loss
        patience_counter = 0
        torch.save({
            'epoch': epoch,
            'model_state_dict': model.state_dict(),
            'optimizer_state_dict': optimizer.state_dict(),
            'scheduler_state_dict': lr_scheduler.state_dict(),
            'best_val_loss': best_val_loss,
        }, model_fp)
        print(f"  Saved best model at epoch {epoch + 1}")
    else:
        patience_counter += 1
        if patience_counter >= early_stopping_patience:
            print(f"  Early stopping at epoch {epoch + 1}")
            return best_val_loss, patience_counter, True
    return best_val_loss, patience_counter, False
