import h5py
import torch
from torch.utils.data import Dataset, DataLoader
import os
import json
import pickle
import time

class HDF5Dataset(Dataset):
    def __init__(self, file_path, dataset_type='train'):
        """
        Args:
            file_path (str): Path to the HDF5 file.
            dataset_type (str): 'train', 'val', or 'test'.
        """
        self.file_path = file_path
        self.dataset_type = dataset_type
        self.hdf5_filename = os.path.basename(self.file_path)

        no_change_idx = [48, 50, 53, 56, 57, 58, 59, 62, 63, 65]
        all_indices = set(range(66))
        self.change_indices = sorted(list(all_indices - set(no_change_idx)))

        self.data = self._load_file()

    def _load_file(self):

        if not os.path.isfile(self.file_path):
            raise FileNotFoundError(f"HDF5 file not found: {self.file_path}")

        start_time = time.time()
        with h5py.File(self.file_path, 'r', libver='latest', swmr=True) as f:
            if 'data' in f:
                data = f['data'][:] 
            else:
                raise KeyError(f"No 'data' dataset found in file {self.hdf5_filename}")
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"  Loaded 【{len(data)}】 samples from 【{self.hdf5_filename}】 takes 【{elapsed_time:.2f}】 seconds.")
        return data

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        """
        Fetch a single data sample.
        """
        sample = self.data[idx]

        init_conc = torch.tensor(sample[:56], dtype=torch.float32)
        final_conc = torch.tensor(sample[56:122], dtype=torch.float32)
        params = torch.tensor(sample[122:148], dtype=torch.float32)
        if_light = torch.tensor(sample[148:149], dtype=torch.float32)
        min_step = torch.tensor(sample[149:150], dtype=torch.float32)



        all_init_conc = torch.full((66,), fill_value=1e-19, dtype=torch.float32)
        all_init_conc[self.change_indices] = init_conc 


        return all_init_conc, final_conc, init_conc, params, if_light, min_step


class HDF5FileIterator:
    """
    Iterator to handle multiple HDF5 files for training, validation, or testing.
    """
    def __init__(self, filepaths, dataset_type, batch_size=256):
        """
        Args:
            filepaths (list): List of HDF5 file paths.
            dataset_type (str): 'train', 'val', or 'test'.
            batch_size (int): Batch size for the DataLoader.
        """
        self.filepaths = [fp for fp in filepaths if os.path.isfile(fp)]
        self.dataset_type = dataset_type
        self.batch_size = batch_size
        self.current_file_index = 0

        if not self.filepaths:
            raise ValueError("No valid HDF5 files provided!")

    def __iter__(self):
        return self

    def __next__(self):
        if self.current_file_index >= len(self.filepaths):
            raise StopIteration

        filepath = self.filepaths[self.current_file_index]
        self.current_file_index += 1

        dataset = HDF5Dataset(filepath, dataset_type=self.dataset_type)

        dataloader = DataLoader(
            dataset,
            batch_size=self.batch_size,
            shuffle=(self.dataset_type == 'train'),
            num_workers=16,
            pin_memory=True
        )

        return dataloader, dataset
