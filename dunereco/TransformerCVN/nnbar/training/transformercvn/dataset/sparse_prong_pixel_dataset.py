from typing import Optional, Tuple

import h5py
import numpy as np
import torch
from numpy import ndarray as Array
from torch import Tensor
from torch.utils.data import Dataset


class CompressedCOOTensor:
    def __init__(self, indices: Array, values: Array, shape: Array):
        self.compressed_index = self.compress_first_index(indices, values, shape)

        self.indices = indices
        self.values = values
        self.shape = shape

    @staticmethod
    def compress_first_index(indices: Array, values: Array, shape: Array) -> Array:
        assert np.all(np.diff(indices[0]) >= 0)

        starts = np.searchsorted(indices[0], np.arange(shape[0]))
        ends = np.concatenate((starts[1:], [indices.shape[1]]))

        return np.stack([starts, ends]).T

    def limit_index(self, lower_limit: int, upper_limit: int):
        lower_mask = self.indices[0] >= lower_limit
        upper_mask = self.indices[0] < upper_limit
        complete_mask = lower_mask & upper_mask

        self.values = self.values[complete_mask]
        self.indices = np.ascontiguousarray(self.indices[:, complete_mask])
        self.indices[0] -= lower_limit

        self.shape = np.concatenate(([upper_limit - lower_limit], self.shape[1:]))
        self.compressed_index = self.compress_first_index(self.indices, self.values, self.shape)

    def dense(self):
        result = np.zeros(self.shape, dtype=np.uint8)
        result[self.indices[0], self.indices[1], self.indices[2]] = self.values
        return result

    def get_sparse_data(self, index):
        target_slice = slice(*self.compressed_index[index])
        new_indices = self.indices[:, target_slice][1:]
        new_values = self.values[target_slice]
        new_shape = self.shape[1:]

        return new_values.view(dtype=np.uint8), new_indices, new_shape

    def get_compressed(self, index, num_prongs):
        target_slice = slice(*self.compressed_index[index])
        new_indices = self.indices[:, target_slice][1:]
        new_values = self.values[target_slice]

        new_shape = self.shape[1:].copy()
        new_shape[0] = num_prongs

        result = np.zeros(new_shape, dtype=np.uint8)
        result[new_indices[0], new_indices[1]] = new_values

        return torch.from_numpy(result)

    def __getitem__(self, index):
        target_slice = slice(*self.compressed_index[index])
        new_indices = self.indices[:, target_slice][1:]
        new_values = self.values[target_slice]
        new_shape = self.shape[1:]

        result = np.zeros(new_shape, dtype=np.uint8)
        result[new_indices[0], new_indices[1]] = new_values

        return torch.from_numpy(result)


class SparseProngPixelDataset(Dataset):
    def __init__(self, data_file: str, limit_index=1.0):
        super(SparseProngPixelDataset, self).__init__()

        self.mean = 0
        self.std = 1

        self.extra_mean = 0
        self.extra_std = 1

        # Needed to make compatible with other datasets
        self.pixel_mean = 0
        self.pixel_std = 1
        self.pixels = None

        with h5py.File(data_file, "r") as file:
            self.num_events = file["data"].shape[0]

            # ---------------------------------------------------------
            # Find the Train / Validation limits for the given dataset.
            # ---------------------------------------------------------
            limit_index = self.compute_limit_index(limit_index)
            min_limit = limit_index.min()
            max_limit = limit_index.max()

            self.min_limit = min_limit
            self.max_limit = max_limit

            # ----------------------------------------------
            # Load data from the HDF5 file and into PyTorch.
            # ----------------------------------------------
            self.data = torch.from_numpy(file["data"][min_limit:max_limit])
            self.mask = torch.from_numpy(file["mask"][min_limit:max_limit])
            self.extra = torch.from_numpy(file["extra"][min_limit:max_limit])
            self.targets = torch.from_numpy(file["target"][min_limit:max_limit])

            prong_pixel_indices = file["prong_pixels_indices"][:]
            prong_pixel_values = file["prong_pixels_values"][:]
            prong_pixel_shape = file["prong_pixels_shape"][:]
            self.prong_pixels = CompressedCOOTensor(
                prong_pixel_indices, prong_pixel_values, prong_pixel_shape
            )
            self.prong_pixels.limit_index(min_limit, max_limit)

            self.mask = self.mask.bool()

            full_pixel_shape = None
            if "full_pixels_shape" in file:
                full_pixel_shape = file["full_pixels_shape"][:]

        self.num_events, self.max_particles, self.num_features = self.data.shape
        self.num_extra = self.extra.shape[1]
        self.num_classes = self.targets.max().item() + 1
        self.target_count = torch.bincount(self.targets)

        self.pixel_features = 2
        self.pixel_shape = int(np.sqrt(self.prong_pixels.shape[2] / 2))
        self.pixel_shape = (self.pixel_shape, self.pixel_shape)

        if full_pixel_shape is not None:
            self.pixel_features, *self.pixel_shape = full_pixel_shape.tolist()
            self.pixel_shape = tuple(self.pixel_shape)

    def compute_limit_index(self, limit_index):
        """Take subsection of the data for training / validation

        Parameters
        ----------
        limit_index : float in [-1, 1], tuple of floats, or array-like
            If a positive float - limit the dataset to the first limit_index percent of the data
            If a negative float - limit the dataset to the last |limit_index| percent of the data
            If a tuple - limit the dataset to [limit_index[0], limit_index[1]] percent of the data
            If array-like or tensor - limit the dataset to the specified indices.

        Returns
        -------
        np.ndarray or torch.Tensor
        """
        # In the float case, we just generate the list with the appropriate bounds
        if isinstance(limit_index, float):
            limit_index = (0.0, limit_index) if limit_index > 0 else (1.0 + limit_index, 1.0)

        # In the list / tuple case, we want a contiguous range
        if isinstance(limit_index, (list, tuple)):
            lower_index = int(round(limit_index[0] * self.num_events))
            upper_index = int(round(limit_index[1] * self.num_events))
            limit_index = np.arange(lower_index, upper_index)

        # Convert to numpy array for simplicity
        if isinstance(limit_index, Tensor):
            limit_index = limit_index.numpy()

        # Make sure the resulting index array is sorted for faster loading.
        return np.sort(limit_index)

    def compute_statistics(
        self,
        mean: Optional[Tensor] = None,
        std: Optional[Tensor] = None,
        extra_mean: Optional[Tensor] = None,
        extra_std: Optional[Tensor] = None,
    ) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor, Tensor]:
        if mean is None:
            masked_data = self.data[self.mask]
            mean = masked_data.mean(0)
            std = masked_data.std(0)

            std[std < 1e-5] = 1

        if extra_mean is None:
            extra_mean = self.extra.mean()
            extra_std = self.extra.std()

        self.mean = mean
        self.std = std

        self.extra_mean = extra_mean
        self.extra_std = extra_std

        return mean, std, extra_mean, extra_std, None, None

    def __len__(self) -> int:
        return self.num_events

    def __getitem__(self, item) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:
        return (
            self.data[item],
            self.prong_pixels[item],
            self.extra[item],
            self.mask[item],
            self.targets[item],
        )
