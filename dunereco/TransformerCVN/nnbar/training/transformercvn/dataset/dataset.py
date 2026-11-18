from typing import Optional, Tuple

import h5py
import numpy as np
import torch
from torch import Tensor
from torch.utils.data import Dataset


class NeutrinoDataset(Dataset):
    def __init__(self, data_file: str, limit_index=1.0):
        super(NeutrinoDataset, self).__init__()

        self.mean = 0
        self.std = 1

        self.extra_mean = 0
        self.extra_std = 1

        self.pixel_mean = 0
        self.pixel_std = 1

        with h5py.File(data_file, "r") as file:
            self.num_events = file["data"].shape[0]

            # ---------------------------------------------------------
            # Find the Train / Validation limits for the given dataset.
            # ---------------------------------------------------------
            limit_index = self.compute_limit_index(limit_index)

            min_limit = limit_index.min()
            max_limit = limit_index.max()

            # ----------------------------------------------
            # Load data from the HDF5 file and into PyTorch.
            # ----------------------------------------------
            self.data = torch.from_numpy(file["data"][min_limit:max_limit])
            self.mask = torch.from_numpy(file["mask"][min_limit:max_limit])
            self.extra = torch.from_numpy(file["extra"][min_limit:max_limit])
            self.targets = torch.from_numpy(file["target"][min_limit:max_limit])

            self.pixels = None
            if "pixels" in file:
                self.pixels = torch.from_numpy(file["pixels"][min_limit:max_limit])

            self.data = self.data.permute(0, 2, 1)
            self.data = self.data.contiguous()

            self.mask = self.mask.bool()

            # ---------------------------------------------------
            # Remove any events that dont have any usable prongs.
            # ---------------------------------------------------
            good_events = self.mask.sum(1) > 0
            self.data = self.data[good_events].contiguous()
            self.mask = self.mask[good_events].contiguous()
            self.extra = self.extra[good_events].contiguous()
            self.targets = self.targets[good_events].contiguous()
            if "pixels" in file:
                self.pixels = self.pixels[good_events].contiguous()

        self.num_events, self.max_particles, self.num_features = self.data.shape
        self.num_classes = self.targets.max().item() + 1

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
        pixel_mean: Optional[Tensor] = None,
        pixel_std: Optional[Tensor] = None,
    ) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor, Tensor]:
        """Compute the mean and standard deviation of layers with normalization enabled in the event file.

        Parameters
        ----------
        mean: Tensor, optional
        std: Tensor, optional
            Give existing values for mean and standard deviation to set this value
            dataset's statistics to those values. This is especially useful for
            normalizing the validation and testing datasets with training statistics.

        Returns
        -------
        (Tensor, Tensor)
            The new mean and standard deviation for this dataset.
        """
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

        if self.pixels is not None:
            if pixel_mean is None:
                pixel_mean = torch.tensor([1.2237777, 0.9603817]).reshape(1, 2, 1, 1)
                pixel_std = torch.tensor([8.210588, 6.342488]).reshape(1, 2, 1, 1)

            self.pixel_mean = pixel_mean
            self.pixel_std = pixel_std

        return mean, std, extra_mean, extra_std, self.pixel_mean, self.pixel_std

    def __len__(self) -> int:
        return self.num_events

    def __getitem__(self, item) -> Tuple[Tensor, Tensor, Tensor, Tensor]:
        return self.data[item], self.extra[item], self.mask[item], self.targets[item]
