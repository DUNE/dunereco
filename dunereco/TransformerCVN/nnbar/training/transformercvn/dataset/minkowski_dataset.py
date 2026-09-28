from typing import List, Optional, Tuple

import h5py
import numpy as np
import torch
from torch import Tensor
from torch.utils.data import Dataset


class MinkowskiCollection:
    @staticmethod
    def collate_sparse(
        coordinates: List[Tensor], values: List[Tensor], masks: List[Tensor]
    ) -> Tuple[Tensor, Tensor]:
        offset = torch.zeros(
            1, coordinates[0].shape[1], dtype=coordinates[0].dtype, device=coordinates[0].device
        )
        offset[0, 0] = 1

        total_prongs = 0
        collated_coordinates = []
        for coord, num_prongs in zip(coordinates, masks):
            collated_coordinates.append(coord + offset * total_prongs)
            total_prongs += num_prongs.sum().item()

        collated_coordinates = torch.cat(collated_coordinates)
        collated_values = torch.cat(values)

        return collated_coordinates, collated_values

    def __call__(self, list_data):
        (
            features,
            extra,
            event_coordinates,
            event_values,
            event_masks,
            prong_coordinates,
            prong_values,
            prong_masks,
            event_targets,
            prong_targets,
        ) = list(zip(*list_data))

        features = torch.stack(features)
        extra = torch.stack(extra)

        event_coordinates, event_values = self.collate_sparse(
            event_coordinates, event_values, event_masks
        )
        prong_coordinates, prong_values = self.collate_sparse(
            prong_coordinates, prong_values, prong_masks
        )

        event_masks = torch.stack(event_masks)
        prong_masks = torch.stack(prong_masks)

        event_targets = torch.stack(event_targets)
        prong_targets = torch.stack(prong_targets)

        return (
            features,
            extra,
            event_coordinates,
            event_values,
            event_masks,
            prong_coordinates,
            prong_values,
            prong_masks,
            event_targets,
            prong_targets,
        )


class MinkowskiDataset(Dataset):
    def __init__(
        self,
        data_file: str,
        limit_index=1.0,
        event_current_targets: bool = False,
        load_full_dataset: bool = False,
    ):
        super(MinkowskiDataset, self).__init__()

        self.load_full_dataset = load_full_dataset

        self.mean = 0
        self.std = 1

        self.extra_mean = 0
        self.extra_std = 1

        # Needed to make compatible with other datasets
        self.pixel_mean = 0
        self.pixel_std = 1
        self.pixels = None

        file = h5py.File(data_file, "r")
        self.num_events = file["features"].shape[0]

        # ---------------------------------------------------------
        # Find the Train / Validation limits for the given dataset.
        # ---------------------------------------------------------
        limit_index = self.compute_limit_index(limit_index)
        self.min_limit = limit_index.min()
        self.max_limit = limit_index.max()

        # ----------------------------------------------------
        # Load basic data from the HDF5 file and into PyTorch.
        # ----------------------------------------------------
        self.features = torch.from_numpy(file["features"][self.min_limit : self.max_limit])
        self.extra = torch.from_numpy(file["extra"][self.min_limit : self.max_limit])

        self.prong_mask = torch.from_numpy(file["prong_mask"][self.min_limit : self.max_limit])
        self.event_targets = torch.from_numpy(file["event_target"][self.min_limit : self.max_limit])
        self.prong_targets = torch.from_numpy(file["prong_target"][self.min_limit : self.max_limit])

        # 0,1 mu; 4,5, e; 8,9, tau; 13, nc; 1000180416, nnbar
        if event_current_targets:
            current_targets = np.zeros_like(self.event_targets)
            current_targets[(self.event_targets <= 1)] = 1
            current_targets[(self.event_targets <= 5) & (self.event_targets > 1)] = 2
            current_targets[self.event_targets > 20] = 3

            self.event_targets = torch.from_numpy(current_targets)

        # -----------------------------
        # Load sparse pixel data.
        # -----------------------------
        self.prong_compressed_index = file["prong_compressed_index"][
            self.min_limit : self.max_limit
        ]
        self.min_prong_index = self.prong_compressed_index[0, 0]
        self.max_prong_index = self.prong_compressed_index[-1, -1]
        self.prong_pixels_shape = file["prong_pixels_shape"][:]

        self.event_compressed_index = file["event_compressed_index"][
            self.min_limit : self.max_limit
        ]
        self.min_event_index = self.event_compressed_index[0, 0]
        self.max_event_index = self.event_compressed_index[-1, -1]
        self.event_pixels_shape = file["event_pixels_shape"][:]

        if load_full_dataset:
            self.prong_pixels_coordinates = torch.from_numpy(
                file["prong_pixels_coordinates"][self.min_prong_index : self.max_prong_index]
            )
            self.prong_pixels_values = torch.from_numpy(
                file["prong_pixels_values"][self.min_prong_index : self.max_prong_index]
            )

            self.event_pixels_coordinates = torch.from_numpy(
                file["event_pixels_coordinates"][self.min_event_index : self.max_event_index]
            )
            self.event_pixels_values = torch.from_numpy(
                file["event_pixels_values"][self.min_event_index : self.max_event_index]
            )

        else:
            self.prong_pixels_coordinates = np.memmap(
                data_file,
                mode="r",
                shape=file["prong_pixels_coordinates"].shape,
                offset=file["prong_pixels_coordinates"].id.get_offset(),
                dtype=file["prong_pixels_coordinates"].dtype,
            )
            self.prong_pixels_values = np.memmap(
                data_file,
                mode="r",
                shape=file["prong_pixels_values"].shape,
                offset=file["prong_pixels_values"].id.get_offset(),
                dtype=file["prong_pixels_values"].dtype,
            )
            self.event_pixels_coordinates = np.memmap(
                data_file,
                mode="r",
                shape=file["event_pixels_coordinates"].shape,
                offset=file["event_pixels_coordinates"].id.get_offset(),
                dtype=file["event_pixels_coordinates"].dtype,
            )
            self.event_pixels_values = np.memmap(
                data_file,
                mode="r",
                shape=file["event_pixels_values"].shape,
                offset=file["event_pixels_values"].id.get_offset(),
                dtype=file["event_pixels_values"].dtype,
            )

        self.full_pixel_shape = file["full_pixels_shape"][:]

        self.num_events, self.max_particles, self.num_features = self.features.shape
        self.num_extra = self.extra.shape[1]

        self.num_event_classes = self.event_targets.max().item() + 1
        self.num_prong_classes = self.prong_targets.max().item() + 1

        self.pixel_features, *self.pixel_shape = self.full_pixel_shape.tolist()
        self.pixel_shape = tuple(self.pixel_shape)

        self.prong_mask = self.prong_mask.bool()
        self.prong_mask[:, 0] = True
        self.event_mask = torch.ones(
            self.num_events, 1, dtype=self.prong_mask.dtype, device=self.prong_mask.device
        )

        self.event_compressed_index -= self.min_event_index
        self.prong_compressed_index -= self.min_prong_index

    def compute_limit_index(self, limit_index) -> np.ndarray:
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
    ) -> Tuple[Tensor, Tensor, Tensor, Tensor, Optional[Tensor], Optional[Tensor]]:
        if mean is None:
            masked_data = self.features[self.prong_mask]
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

    def __getitem__(self, item) -> Tuple[Tensor, ...]:
        event_lower, event_upper = self.event_compressed_index[item]
        prong_lower, prong_upper = self.prong_compressed_index[item]

        if self.load_full_dataset:
            event_coordinates = self.event_pixels_coordinates[event_lower:event_upper]
            event_values = self.event_pixels_values[event_lower:event_upper]

            prong_coordinates = self.prong_pixels_coordinates[prong_lower:prong_upper]
            prong_values = self.prong_pixels_values[prong_lower:prong_upper]

        else:
            event_lower += self.min_event_index
            event_upper += self.min_event_index
            prong_lower += self.min_prong_index
            prong_upper += self.min_prong_index

            event_coordinates = torch.from_numpy(
                np.array(self.event_pixels_coordinates[event_lower:event_upper])
            )
            event_values = torch.from_numpy(
                np.array(self.event_pixels_values[event_lower:event_upper])
            )

            prong_coordinates = torch.from_numpy(
                np.array(self.prong_pixels_coordinates[prong_lower:prong_upper])
            )
            prong_values = torch.from_numpy(
                np.array(self.prong_pixels_values[prong_lower:prong_upper])
            )

        return (
            self.features[item],
            self.extra[item],
            event_coordinates,
            event_values,
            self.event_mask[item],
            prong_coordinates,
            prong_values,
            self.prong_mask[item],
            self.event_targets[item],
            self.prong_targets[item],
        )
