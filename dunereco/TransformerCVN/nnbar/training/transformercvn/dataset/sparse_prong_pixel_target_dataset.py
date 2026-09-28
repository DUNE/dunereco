from typing import Tuple

import h5py
import torch
from torch import Tensor

from transformercvn.dataset.sparse_prong_pixel_dataset import SparseProngPixelDataset


class SparseProngPixelTargetDataset(SparseProngPixelDataset):
    def __init__(self, data_file: str, limit_index=1.0):
        super(SparseProngPixelTargetDataset, self).__init__(data_file, limit_index)

        with h5py.File(data_file, "r") as file:
            self.prong_targets = torch.from_numpy(
                file["prong_target"][self.min_limit : self.max_limit]
            )
            self.num_prong_classes = self.prong_targets.max().item() + 1

            self.prong_target_count = torch.bincount(self.prong_targets[self.prong_targets >= 0])
            self.prong_target_count = torch.clip(self.prong_target_count, 1, None)

    def __getitem__(self, item) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor, Tensor]:
        return (
            self.data[item],
            self.prong_pixels[item],
            self.extra[item],
            self.mask[item],
            self.targets[item],
            self.prong_targets[item],
        )
