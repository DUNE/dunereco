from typing import Tuple

import numpy as np
import torch
from MinkowskiEngine import SparseTensor
from torch import Tensor
from torch.nn import functional as F

from transformercvn.network.networks.neutrino_full_sparse_network import NeutrinoSparseNetwork
from transformercvn.network.trainers.neutrino_full_base_trainer import NeutrinoFullBaseTrainer
from transformercvn.options import Options

TArray = np.ndarray


class NeutrinoFullSparseTrainer(NeutrinoFullBaseTrainer):
    def create_network(
        self,
        options: Options,
        features_dim: int,
        extra_dim: int,
        pixel_dim: int,
        num_prong_classes: int,
        num_event_classes: int,
    ) -> NeutrinoSparseNetwork:
        return NeutrinoSparseNetwork(
            options, features_dim, extra_dim, pixel_dim, num_prong_classes, num_event_classes
        )

    def preprocess_pixels(
        self, pixel_coords: Tensor, pixel_values: Tensor, image_size: Tuple[int, ...]
    ) -> SparseTensor:
        if self.options.one_hot_pixels:
            pixel_type = pixel_values.dtype
            sparse_size, pixel_features = pixel_values.shape
            pixel_values = F.one_hot(pixel_values.long(), 256)
            pixel_values = pixel_values.reshape(sparse_size, 256 * pixel_features)
            pixel_values = pixel_values.to(pixel_type)

        else:
            if self.options.log_pixels:
                pixel_values = pixel_values + 1
                pixel_values = torch.log(pixel_values)

            else:
                pixel_values = pixel_values / 255.0

            if self.training:
                noise = torch.randn_like(pixel_values)
                noise = 1 + noise * self.options.pixel_noise_std
                pixel_values = pixel_values * noise

        return SparseTensor(pixel_values, pixel_coords)
