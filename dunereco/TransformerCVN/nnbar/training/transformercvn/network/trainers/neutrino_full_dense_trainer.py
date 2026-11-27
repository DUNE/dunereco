from typing import Tuple

import numpy as np
import torch
from torch import Tensor
from torch.nn import functional as F

from transformercvn.network.networks.neutrino_full_dense_network import NeutrinoDenseNetwork
from transformercvn.network.trainers.neutrino_full_base_trainer import NeutrinoFullBaseTrainer
from transformercvn.options import Options

TArray = np.ndarray


def sparse_to_dense(features: Tensor, coordinates: Tensor, image_size: Tuple[int, ...]):
    num_features = features.shape[1]

    coordinates = coordinates.T.long()
    batch_size = coordinates[0, -1].item() + 1

    output = torch.zeros(
        batch_size, *image_size, num_features, dtype=features.dtype, device=features.device
    )
    output[coordinates[0], coordinates[1], coordinates[2]] += features

    return output.permute(0, 3, 1, 2).contiguous()


class NeutrinoFullDenseTrainer(NeutrinoFullBaseTrainer):
    def create_network(
        self,
        options: Options,
        features_dim: int,
        extra_dim: int,
        pixel_dim: int,
        num_prong_classes: int,
        num_event_classes: int,
    ) -> NeutrinoDenseNetwork:
        return NeutrinoDenseNetwork(
            options, features_dim, extra_dim, pixel_dim, num_prong_classes, num_event_classes
        )

    def preprocess_pixels(
        self, pixel_coords: Tensor, pixel_values: Tensor, image_size: Tuple[int, ...]
    ) -> Tensor:
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

        return sparse_to_dense(pixel_values, pixel_coords, image_size)
