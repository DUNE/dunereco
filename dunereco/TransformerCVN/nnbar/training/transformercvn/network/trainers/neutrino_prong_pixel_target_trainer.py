import numpy as np
import torch
import torchmetrics as metrics
from torch import Tensor, jit
from torch.nn import functional as F

from transformercvn.dataset.sparse_prong_pixel_target_dataset import SparseProngPixelTargetDataset
from transformercvn.network.networks.neutrino_prong_pixel_target_network import (
    NeutrinoProngTargetNetwork,
)
from transformercvn.network.trainers.neutrino_base import NeutrinoBase
from transformercvn.options import Options

TArray = np.ndarray


class NeutrinoProngPixelTargetTrainer(NeutrinoBase):
    def __init__(self, options: Options):
        """

        Parameters
        ----------
        options: Options
            Global options for the entire network.
            See network.options.Options
        """
        super(NeutrinoProngPixelTargetTrainer, self).__init__(options)

        self.hidden_dim = options.hidden_dim
        self.network = jit.script(
            NeutrinoProngTargetNetwork(
                options,
                self.training_dataset.num_features,
                self.training_dataset.pixel_features,
                self.training_dataset.pixel_shape,
                self.training_dataset.num_prong_classes,
            )
        )

        self.accuracy = metrics.Accuracy()
        self.beta = 1 - 1 / len(self.training_dataset)

    @property
    def dataset(self):
        return SparseProngPixelTargetDataset

    def forward(self, features: Tensor, pixels: Tensor, extra: Tensor, mask: Tensor) -> Tensor:
        batch_size, timesteps, _ = features.shape

        # Normalize the high level layers
        features = features.clone()
        features[mask] -= self.mean
        features[mask] /= self.std

        #  Completely disable extra features
        # features *= 0

        # Normalize the event level layers
        extra = extra.clone()
        extra -= self.extra_mean
        extra /= self.extra_std

        # Convert the prong_pixels into a standardized floating point range
        pixels = pixels.float()
        pixels = pixels.reshape(
            batch_size,
            timesteps,
            self.training_dataset.pixel_features,
            *self.training_dataset.pixel_shape,
        )
        pixels /= 255

        if self.training:
            non_zero_pixels = pixels > 0
            noise = torch.randn_like(pixels[non_zero_pixels])
            noise = 1 + noise * self.options.pixel_noise_std
            pixels[non_zero_pixels] *= noise

        return self.network(features, pixels, extra, mask)

    def training_step(self, batch, batch_idx):
        data, pixels, extra, mask, targets, prong_targets = batch

        predictions = self.forward(data, pixels, extra, mask)

        prong_mask = prong_targets >= 0
        masked_predictions = torch.masked_select(predictions, prong_mask.unsqueeze(-1))
        masked_predictions = masked_predictions.reshape(-1, predictions.shape[-1])
        masked_targets = torch.masked_select(prong_targets.long(), prong_mask)

        loss = F.cross_entropy(masked_predictions, masked_targets)

        self.log("train_loss", loss)

        return loss

    def validation_step(self, batch, batch_idx):
        data, pixels, extra, mask, targets, prong_targets = batch

        predictions = self.forward(data, pixels, extra, mask)

        prong_mask = prong_targets >= 0
        masked_predictions = torch.masked_select(predictions, prong_mask.unsqueeze(-1))
        masked_predictions = masked_predictions.reshape(-1, predictions.shape[-1])
        masked_targets = torch.masked_select(prong_targets.long(), prong_mask)

        masked_predictions = torch.softmax(masked_predictions, dim=1)
        accuracy = self.accuracy(masked_predictions, masked_targets)
        self.log("val_accuracy", accuracy)
        self.log("val_epoch_accuracy", accuracy)
