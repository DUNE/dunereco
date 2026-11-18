from typing import Tuple

import numpy as np
import torch
from torch import Tensor, jit
from torch.nn import functional as F
from torchmetrics import Accuracy

from transformercvn.dataset.sparse_prong_pixel_target_dataset import SparseProngPixelTargetDataset
from transformercvn.network.networks.neutrino_combined_network import NeutrinoCombinedNetwork
from transformercvn.network.trainers.neutrino_base import NeutrinoBase
from transformercvn.options import Options

TArray = np.ndarray


class NeutrinoCombinedTrainer(NeutrinoBase):
    def __init__(self, options: Options):
        """

        Parameters
        ----------
        options: Options
            Global options for the entire network.
            See network.options.Options
        """
        super(NeutrinoCombinedTrainer, self).__init__(options)

        self.hidden_dim = options.hidden_dim

        self.network = jit.script(
            NeutrinoCombinedNetwork(
                options,
                self.training_dataset.num_features,
                self.training_dataset.pixel_features,
                self.training_dataset.pixel_shape,
                self.training_dataset.num_prong_classes,
                self.training_dataset.num_classes,
            )
        )

        self.prong_accuracy = Accuracy()
        self.event_accuracy = Accuracy()
        self.beta = 1 - 1 / len(self.training_dataset)

    @property
    def dataset(self):
        return SparseProngPixelTargetDataset

    def forward(
        self, features: Tensor, pixels: Tensor, extra: Tensor, mask: Tensor
    ) -> Tuple[Tensor, Tensor]:
        batch_size, timesteps, _ = features.shape

        # Normalize the high level layers
        features = features.clone()
        features[mask] -= self.mean
        features[mask] /= self.std

        # Normalize the event level layers
        extra = extra.clone()
        extra -= self.extra_mean
        extra /= self.extra_std

        pixels = pixels.reshape(
            batch_size,
            timesteps,
            self.training_dataset.pixel_features,
            *self.training_dataset.pixel_shape,
        )

        if not self.options.one_hot_pixels:
            # Convert the prong_pixels into a standardized floating point range
            # pixel_non_zero = pixels > 0
            pixels = pixels.float()

            if self.options.log_pixels:
                pixels += 1
                torch.log_(pixels)

            else:
                pixels /= 255

            if self.training:
                non_zero_pixels = pixels > 0
                noise = torch.randn_like(pixels[non_zero_pixels])
                noise = 1 + noise * self.options.pixel_noise_std
                pixels[non_zero_pixels] *= noise

        return self.network(features, pixels, extra, mask)

    def training_step(self, batch, batch_idx):
        data, pixels, extra, mask, targets, prong_targets = batch

        event_predictions, prong_predictions = self.forward(data, pixels, extra, mask)

        # Compute Event loss from simple prediction
        event_loss = F.cross_entropy(event_predictions, targets)

        # Compute Prong Loss after masking away padding prongs
        prong_mask = prong_targets >= 0
        masked_predictions = torch.masked_select(prong_predictions, prong_mask.unsqueeze(-1))
        masked_predictions = masked_predictions.reshape(-1, prong_predictions.shape[-1])
        masked_targets = torch.masked_select(prong_targets.long(), prong_mask)

        prong_loss = F.cross_entropy(masked_predictions, masked_targets)

        # Total Loss is just sum
        total_loss = event_loss + prong_loss

        self.log("prong_loss", prong_loss)
        self.log("event_loss", event_loss)
        self.log("train_loss", total_loss)

        return total_loss

    def validation_step(self, batch, batch_idx):
        data, pixels, extra, mask, targets, prong_targets = batch

        event_predictions, prong_predictions = self.forward(data, pixels, extra, mask)

        # Prong accuracy
        prong_mask = prong_targets >= 0
        masked_predictions = torch.masked_select(prong_predictions, prong_mask.unsqueeze(-1))
        masked_predictions = masked_predictions.reshape(-1, prong_predictions.shape[-1])
        masked_targets = torch.masked_select(prong_targets.long(), prong_mask)

        masked_predictions = torch.softmax(masked_predictions, dim=1)
        prong_accuracy = self.prong_accuracy(masked_predictions, masked_targets)

        # Event accuracy
        event_predictions = torch.softmax(event_predictions, dim=1)
        event_accuracy = self.event_accuracy(event_predictions, targets)

        self.log("prong_accuracy", prong_accuracy)
        self.log("event_accuracy", event_accuracy)
        self.log("val_accuracy", (prong_accuracy + event_accuracy) / 2)

    def validation_epoch_end(self, outputs) -> None:
        # =========================================================================================
        # Basic accuracy based metrics
        # =========================================================================================
        average_prong_accuracy = self.prong_accuracy.compute()
        average_event_accuracy = self.event_accuracy.compute()

        self.log("val_epoch_accuracy", (average_prong_accuracy + average_event_accuracy) / 2)
        self.log("event_epoch_accuracy", average_event_accuracy)
        self.log("prong_epoch_accuracy", average_prong_accuracy)
