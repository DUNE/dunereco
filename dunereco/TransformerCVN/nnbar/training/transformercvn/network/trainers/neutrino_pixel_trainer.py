from typing import Tuple

import numpy as np
import torch
from matplotlib import pyplot as plt
from sklearn.metrics import ConfusionMatrixDisplay
from torch import Tensor, nn
from torchmetrics import Accuracy, ConfusionMatrix

from transformercvn.dataset.pixel_dataset import PixelNeutrinoDataset
from transformercvn.focal_loss import CB_loss
from transformercvn.network.layers.combiner import Combiner
from transformercvn.network.layers.decoder import Decoder
from transformercvn.network.layers.encoder import Encoder
from transformercvn.network.layers.pixel_encoder import PixelEncoder
from transformercvn.network.trainers.neutrino_base import NeutrinoBase
from transformercvn.options import Options

TArray = np.ndarray


class NeutrinoPixelTrainer(NeutrinoBase):
    def __init__(self, options: Options):
        """Base class defining the SPANet architecture.

        Parameters
        ----------
        options: Options
            Global options for the entire network.
            See network.options.Options
        """
        super(NeutrinoPixelTrainer, self).__init__(options)

        self.hidden_dim = options.hidden_dim

        # Shared options for all transformer layers
        transformer_options = (
            options.hidden_dim,
            options.num_attention_heads,
            options.hidden_dim,
            options.dropout,
            options.transformer_activation,
        )

        self.encoder = Encoder(options, self.training_dataset.num_features, transformer_options)
        self.combiner = Combiner(options)

        self.pixel_encoder = PixelEncoder(
            options,
            self.training_dataset.pixel_features,
            self.training_dataset.pixel_shape,
            self.options.resnet_block_depth,
        )

        self.current_decoder = Decoder(
            options, self.training_dataset.num_current_classes, hidden_dim_factor=2
        )
        self.generation_decoder = Decoder(
            options, self.training_dataset.num_generation_classes, hidden_dim_factor=2
        )

        self.current_loss = nn.CrossEntropyLoss()
        self.generation_loss = nn.CrossEntropyLoss(reduction="none")

        self.current_accuracy = Accuracy()
        self.current_confusion = ConfusionMatrix(
            self.training_dataset.num_current_classes, normalize="true"
        )
        self.current_confusion_pred = ConfusionMatrix(
            self.training_dataset.num_current_classes, normalize="pred"
        )

        self.generation_accuracy = Accuracy()
        self.generation_confusion = ConfusionMatrix(
            self.training_dataset.num_current_classes, normalize="true"
        )

        self.current_class_count = self.training_dataset.current_target_count
        self.generation_class_count = self.training_dataset.generation_target_count

        self.beta = self.options.loss_beta
        if self.options.loss_beta < 0.01:
            self.beta = 1 - 1 / len(self.training_dataset)

    @property
    def dataset(self):
        return PixelNeutrinoDataset

    def forward(
        self, data: Tensor, extra: Tensor, pixels: Tensor, mask: Tensor
    ) -> Tuple[Tensor, Tensor]:
        # Normalize incoming data
        # This operation is gradient-free, so we can use inplace operations.
        data = data.clone()
        data[mask] -= self.mean
        data[mask] /= self.std

        pixels = pixels.float().clone()
        pixels -= self.pixel_mean
        pixels /= self.pixel_std

        if self.training:
            noise = torch.randn_like(pixels)
            noise = 1 + noise * self.options.pixel_noise_std
            pixels = pixels * noise

        extra = extra - self.extra_mean
        extra = extra / self.extra_std

        # Extract layers from data using transformer
        hidden, padding_mask, sequence_mask = self.encoder(data, extra, mask)
        hidden = self.combiner(hidden, sequence_mask)

        pixels = self.pixel_encoder(pixels)
        pixels = pixels
        hidden = torch.cat([hidden, pixels], dim=-1)

        return self.current_decoder(hidden), self.generation_decoder(hidden)

    def training_step(self, batch, batch_idx):
        data, extra, pixels, mask, current_target, generation_target = batch

        # Only compute generation loss when there is a generation to predict
        # generation_valid = current_target < 2

        current_prediction, generation_prediction = self.forward(data, extra, pixels, mask)

        current_loss = CB_loss(
            current_target,
            current_prediction,
            self.current_class_count,
            self.training_dataset.num_current_classes,
            "focal",
            self.beta,
            self.options.loss_gamma,
        )

        generation_loss = CB_loss(
            generation_target,
            generation_prediction,
            self.generation_class_count,
            self.training_dataset.num_generation_classes,
            "focal",
            self.beta,
            self.options.loss_gamma,
        )

        total_loss = current_loss + self.options.event_prong_loss_proportion * generation_loss

        self.log("generation_loss", generation_loss)
        self.log("current_loss", current_loss)
        self.log("train_loss", total_loss)

        return total_loss

    def validation_step(self, batch, batch_idx):
        data, extra, pixels, mask, current_target, generation_target = batch
        generation_valid = current_target < 2

        current_prediction, generation_prediction = self.forward(data, extra, pixels, mask)

        generation_target = generation_target[generation_valid]
        generation_prediction = generation_prediction[generation_valid]

        current_prediction = torch.softmax(current_prediction, dim=-1)
        generation_prediction = torch.softmax(generation_prediction, dim=-1)

        # Mutate the current targets a bit to remove any background samples and boost the weight of neutral currents
        non_background = current_target != 3
        current_prediction = current_prediction[non_background]
        current_target = current_target[non_background]

        neutral_current = current_target == 2
        current_prediction = torch.cat([current_prediction, current_prediction[neutral_current]])
        current_target = torch.cat([current_target, current_target[neutral_current]])

        current_accuracy = self.current_accuracy(current_prediction, current_target)
        generation_accuracy = self.generation_accuracy(generation_prediction, generation_target)

        val_accuracy = (
            current_accuracy + self.options.event_prong_loss_proportion * generation_accuracy
        )
        val_accuracy = val_accuracy / (1 + self.options.event_prong_loss_proportion)

        metrics = {
            "current_accuracy": current_accuracy,
            "generation_accuracy": generation_accuracy,
            "val_accuracy": val_accuracy,
        }

        for key, val in metrics.items():
            self.log(key, val)

        # Add to the confusion matrix
        self.current_confusion(current_prediction, current_target)
        self.current_confusion_pred(current_prediction, current_target)
        self.generation_confusion(generation_prediction, generation_target)

        return metrics

    def validation_epoch_end(self, outputs) -> None:
        tensorboard = self.logger.experiment

        # =========================================================================================
        # Basic accuracy based metrics
        # =========================================================================================
        average_current_accuracy = self.current_accuracy.compute()
        average_generation_accuracy = self.generation_accuracy.compute()

        average_val_accuracy = (
            average_current_accuracy
            + self.options.event_prong_loss_proportion * average_generation_accuracy
        )
        average_val_accuracy = average_val_accuracy / (1 + self.options.event_prong_loss_proportion)

        self.log("val_epoch_accuracy", average_val_accuracy)
        self.log("current_epoch_accuracy", average_current_accuracy)
        self.log("generation_epoch_accuracy", average_generation_accuracy)

        # =========================================================================================
        # Target Normalized Current Confusion Matrix
        # =========================================================================================
        figure, ax = plt.subplots(figsize=(6, 6))

        confusion_matrix = self.current_confusion.compute().cpu().numpy()

        display = ConfusionMatrixDisplay(confusion_matrix * 100)
        display.plot(cmap="inferno", ax=ax, values_format=".3g")
        plt.title("Transformer Confusion Matrix \n Normalized by Truth", y=1.1)
        plt.tight_layout()

        tensorboard.add_figure("Current Confusion Matrix", figure, self.global_step)
        plt.close(figure)

        # =========================================================================================
        # Prediction Normalized Current Confusion Matrix
        # =========================================================================================
        figure, ax = plt.subplots(figsize=(6, 6))

        confusion_matrix = self.current_confusion_pred.compute().cpu().numpy()

        display = ConfusionMatrixDisplay(confusion_matrix * 100)
        display.plot(cmap="inferno", ax=ax, values_format=".3g")
        plt.title("Transformer Confusion Matrix \n Normalized by Pred", y=1.1)
        plt.tight_layout()

        tensorboard.add_figure("Current Confusion Matrix Pred", figure, self.global_step)
        plt.close(figure)

        # =========================================================================================
        # Target Normalized Generation Confusion Matrix
        # =========================================================================================
        # figure, ax = plt.subplots(figsize=(6, 6))
        #
        # confusion_matrix = self.generation_confusion.compute().cpu().numpy()
        #
        # display = ConfusionMatrixDisplay(confusion_matrix * 100)
        # display.plot(cmap='inferno', ax=ax, values_format=".3g")
        # plt.title("Transformer Confusion Matrix \n Normalized by Truth", y=1.1)
        # plt.tight_layout()
        #
        # tensorboard.add_figure("Current Confusion Matrix", figure, self.global_step)
        # plt.close(figure)

        # =========================================================================================
        # Reset all epoch-based metrics
        # =========================================================================================
        self.current_accuracy.reset()
        self.generation_accuracy.reset()
        self.current_confusion.reset()
        self.current_confusion_pred.reset()
        self.generation_confusion.reset()

        return average_val_accuracy
