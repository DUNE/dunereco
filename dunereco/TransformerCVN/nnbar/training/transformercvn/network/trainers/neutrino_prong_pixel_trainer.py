import numpy as np
import torch
import torchmetrics as metrics
from matplotlib import pyplot as plt
from sklearn.metrics import ConfusionMatrixDisplay
from torch import Tensor, jit

from transformercvn.dataset.sparse_prong_pixel_dataset import SparseProngPixelDataset
from transformercvn.focal_loss import CB_loss
from transformercvn.network.networks.neutrino_prong_pixel_network import NeutrinoProngPixelNetwork
from transformercvn.network.trainers.neutrino_base import NeutrinoBase
from transformercvn.options import Options

TArray = np.ndarray


class NeutrinoProngPixelTrainer(NeutrinoBase):
    def __init__(self, options: Options):
        """

        Parameters
        ----------
        options: Options
            Global options for the entire network.
            See network.options.Options
        """
        super(NeutrinoProngPixelTrainer, self).__init__(options)

        self.hidden_dim = options.hidden_dim
        self.network = jit.script(
            NeutrinoProngPixelNetwork(
                options,
                self.training_dataset.num_features,
                self.training_dataset.pixel_features,
                self.training_dataset.num_classes,
            )
        )

        self.accuracy = metrics.Accuracy()
        self.confusion_target = metrics.ConfusionMatrix(
            self.training_dataset.num_classes, normalize="true"
        )
        self.confusion_prediction = metrics.ConfusionMatrix(
            self.training_dataset.num_classes, normalize="pred"
        )

        self.beta = self.options.loss_beta
        if self.options.loss_beta < 0.01:
            self.beta = 1 - 1 / len(self.training_dataset)

    @property
    def dataset(self):
        return SparseProngPixelDataset

    def forward(self, features: Tensor, pixels: Tensor, extra: Tensor, mask: Tensor) -> Tensor:
        batch_size, timesteps, _ = features.shape

        # Normalize the high level layers
        features = features.clone()
        features[mask] -= self.mean
        features[mask] /= self.std

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
        data, pixels, extra, mask, targets = batch

        predictions = self.forward(data, pixels, extra, mask)

        loss = CB_loss(
            targets,
            predictions,
            self.training_dataset.target_count,
            self.training_dataset.num_classes,
            "focal",
            self.beta,
            self.options.loss_gamma,
        )

        self.log("train_loss", loss)

        return loss

    def validation_step(self, batch, batch_idx):
        data, pixels, extra, mask, targets = batch

        predictions = self.forward(data, pixels, extra, mask)
        predictions = torch.softmax(predictions, dim=-1)

        neutral_current_mask = targets == 2
        accuracy_predictions = torch.cat([predictions, predictions[neutral_current_mask]])
        accuracy_targets = torch.cat([targets, targets[neutral_current_mask]])

        neutral_current_mask = accuracy_targets == 1
        accuracy_predictions = torch.cat(
            [accuracy_predictions, accuracy_predictions[neutral_current_mask]]
        )
        accuracy_targets = torch.cat([accuracy_targets, accuracy_targets[neutral_current_mask]])

        # Mutate the current targets a bit to remove any background samples and boost the weight of neutral currents
        non_background_mask = accuracy_targets != 3
        accuracy_predictions = accuracy_predictions[non_background_mask]
        accuracy_targets = accuracy_targets[non_background_mask]

        accuracy = self.accuracy(accuracy_predictions, accuracy_targets)
        self.log("val_accuracy", accuracy)

        # Add to the confusion matrix
        self.confusion_target(predictions, targets)
        self.confusion_prediction(predictions, targets)

        return metrics

    def validation_epoch_end(self, outputs) -> None:
        tensorboard = self.logger.experiment

        # =========================================================================================
        # Basic accuracy based metrics
        # =========================================================================================
        average_accuracy = self.accuracy.compute()

        self.log("val_epoch_accuracy", average_accuracy)

        # =========================================================================================
        # Target Normalized Current Confusion Matrix
        # =========================================================================================
        figure, ax = plt.subplots(figsize=(9, 9))

        confusion_matrix = self.confusion_target.compute().cpu().numpy()

        display = ConfusionMatrixDisplay(confusion_matrix * 100)
        display.plot(cmap="inferno", ax=ax, values_format=".3g")
        plt.title("Transformer Confusion Matrix \n Normalized by Truth", y=1.1)
        plt.tight_layout()

        tensorboard.add_figure("Current Confusion Matrix", figure, self.global_step)
        plt.close(figure)

        # =========================================================================================
        # Prediction Normalized Current Confusion Matrix
        # =========================================================================================
        figure, ax = plt.subplots(figsize=(9, 9))

        confusion_matrix = self.confusion_prediction.compute().cpu().numpy()

        display = ConfusionMatrixDisplay(confusion_matrix * 100)
        display.plot(cmap="inferno", ax=ax, values_format=".3g")
        plt.title("Transformer Confusion Matrix \n Normalized by Pred", y=1.1)
        plt.tight_layout()

        tensorboard.add_figure("Current Confusion Matrix Pred", figure, self.global_step)
        plt.close(figure)

        # =========================================================================================
        # Reset all epoch-based metrics
        # =========================================================================================
        self.accuracy.reset()
        self.confusion_target.reset()
        self.confusion_prediction.reset()

        return average_accuracy
