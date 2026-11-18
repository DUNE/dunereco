from typing import Tuple

import numpy as np
import torchmetrics as metrics
from matplotlib import pyplot as plt
from torch import Tensor, nn

from transformercvn.dataset.split_dataset import SplitNeutrinoDataset
from transformercvn.focal_loss import CB_loss
from transformercvn.network.layers.combiner import Combiner
from transformercvn.network.layers.decoder import Decoder
from transformercvn.network.layers.encoder import Encoder
from transformercvn.network.trainers.neutrino_base import NeutrinoBase
from transformercvn.options import Options

TArray = np.ndarray


class NeutrinoSplitTrainer(NeutrinoBase):
    def __init__(self, options: Options):
        """Base class defining the SPANet architecture.

        Parameters
        ----------
        options: Options
            Global options for the entire network.
            See network.options.Options
        """
        super(NeutrinoSplitTrainer, self).__init__(options)

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

        self.current_decoder = Decoder(options, self.training_dataset.num_current_classes)
        self.generation_decoder = Decoder(options, self.training_dataset.num_generation_classes)

        self.current_loss = nn.CrossEntropyLoss()
        self.generation_loss = nn.CrossEntropyLoss(reduction="none")

        self.current_accuracy = metrics.Accuracy()
        self.current_confusion = metrics.ConfusionMatrix(
            self.training_dataset.num_current_classes, normalize="true"
        )

        self.generation_accuracy = metrics.Accuracy()
        self.generation_confusion = metrics.ConfusionMatrix(
            self.training_dataset.num_current_classes, normalize="true"
        )

        self.current_class_count = self.training_dataset.current_target_count
        self.generation_class_count = self.training_dataset.generation_target_count

    @property
    def dataset(self):
        return SplitNeutrinoDataset

    def forward(self, data: Tensor, extra: Tensor, mask: Tensor) -> Tuple[Tensor, Tensor]:
        # Normalize incoming data
        # This operation is gradient-free, so we can use inplace operations.
        data = data.clone()
        data[mask] -= self.mean
        data[mask] /= self.std

        extra = extra - self.extra_mean
        extra = extra / self.extra_std

        # Extract layers from data using transformer
        hidden, padding_mask, sequence_mask = self.encoder(data, extra, mask)
        hidden = self.combiner(hidden, sequence_mask)
        return self.current_decoder(hidden), self.generation_decoder(hidden)

    def training_step(self, batch, batch_idx):
        data, extra, mask, current_target, generation_target = batch

        # Only compute generation loss when there is a generation to predict
        generation_valid = current_target < 2

        current_prediction, generation_prediction = self.forward(data, extra, mask)

        generation_target = generation_target[generation_valid]
        generation_prediction = generation_prediction[generation_valid]

        current_loss = CB_loss(
            current_target,
            current_prediction,
            self.current_class_count,
            self.training_dataset.num_current_classes,
            "focal",
            self.options.loss_beta,
            self.options.loss_gamma,
        )

        generation_loss = CB_loss(
            generation_target,
            generation_prediction,
            self.generation_class_count,
            self.training_dataset.num_generation_classes,
            "focal",
            self.options.loss_beta,
            self.options.loss_gamma,
        )

        total_loss = current_loss + self.options.event_prong_loss_proportion * generation_loss

        self.log("generation_loss", generation_loss)
        self.log("current_loss", current_loss)
        self.log("train_loss", total_loss)

        return total_loss

    def validation_step(self, batch, batch_idx):
        data, extra, mask, current_target, generation_target = batch
        generation_valid = current_target < 2

        current_prediction, generation_prediction = self.forward(data, extra, mask)

        generation_target = generation_target[generation_valid]
        generation_prediction = generation_prediction[generation_valid]

        current_accuracy = self.current_accuracy(current_prediction, current_target)
        generation_accuracy = self.generation_accuracy(generation_prediction, generation_target)

        val_accuracy = (
            current_accuracy + self.options.event_prong_loss_proportion * current_accuracy
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
        self.generation_confusion(generation_prediction, generation_target)

        return metrics

    def validation_epoch_end(self, outputs) -> None:
        tensorboard = self.logger.experiment

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

        figure = plt.figure(figsize=(8, 6))

        confusion_matrix = self.current_confusion.compute().cpu().numpy()
        plt.imshow(confusion_matrix, vmin=0.0, vmax=1.0)
        plt.xlabel("Prediction")
        plt.ylabel("Target")
        plt.colorbar()

        tensorboard.add_figure("Current Confusion Matrix", figure, self.global_step)

        figure = plt.figure(figsize=(8, 6))

        confusion_matrix = self.generation_confusion.compute().cpu().numpy()
        plt.imshow(confusion_matrix, vmin=0.0, vmax=1.0)
        plt.xlabel("Prediction")
        plt.ylabel("Target")
        plt.colorbar()

        tensorboard.add_figure("Generation Confusion Matrix", figure, self.global_step)

        return average_val_accuracy
