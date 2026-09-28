import numpy as np
import torchmetrics as metrics
from matplotlib import pyplot as plt
from torch import Tensor, nn

from transformercvn.network.layers.combiner import Combiner
from transformercvn.network.layers.decoder import Decoder
from transformercvn.network.layers.encoder import Encoder
from transformercvn.network.trainers.neutrino_base import NeutrinoBase
from transformercvn.options import Options

TArray = np.ndarray


class NeutrinoTrainer(NeutrinoBase):
    def __init__(self, options: Options):
        """Base class defining the SPANet architecture.

        Parameters
        ----------
        options: Options
            Global options for the entire network.
            See network.options.Options
        """
        super(NeutrinoTrainer, self).__init__(options)

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
        self.decoder = Decoder(options, self.training_dataset.num_classes)
        self.combiner = Combiner(options)

        self.loss = nn.CrossEntropyLoss()

        self.metric_acc = metrics.Accuracy()
        self.metric_pre = metrics.Precision(self.training_dataset.num_classes)
        self.metric_rec = metrics.Recall(self.training_dataset.num_classes)
        self.metric_f1s = metrics.F1(self.training_dataset.num_classes)

        self.metric_confusion = metrics.ConfusionMatrix(
            self.training_dataset.num_classes, normalize="true"
        )
        # An example images for generating the network's graph, batch size of 2
        # self.example_input_array = tuple(images.contiguous() for images in self.training_dataset[:2][0])

    def forward(self, data: Tensor, extra: Tensor, mask: Tensor) -> Tensor:
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
        return self.decoder(hidden)

    def training_step(self, batch, batch_idx):
        data, extra, mask, targets = batch

        predictions = self.forward(data, extra, mask)
        loss = self.loss(predictions, targets)

        self.log("train_loss", loss)
        return loss

    def validation_step(self, batch, batch_idx):
        data, extra, mask, targets = batch

        predictions = self.forward(data, extra, mask)

        metrics = {
            "val_accuracy": self.metric_acc(predictions, targets),
            "val_precision": self.metric_pre(predictions, targets),
            "val_recall": self.metric_rec(predictions, targets),
            "val_f1_score": self.metric_f1s(predictions, targets),
        }

        for key, val in metrics.items():
            self.log(key, val)

        # Add to the confusion matrix
        self.metric_confusion(predictions, targets)

        return metrics

    def validation_epoch_end(self, outputs) -> None:
        tensorboard = self.logger.experiment
        average_accuracy = self.metric_acc.compute()

        self.log("val_epoch_accuracy", average_accuracy)

        figure = plt.figure(figsize=(8, 6))

        confusion_matrix = self.metric_confusion.compute().cpu().numpy()
        plt.imshow(confusion_matrix, vmin=0.0, vmax=1.0)
        plt.xlabel("Prediction")
        plt.ylabel("Target")
        plt.colorbar()

        tensorboard.add_figure("Confusion Matrix", figure, self.global_step)

        return average_accuracy
