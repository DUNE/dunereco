from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, List, Tuple, Union

import numpy as np
import torch
import torch.distributed as dist
from sklearn.metrics import roc_curve
from torch import Tensor
from torch.nn import functional as F
from torchmetrics import AUROC, Accuracy

if TYPE_CHECKING:
    from MinkowskiEngine import SparseTensor

from transformercvn.dataset.minkowski_dataset import MinkowskiCollection, MinkowskiDataset
from transformercvn.network.networks.neutrino_full_base_network import NeutrinoBaseNetwork
from transformercvn.network.trainers.neutrino_base import NeutrinoBase
from transformercvn.options import Options

TArray = np.ndarray

# Historical NNbar checkpoint metric: one-vs-rest ROC for the last event class.
NNBAR_OVR_TARGET_FPR = 0.00054


def tpr_at_fpr_ovr_last_class(
    targets_np: np.ndarray, probs_np: np.ndarray, target_fpr: float
) -> float:
    """Return TPR at the first ROC point with FPR >= target_fpr for the last class."""
    n_classes = probs_np.shape[1]
    if n_classes < 2:
        return 0.0
    i = n_classes - 1
    y_score = probs_np[:, i]
    if np.sum(targets_np == i) == 0 or np.sum(targets_np != i) == 0:
        return 0.0
    fpr, tpr, _ = roc_curve(targets_np, y_score, pos_label=i, drop_intermediate=False)
    idx = int(np.searchsorted(fpr, target_fpr))
    idx = int(np.clip(idx, 0, len(fpr) - 1))
    return float(tpr[idx])


class NeutrinoFullBaseTrainer(NeutrinoBase, ABC):
    @abstractmethod
    def create_network(
        self,
        options: Options,
        features_dim: int,
        extra_dim: int,
        pixel_dim: int,
        num_prong_classes: int,
        num_event_classes: int,
    ) -> NeutrinoBaseNetwork:
        raise NotImplementedError()

    @abstractmethod
    def preprocess_pixels(
        self, pixel_coords: Tensor, pixel_values: Tensor, image_size: Tuple[int, ...]
    ) -> Union[Tensor, SparseTensor]:
        raise NotImplementedError()

    def __init__(self, options: Options):
        """

        Parameters
        ----------
        options: Options
            Global options for the entire network.
            See network.options.Options
        """
        super(NeutrinoFullBaseTrainer, self).__init__(options)

        self.hidden_dim = options.hidden_dim

        self.network = self.create_network(
            options,
            self.training_dataset.num_features,
            self.training_dataset.num_extra,
            self.training_dataset.pixel_features,
            self.training_dataset.num_prong_classes,
            self.training_dataset.num_event_classes,
        )

        self.beta = 1 - 1 / len(self.training_dataset)
        self.gamma = options.loss_gamma

        self.event_loss_scale = options.event_prong_loss_proportion
        self.prong_loss_scale = 1.0 - options.event_prong_loss_proportion

        self.event_accuracy = Accuracy(
            task="multiclass", num_classes=self.training_dataset.num_event_classes
        )
        self.prong_accuracy = Accuracy(
            task="multiclass", num_classes=self.training_dataset.num_prong_classes
        )

        self.event_auc = AUROC(
            task="multiclass", num_classes=self.training_dataset.num_event_classes
        )
        self.prong_auc = AUROC(
            task="multiclass", num_classes=self.training_dataset.num_prong_classes
        )

        self._val_epoch_event_targets: List[Tensor] = []
        self._val_epoch_event_probs: List[Tensor] = []

    @property
    def dataset(self):
        return MinkowskiDataset

    @property
    def dataloader_options(self):
        return {
            "drop_last": True,
            "batch_size": self.options.batch_size,
            # "pin_memory": self.options.num_gpu > 0,
            "num_workers": self.options.num_dataloader_workers,
            "collate_fn": MinkowskiCollection(),
        }

    def forward(
        self,
        features: Tensor,
        extra: Tensor,
        event_coords: Tensor,
        event_values: Tensor,
        event_mask: Tensor,
        prong_coords: Tensor,
        prong_values: Tensor,
        prong_mask: Tensor,
    ) -> Tuple[Tensor, Tensor]:
        batch_size, timesteps, _ = features.shape

        # Normalize the high level layers
        features = features.clone()
        features[prong_mask] -= self.mean
        features[prong_mask] /= self.std

        # Normalize the event level layers
        extra = extra.clone()
        extra -= self.extra_mean
        extra /= self.extra_std

        event_pixels = self.preprocess_pixels(
            event_coords, event_values, self.training_dataset.pixel_shape
        )
        prong_pixels = self.preprocess_pixels(
            prong_coords, prong_values, self.training_dataset.pixel_shape
        )

        return self.network(features, extra, event_pixels, event_mask, prong_pixels, prong_mask)

    def shared_step(self, batch):
        (
            features,
            extra,
            event_coordinates,
            event_pixel_values,
            event_masks,
            prong_coordinates,
            prong_pixel_values,
            prong_masks,
            event_targets,
            prong_targets,
        ) = batch

        max_prongs_in_batch = prong_masks.sum(1).max()
        features = features[:, :max_prongs_in_batch].contiguous()
        prong_masks = prong_masks[:, :max_prongs_in_batch].contiguous()
        prong_targets = prong_targets[:, :max_prongs_in_batch].contiguous()

        return (
            event_targets,
            prong_targets,
            *self.forward(
                features,
                extra,
                event_coordinates,
                event_pixel_values,
                event_masks,
                prong_coordinates,
                prong_pixel_values,
                prong_masks,
            ),
        )

    def loss(self, logits, targets):
        if self.gamma == 0:
            return F.cross_entropy(logits, targets)

        one_hot_targets = F.one_hot(targets.to(torch.int64), logits.shape[1]) > 0.5
        log_probabilities = torch.log_softmax(logits, dim=-1)
        probabilities = torch.softmax(logits, dim=-1)

        log_probabilities = torch.masked_select(log_probabilities, one_hot_targets)
        probabilities = torch.masked_select(probabilities, one_hot_targets)

        loss = -log_probabilities * (1 - probabilities) ** self.gamma
        return loss.mean()

    def training_step(self, batch, batch_idx):
        event_targets, prong_targets, event_logits, prong_logits = self.shared_step(batch)

        # Compute Event loss from simple prediction
        event_loss = self.loss(event_logits, event_targets)

        # Compute Prong Loss after masking away padding prongs
        prong_mask = prong_targets >= 0
        masked_prong_logits = torch.masked_select(prong_logits, prong_mask.unsqueeze(-1))
        masked_prong_logits = masked_prong_logits.reshape(-1, prong_logits.shape[-1])
        masked_prong_targets = torch.masked_select(prong_targets.long(), prong_mask)

        prong_loss = self.loss(masked_prong_logits, masked_prong_targets)

        # Total Loss is just sum
        total_loss = self.event_loss_scale * event_loss + self.prong_loss_scale * prong_loss

        # Prong Accuracy
        prong_accuracy = (masked_prong_logits.argmax(1) == masked_prong_targets).float().mean()

        # Event accuracy
        event_accuracy = (event_logits.argmax(1) == event_targets).float().mean()

        self.log("prong_loss", prong_loss)
        self.log("event_loss", event_loss)
        self.log("train_loss", total_loss)

        self.log("train_event_accuracy", event_accuracy)
        self.log("train_prong_accuracy", prong_accuracy)

        return total_loss

    def validation_step(self, batch, batch_idx):
        event_targets, prong_targets, event_logits, prong_logits = self.shared_step(batch)

        prong_mask = prong_targets >= 0
        masked_prong_logits = torch.masked_select(prong_logits, prong_mask.unsqueeze(-1))
        masked_prong_logits = masked_prong_logits.reshape(-1, prong_logits.shape[-1])
        masked_prong_targets = torch.masked_select(prong_targets.long(), prong_mask)

        event_probabilities = torch.softmax(event_logits, dim=-1)
        masked_prong_probabilities = torch.softmax(masked_prong_logits, dim=-1)

        self.prong_accuracy.update(masked_prong_probabilities, masked_prong_targets)
        self.event_accuracy.update(event_probabilities, event_targets)

        self.prong_auc.update(masked_prong_probabilities, masked_prong_targets)
        self.event_auc.update(event_probabilities, event_targets)

        self._val_epoch_event_targets.append(event_targets.detach().cpu())
        self._val_epoch_event_probs.append(event_probabilities.detach().cpu())

    def on_validation_epoch_start(self) -> None:
        self._val_epoch_event_targets.clear()
        self._val_epoch_event_probs.clear()

    def _gather_validation_event_predictions(
        self, targets: Tensor, probs: Tensor
    ) -> Tuple[Tensor, Tensor]:
        if not dist.is_available() or not dist.is_initialized() or dist.get_world_size() <= 1:
            return targets, probs
        t_np = targets.numpy()
        p_np = probs.numpy()
        gathered: List = [None] * dist.get_world_size()
        dist.all_gather_object(gathered, (t_np, p_np))
        t_cat = np.concatenate([g[0] for g in gathered])
        p_cat = np.concatenate([g[1] for g in gathered])
        return torch.from_numpy(t_cat), torch.from_numpy(p_cat)

    def validation_epoch_end(self, outputs) -> None:
        event_accuracy = self.event_accuracy.compute()
        prong_accuracy = self.prong_accuracy.compute()

        event_auc = self.event_auc.compute()
        prong_auc = self.prong_auc.compute()

        self.log("val_epoch_accuracy", (prong_accuracy + event_accuracy) / 2, sync_dist=True)
        self.log("event_epoch_accuracy", event_accuracy, sync_dist=True)
        self.log("prong_epoch_accuracy", prong_accuracy, sync_dist=True)

        self.log("val_epoch_AUC", (event_auc + prong_auc) / 2, sync_dist=True)
        self.log("event_epoch_AUC", event_auc, sync_dist=True)
        self.log("prong_epoch_AUC", prong_auc, sync_dist=True)

        if self._val_epoch_event_targets:
            targets_cpu = torch.cat(self._val_epoch_event_targets)
            probs_cpu = torch.cat(self._val_epoch_event_probs)
            targets_cpu, probs_cpu = self._gather_validation_event_predictions(
                targets_cpu, probs_cpu
            )
            nnbar_ovr_tpr = tpr_at_fpr_ovr_last_class(
                targets_cpu.numpy(),
                probs_cpu.numpy(),
                NNBAR_OVR_TARGET_FPR,
            )
        else:
            nnbar_ovr_tpr = 0.0
        self.log("val_epoch_nnbar_ovr_tpr", nnbar_ovr_tpr, sync_dist=True)

        self.event_accuracy.reset()
        self.prong_accuracy.reset()

        self.event_auc.reset()
        self.prong_auc.reset()
