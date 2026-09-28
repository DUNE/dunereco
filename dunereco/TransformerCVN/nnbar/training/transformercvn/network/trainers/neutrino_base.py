import pytorch_lightning as pl
import torch
from torch.utils.data import DataLoader

from transformercvn.dataset import NeutrinoDataset
from transformercvn.network.networks.learning_rate_schedules import (
    get_cosine_with_hard_restarts_schedule_with_warmup,
    get_linear_schedule_with_warmup,
)
from transformercvn.options import Options


class NeutrinoBase(pl.LightningModule):
    def __init__(self, options: Options):
        super(NeutrinoBase, self).__init__()

        self.options = options

        self.training_dataset, self.validation_dataset, self.testing_dataset = (
            self.create_datasets()
        )

        self.mean = 0
        self.std = 1

        self.extra_mean = 0
        self.extra_std = 1

        self.pixel_mean = 0
        self.pixel_std = 1

        # Normalize datasets using training dataset statistics
        if self.options.normalize_features:
            (
                self.mean,
                self.std,
                self.extra_mean,
                self.extra_std,
                self.pixel_mean,
                self.pixel_std,
            ) = self.training_dataset.compute_statistics()

            self.mean = torch.nn.Parameter(self.mean, requires_grad=False)
            self.std = torch.nn.Parameter(self.std, requires_grad=False)

            self.extra_mean = torch.nn.Parameter(self.extra_mean, requires_grad=False)
            self.extra_std = torch.nn.Parameter(self.extra_std, requires_grad=False)

            if self.training_dataset.pixels is not None:
                self.pixel_mean = torch.nn.Parameter(self.pixel_mean, requires_grad=False)
                self.pixel_std = torch.nn.Parameter(self.pixel_std, requires_grad=False)

        self.steps_per_epoch = len(self.training_dataset) // (
            self.options.batch_size * max(1, self.options.num_gpu)
        )
        self.total_steps = self.steps_per_epoch * self.options.epochs
        self.warmup_steps = int(
            round(self.steps_per_epoch * self.options.learning_rate_warmup_epochs)
        )

    @property
    def dataset(self):
        return NeutrinoDataset

    @property
    def dataloader(self):
        return DataLoader

    @property
    def dataloader_options(self):
        return {
            "drop_last": True,
            "batch_size": self.options.batch_size,
            "pin_memory": self.options.num_gpu > 0,
            "num_workers": self.options.num_dataloader_workers,
        }

    def create_datasets(self):
        if len(self.options.validation_file) > 0:
            training_dataset = self.dataset(
                self.options.training_file,
                event_current_targets=self.options.event_current_targets,
                load_full_dataset=self.options.load_full_dataset,
            )
            validation_dataset = self.dataset(
                self.options.validation_file,
                event_current_targets=self.options.event_current_targets,
                load_full_dataset=self.options.load_full_dataset,
            )
        else:
            # Compute the training / validation ranges based on the data-split and the limiting percentage.
            train_validation_split = (
                self.options.dataset_limit * self.options.train_validation_split
            )
            train_range = (0.0, train_validation_split)
            validation_range = (train_validation_split, self.options.dataset_limit)

            training_dataset = self.dataset(
                self.options.training_file,
                train_range,
                event_current_targets=self.options.event_current_targets,
                load_full_dataset=self.options.load_full_dataset,
            )

            validation_dataset = self.dataset(
                self.options.training_file,
                validation_range,
                event_current_targets=self.options.event_current_targets,
                load_full_dataset=self.options.load_full_dataset,
            )

        testing_dataset = None
        if len(self.options.testing_file) > 0:
            testing_dataset = self.dataset(
                self.options.testing_file,
                event_current_targets=self.options.event_current_targets,
                load_full_dataset=self.options.load_full_dataset,
            )

        return training_dataset, validation_dataset, testing_dataset

    def configure_optimizers(self):
        optimizer = None

        if "apex" in self.options.optimizer:
            try:
                import apex.optimizers

                if self.options.optimizer == "apex_adam":
                    optimizer = apex.optimizers.FusedAdam

                elif self.options.optimizer == "apex_lamb":
                    optimizer = apex.optimizers.FusedLAMB

                else:
                    optimizer = apex.optimizers.FusedSGD

            except ImportError:
                pass

        else:
            optimizer = getattr(torch.optim, self.options.optimizer)

        if optimizer is None:
            print(f"Unable to load desired optimizer: {self.options.optimizer}.")
            print("Using pytorch AdamW as a default.")
            optimizer = torch.optim.AdamW

        decay_mask = ["bias", "LayerNorm.weight"]
        optimizer_grouped_parameters = [
            {
                "params": [
                    param
                    for name, param in self.named_parameters()
                    if not any(no_decay in name for no_decay in decay_mask)
                ],
                "weight_decay": self.options.l2_penalty,
            },
            {
                "params": [
                    param
                    for name, param in self.named_parameters()
                    if any(no_decay in name for no_decay in decay_mask)
                ],
                "weight_decay": 0.0,
            },
        ]

        optimizer = optimizer(optimizer_grouped_parameters, lr=self.options.learning_rate)

        if self.options.learning_rate_cycles < 1:
            scheduler = get_linear_schedule_with_warmup(
                optimizer, num_warmup_steps=self.warmup_steps, num_training_steps=self.total_steps
            )
        else:
            scheduler = get_cosine_with_hard_restarts_schedule_with_warmup(
                optimizer,
                num_warmup_steps=self.warmup_steps,
                num_training_steps=self.total_steps,
                num_cycles=self.options.learning_rate_cycles,
            )

        scheduler = {"scheduler": scheduler, "interval": "step", "frequency": 1}

        return [optimizer], [scheduler]

    def train_dataloader(self) -> DataLoader:
        return self.dataloader(self.training_dataset, shuffle=True, **self.dataloader_options)

    def val_dataloader(self) -> DataLoader:
        return self.dataloader(self.validation_dataset, **self.dataloader_options)

    def test_dataloader(self) -> DataLoader:
        if self.testing_dataset is None:
            raise ValueError("Testing dataset not provided.")

        return self.dataloader(self.testing_dataset, **self.dataloader_options)
