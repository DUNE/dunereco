import json
from argparse import Namespace
from typing import List, Optional


class Options(Namespace):
    def __init__(self, training_file: str = "", testing_file: str = "", validation_file: str = ""):
        super(Options, self).__init__()

        # =========================================================================================
        # Network Architecture
        # =========================================================================================

        # Dimensions used internally by all hidden layers / transformers (should be a multiple of 8 for efficiency).
        self.hidden_dim: int = 128

        # Hidden dimensionality of the first embedding layers after input (should be multiples of 8 for efficiency).
        self.initial_feature_dim: int = 32
        self.initial_pixel_dim: int = 16

        # Embedding Dimension sizes (should be multiples of 8 for efficiency).
        self.feature_embedding_dim: int = 8
        self.pixel_embedding_dim: int = 512
        self.position_embedding_dim: int = 16

        # Smallest layer size to drop to during decoder layers (should be a multiple of 8 for efficiency).
        self.final_decoder_dim: int = 16

        # Maximum Number of double-sized embedding layers to add between the layers and the encoder.
        # The size of the embedding dimension will be capped at the hidden_dim,
        # So setting this option to a very large integer will just keep embedding up to the hidden_dim.
        self.num_embedding_layers: int = 100

        # Number of encoder layers for the central shared transformer.
        self.num_encoder_layers: int = 5

        # Number of hidden layers for the classification decoder.
        self.num_decoder_layers: int = 100

        self.num_prong_decoder_layers: int = 4

        # Number of heads for multi-head attention, used in all transformer layers.
        self.num_attention_heads: int = 8

        # Activation function for all transformer layers, 'relu' or 'gelu'.
        self.transformer_activation: str = "gelu"

        # Use the more common variant of transformer with the norm in the front
        self.transformer_norm_first: bool = False

        # Whether or not to use PreLU activation on linear / embedding layers,
        # Otherwise a regular relu will be used.
        self.linear_prelu_activation: bool = True

        # Whether or not to apply batch norm on linear / embedding layers.
        self.linear_batch_norm: bool = True

        # Whether or not to zero out reconstructed variable inputs
        self.disable_smart_features: bool = False

        # Normalize inputs using mean and std of datasets
        self.normalize_features: bool = True

        self.one_hot_pixels: bool = False
        self.log_pixels: bool = False

        self.mobilenet_structure: Optional[List[List[int]]] = None

        self.densenet_structure: List[int] = [6, 12, 24, 16]
        self.densenet_growth_rate: int = 16
        self.densenet_batch_norm_size = 4

        # =========================================================================================
        # Dataset Options
        # =========================================================================================

        self.training_file: str = training_file
        self.testing_file: str = testing_file
        self.validation_file: str = validation_file

        # Limit the dataset to the first images% of the data.
        self.dataset_limit: float = 1.0

        # Percent of data to use for training vs. validation.
        self.train_validation_split: float = 0.95

        # Training batch size.
        self.batch_size: int = 2048

        # Number of processes to spawn for data collection.
        self.num_dataloader_workers: int = 8

        # Load full dataset into RAM at start.
        self.load_full_dataset: bool = False

        # Use basic current targets instead of detailed.
        self.event_current_targets: bool = False

        # =========================================================================================
        # Training Options
        # =========================================================================================

        # The optimizer to use for training the network.
        # This must be a valid class in torch.optim or nvidia apex with 'apex' prefix.
        self.optimizer: str = "AdamW"

        # Optimizer learning rate.
        self.learning_rate: float = 0.0001

        # Optimizer l2 penalty based on weight values.
        self.l2_penalty: float = 0.015

        # Clip the L2 norm of the gradient. Set to 0.0 to disable.
        self.gradient_clip: float = 90.0

        # Dropout added to all layers.
        self.dropout: float = 0.0

        # Number of epochs to train for.
        self.epochs: int = 25

        # Number of epochs to ramp up the learning rate up to the given value. Can be fractional.
        self.learning_rate_warmup_epochs: float = 1.0

        # Number of times to cycles the learning rate through cosine annealing with hard resets.
        # Set to 0 to disable cosine annealing and just use a decaying learning rate.
        self.learning_rate_cycles: int = 1

        # Total number of GPUs to use.
        self.num_gpu: int = 1

        self.event_prong_loss_proportion: float = 0.5

        # Not implemented. Previously used to upweight the neutral current class in CB_loss.
        self.loss_beta: float = 2.5

        # Exponent in the focal loss term (1 - p). If 0, simple cross-entropy is used.
        self.loss_gamma: float = 0.0

        # Standard deviation of the noise to add to pixel-maps
        self.pixel_noise_std = 0.01

        # =========================================================================================
        # Miscellaneous Options
        # =========================================================================================

        # Whether or not to print additional information during training and log extra metrics.
        self.verbose_output: bool = True

        # Misc parameters used by sherpa to delegate GPUs and output directories.
        # These should not be set manually.
        self.usable_gpus: str = ""

        self.trial_time: str = ""

        self.trial_output_dir: str = "./test_output"

    def update_options(self, new_options):
        """Apply overrides, preserving boolean types and historical extension keys."""
        integer_options = {
            key
            for key, val in self.__dict__.items()
            if isinstance(val, int) and not isinstance(val, bool)
        }
        boolean_options = {key for key, val in self.__dict__.items() if isinstance(val, bool)}
        for key, value in new_options.items():
            if key in integer_options:
                setattr(self, key, int(value))
            elif key in boolean_options:
                if not isinstance(value, (bool, int)) or value not in (False, True, 0, 1):
                    raise ValueError(f"{key} must be a boolean or 0/1, got {value!r}")
                setattr(self, key, bool(value))
            else:
                setattr(self, key, value)

    @classmethod
    def load(cls, filepath: str):
        options = cls()
        with open(filepath, "r") as json_file:
            options.update_options(json.load(json_file))
        return options

    def display(self):
        print("=" * 70)
        print("Options")
        print("-" * 70)
        for key, val in sorted(vars(self).items()):
            print(f"{key:32}: {val}")
        print("=" * 70)
