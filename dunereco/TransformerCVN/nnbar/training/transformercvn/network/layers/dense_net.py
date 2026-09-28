from collections import OrderedDict
from typing import Tuple

import torch
from torch import Tensor, nn


class Bottleneck(nn.Module):
    def __init__(
        self, input_features: int, growth_rate: int, batch_norm_size: int, dropout: float
    ) -> None:
        super(Bottleneck, self).__init__()

        self.bottleneck_block = nn.Sequential(
            OrderedDict(
                [
                    ("norm1", nn.BatchNorm2d(input_features)),
                    ("relu1", nn.PReLU(input_features)),
                    (
                        "conv1",
                        nn.Conv2d(
                            input_features,
                            batch_norm_size * growth_rate,
                            kernel_size=1,
                            stride=1,
                        ),
                    ),
                ]
            )
        )

        self.output_block = nn.Sequential(
            OrderedDict(
                [
                    ("norm2", nn.BatchNorm2d(batch_norm_size * growth_rate)),
                    ("relu2", nn.PReLU(batch_norm_size * growth_rate)),
                    (
                        "conv2",
                        nn.Conv2d(
                            batch_norm_size * growth_rate,
                            growth_rate,
                            kernel_size=3,
                            padding=1,
                            stride=1,
                        ),
                    ),
                    ("dropout", nn.Dropout(dropout)),
                ]
            )
        )

    def forward(self, features: Tensor) -> Tensor:
        bottleneck = self.bottleneck_block(features)
        output = self.output_block(bottleneck)
        return torch.cat((features, output), dim=1)


class DenseBlock(nn.Module):
    def __init__(
        self,
        num_layers: int,
        input_features: int,
        batch_norm_size: int,
        growth_rate: int,
        dropout: float,
    ) -> None:
        super(DenseBlock, self).__init__()

        self.layers = nn.ModuleList(
            [
                Bottleneck(input_features + i * growth_rate, growth_rate, batch_norm_size, dropout)
                for i in range(num_layers)
            ]
        )

    def forward(self, init_features: Tensor) -> Tensor:
        features = init_features

        for layer in self.layers:
            features = layer(features)

        return features


class Transition(nn.Sequential):
    def __init__(self, input_features: int, output_features: int) -> None:
        super(Transition, self).__init__(
            OrderedDict(
                [
                    ("norm", nn.BatchNorm2d(input_features)),
                    ("relu", nn.PReLU(input_features)),
                    ("conv", nn.Conv2d(input_features, output_features, kernel_size=1, stride=1)),
                    ("pooling", nn.AvgPool2d(kernel_size=2, stride=2)),
                ]
            )
        )


class DenseNet(nn.Module):
    def __init__(
        self,
        input_features: int,
        output_features: int,
        initial_latent_features: int = 64,
        growth_rate: int = 32,
        batch_norm_size: int = 4,
        block_config: Tuple[int, ...] = (6, 12, 24, 16),
        dropout: float = 0.0,
    ) -> None:
        super().__init__()

        # Construct first convolution to get input input latent space.
        self.features = nn.Sequential(
            OrderedDict(
                [
                    (
                        "conv0",
                        nn.Conv2d(
                            input_features,
                            initial_latent_features,
                            kernel_size=7,
                            padding=3,
                            stride=2,
                        ),
                    ),
                    ("norm0", nn.BatchNorm2d(initial_latent_features)),
                    ("relu0", nn.PReLU(initial_latent_features)),
                    ("pooling0", nn.AvgPool2d(kernel_size=3, stride=2)),
                ]
            )
        )

        # Construct each dense block in increasing latent dimensions.
        num_features = initial_latent_features
        for i, num_layers in enumerate(block_config):
            self.features.add_module(
                f"dense{(i + 1)}",
                DenseBlock(num_layers, num_features, batch_norm_size, growth_rate, dropout),
            )

            # Next layer will have slightly more features
            num_features = num_features + num_layers * growth_rate

            if i != len(block_config) - 1:
                self.features.add_module(
                    f"transition{(i + 1)}",
                    Transition(input_features=num_features, output_features=num_features // 2),
                )

                num_features = num_features // 2

        # Final non-linearity and normalization layers.
        self.features.add_module("final_norm", nn.BatchNorm2d(num_features))
        self.features.add_module("final_relu", nn.PReLU(num_features))

        # Final pooling layer to convert whatever is left of the image into a dense vector.
        self.condense = nn.Sequential(nn.AdaptiveAvgPool2d((1, 1)), nn.Flatten())

        # Output linear layer to convert whatever latent dimension we end up on to the desired output dimension.
        self.output_block = nn.Sequential(
            OrderedDict(
                [
                    ("linear", nn.Linear(num_features, output_features, bias=False)),
                    ("norm", nn.BatchNorm1d(output_features)),
                    ("relu", nn.PReLU(output_features)),
                    ("dropout", nn.Dropout(dropout)),
                ]
            )
        )

    def forward(self, x: Tensor) -> Tensor:
        features = self.features(x)
        output = self.condense(features)
        return self.output_block(output)
