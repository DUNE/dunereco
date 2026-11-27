from collections import OrderedDict
from typing import Tuple

import torch
from torch import Tensor, nn

from transformercvn.network.layers.sparse_layers import (
    SparseAvgPooling2D,
    SparseBatchNorm2D,
    SparseConcatenation,
    SparseCondense,
    SparseConv2D,
    SparseDropout,
    SparsePReLU,
    SparseTensor,
)


class SparseDenseLayer(nn.Module):
    def __init__(
        self,
        input_features: int,
        growth_rate: int,
        batch_norm_size: int,
        dropout: float,
        expand_coordinates: bool = True,
    ) -> None:
        super(SparseDenseLayer, self).__init__()

        self.bottleneck_block = nn.Sequential(
            OrderedDict(
                [
                    ("norm1", SparseBatchNorm2D(input_features)),
                    ("relu1", SparsePReLU(input_features)),
                    (
                        "conv1",
                        SparseConv2D(
                            input_features,
                            batch_norm_size * growth_rate,
                            kernel_size=1,
                            stride=1,
                            expand_coordinates=expand_coordinates,
                        ),
                    ),
                ]
            )
        )

        self.output_block = nn.Sequential(
            OrderedDict(
                [
                    ("norm2", SparseBatchNorm2D(batch_norm_size * growth_rate)),
                    ("relu2", SparsePReLU(batch_norm_size * growth_rate)),
                    (
                        "conv2",
                        SparseConv2D(
                            batch_norm_size * growth_rate,
                            growth_rate,
                            kernel_size=3,
                            stride=1,
                            expand_coordinates=expand_coordinates,
                        ),
                    ),
                    ("dropout", SparseDropout(dropout)),
                ]
            )
        )

        self.concatenate = SparseConcatenation()

        self.expand_coordinates = expand_coordinates
        if expand_coordinates:
            self.skip_connection = SparseConv2D(
                input_features, input_features, kernel_size=3, expand_coordinates=expand_coordinates
            )
            self.skip_connection.kernel.requires_grad_(False)
            self.skip_connection.kernel.zero_()
            self.skip_connection.kernel.view(3, 3, input_features, input_features)[1, 1] += (
                torch.eye(input_features)
            )

    def forward(self, features: SparseTensor) -> Tensor:
        bottleneck = self.bottleneck_block(features)
        output = self.output_block(bottleneck)

        if self.expand_coordinates:
            features = self.skip_connection(features, output.coordinate_map_key)

        return self.concatenate(features, output)


class SparseDenseBlock(nn.Module):
    def __init__(
        self,
        num_layers: int,
        input_features: int,
        batch_norm_size: int,
        growth_rate: int,
        dropout: float,
        expand_coordinates: bool = True,
    ) -> None:
        super(SparseDenseBlock, self).__init__()

        self.layers = nn.ModuleList(
            [
                SparseDenseLayer(
                    input_features + i * growth_rate,
                    growth_rate,
                    batch_norm_size,
                    dropout,
                    expand_coordinates=expand_coordinates,
                )
                for i in range(num_layers)
            ]
        )

    def forward(self, init_features: SparseTensor) -> SparseTensor:
        features = init_features

        for layer in self.layers:
            features = layer(features)

        return features


class SparseTransition(nn.Sequential):
    def __init__(
        self, input_features: int, output_features: int, expand_coordinates: bool = True
    ) -> None:
        super(SparseTransition, self).__init__(
            OrderedDict(
                [
                    ("norm", SparseBatchNorm2D(input_features)),
                    ("relu", SparsePReLU(input_features)),
                    (
                        "conv",
                        SparseConv2D(
                            input_features,
                            output_features,
                            kernel_size=1,
                            stride=1,
                            expand_coordinates=expand_coordinates,
                        ),
                    ),
                    ("pooling", SparseAvgPooling2D(kernel_size=2, stride=2)),
                ]
            )
        )


class SparseDenseNet(nn.Module):
    def __init__(
        self,
        input_features: int,
        output_features: int,
        initial_latent_features: int = 64,
        growth_rate: int = 32,
        batch_norm_size: int = 4,
        block_config: Tuple[int, ...] = (6, 12, 24, 16),
        dropout: float = 0.0,
        expand_coordinates: bool = True,
    ) -> None:
        super().__init__()

        # Construct first convolution to get input input latent space.
        self.features = nn.Sequential(
            OrderedDict(
                [
                    (
                        "conv0",
                        SparseConv2D(
                            input_features,
                            initial_latent_features,
                            kernel_size=7,
                            stride=2,
                            expand_coordinates=expand_coordinates,
                        ),
                    ),
                    ("norm0", SparseBatchNorm2D(initial_latent_features)),
                    ("relu0", SparsePReLU(initial_latent_features)),
                    ("pooling0", SparseAvgPooling2D(kernel_size=3, stride=2)),
                ]
            )
        )

        # Construct each dense block in increasing latent dimensions.
        num_features = initial_latent_features
        for i, num_layers in enumerate(block_config):
            self.features.add_module(
                f"dense{(i + 1)}",
                SparseDenseBlock(
                    num_layers,
                    num_features,
                    batch_norm_size,
                    growth_rate,
                    dropout,
                    expand_coordinates=expand_coordinates,
                ),
            )

            # Next layer will have slightly more features
            num_features = num_features + num_layers * growth_rate

            if i != len(block_config) - 1:
                self.features.add_module(
                    f"transition{(i + 1)}",
                    SparseTransition(
                        input_features=num_features,
                        output_features=num_features // 2,
                        expand_coordinates=expand_coordinates,
                    ),
                )

                num_features = num_features // 2

        # Final non-linearity and normalization layers.
        self.features.add_module("final_norm", SparseBatchNorm2D(num_features))
        self.features.add_module("final_relu", SparsePReLU(num_features))

        # Final pooling layer to convert whatever is left of the image into a dense vector.
        self.condense = SparseCondense()

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
