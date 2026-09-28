# Copyright (c) Meta Platforms, Inc. and affiliates.

# All rights reserved.

# This source code is licensed under the license found in the
# LICENSE file in the root directory of this source tree.

from collections import OrderedDict
from typing import Tuple

import torch
import torch.nn as nn

from transformercvn.network.layers.sparse_layers import (
    SparseCondense,
    SparseConv2D,
    SparseDropPath,
    SparseGeLU,
    SparseLayerNorm,
    SparseTensor,
)


class ConvNextBlock(nn.Module):
    r"""ConvNeXt Block. There are two equivalent implementations:
    (1) DwConv -> LayerNorm (channels_first) -> 1x1 Conv -> GELU -> 1x1 Conv; all in (N, C, H, W)
    (2) DwConv -> Permute to (N, H, W, C); LayerNorm (channels_last) -> Linear -> GELU -> Linear; Permute back
    We use (2) as we find it slightly faster in PyTorch

    Args:
        hidden_dim (int): Number of input channels.
        drop_path (float): Stochastic depth rate. Default: 0.0
        layer_scale_init_value (float): Init value for Layer Scale. Default: 1e-6.
    """

    def __init__(
        self,
        hidden_dim,
        kernel_size: int = 3,
        expansion: int = 4,
        drop_path: float = 0.0,
        layer_scale_init_value: float = 1e-6,
        learned_skip_connection: bool = False,
    ):
        super().__init__()
        self.conv = SparseConv2D(
            hidden_dim, hidden_dim, kernel_size=kernel_size
        )  # depthwise conv TODO
        nn.init.trunc_normal_(self.conv.kernel, std=0.02)

        self.expand_conv = SparseConv2D(hidden_dim, expansion * hidden_dim, kernel_size=1)
        nn.init.trunc_normal_(self.expand_conv.kernel, std=0.02)

        self.contract_conv = SparseConv2D(expansion * hidden_dim, hidden_dim, kernel_size=1)
        nn.init.trunc_normal_(self.contract_conv.kernel, std=0.02)

        self.norm = SparseLayerNorm(hidden_dim, eps=1e-6)
        self.gelu = SparseGeLU()

        self.skip_connection = SparseConv2D(hidden_dim, hidden_dim, kernel_size=kernel_size)
        if not learned_skip_connection:
            self.skip_connection.kernel.requires_grad_(False)
            self.skip_connection.kernel.zero_()
            self.skip_connection.kernel.view(kernel_size, kernel_size, hidden_dim, hidden_dim)[
                kernel_size // 2, kernel_size // 2
            ] += torch.eye(hidden_dim)

        if layer_scale_init_value > 0:
            self.gamma = nn.Parameter(
                layer_scale_init_value * torch.ones(hidden_dim), requires_grad=True
            )
        else:
            self.gamma = 1.0

        if drop_path > 0:
            self.drop_path = SparseDropPath(drop_path)
        else:
            self.drop_path = nn.Identity()

    def forward(self, x):
        input = x

        x = self.conv(x)
        x = self.norm(x)
        x = self.expand_conv(x)
        x = self.gelu(x)
        x = self.contract_conv(x)

        x = SparseTensor(
            self.gamma * x.F,
            coordinate_map_key=x.coordinate_map_key,
            coordinate_manager=x.coordinate_manager,
        )

        x = self.skip_connection(input) + self.drop_path(x)

        return x


class SparseConvNeXt(nn.Module):
    def __init__(
        self,
        input_features: int,
        output_features: int,
        kernel_size: int = 3,
        hidden_features: Tuple[int, ...] = (96, 192, 384, 768),
        hidden_depths: Tuple[int, ...] = (3, 3, 3, 3),
        drop_path_rate: float = 0.0,
        layer_scale_init_value=1e-6,
    ) -> None:
        super().__init__()

        self.num_blocks = len(hidden_depths)
        self.total_layers = sum(hidden_depths)

        features = nn.Sequential(
            SparseConv2D(input_features, hidden_features[0], kernel_size=4, stride=4),
            SparseLayerNorm(hidden_features[0], eps=1e-6),
        )

        self.downsample_layers = nn.ModuleList(
            [features]
            + [
                nn.Sequential(
                    SparseLayerNorm(hidden_features[i], eps=1e-6),
                    SparseConv2D(
                        hidden_features[i], hidden_features[i + 1], kernel_size=2, stride=2
                    ),
                )
                for i in range(self.num_blocks - 1)
            ]
        )

        self.stages = nn.ModuleList()

        drop_rates = [x.item() for x in torch.linspace(0, drop_path_rate, self.total_layers)]
        current_layer = 0

        for i in range(self.num_blocks):
            stage = nn.Sequential(
                *(
                    ConvNextBlock(
                        hidden_dim=hidden_features[i],
                        kernel_size=kernel_size,
                        drop_path=drop_rates[current_layer + j],
                        layer_scale_init_value=layer_scale_init_value,
                    )
                    for j in range(hidden_depths[i])
                )
            )

            self.stages.append(stage)
            current_layer += hidden_depths[i]

        # Final pooling layer to convert whatever is left of the image into a dense vector.
        self.condense = SparseCondense()

        # Output linear layer to convert whatever latent dimension we end up on to the desired output dimension.
        self.output_block = nn.Sequential(
            OrderedDict(
                [
                    ("norm", nn.LayerNorm(hidden_features[-1])),
                    ("linear", nn.Linear(hidden_features[-1], output_features, bias=False)),
                    ("gelu", nn.GELU()),
                ]
            )
        )

    def forward(self, x):
        for i in range(4):
            x = self.downsample_layers[i](x)
            x = self.stages[i](x)

        x = self.condense(x)
        return self.output_block(x)
