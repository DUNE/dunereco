from typing import Tuple

import numpy as np
from torch import nn

from transformercvn.network.layers.resnet import ResNetStack
from transformercvn.options import Options


class PixelEncoder(nn.Module):
    def __init__(self, options: Options, input_dim: int, input_shape: Tuple[int, int], block_depth):
        super(PixelEncoder, self).__init__()

        self.options = options

        current_hidden_dim = self.options.initial_pixel_dim
        current_blocks = 0
        max_blocks = round(int(np.log2(max(input_shape))))

        block_sizes = [current_hidden_dim]
        block_depths = [block_depth]

        while (current_hidden_dim < options.hidden_dim) and (current_blocks < max_blocks):
            current_hidden_dim = 2 * current_hidden_dim
            current_blocks += 1

            block_sizes.append(current_hidden_dim)
            block_depths.append(block_depth)

        self.embedding = nn.Sequential(
            nn.Conv2d(input_dim, block_sizes[0], kernel_size=3, stride=1, padding=1, bias=False),
            nn.BatchNorm2d(block_sizes[0]),
            nn.PReLU(block_sizes[0]),
        )

        self.block = ResNetStack(block_sizes, block_depths)

        self.input_shape = np.array(input_shape)
        self.output_shape = self.input_shape // (2 ** (len(block_sizes) - 1))

        self.output = nn.Identity()
        if (self.output_shape > 1).any():
            self.output = nn.Sequential(
                nn.Conv2d(
                    block_sizes[-1], options.hidden_dim, kernel_size=tuple(self.output_shape)
                ),
                nn.BatchNorm2d(options.hidden_dim),
                nn.PReLU(options.hidden_dim),
            )

    def forward(self, pixels):
        y = self.embedding(pixels)
        y = self.block(y)
        y = self.output(y)
        return y.reshape(-1, self.options.hidden_dim)
