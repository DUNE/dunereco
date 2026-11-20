from typing import Tuple

import numpy as np
from torch import Tensor, nn

from transformercvn.network.layers.resnet import ResNetStack


class ProngPixelEmbedding(nn.Module):
    __constants__ = ["input_dim", "hidden_dim", "block_depth"]

    def __init__(
        self,
        input_dim: int,
        input_shape: Tuple[int, int],
        hidden_dim: int,
        block_depth: int = 1,
        final_size: int = 1,
    ):
        super(ProngPixelEmbedding, self).__init__()

        self.input_dim = input_dim
        self.hidden_dim = hidden_dim
        self.block_depth = block_depth

        max_blocks = int(np.log2(min(input_shape) / final_size))
        initial_cnn_dim = hidden_dim // (2**max_blocks)

        current_block_dim = initial_cnn_dim
        current_block_size = final_size
        current_block_count = 0

        block_dims = [current_block_dim]
        block_depths = [block_depth]

        while (current_block_dim < hidden_dim) and (current_block_count < max_blocks):
            current_block_dim *= 2
            current_block_size *= 2
            current_block_count += 1

            block_dims.append(current_block_dim)
            block_depths.append(block_depth)

        initial_kernel_size = np.array(input_shape) - np.array(
            [current_block_size, current_block_size]
        )
        initial_kernel_size = initial_kernel_size + 3

        self.initial_cnn = nn.Sequential(
            nn.Conv2d(
                input_dim,
                initial_cnn_dim,
                kernel_size=tuple(map(int, initial_kernel_size)),
                stride=1,
                padding=1,
                bias=False,
            ),
            nn.BatchNorm2d(block_dims[0]),
            nn.PReLU(block_dims[0]),
        )

        self.resnet = ResNetStack(block_dims, block_depths)

        input_shape = np.array([current_block_size, current_block_size])
        output_shape = input_shape // (2 ** (len(block_dims) - 1))

        self.output = nn.Identity()
        if (output_shape > 1).any():
            self.output = nn.Sequential(
                nn.Conv2d(block_dims[-1], hidden_dim, kernel_size=tuple(map(int, output_shape))),
                nn.BatchNorm2d(hidden_dim),
                nn.PReLU(hidden_dim),
            )

    def forward(self, pixels: Tensor, mask: Tensor) -> Tensor:
        batch_size, max_particles, channels, height, width = pixels.shape

        pixel_mask = mask.view(batch_size, max_particles, 1)

        # Apply the CNN to every timestep independently
        pixels = pixels.view(batch_size * max_particles, channels, height, width)

        y = self.initial_cnn(pixels)
        y = self.resnet(y)
        y = self.output(y)
        y = y.reshape(batch_size, max_particles, self.hidden_dim)

        return y * pixel_mask
