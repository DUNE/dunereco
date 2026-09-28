from collections import OrderedDict
from functools import partial
from typing import Tuple

from torch import nn


class Conv2dAuto(nn.Conv2d):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # dynamic add padding based on the kernel_size
        self.padding = (self.kernel_size[0] // 2, self.kernel_size[1] // 2)


conv3x3 = partial(Conv2dAuto, kernel_size=3, bias=False)


class ResidualBlock(nn.Module):
    def __init__(self, in_channels, out_channels):
        super().__init__()
        self.in_channels, self.out_channels = in_channels, out_channels
        self.blocks = nn.Identity()
        self.shortcut = nn.Identity()

    def forward(self, x):
        residual = x
        if self.should_apply_shortcut:
            residual = self.shortcut(x)
        x = self.blocks(x)
        x += residual
        return x

    @property
    def should_apply_shortcut(self):
        return self.in_channels != self.out_channels


class ResNetResidualBlock(ResidualBlock):
    def __init__(
        self, in_channels, out_channels, expansion=1, down_sampling=1, conv=conv3x3, *args, **kwargs
    ):
        super().__init__(in_channels, out_channels)
        self.expansion, self.down_sampling, self.conv = expansion, down_sampling, conv
        self.shortcut = (
            nn.Sequential(
                OrderedDict(
                    {
                        "conv": nn.Conv2d(
                            self.in_channels,
                            self.expanded_channels,
                            kernel_size=1,
                            stride=self.down_sampling,
                            bias=False,
                        ),
                        "bn": nn.BatchNorm2d(self.expanded_channels),
                    }
                )
            )
            if self.should_apply_shortcut
            else None
        )

    @property
    def expanded_channels(self):
        return self.out_channels * self.expansion

    @property
    def should_apply_shortcut(self):
        return self.in_channels != self.expanded_channels


def conv_bn(in_channels, out_channels, conv, *args, **kwargs):
    return nn.Sequential(
        OrderedDict(
            {
                "conv": conv(in_channels, out_channels, *args, **kwargs),
                "bn": nn.BatchNorm2d(out_channels),
            }
        )
    )


class ResNetBasicBlock(ResNetResidualBlock):
    expansion = 1

    def __init__(self, in_channels, out_channels, *args, **kwargs):
        super().__init__(in_channels, out_channels, *args, **kwargs)
        self.blocks = nn.Sequential(
            conv_bn(
                self.in_channels,
                self.out_channels,
                conv=self.conv,
                bias=False,
                stride=self.down_sampling,
            ),
            nn.PReLU(self.out_channels),
            conv_bn(self.out_channels, self.expanded_channels, conv=self.conv, bias=False),
        )


class ResNetBottleNeckBlock(ResNetResidualBlock):
    expansion = 4

    def __init__(self, in_channels, out_channels, *args, **kwargs):
        super(ResNetBottleNeckBlock, self).__init__(
            in_channels, out_channels, expansion=4, *args, **kwargs
        )

        self.blocks = nn.Sequential(
            conv_bn(self.in_channels, self.out_channels, self.conv, kernel_size=1),
            nn.PReLU(self.out_channels),
            conv_bn(
                self.out_channels,
                self.out_channels,
                self.conv,
                kernel_size=3,
                stride=self.down_sampling,
            ),
            nn.PReLU(self.out_channels),
            conv_bn(self.out_channels, self.expanded_channels, self.conv, kernel_size=1),
        )


class ResNetLayer(nn.Module):
    def __init__(self, in_channels, out_channels, block=ResNetBasicBlock, n=1, *args, **kwargs):
        super(ResNetLayer, self).__init__()

        down_sampling = 2 if in_channels != out_channels else 1

        self.blocks = nn.Sequential(
            block(in_channels, out_channels, *args, **kwargs, down_sampling=down_sampling),
            *[
                block(
                    out_channels * block.expansion, out_channels, down_sampling=1, *args, **kwargs
                )
                for _ in range(n - 1)
            ],
        )

    def forward(self, x):
        x = self.blocks(x)
        return x


class ResNetStack(nn.Module):
    def __init__(
        self,
        blocks_sizes: Tuple[int] = (64, 128, 256, 512),
        depths: Tuple[int] = (2, 2, 2, 2),
        block=ResNetBasicBlock,
        *args,
        **kwargs,
    ):
        super(ResNetStack, self).__init__()

        if len(blocks_sizes) != len(depths):
            raise ValueError("Length of ResNet block depths and sizes do not match.")

        self.in_out_block_sizes = list(zip(blocks_sizes, blocks_sizes[1:]))
        self.blocks_sizes = blocks_sizes

        self.blocks = nn.ModuleList(
            [
                ResNetLayer(
                    blocks_sizes[0], blocks_sizes[0], n=depths[0], block=block, *args, **kwargs
                ),
                *[
                    ResNetLayer(
                        in_channels * block.expansion,
                        out_channels,
                        n=n,
                        block=block,
                        *args,
                        **kwargs,
                    )
                    for (in_channels, out_channels), n in zip(self.in_out_block_sizes, depths[1:])
                ],
            ]
        )

    def forward(self, x):
        for block in self.blocks:
            x = block(x)
        return x
