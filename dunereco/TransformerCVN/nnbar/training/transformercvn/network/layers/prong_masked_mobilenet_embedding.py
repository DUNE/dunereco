from typing import List, Optional, Tuple

import torch
from torch import Tensor, nn
from torch.nn import functional as F


def make_divisible_channel_count(v: float, divisor: int, min_value: Optional[int] = None) -> int:
    """
    This function is taken from the original tf repo.
    It ensures that all layers have a channel number that is divisible by 8
    It can be seen here:
    https://github.com/tensorflow/models/blob/master/research/slim/nets/mobilenet/mobilenet.py
    """
    if min_value is None:
        min_value = divisor
    new_v = max(min_value, int(v + divisor / 2) // divisor * divisor)
    # Make sure that round down does not go down by more than 10%.
    if new_v < 0.9 * v:
        new_v += divisor
    return new_v


class ConvolutionBlock(nn.Module):
    def __init__(
        self,
        in_planes: int,
        out_planes: int,
        kernel_size: int = 3,
        stride: int = 1,
        groups: int = 1,
        dilation: int = 1,
        dropout: float = 0.1,
    ):
        super(ConvolutionBlock, self).__init__()

        if isinstance(kernel_size, (tuple, list)):
            padding = ((kernel_size[0] - 1) // 2 * dilation, (kernel_size[1] - 1) // 2 * dilation)
        else:
            padding = (kernel_size - 1) // 2 * dilation

        self.out_channels = out_planes

        self.conv = nn.Conv2d(
            in_planes,
            out_planes,
            kernel_size,
            stride,
            padding,
            dilation=dilation,
            groups=groups,
            bias=False,
        )

        self.norm = nn.BatchNorm2d(out_planes)

        self.activation = nn.SiLU()

        self.dropout = nn.Dropout(dropout)

    def forward(self, images: Tensor) -> Tensor:
        features = self.conv(images)
        features = self.norm(features)
        features = self.activation(features)
        features = self.dropout(features)

        return features


class SqueezeAndExcitation(nn.Module):
    """
    Re-implementation of Squeeze-and-Excitation (SE) block described in:
        *Hu et al., Squeeze-and-Excitation Networks, arXiv:1709.01507*
    """

    __constants__ = ["reduction_ratio"]

    def __init__(self, num_channels: int, reduction_ratio: int = 2):
        super(SqueezeAndExcitation, self).__init__()

        num_channels_reduced = num_channels // reduction_ratio
        self.reduction_ratio = reduction_ratio
        self.fc1 = nn.Linear(num_channels, num_channels_reduced, bias=True)
        self.fc2 = nn.Linear(num_channels_reduced, num_channels, bias=True)

        self.relu = nn.SiLU()

        self.sigmoid = nn.Sigmoid()

    def forward(self, images: Tensor) -> Tensor:
        B, C, H, W = images.shape

        # Average along each channel
        squeeze_tensor = images.view(B, C, -1).mean(dim=2)

        # channel excitation
        fc_out_1 = self.relu(self.fc1(squeeze_tensor))
        fc_out_2 = self.sigmoid(self.fc2(fc_out_1))

        a, b = squeeze_tensor.size()
        output_tensor = torch.mul(images, fc_out_2.view(a, b, 1, 1))

        return output_tensor


class InvertedResidual(nn.Module):
    __constants__ = ["use_res_connect", "out_channels", "_is_cn"]

    def __init__(
        self, input_dim: int, output_dim: int, stride: int, expand_ratio: int, dropout: float = 0.1
    ):
        super(InvertedResidual, self).__init__()

        self.use_res_connect = stride == 1 and input_dim == output_dim
        self.out_channels = output_dim
        self._is_cn = stride > 1
        self.stride = stride

        assert stride in [1, 2]

        hidden_dim = int(round(input_dim * expand_ratio))

        layers: List[nn.Module] = []
        if expand_ratio != 1:
            # Down-sample resnet block
            layers.append(ConvolutionBlock(input_dim, hidden_dim, kernel_size=1, dropout=dropout))

        layers.extend(
            [
                # Base resnet block
                ConvolutionBlock(
                    hidden_dim, hidden_dim, stride=stride, groups=hidden_dim, dropout=dropout
                ),
                # Squeeze and Excitation layer to extend non-linearity
                SqueezeAndExcitation(hidden_dim),
                # Final convolution semi-block to project results back to output dim
                nn.Conv2d(hidden_dim, output_dim, 1, 1, 0, bias=False),
                nn.BatchNorm2d(output_dim),
                nn.Dropout(dropout),
            ]
        )

        self.convolutions = nn.Sequential(*layers)

    def forward(self, images: Tensor) -> Tensor:
        B, C, H, W = images.shape

        if self.use_res_connect:
            return images + self.convolutions(images)
        else:
            return self.convolutions(images)


class MaskedProngMobileNetEmbedding(nn.Module):
    def __init__(
        self,
        input_shape: Tuple[int, int],
        input_dim: int,
        hidden_dim: int,
        width_multiplier: float = 1.0,
        dropout: float = 0.1,
        initial_dimension: int = 32,
        inverted_residual_setting: Optional[List[List[int]]] = None,
        round_nearest: int = 8,
    ):
        """
        MobileNet V2 main class
        Args:
            width_multiplier (float): Width multiplier - adjusts number of channels in each layer by this amount
            inverted_residual_setting: Network structure
            round_nearest (int): Round the number of channels in each layer to be a multiple of this number
            Set to 1 to turn off rounding
        """
        super(MaskedProngMobileNetEmbedding, self).__init__()

        input_channel = initial_dimension
        last_channel = hidden_dim

        if inverted_residual_setting is None:
            inverted_residual_setting = [
                # t, c, n, s
                # expand_ratio, channels, num_layers, stride
                [1, 8, 1, 1],
                [6, 16, 2, 2],
                [6, 32, 3, 2],
                [6, 64, 4, 2],
                [6, 128, 3, 1],
            ]

        # only check the first element, assuming user knows t,c,n,s are required
        if len(inverted_residual_setting) == 0 or len(inverted_residual_setting[0]) != 4:
            raise ValueError(
                "inverted_residual_setting should be non-empty or a 4-element list, got {}".format(
                    inverted_residual_setting
                )
            )

        # building first layer
        self.last_channel = make_divisible_channel_count(
            last_channel * max(1.0, width_multiplier), round_nearest
        )
        input_channel = make_divisible_channel_count(
            input_channel * width_multiplier, round_nearest
        )
        initial_kernel_size = 3
        if input_shape is not None:
            delta = max(input_shape) - min(input_shape)
            initial_kernel_size = (
                (3, 3 + delta) if input_shape[1] > input_shape[0] else (3 + delta, 3)
            )

        layers: List[nn.Module] = [
            ConvolutionBlock(
                input_dim, input_channel, initial_kernel_size, stride=2, dropout=dropout
            )
        ]

        # building inverted residual blocks
        for t, c, n, s in inverted_residual_setting:
            output_channel = make_divisible_channel_count(c * width_multiplier, round_nearest)

            for i in range(n):
                stride = s if i == 0 else 1
                layers.append(
                    InvertedResidual(
                        input_channel, output_channel, stride, expand_ratio=t, dropout=dropout
                    )
                )
                input_channel = output_channel

        # building last several layers
        layers.append(
            ConvolutionBlock(input_channel, self.last_channel, kernel_size=1, dropout=dropout)
        )

        # make it nn.Sequential
        self.resnet = nn.Sequential(*layers)

        # weight initialization
        for m in self.modules():
            if isinstance(m, nn.Conv2d):
                nn.init.kaiming_normal_(m.weight, mode="fan_out")
                if m.bias is not None:
                    nn.init.zeros_(m.bias)
            elif isinstance(m, (nn.BatchNorm2d, nn.GroupNorm)):
                nn.init.ones_(m.weight)
                nn.init.zeros_(m.bias)
            elif isinstance(m, nn.Linear):
                nn.init.normal_(m.weight, 0, 0.01)
                nn.init.zeros_(m.bias)

    def _forward_impl(self, images: Tensor) -> Tensor:
        B, C, H, W = images.shape

        # Sequential layer
        features = self.resnet(images)

        # Cannot use "squeeze" as batch-size can be 1
        features = F.adaptive_avg_pool2d(features, (1, 1))
        features = torch.flatten(features, 1)

        return features

    def forward(self, pixels: Tensor, mask: Tensor) -> Tensor:

        hidden = pixels
        hidden = self._forward_impl(hidden)
        return hidden
