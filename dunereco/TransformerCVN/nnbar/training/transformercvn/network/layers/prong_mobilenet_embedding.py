from typing import Callable, List, Optional

import torch
from torch import Tensor, nn


def _make_divisible(v: float, divisor: int, min_value: Optional[int] = None) -> int:
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


class ConvBNReLU(nn.Sequential):
    def __init__(
        self,
        in_planes: int,
        out_planes: int,
        kernel_size: int = 3,
        stride: int = 1,
        groups: int = 1,
        norm_layer: Optional[Callable[..., nn.Module]] = None,
        activation_layer: Optional[Callable[..., nn.Module]] = None,
        dilation: int = 1,
    ):
        padding = (kernel_size - 1) // 2 * dilation
        if norm_layer is None:
            norm_layer = nn.BatchNorm2d
        if activation_layer is None:
            activation_layer = nn.ReLU6
        super(ConvBNReLU, self).__init__(
            nn.Conv2d(
                in_planes,
                out_planes,
                kernel_size,
                stride,
                padding,
                dilation=dilation,
                groups=groups,
                bias=False,
            ),
            norm_layer(out_planes),
            activation_layer(inplace=True),
        )
        self.out_channels = out_planes


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
        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()

    def forward(self, input_tensor):
        """
        :param input_tensor: X, shape = (batch_size, num_channels, H, W)
        :return: output tensor
        """
        batch_size, num_channels, H, W = input_tensor.size()
        # Average along each channel
        squeeze_tensor = input_tensor.view(batch_size, num_channels, -1).mean(dim=2)

        # channel excitation
        fc_out_1 = self.relu(self.fc1(squeeze_tensor))
        fc_out_2 = self.sigmoid(self.fc2(fc_out_1))

        a, b = squeeze_tensor.size()
        output_tensor = torch.mul(input_tensor, fc_out_2.view(a, b, 1, 1))
        return output_tensor


class InvertedResidual(nn.Module):
    __constants__ = ["use_res_connect", "out_channels", "_is_cn"]

    def __init__(
        self,
        input_dim: int,
        output_dim: int,
        stride: int,
        expand_ratio: int,
        norm_layer: Optional[Callable[..., nn.Module]] = None,
    ):
        super(InvertedResidual, self).__init__()

        self.stride = stride
        assert stride in [1, 2]

        if norm_layer is None:
            norm_layer = nn.BatchNorm2d

        hidden_dim = int(round(input_dim * expand_ratio))
        self.use_res_connect = self.stride == 1 and input_dim == output_dim

        layers: List[nn.Module] = []
        if expand_ratio != 1:
            layers.append(ConvBNReLU(input_dim, hidden_dim, kernel_size=1, norm_layer=norm_layer))

        layers.extend(
            [
                ConvBNReLU(
                    hidden_dim, hidden_dim, stride=stride, groups=hidden_dim, norm_layer=norm_layer
                ),
                nn.Conv2d(hidden_dim, output_dim, 1, 1, 0, bias=False),
                SqueezeAndExcitation(output_dim),
                norm_layer(output_dim),
            ]
        )

        self.conv = nn.Sequential(*layers)
        self.out_channels = output_dim
        self._is_cn = stride > 1

    def forward(self, x: Tensor) -> Tensor:
        if self.use_res_connect:
            return x + self.conv(x)
        else:
            return self.conv(x)


class ProngMobileNetEmbedding(nn.Module):
    def __init__(
        self,
        input_dim: int,
        hidden_dim: int,
        width_multiplier: float = 1.0,
        inverted_residual_setting: Optional[List[List[int]]] = None,
        round_nearest: int = 8,
        block: Optional[Callable[..., nn.Module]] = None,
        norm_layer: Optional[Callable[..., nn.Module]] = None,
    ):
        """
        MobileNet V2 main class
        Args:
            width_multiplier (float): Width multiplier - adjusts number of channels in each layer by this amount
            inverted_residual_setting: Network structure
            round_nearest (int): Round the number of channels in each layer to be a multiple of this number
            Set to 1 to turn off rounding
            block: Module specifying inverted residual building block for mobilenet
            norm_layer: Module specifying the normalization layer to use
        """
        super(ProngMobileNetEmbedding, self).__init__()

        if block is None:
            block = InvertedResidual

        if norm_layer is None:
            norm_layer = nn.BatchNorm2d

        input_channel = 32
        last_channel = hidden_dim

        if inverted_residual_setting is None:
            inverted_residual_setting = [
                # t, c, n, s
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
        input_channel = _make_divisible(input_channel * width_multiplier, round_nearest)
        self.last_channel = _make_divisible(
            last_channel * max(1.0, width_multiplier), round_nearest
        )
        features: List[nn.Module] = [
            ConvBNReLU(input_dim, input_channel, stride=2, norm_layer=norm_layer)
        ]
        # building inverted residual blocks
        for t, c, n, s in inverted_residual_setting:
            output_channel = _make_divisible(c * width_multiplier, round_nearest)
            for i in range(n):
                stride = s if i == 0 else 1
                features.append(
                    block(
                        input_channel, output_channel, stride, expand_ratio=t, norm_layer=norm_layer
                    )
                )
                input_channel = output_channel
        # building last several layers
        features.append(
            ConvBNReLU(input_channel, self.last_channel, kernel_size=1, norm_layer=norm_layer)
        )

        # make it nn.Sequential
        self.features = nn.Sequential(*features)

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

    def _forward_impl(self, x: Tensor) -> Tensor:
        # This exists since TorchScript doesn't support inheritance, so the superclass method
        # (this one) needs to have a name other than `forward` that can be accessed in a subclass
        x = self.features(x)

        # Cannot use "squeeze" as batch-size can be 1
        x = nn.functional.adaptive_avg_pool2d(x, (1, 1))
        x = torch.flatten(x, 1)

        return x

    def forward(self, pixels: Tensor, mask: Tensor) -> Tensor:
        batch_size, max_particles, channels, height, width = pixels.shape

        pixel_mask = mask.view(batch_size, max_particles, 1)

        # Apply the CNN to every timestep independently
        pixels = pixels.view(batch_size * max_particles, channels, height, width)

        hidden = self._forward_impl(pixels)
        hidden = hidden.reshape(batch_size, max_particles, self.last_channel)

        return hidden * pixel_mask
