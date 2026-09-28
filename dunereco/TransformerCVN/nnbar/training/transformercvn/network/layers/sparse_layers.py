from typing import List, Tuple, Union

import MinkowskiEngine as ME
import torch
from torch import Tensor, nn

SparseTensor = ME.SparseTensor


def SparseCat(tensors: Union[List[Tensor], Tuple[Tensor]]) -> Tensor:
    if len(tensors) > 1:
        return ME.cat(*tensors)
    else:
        return tensors[0]


class SparseConcatenation(nn.Module):
    def forward(self, *tensors):
        return ME.cat(*tensors)


def SparseConv2D(
    in_channels,
    out_channels,
    kernel_size=-1,
    stride=1,
    dilation=1,
    bias=False,
    expand_coordinates=True,
):
    return ME.MinkowskiConvolution(
        in_channels=in_channels,
        out_channels=out_channels,
        kernel_size=kernel_size,
        stride=stride,
        dilation=dilation,
        bias=bias,
        dimension=2,
        expand_coordinates=expand_coordinates,
    )


def SparseChannelConv2D(
    in_channels, kernel_size=-1, stride=1, dilation=1, bias=False, expand_coordinates=True
):
    return ME.MinkowskiChannelwiseConvolution(
        in_channels=in_channels,
        kernel_size=kernel_size,
        stride=stride,
        dilation=dilation,
        bias=bias,
        dimension=2,
        expand_coordinates=expand_coordinates,
    )


def SparseBatchNorm2D(num_features, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True):
    return ME.MinkowskiBatchNorm(
        num_features=num_features,
        eps=eps,
        momentum=momentum,
        affine=affine,
        track_running_stats=track_running_stats,
    )


class SparseLayerNorm(nn.Module):
    r"""A batch normalization layer for a sparse tensor.
    See the pytorch :attr:`torch.nn.BatchNorm1d` for more details.
    """

    def __init__(self, num_features, eps=1e-5, elementwise_affine: bool = True):
        super(SparseLayerNorm, self).__init__()
        self.ln = torch.nn.LayerNorm(num_features, eps, elementwise_affine)

    def forward(self, x):
        return SparseTensor(
            self.ln(x.F),
            coordinate_map_key=x.coordinate_map_key,
            coordinate_manager=x.coordinate_manager,
        )

    def __repr__(self):
        s = "({}, eps={}, elementwise_affine={})".format(
            self.ln.num_features, self.ln.eps, self.ln.elementwise_affine
        )

        return self.__class__.__name__ + s


def SparsePReLU(num_parameters=1, init=0.25):
    return ME.MinkowskiPReLU(num_parameters=num_parameters, init=init)


def SparseGeLU():
    return ME.MinkowskiGELU()


def SparseDropout(p=0.5, inplace=False):
    return ME.MinkowskiDropout(p=p, inplace=inplace)


def SparseAvgPooling2D(kernel_size=-1, stride=1, dilation=1):
    return ME.MinkowskiAvgPooling(
        kernel_size=kernel_size, stride=stride, dilation=dilation, dimension=2
    )


def SparseGlobalAvgPooling2D():
    return ME.MinkowskiGlobalAvgPooling()


class SparseCondense(nn.Module):
    def __init__(self):
        super(SparseCondense, self).__init__()

        self.pooling = ME.MinkowskiGlobalAvgPooling()

    def forward(self, features):
        output = self.pooling(features)
        return output.F[torch.argsort(output.C[:, 0])].contiguous()


def drop_path(x, drop_prob: float = 0.0, training: bool = False, scale_by_keep: bool = True):
    if drop_prob == 0.0 or not training:
        return x
    keep_prob = 1 - drop_prob
    shape = (x.shape[0],) + (1,) * (x.ndim - 1)  # work with diff dim tensors, not just 2D ConvNets
    random_tensor = x.new_empty(shape).bernoulli_(keep_prob)
    if keep_prob > 0.0 and scale_by_keep:
        random_tensor.div_(keep_prob)
    return x * random_tensor


class DropPath(nn.Module):
    def __init__(self, drop_prob=None, scale_by_keep=True):
        super(DropPath, self).__init__()
        self.drop_prob = drop_prob
        self.scale_by_keep = scale_by_keep

    def forward(self, x):
        return drop_path(x, self.drop_prob, self.training, self.scale_by_keep)


class SparseDropPath(nn.Module):
    r"""A batch normalization layer for a sparse tensor.
    See the pytorch :attr:`torch.nn.BatchNorm1d` for more details.
    """

    def __init__(self, drop_prob=None, scale_by_keep=True):
        super(SparseDropPath, self).__init__()
        self.drop_path = DropPath(drop_prob, scale_by_keep)

    def forward(self, x):
        return SparseTensor(
            self.drop_path(x.F),
            coordinate_map_key=x.coordinate_map_key,
            coordinate_manager=x.coordinate_manager,
        )
