from typing import Optional

import torch
from torch import Tensor, nn
from torch.nn import init


class MaskedBatchNorm0D(nn.Module):
    """A masked version of nn.BatchNorm1d. Only tested for 3D inputs.

    Args:
        num_features: :math:`C` from an expected images of size
            :math:`(N, C, L)`
        eps: a value added to the denominator for numerical stability.
            Default: 1e-5
        momentum: the value used for the running_mean and running_var
            computation. Can be set to ``None`` for cumulative moving average
            (i.e. simple average). Default: 0.1
        affine: a boolean value that when set to ``True``, this module has
            learnable affine parameters. Default: ``True``
        track_running_stats: a boolean value that when set to ``True``, this
            module tracks the running mean and variance, and when set to ``False``,
            this module does not track such statistics and always uses batch
            statistics in both training and eval modes. Default: ``True``

    Shape:
        - Input: :math:`(N, C)`
        - mask: (N) tensor of ones and zeros, where the zeros indicate locations not to use.
        - Output: :math:`(N, C)`
    """

    __constants__ = ["num_features", "eps", "momentum", "affine", "track_running_stats"]

    def __init__(
        self,
        num_features: int,
        eps: float = 1e-5,
        momentum: float = 0.1,
        affine: bool = True,
        track_running_stats: bool = True,
    ):
        super(MaskedBatchNorm0D, self).__init__()

        self.track_running_stats = track_running_stats
        self.num_features = num_features
        self.momentum = momentum
        self.affine = affine
        self.eps = eps

        # Register affine transform learnable parameters
        if affine:
            self.weight = nn.Parameter(torch.Tensor(num_features))
            self.bias = nn.Parameter(torch.Tensor(num_features))
        else:
            self.register_parameter("weight", None)
            self.register_parameter("bias", None)

        # Register moving average storable parameters
        if self.track_running_stats:
            self.register_buffer("running_mean", torch.zeros(num_features))
            self.register_buffer("running_var", torch.ones(num_features))
            self.register_buffer("num_batches_tracked", torch.tensor(0, dtype=torch.long))
        else:
            self.register_parameter("running_mean", None)
            self.register_parameter("running_var", None)
            self.register_parameter("num_batches_tracked", None)

        self.reset_parameters()

    def reset_running_stats(self):
        if self.track_running_stats:
            self.running_mean.zero_()
            self.running_var.fill_(1)
            self.num_batches_tracked.zero_()

    def reset_parameters(self):
        self.reset_running_stats()
        if self.affine:
            init.ones_(self.weight)
            init.zeros_(self.bias)

    def forward(self, images: Tensor, mask: Optional[Tensor] = None) -> Tensor:
        # Calculate the masked mean and variance
        B, C = images.shape

        if mask is not None and mask.shape != (B,):
            raise ValueError("Mask should have shape (B, ).")

        if C != self.num_features:
            raise ValueError(
                "Expected %d channels but images has %d channels" % (self.num_features, C)
            )

        if mask is not None:
            masked_images = images * mask.view(B, 1)
            n = mask.sum()
        else:
            masked_images = images
            n = torch.scalar_tensor(B, dtype=torch.int64)

        # Find the masked sum of the images
        masked_sum = masked_images.sum(dim=0, keepdim=True)

        # Compute masked image statistics
        current_mean = masked_sum / n
        current_var = ((masked_images - current_mean) ** 2).sum(dim=0, keepdim=True) / n

        # Update running stats
        if self.track_running_stats and self.training:
            if self.num_batches_tracked == 0:
                self.running_mean = current_mean
                self.running_var = current_var
            else:
                self.running_mean = (
                    1 - self.momentum
                ) * self.running_mean + self.momentum * current_mean.detach()
                self.running_var = (
                    1 - self.momentum
                ) * self.running_var + self.momentum * current_var.detach()
            self.num_batches_tracked += 1

        # Norm the images
        if self.track_running_stats and not self.training:
            normed_images = (masked_images - self.running_mean) / (
                torch.sqrt(self.running_var + self.eps)
            )
        else:
            normed_images = (masked_images - current_mean) / (torch.sqrt(current_var + self.eps))

        # Apply affine parameters
        if self.affine:
            normed_images = normed_images * self.weight + self.bias

        return normed_images
