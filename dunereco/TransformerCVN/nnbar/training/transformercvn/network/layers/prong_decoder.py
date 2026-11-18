import torch
from torch import Tensor, nn

from transformercvn.options import Options


class ProngDecoder(nn.Module):
    def __init__(self, options: Options, output_dim: int, hidden_dim_factor: int = 1):
        super(ProngDecoder, self).__init__()

        self.options = options
        self.hidden_dim_factor = hidden_dim_factor
        self.hidden_layer = nn.Linear(hidden_dim_factor * options.hidden_dim, output_dim)

    def forward(self, hidden: torch.Tensor) -> Tensor:
        return self.hidden_layer(hidden)
