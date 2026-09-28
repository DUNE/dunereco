import torch
from torch import Tensor, jit, nn

from transformercvn.network.layers.encoder import create_linear_block
from transformercvn.options import Options


class Decoder(nn.Module):
    def __init__(self, options: Options, output_dim: int, hidden_dim_factor: int = 1):
        super(Decoder, self).__init__()

        self.options = options
        self.hidden_dim_factor = hidden_dim_factor
        self.hidden_layers = jit.script(self.create_decoder_layers(output_dim))

    def create_decoder_layers(self, output_dim: int):
        current_hidden_dim = self.hidden_dim_factor * self.options.hidden_dim
        hidden_layers = []

        for i in range(self.options.num_decoder_layers):
            next_hidden_dim = current_hidden_dim // 2
            if next_hidden_dim < self.options.final_decoder_dim:
                break

            hidden_layers.extend(
                create_linear_block(current_hidden_dim, next_hidden_dim, self.options)
            )
            current_hidden_dim = next_hidden_dim

        hidden_layers.append(nn.Linear(current_hidden_dim, output_dim))
        return nn.Sequential(*hidden_layers)

    def forward(self, hidden: torch.Tensor) -> Tensor:
        return self.hidden_layers(hidden)
