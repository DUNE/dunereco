import torch
from torch import Tensor, nn

from transformercvn.network.layers.encoder import create_linear_block
from transformercvn.options import Options


class ProngTargetDecoder(nn.Module):
    __constants__ = ["output_dim"]

    def __init__(self, options: Options, num_hidden: int, output_dim: int):
        super(ProngTargetDecoder, self).__init__()

        self.output_dim = output_dim

        self.hidden_layers, final_dimension = self.create_decoder_layers(options, num_hidden)
        self.output_layer = nn.Linear(final_dimension, output_dim)

    @staticmethod
    def create_decoder_layers(options: Options, num_layers: int):
        current_hidden_dim = options.hidden_dim
        hidden_layers = []

        for i in range(num_layers):
            next_hidden_dim = current_hidden_dim // 2
            if next_hidden_dim < 8:
                break

            hidden_layers.extend(create_linear_block(current_hidden_dim, next_hidden_dim, options))
            current_hidden_dim = next_hidden_dim

        return nn.Sequential(*hidden_layers), next_hidden_dim

    def forward(self, hidden: torch.Tensor) -> Tensor:
        timesteps, batch_size, hidden_dim = hidden.shape

        hidden = hidden.reshape(timesteps * batch_size, hidden_dim)
        hidden = self.hidden_layers(hidden)
        hidden = self.output_layer(hidden)

        return hidden.reshape(timesteps, batch_size, self.output_dim)
