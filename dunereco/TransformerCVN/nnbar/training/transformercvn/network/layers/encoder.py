from typing import Tuple

import torch
from torch import Tensor, jit, nn

from transformercvn.options import Options


def create_linear_block(input_dim, output_dim, options: Options):
    layers = [nn.Linear(input_dim, output_dim)]

    if options.linear_batch_norm:
        layers.append(nn.BatchNorm1d(output_dim))

    if options.linear_prelu_activation:
        layers.append(nn.PReLU(output_dim))
    else:
        layers.append(nn.ReLU())

    if options.dropout > 0.0:
        layers.append(nn.Dropout(options.dropout))

    return layers


class Encoder(nn.Module):
    def __init__(
        self,
        options: Options,
        input_dim: int,
        transformer_options: Tuple[int, int, int, float, str],
    ):
        super(Encoder, self).__init__()

        self.options = options
        self.embedding = jit.script(self.create_embedding_layers(options, input_dim))
        self.embedding_dim = self.options.hidden_dim

        encoder = nn.TransformerEncoder(
            self.encoder_layer(*transformer_options), options.num_encoder_layers
        )
        self.encoder = jit.script(encoder)

    @property
    def encoder_layer(self):
        return nn.TransformerEncoderLayer

    @staticmethod
    def create_embedding_layers(options, input_dim):
        embedding_layers = create_linear_block(input_dim, options.initial_dimension, options)
        current_embedding_dim = options.initial_dimension

        for i in range(options.num_embedding_layers):
            next_embedding_dim = 2 * current_embedding_dim
            if next_embedding_dim >= options.hidden_dim:
                break

            embedding_layers.extend(
                create_linear_block(current_embedding_dim, next_embedding_dim, options)
            )
            current_embedding_dim = next_embedding_dim

        embedding_layers.extend(
            create_linear_block(current_embedding_dim, options.hidden_dim - 1, options)
        )

        return nn.Sequential(*embedding_layers)

    def forward(self, data: Tensor, extra: Tensor, mask: Tensor) -> Tuple[Tensor, Tensor, Tensor]:
        batch_size, max_particles, input_dim = data.shape

        # Setup variants of mask
        padding_mask = ~mask

        sequence_mask = mask.view(batch_size, max_particles, 1).transpose(0, 1).contiguous()

        # Perform embedding on all of the vectors uniformly
        hidden = self.embedding(data.view(-1, input_dim))
        hidden = hidden.view(batch_size, max_particles, self.embedding_dim - 1)

        # Append the final extra dimension to the hidden layer
        extra = extra.reshape(batch_size, 1, 1)
        extra = extra.repeat(1, max_particles, 1)
        hidden = torch.cat((hidden, extra), dim=-1)

        # Reshape vector to have time axis first for transformer images
        hidden = hidden.transpose(0, 1) * sequence_mask

        # Primary transformer encoder
        hidden = self.encoder(hidden, src_key_padding_mask=padding_mask) * sequence_mask

        return hidden, padding_mask, sequence_mask
