from typing import Tuple

import torch
from torch import Tensor, jit, nn

from transformercvn.options import Options


class ProngEncoder(nn.Module):
    def __init__(
        self, options: Options, hidden_dim: int, num_heads: int, dropout: float, activation: str
    ):
        super(ProngEncoder, self).__init__()

        self.options = options

        encoder_layer = self.encoder_layer(hidden_dim, num_heads, hidden_dim, dropout, activation)
        encoder = nn.TransformerEncoder(encoder_layer, options.num_encoder_layers)
        self.encoder = jit.script(encoder)

    @property
    def encoder_layer(self):
        return nn.TransformerEncoderLayer

    def forward(
        self, hidden_features: Tensor, hidden_pixels: Tensor, mask: Tensor
    ) -> Tuple[Tensor, Tensor, Tensor]:
        batch_size, max_particles, _ = hidden_features.shape

        # Setup variants of mask
        padding_mask = ~mask
        sequence_mask = mask.view(batch_size, max_particles, 1).transpose(0, 1).contiguous()

        # Combine the different layers into a single vector
        hidden_sequence = torch.cat([hidden_features, hidden_pixels], dim=2)

        # Reshape vector to have time axis first for transformer images
        hidden_sequence = hidden_sequence.transpose(0, 1).contiguous() * sequence_mask

        # Primary transformer encoder
        hidden_sequence = (
            self.encoder(hidden_sequence, src_key_padding_mask=padding_mask) * sequence_mask
        )

        return hidden_sequence, padding_mask, sequence_mask
