from typing import Optional, Tuple

import torch
from torch import Tensor, nn

from transformercvn.options import Options


class InducedSetAttentionBlock(nn.Module):
    def __init__(
        self,
        hidden_dim: int,
        num_heads: int,
        num_indices: int,
        dropout: float = 0.1,
        activation: str = "gelu",
    ):
        super(InducedSetAttentionBlock, self).__init__()

        self.induction_weights = nn.Parameter(torch.Tensor(num_indices, 1, hidden_dim))
        self.attention_1 = nn.TransformerDecoderLayer(
            hidden_dim, num_heads, hidden_dim, dropout, activation
        )
        self.attention_2 = nn.TransformerDecoderLayer(
            hidden_dim, num_heads, hidden_dim, dropout, activation
        )

        nn.init.xavier_uniform_(self.induction_weights)

    def forward(
        self,
        src: Tensor,
        src_mask: Optional[Tensor] = None,
        src_key_padding_mask: Optional[Tensor] = None,
    ) -> Tensor:

        weights = self.induction_weights.repeat(1, src.size(1), 1)
        hidden = self.attention_1(
            weights, src, memory_mask=src_mask, memory_key_padding_mask=src_key_padding_mask
        )
        return self.attention_2(
            src, hidden, tgt_mask=src_mask, tgt_key_padding_mask=src_key_padding_mask
        )


class ProngCustomBertEncoder(nn.Module):
    def __init__(
        self,
        options: Options,
        hidden_dim: int,
        num_heads: int,
        dropout: float,
        activation: str,
        norm_first: bool,
    ):
        super(ProngCustomBertEncoder, self).__init__()

        self.options = options

        encoder_layer = nn.TransformerEncoderLayer(
            hidden_dim, num_heads, hidden_dim, dropout, activation, norm_first=norm_first
        )

        self.encoder = nn.TransformerEncoder(encoder_layer, options.num_encoder_layers)

    def forward(self, embeddings: Tensor, mask: Tensor) -> Tuple[Tensor, Tensor, Tensor]:
        batch_size, max_particles, _ = embeddings.shape

        # Setup variants of mask
        padding_mask = ~mask
        sequence_mask = mask.view(batch_size, max_particles, 1).transpose(0, 1).contiguous()

        # Combine the different layers into a single vector

        # Reshape vector to have time axis first for transformer images
        hidden_sequence = embeddings.transpose(0, 1).contiguous() * sequence_mask

        # Primary transformer encoder
        hidden_sequence = (
            self.encoder(hidden_sequence, src_key_padding_mask=padding_mask) * sequence_mask
        )

        return hidden_sequence, padding_mask, sequence_mask
