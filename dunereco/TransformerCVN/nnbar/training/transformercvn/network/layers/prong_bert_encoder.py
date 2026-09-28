from typing import Tuple

import torch
from torch import Tensor, nn
from transformers.models.bert import BertConfig, BertModel

from transformercvn.network.layers.bert_embeddings_override import override_bert_embeddings
from transformercvn.options import Options

override_bert_embeddings()


class ProngBertEncoder(nn.Module):
    def __init__(
        self, options: Options, hidden_dim: int, num_heads: int, dropout: float, activation: str
    ):
        super(ProngBertEncoder, self).__init__()

        self.options = options
        self.encoder_config = BertConfig(
            vocab_size=1,
            hidden_size=hidden_dim,
            num_hidden_layers=options.num_encoder_layers,
            num_attention_heads=num_heads,
            intermediate_size=hidden_dim,
            hidden_act=activation,
            hidden_dropout_prob=dropout,
            attention_probs_dropout_prob=dropout,
            position_embedding_type="none",
        )

        self.encoder = BertModel(self.encoder_config, add_pooling_layer=False)

    def forward(
        self, hidden_features: Tensor, hidden_pixels: Tensor, mask: Tensor
    ) -> Tuple[Tensor, Tensor, Tensor]:
        batch_size, max_particles, hidden_dim = hidden_features.shape

        # Combine the different layers into a single vector
        hidden_sequence = torch.cat([hidden_features, hidden_pixels], dim=2)

        # Primary transformer encoder
        hidden_sequence = self.encoder(
            inputs_embeds=hidden_sequence, attention_mask=mask, return_dict=False
        )[0]

        return hidden_sequence, mask, mask
