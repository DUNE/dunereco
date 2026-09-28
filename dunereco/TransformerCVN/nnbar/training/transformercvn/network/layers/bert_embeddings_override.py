from typing import Optional

import torch
import transformers
from torch import Tensor, nn


class BertEmbeddings(nn.Module):
    """A simplified BertEmbeddings that avoids making extra embedding layers."""

    def __init__(self, config):
        super().__init__()
        # self.LayerNorm is not snake-cased to stick with TensorFlow model variable name and be able to load
        # any TensorFlow checkpoint file
        self.LayerNorm = nn.LayerNorm(config.hidden_size, eps=config.layer_norm_eps)
        self.dropout = nn.Dropout(config.hidden_dropout_prob)

        # position_ids (1, len position emb) is contiguous in memory and exported when serialized
        self.register_buffer(
            "position_ids", torch.arange(config.max_position_embeddings).expand((1, -1))
        )
        self.position_embedding_type = getattr(config, "position_embedding_type", "absolute")

        self.hidden_size = config.hidden_size

    def forward(
        self,
        input_ids: Optional[Tensor] = None,
        token_type_ids: Optional[Tensor] = None,
        position_ids: Optional[Tensor] = None,
        inputs_embeds: Optional[Tensor] = None,
        past_key_values_length: int = 0,
    ):
        embeddings = inputs_embeds
        embeddings = self.LayerNorm(embeddings)
        embeddings = self.dropout(embeddings)
        return embeddings


def override_bert_embeddings():
    transformers.models.bert.modeling_bert.BertEmbeddings = BertEmbeddings
