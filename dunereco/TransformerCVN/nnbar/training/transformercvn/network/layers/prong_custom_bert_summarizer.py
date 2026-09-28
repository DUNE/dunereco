import torch
from torch import Tensor, nn

from transformercvn.options import Options


class MultiHeadPooling(nn.Module):
    def __init__(
        self,
        hidden_dim: int,
        num_heads: int,
        dropout: float = 0.1,
        activation: str = "gelu",
        num_outputs: int = 1,
    ):

        super(MultiHeadPooling, self).__init__()
        self.attention = nn.TransformerDecoderLayer(
            hidden_dim, num_heads, hidden_dim, dropout, activation
        )
        self.pooling_weights = nn.Parameter(torch.Tensor(num_outputs, 1, hidden_dim))
        nn.init.xavier_uniform_(self.pooling_weights)

    def forward(self, hidden: Tensor, padding_mask: Tensor) -> Tensor:
        pooling_weights = self.pooling_weights.repeat(1, hidden.size(1), 1)
        return self.attention(pooling_weights, hidden, memory_key_padding_mask=padding_mask)


class ProngCustomBertSummarizer(nn.Module):
    def __init__(self, options: Options):
        super(ProngCustomBertSummarizer, self).__init__()

        self.pooling = MultiHeadPooling(
            options.hidden_dim,
            options.num_attention_heads,
            options.dropout,
            options.transformer_activation,
        )

        self.dense = nn.Linear(options.hidden_dim, options.hidden_dim)
        self.activation = nn.PReLU(options.hidden_dim)

    def forward(self, hidden_states: Tensor, padding_mask: Tensor, sequence_mask: Tensor) -> Tensor:
        pooled_output = self.pooling(hidden_states, padding_mask)[0]
        pooled_output = self.dense(pooled_output)
        pooled_output = self.activation(pooled_output)
        return pooled_output
