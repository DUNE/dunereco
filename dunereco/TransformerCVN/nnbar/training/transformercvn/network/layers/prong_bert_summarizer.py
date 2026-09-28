from torch import Tensor, nn

from transformercvn.options import Options


class ProngBertSummarizer(nn.Module):
    def __init__(self, options: Options):
        super(ProngBertSummarizer, self).__init__()

        self.dense = nn.Linear(options.hidden_dim, options.hidden_dim)
        self.activation = nn.Tanh()

    def forward(self, hidden_states: Tensor, sequence_mask: Tensor) -> Tensor:
        # We "pool" the model by simply taking the hidden state corresponding
        # to the first token.
        first_token_tensor = hidden_states[:, 0]
        pooled_output = self.dense(first_token_tensor)
        pooled_output = self.activation(pooled_output)
        return pooled_output
