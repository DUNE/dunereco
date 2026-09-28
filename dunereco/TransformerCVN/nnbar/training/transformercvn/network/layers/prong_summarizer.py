from torch import Tensor, nn

from transformercvn.options import Options


class ProngSummarizer(nn.Module):
    def __init__(self, options: Options):
        super(ProngSummarizer, self).__init__()

    def forward(self, hidden_sequence: Tensor, sequence_mask: Tensor) -> Tensor:
        return hidden_sequence[0]
