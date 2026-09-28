import torch
from torch import Tensor, nn, tensor

from transformercvn.options import Options


class Combiner(nn.Module):
    def __init__(self, options: Options):
        super(Combiner, self).__init__()

        self.hidden_dim = options.hidden_dim
        self.weight_network = nn.Linear(options.hidden_dim, 1)

    @staticmethod
    def masked_softmax(
        x: Tensor, mask: Tensor, dim: int = 1, eps: Tensor = tensor(1e-6, dtype=torch.float)
    ) -> Tensor:

        offset = x.max(dim, keepdim=True).values
        output = torch.exp(x - offset) * mask

        normalizing_sum = output.sum(dim, keepdim=True) + eps
        return output / normalizing_sum

    def forward(self, hidden: Tensor, sequence_mask: Tensor) -> Tensor:
        max_particles, batch_size, hidden_dim = hidden.shape
        sequence_mask = sequence_mask.view(max_particles, batch_size)

        weights = hidden.reshape(-1, hidden_dim)
        weights = self.weight_network(weights).view(max_particles, batch_size)
        weights = self.masked_softmax(weights, sequence_mask, dim=0)
        weights = weights.view(max_particles, batch_size, 1)

        return (weights * hidden).sum(0)
