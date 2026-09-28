from typing import Optional

from torch import Tensor, nn


class MaskedSequential(nn.Module):
    def __init__(self, *layers):
        super(MaskedSequential, self).__init__()

        self.layers = nn.ModuleList(layers)

    def forward(self, data: Tensor, mask: Optional[Tensor] = None) -> Tensor:
        output = data
        for layer in self.layers:
            output = layer(output, mask)

        return output
