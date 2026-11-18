from typing import Tuple

from torch import Tensor, nn

from transformercvn.network.layers.prong_target_decoder import ProngTargetDecoder
from transformercvn.network.networks.neutrino_prong_pixel_network import NeutrinoProngPixelEncoder
from transformercvn.options import Options


class NeutrinoProngTargetNetwork(nn.Module):
    DECODER = ProngTargetDecoder

    def __init__(
        self,
        options: Options,
        features_dim: int,
        pixel_dim: int,
        pixel_shape: Tuple[int, int],
        num_prong_classes: int,
    ):
        super(NeutrinoProngTargetNetwork, self).__init__()

        self.hidden_dim = options.hidden_dim
        self.encoder = NeutrinoProngPixelEncoder(options, features_dim, pixel_dim, pixel_shape)

        self.decoder = self.DECODER(options, options.num_prong_decoder_layers, num_prong_classes)

    def forward(self, features: Tensor, pixels: Tensor, extra: Tensor, mask: Tensor) -> Tensor:
        batch_size, max_prongs, _ = features.shape
        encoded_hidden, padding_mask, sequence_mask = self.encoder(features, pixels, extra, mask)
        return self.decoder(encoded_hidden).transpose(0, 1)
