from transformercvn.network.layers.sdxl_net import SDXLNet
from transformercvn.network.networks.neutrino_full_base_network import (
    BaseProngEmbedding,
    NeutrinoBaseNetwork,
)
from transformercvn.options import Options


class SDXLProngEmbedding(BaseProngEmbedding):
    def create_pixel_embedding(self, options: Options, pixel_dim: int, output_dim: int):
        return SDXLNet(
            input_features=pixel_dim * 256 if self.one_hot_pixels else pixel_dim,
            output_features=output_dim,
            init_block_dim=options.initial_pixel_dim,
            repeat_block_dim=2,
            num_blocks=4,
            norm_num_groups=1,
        )


class NeutrinoSDXLNetwork(NeutrinoBaseNetwork):
    def create_prong_embedding(
        self, options: Options, features_dim: int, extra_dim: int, pixel_dim: int
    ):
        return SDXLProngEmbedding(options, features_dim, extra_dim, pixel_dim)
