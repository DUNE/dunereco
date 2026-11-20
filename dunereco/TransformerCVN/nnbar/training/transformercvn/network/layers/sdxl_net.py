from diffusers.models.vae import Encoder
from torch import nn


class SDXLNet(nn.Module):
    def __init__(
        self,
        input_features: int,
        output_features: int,
        init_block_dim: int,
        repeat_block_dim: int,
        num_blocks: int,
        norm_num_groups: int = 8,
    ):
        super().__init__()

        block_out_channels = []
        block_dim = init_block_dim
        for _ in range(num_blocks):
            for _ in range(repeat_block_dim):
                block_out_channels.append(block_dim)
            block_dim = block_dim * 2
        block_out_channels.append(output_features)

        self.encoder = Encoder(
            in_channels=input_features,
            out_channels=output_features,
            down_block_types=("DownEncoderBlock2D",) * len(block_out_channels),
            block_out_channels=block_out_channels,
            norm_num_groups=norm_num_groups,
            double_z=False,
        )

        self.output_layer = nn.Sequential(
            nn.Flatten(),
            nn.Linear(output_features, output_features),
        )

    def forward(self, x):
        return self.output_layer(self.encoder(x))
