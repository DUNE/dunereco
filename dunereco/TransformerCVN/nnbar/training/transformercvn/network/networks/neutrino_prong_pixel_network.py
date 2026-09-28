from typing import Optional, Tuple

import torch
from torch import Tensor, nn

from transformercvn.network.layers.packed_data import (
    masked_pack_1d_precomputed,
    masked_pack_3d_precomputed,
    masked_pad_1d_precomputed,
)
from transformercvn.network.layers.prong_custom_bert_encoder import ProngCustomBertEncoder
from transformercvn.network.layers.prong_custom_bert_summarizer import ProngCustomBertSummarizer
from transformercvn.network.layers.prong_decoder import ProngDecoder
from transformercvn.network.layers.prong_feature_embedding import ProngFeatureEmbedding
from transformercvn.network.layers.prong_masked_mobilenet_embedding import (
    MaskedProngMobileNetEmbedding,
    make_divisible_channel_count,
)
from transformercvn.options import Options


class NeutrinoProngPixelEncoder(nn.Module):
    FEATURE_EMBEDDING = ProngFeatureEmbedding

    PIXEL_EMBEDDING = MaskedProngMobileNetEmbedding

    ENCODER = ProngCustomBertEncoder

    __constants__ = ["hidden_dim", "cnn_hidden_dim", "feature_hidden_dim"]

    def __init__(
        self,
        options: Options,
        features_dim: int,
        pixel_dim: int,
        pixel_shape: Optional[Tuple[int, int]] = None,
    ):
        super(NeutrinoProngPixelEncoder, self).__init__()

        self.hidden_dim = options.hidden_dim

        # Split the hidden representation accordingly
        self.cnn_hidden_dim = make_divisible_channel_count(
            options.hidden_dim * options.cnn_embedding_proportion, 8
        )
        self.feature_hidden_dim = options.hidden_dim - self.cnn_hidden_dim

        self.feature_embedding = self.FEATURE_EMBEDDING(
            options=options,
            sequence_dim=features_dim,
            extra_dim=1,
            output_dim=self.feature_hidden_dim,
        )

        self.prong_pixel_embedding = self.PIXEL_EMBEDDING(
            input_shape=pixel_shape,
            input_dim=pixel_dim,
            hidden_dim=self.cnn_hidden_dim,
            dropout=options.dropout,
            inverted_residual_setting=options.mobilenet_structure,
            initial_dimension=options.initial_pixel_dim,
        )

        self.prong_encoder = self.ENCODER(
            options,
            options.hidden_dim,
            options.num_attention_heads,
            options.dropout,
            options.transformer_activation,
        )

    def forward(
        self, features: Tensor, pixels: Tensor, extra: Tensor, mask: Tensor
    ) -> Tuple[Tensor, Tensor, Tensor]:
        batch_size, max_prongs, _ = features.shape

        features, I1, I2 = masked_pack_1d_precomputed(features, mask)
        pixels, I13, I23 = masked_pack_3d_precomputed(pixels, mask)
        extra = extra[I1]

        embedded_features = self.feature_embedding(features, extra, mask)
        embedded_pixels = self.prong_pixel_embedding(pixels, mask)

        embeddings = torch.cat([embedded_features, embedded_pixels], dim=1)
        embeddings = masked_pad_1d_precomputed(embeddings, I1, I2, batch_size, max_prongs)
        encoded_hidden, padding_mask, sequence_mask = self.prong_encoder(embeddings, mask)
        return encoded_hidden, padding_mask, sequence_mask


class NeutrinoProngPixelNetwork(nn.Module):
    SUMMARIZER = ProngCustomBertSummarizer

    DECODER = ProngDecoder

    __constants__ = ["hidden_dim", "cnn_hidden_dim", "feature_hidden_dim"]

    def __init__(self, options: Options, features_dim: int, pixel_dim: int, num_classes: int):
        super(NeutrinoProngPixelNetwork, self).__init__()

        self.hidden_dim = options.hidden_dim
        self.encoder = NeutrinoProngPixelEncoder(options, features_dim, pixel_dim)

        self.summarizer = self.SUMMARIZER(options)

        self.decoder = self.DECODER(options, num_classes)

    def forward(self, features: Tensor, pixels: Tensor, extra: Tensor, mask: Tensor) -> Tensor:
        batch_size, max_prongs, _ = features.shape
        encoded_hidden, padding_mask, sequence_mask = self.encoder(features, pixels, extra, mask)
        summarized_hidden = self.summarizer(encoded_hidden, padding_mask, sequence_mask)
        return self.decoder(summarized_hidden)
