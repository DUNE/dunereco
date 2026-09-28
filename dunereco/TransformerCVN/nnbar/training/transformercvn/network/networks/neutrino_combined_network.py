from typing import Optional, Tuple

import torch
from torch import Tensor, nn
from torch.nn import functional as F

from transformercvn.network.layers.packed_data import (
    masked_pack_1d_precomputed,
    masked_pack_3d_precomputed,
    masked_pad_1d_precomputed,
)
from transformercvn.network.layers.prong_custom_bert_encoder import ProngCustomBertEncoder
from transformercvn.network.layers.prong_decoder import ProngDecoder
from transformercvn.network.layers.prong_feature_embedding import LinearBlock, ProngFeatureEmbedding
from transformercvn.network.layers.prong_masked_mobilenet_embedding import (
    MaskedProngMobileNetEmbedding,
    make_divisible_channel_count,
)
from transformercvn.network.layers.prong_target_decoder import ProngTargetDecoder
from transformercvn.options import Options


class NeutrinoProngEmbedding(nn.Module):
    __constants__ = [
        "hidden_dim",
        "cnn_hidden_dim",
        "feature_hidden_dim",
        "disable_smart_features",
        "position_dim",
        "one_hot_pixels",
    ]

    def __init__(
        self,
        options: Options,
        features_dim: int,
        pixel_dim: int,
        pixel_shape: Optional[Tuple[int, int]] = None,
    ):
        super(NeutrinoProngEmbedding, self).__init__()

        self.one_hot_pixels = options.one_hot_pixels
        self.hidden_dim = options.hidden_dim
        self.position_dim = options.position_embedding_dim

        # Split the hidden representation accordingly
        self.cnn_hidden_dim = make_divisible_channel_count(
            options.hidden_dim * options.cnn_embedding_proportion, 8
        )
        self.cnn_hidden_dim = min(self.cnn_hidden_dim, ((options.hidden_dim // 8) - 1) * 8)
        self.feature_hidden_dim = options.hidden_dim - self.cnn_hidden_dim
        self.disable_smart_features = options.disable_smart_features

        self.feature_embedding = ProngFeatureEmbedding(
            options=options,
            sequence_dim=features_dim,
            extra_dim=1,
            output_dim=self.feature_hidden_dim,
        )

        self.prong_pixel_embedding = MaskedProngMobileNetEmbedding(
            input_shape=pixel_shape,
            input_dim=pixel_dim * 256 if self.one_hot_pixels else pixel_dim,
            hidden_dim=self.cnn_hidden_dim,
            dropout=options.dropout,
            inverted_residual_setting=options.mobilenet_structure,
            initial_dimension=options.initial_pixel_dim,
        )

        self.contextual_position = nn.Parameter(torch.randn(1, options.position_embedding_dim))

        self.prong_position = nn.Parameter(torch.randn(1, options.position_embedding_dim))

        self.position_embedding = LinearBlock(
            options, options.hidden_dim + options.position_embedding_dim, options.hidden_dim
        )

    def forward(self, features: Tensor, pixels: Tensor, extra: Tensor, mask: Tensor) -> Tensor:
        batch_size, max_prongs, _ = features.shape

        features, I1, I2 = masked_pack_1d_precomputed(features, mask)
        pixels, _, _ = masked_pack_3d_precomputed(pixels, mask)
        extra = extra[I1]

        if self.one_hot_pixels:
            pixels = F.one_hot(pixels.long(), 256)
            pixels = pixels.permute(0, 1, 4, 2, 3)
            pixels = pixels.reshape(pixels.shape[0], -1, pixels.shape[3], pixels.shape[4])
            pixels = pixels.float()

        embedded_features = self.feature_embedding(features, extra)
        embedded_pixels = self.prong_pixel_embedding(pixels, mask)

        if self.disable_smart_features:
            embedded_features = embedded_features * 0.0

        embeddings = torch.cat([embedded_features, embedded_pixels], dim=1)

        # Construct Position encoding vectors on packed data
        contextual_position = (I2 == 0).unsqueeze(1)
        prong_position = ~contextual_position

        contextual_position = self.contextual_position * contextual_position
        prong_position = self.prong_position * prong_position

        position = contextual_position + prong_position

        # Add position information to embeddings and run through another layer
        embeddings = torch.cat((embeddings, position), dim=1)
        embeddings = self.position_embedding(embeddings)

        # Pad back into zero padded array for transformer
        embeddings = masked_pad_1d_precomputed(embeddings, I1, I2, batch_size, max_prongs)

        return embeddings


class ClassifierProng(nn.Module):
    def __init__(self, hidden_dim: int):
        super(ClassifierProng, self).__init__()

        self.classifier_embedding = torch.nn.Parameter(torch.randn(1, 1, hidden_dim))

    def forward(self, embeddings: Tensor, mask: Tensor) -> Tuple[Tensor, Tensor]:
        batch_size, timesteps, hidden_dim = embeddings.shape

        new_embeddings = torch.cat(
            (self.classifier_embedding.expand(batch_size, 1, hidden_dim), embeddings), dim=1
        )
        new_mask = torch.cat(
            (torch.ones(batch_size, 1, dtype=mask.dtype, device=mask.device), mask), dim=1
        )

        return new_embeddings, new_mask


class NeutrinoCombinedNetwork(nn.Module):
    __constants__ = ["disable_smart_features"]

    def __init__(
        self,
        options: Options,
        features_dim: int,
        pixel_dim: int,
        pixel_shape: Tuple[int, int],
        num_prong_classes: int,
        num_event_classes: int,
    ):
        super(NeutrinoCombinedNetwork, self).__init__()

        self.disable_smart_features = options.disable_smart_features

        self.embedding = NeutrinoProngEmbedding(options, features_dim, pixel_dim, pixel_shape)
        self.classifier_embedding = ClassifierProng(options.hidden_dim)

        self.encoder = ProngCustomBertEncoder(
            options,
            options.hidden_dim,
            options.num_attention_heads,
            options.dropout,
            options.transformer_activation,
        )

        self.event_decoder = ProngDecoder(options, num_event_classes)

        self.prong_decoder = ProngTargetDecoder(
            options, options.num_prong_decoder_layers, num_prong_classes
        )

    def forward(
        self, features: Tensor, pixels: Tensor, extra: Tensor, mask: Tensor
    ) -> Tuple[Tensor, Tensor]:
        if self.disable_smart_features:
            features = features * 0.0
            extra = extra * 0.0

        embeddings = self.embedding(features, pixels, extra, mask)
        embeddings, mask = self.classifier_embedding(embeddings, mask)
        hidden_features, padding_mask, sequence_mask = self.encoder(embeddings, mask)

        classification_features, prong_features = hidden_features[0], hidden_features[1:]

        return self.event_decoder(classification_features), self.prong_decoder(
            prong_features
        ).transpose(0, 1)
