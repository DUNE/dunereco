from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Tuple

import torch

if TYPE_CHECKING:
    pass
from torch import Tensor, nn

from transformercvn.network.layers.packed_data import (
    masked_pack_1d_precomputed,
    masked_pad_1d_precomputed,
)
from transformercvn.network.layers.prong_custom_bert_encoder import ProngCustomBertEncoder
from transformercvn.network.layers.prong_decoder import ProngDecoder
from transformercvn.network.layers.prong_feature_embedding import LinearBlock, ProngFeatureEmbedding
from transformercvn.network.layers.prong_masked_mobilenet_embedding import (
    make_divisible_channel_count,
)
from transformercvn.network.layers.prong_target_decoder import ProngTargetDecoder
from transformercvn.options import Options


class BaseProngEmbedding(nn.Module, ABC):
    __constants__ = [
        "hidden_dim",
        "one_hot_pixels",
        "pixel_embedding_dim",
        "feature_embedding_dim",
        "position_embedding_dim",
    ]

    @abstractmethod
    def create_pixel_embedding(self, options: Options, pixel_dim: int, output_dim: int):
        raise NotImplementedError()

    def create_feature_embedding(self, options: Options, features_dim: int, extra_dim: int):
        return ProngFeatureEmbedding(
            options=options,
            sequence_dim=features_dim,
            extra_dim=extra_dim,
            output_dim=self.feature_embedding_dim,
        )

    def __init__(self, options: Options, features_dim: int, extra_dim: int, pixel_dim: int):
        super(BaseProngEmbedding, self).__init__()

        self.hidden_dim = options.hidden_dim
        self.one_hot_pixels = options.one_hot_pixels

        # Split the hidden representation accordingly
        self.pixel_embedding_dim = make_divisible_channel_count(options.pixel_embedding_dim, 8)
        self.feature_embedding_dim = make_divisible_channel_count(options.feature_embedding_dim, 8)
        self.position_embedding_dim = make_divisible_channel_count(
            options.position_embedding_dim, 8
        )

        self.feature_embedding = self.create_feature_embedding(options, features_dim, extra_dim)

        self.prong_pixel_embedding = self.create_pixel_embedding(
            options, pixel_dim, output_dim=self.pixel_embedding_dim
        )

        self.event_pixel_embedding = self.create_pixel_embedding(
            options, pixel_dim, output_dim=self.pixel_embedding_dim + self.feature_embedding_dim
        )

        self.event_position_embedding = nn.Parameter(torch.randn(1, self.position_embedding_dim))

        self.prong_position_embedding = nn.Parameter(torch.randn(1, self.position_embedding_dim))

        self.combined_embedding = LinearBlock(
            options,
            self.feature_embedding_dim + self.pixel_embedding_dim + self.position_embedding_dim,
            options.hidden_dim,
        )

    def forward(
        self,
        features: Tensor,
        extra: Tensor,
        event_pixels: Tensor,
        event_mask: Tensor,
        prong_pixels: Tensor,
        prong_mask: Tensor,
    ) -> Tuple[Tensor, Tensor]:
        batch_size, max_prongs, feature_dim = features.shape

        # Event pixel embeddings run through dedicated CNN and get their own position embedding.
        event_pixel_embeddings = self.event_pixel_embedding(event_pixels)
        event_position_embedding = self.event_position_embedding.expand(batch_size, -1)
        event_embeddings = torch.cat((event_pixel_embeddings, event_position_embedding), dim=1)

        # Prong features and pixels each run through their DNN and add prong position markers to them.
        features, I1, I2 = masked_pack_1d_precomputed(features, prong_mask)
        prong_feature_embeddings = self.feature_embedding(features, extra[I1])
        prong_pixel_embeddings = self.prong_pixel_embedding(prong_pixels)
        prong_position_embedding = self.event_position_embedding.expand(
            prong_pixel_embeddings.shape[0], -1
        )
        prong_embeddings = torch.cat(
            (prong_feature_embeddings, prong_pixel_embeddings, prong_position_embedding), dim=1
        )

        # Combine the two types of embeddings and then through a shared combined embedding.
        combined_embeddings = torch.cat((event_embeddings, prong_embeddings), dim=0)
        combined_embeddings = self.combined_embedding(combined_embeddings)

        # Split the embedding back into their pairs since we need to pad back the prong embeddings.
        event_embeddings, prong_embeddings = (
            combined_embeddings[:batch_size],
            combined_embeddings[batch_size:],
        )
        prong_embeddings = masked_pad_1d_precomputed(
            prong_embeddings, I1, I2, batch_size, max_prongs
        )

        # Add the event embeddings as the first element in the sequence of prongs.
        combined_embeddings = torch.cat(
            (event_embeddings.view(batch_size, 1, -1), prong_embeddings), dim=1
        )
        combined_mask = torch.cat((event_mask, prong_mask), dim=1)

        return combined_embeddings, combined_mask


class NeutrinoBaseNetwork(nn.Module):
    @abstractmethod
    def create_prong_embedding(
        self, options: Options, features_dim: int, extra_dim: int, pixel_dim: int
    ):
        raise NotImplementedError

    def __init__(
        self,
        options: Options,
        features_dim: int,
        extra_dim: int,
        pixel_dim: int,
        num_prong_classes: int,
        num_event_classes: int,
    ):
        super(NeutrinoBaseNetwork, self).__init__()

        self.prong_embedding = self.create_prong_embedding(
            options, features_dim, extra_dim, pixel_dim
        )

        self.encoder = ProngCustomBertEncoder(
            options,
            options.hidden_dim,
            options.num_attention_heads,
            options.dropout,
            options.transformer_activation,
            options.transformer_norm_first,
        )

        self.event_decoder = ProngDecoder(options, num_event_classes)

        self.prong_decoder = ProngTargetDecoder(
            options, options.num_prong_decoder_layers, num_prong_classes
        )

    def forward(
        self,
        features: Tensor,
        extra: Tensor,
        event_pixels: Tensor,
        event_mask: Tensor,
        prong_pixels: Tensor,
        prong_mask: Tensor,
    ) -> Tuple[Tensor, Tensor]:
        combined_embeddings, combined_mask = self.prong_embedding(
            features, extra, event_pixels, event_mask, prong_pixels, prong_mask
        )

        hidden_features, padding_mask, sequence_mask = self.encoder(
            combined_embeddings, combined_mask
        )

        event_features, prong_features = hidden_features[0], hidden_features[1:]

        return self.event_decoder(event_features), self.prong_decoder(prong_features).transpose(
            0, 1
        )
