import torch
from torch import Tensor, nn

from transformercvn.options import Options


class LinearBlock(nn.Module):
    def __init__(self, options: Options, input_dim: int, output_dim: int):
        super(LinearBlock, self).__init__()

        self.linear = nn.Linear(input_dim, output_dim, bias=not options.linear_batch_norm)

        if options.linear_batch_norm:
            self.norm = nn.BatchNorm1d(output_dim)
        else:
            self.norm = nn.Identity()

        if options.linear_prelu_activation:
            self.activation = nn.PReLU(output_dim)
        else:
            self.activation = nn.ReLU()

        self.dropout = nn.Dropout(options.dropout)

    def forward(self, data: Tensor) -> Tensor:
        B, C = data.shape

        features = self.linear(data)
        features = self.norm(features)
        features = self.activation(features)
        features = self.dropout(features)

        return features


class ProngFeatureEmbedding(nn.Module):
    __constants__ = [
        "extra_dim",
        "output_dim",
        "sequence_dim",
        "embedding_dim",
        "disable_smart_features",
    ]

    def __init__(self, options: Options, sequence_dim: int, extra_dim: int, output_dim: int):
        super(ProngFeatureEmbedding, self).__init__()

        self.extra_dim = extra_dim
        self.output_dim = output_dim
        self.sequence_dim = sequence_dim
        self.disable_smart_features = options.disable_smart_features

        self.embedding = self.create_embedding_layers(options, sequence_dim + extra_dim, output_dim)
        self.embedding_dim = options.hidden_dim

    @staticmethod
    def create_embedding_layers(options, input_dim, output_dim):
        embedding_layers = [LinearBlock(options, input_dim, options.initial_feature_dim)]
        current_embedding_dim = options.initial_feature_dim

        for i in range(options.num_embedding_layers):
            next_embedding_dim = 2 * current_embedding_dim
            if next_embedding_dim >= output_dim:
                break

            embedding_layers.append(LinearBlock(options, current_embedding_dim, next_embedding_dim))
            current_embedding_dim = next_embedding_dim

        embedding_layers.append(LinearBlock(options, current_embedding_dim, output_dim))

        return nn.Sequential(*embedding_layers)

    def forward(self, data: Tensor, extra: Tensor) -> Tensor:
        packed_size, feature_dim = data.shape
        if self.disable_smart_features:
            return torch.zeros(packed_size, self.output_dim, dtype=data.dtype, device=data.device)
        else:
            return self.embedding(torch.cat([data, extra], dim=1))
