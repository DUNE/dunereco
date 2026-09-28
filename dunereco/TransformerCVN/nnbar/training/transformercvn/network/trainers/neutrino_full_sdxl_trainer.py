from transformercvn.network.networks.neutrino_full_sdxl_network import NeutrinoSDXLNetwork
from transformercvn.network.trainers.neutrino_full_dense_trainer import NeutrinoFullDenseTrainer
from transformercvn.options import Options


class NeutrinoFullSDXLTrainer(NeutrinoFullDenseTrainer):
    def create_network(
        self,
        options: Options,
        features_dim: int,
        extra_dim: int,
        pixel_dim: int,
        num_prong_classes: int,
        num_event_classes: int,
    ) -> NeutrinoSDXLNetwork:
        return NeutrinoSDXLNetwork(
            options, features_dim, extra_dim, pixel_dim, num_prong_classes, num_event_classes
        )
