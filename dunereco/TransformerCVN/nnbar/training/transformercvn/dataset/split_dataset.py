import torch

from transformercvn.dataset.dataset import NeutrinoDataset


class SplitNeutrinoDataset(NeutrinoDataset):
    def __init__(self, data_file: str, limit_index=1.0):
        super(SplitNeutrinoDataset, self).__init__(data_file, limit_index)

        self.current_target = torch.zeros_like(self.targets)
        self.current_target[(self.targets > 3) & (self.targets <= 7)] = 1
        self.current_target[self.targets == 8] = 2
        self.current_target[self.targets == 9] = 3

        self.generation_target = torch.zeros_like(self.targets)
        self.generation_target[(self.targets == 0) | (self.targets == 4)] = 0
        self.generation_target[(self.targets == 1) | (self.targets == 5)] = 1
        self.generation_target[(self.targets == 2) | (self.targets == 6)] = 2
        self.generation_target[(self.targets == 3) | (self.targets == 7)] = 3

        self.generation_valid = self.targets < 8

        self.num_current_classes = 4
        self.num_generation_classes = 4

        self.current_target_count = torch.bincount(self.current_target)
        self.generation_target_count = torch.bincount(self.generation_target[self.generation_valid])

    def __getitem__(self, item):
        return (
            self.data[item],
            self.extra[item],
            self.mask[item],
            self.current_target[item],
            self.generation_target[item],
        )
