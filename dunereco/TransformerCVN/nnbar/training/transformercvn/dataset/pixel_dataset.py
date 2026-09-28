import torch

from transformercvn.dataset.split_dataset import SplitNeutrinoDataset


class PixelNeutrinoDataset(SplitNeutrinoDataset):
    def __init__(self, data_file: str, limit_index=1.0):
        super(PixelNeutrinoDataset, self).__init__(data_file, limit_index)

        if self.pixels is None:
            raise ValueError(f"{data_file} does not contain pixel information.")

        self.pixel_features = self.pixels.shape[1]
        self.pixel_shape = self.pixels.shape[2:]
        self.current_target_count * torch.tensor([1, 1, 1, 10])

    def __getitem__(self, item):
        return (
            self.data[item],
            self.extra[item],
            self.pixels[item],
            self.mask[item],
            self.current_target[item],
            self.generation_target[item],
        )
