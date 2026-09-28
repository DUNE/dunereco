import torch
from torch.nn import functional as F

from transformercvn.network.trainers.neutrino_prong_pixel_trainer import NeutrinoProngPixelTrainer


class ElectronProngPixelTrainer(NeutrinoProngPixelTrainer):
    def training_step(self, batch, batch_idx):
        data, pixels, extra, mask, targets = batch

        predictions = self.forward(data, pixels, extra, mask)
        one_hot_targets = F.one_hot(targets).float()
        loss = F.binary_cross_entropy_with_logits(predictions, one_hot_targets)

        self.log("train_loss", loss)

        return loss

    def validation_step(self, batch, batch_idx):
        data, pixels, extra, mask, targets = batch

        predictions = self.forward(data, pixels, extra, mask)
        predictions = torch.sigmoid(predictions) > 0.5

        outputs = {}

        for i, name in enumerate(["numu", "nue", "neutral", "cosmic"]):
            particle_target = (targets == i).cpu().numpy()
            particle_prediction = predictions[:, i].cpu().numpy()

            accuracy = (particle_target == particle_prediction).mean()
            precision = (particle_target & particle_prediction).sum() / particle_target.sum()
            outputs[f"{name}_accuracy"] = accuracy
            outputs[f"{name}_precision"] = precision

            self.log(f"{name}_accuracy", accuracy)
            self.log(f"{name}_precision", precision)

        outputs["val_epoch_accuracy"] = outputs["nue_accuracy"]
        self.log("val_epoch_accuracy", outputs["nue_accuracy"])

        return outputs

    def validation_epoch_end(self, outputs) -> None:
        pass
