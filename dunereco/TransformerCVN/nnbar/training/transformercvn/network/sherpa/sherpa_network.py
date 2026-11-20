from typing import Optional

from sherpa import Client, Trial

from transformercvn.network.networks.neutrino_network import NeutrinoNetwork
from transformercvn.options import Options


class SherpaNetwork(NeutrinoNetwork):
    def __init__(
        self, options: Options, client: Optional[Client] = None, trial: Optional[Trial] = None
    ):
        super(SherpaNetwork, self).__init__(options)

        # Sherpa Objects
        self.client = client
        self.trial = trial
        self.sherpa_iteration = 0

        # study => trial
        assert (not self.client) or self.trial

    def commit_sherpa(self, objective, context: Optional[dict] = None):
        if context is None:
            context = {}

        if self.client:
            self.sherpa_iteration += 1
            self.client.send_metrics(
                trial=self.trial,
                iteration=self.sherpa_iteration,
                objective=objective.item(),
                context={key: val.item() for key, val in context.items()},
            )

    def validation_epoch_end(self, output):
        accuracy = super(SherpaNetwork, self).validation_epoch_end(output)
        self.commit_sherpa(accuracy)
