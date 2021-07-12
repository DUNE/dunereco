"""
HolePixelLoss
This loss mimics nividia's pixelwise loss for holes (L1)
used in the infill network

Loss function designed for a network that infills tracks. Separate functions for collection views
(unipolar signals) and induction views (bipolar signal)
Taken from github.com/NuTufts/sparse_infill.
"""

import torch
import torch.nn as nn


# Taken from torch.nn.modules.loss
def _assert_no_grad(variable):
    assert not variable.requires_grad, \
        "nn criterions don't compute the gradient w.r.t. targets - please " \
        "mark these variables as not requiring gradients"


class InfillLossInduction(nn.modules.loss._WeightedLoss):
    def __init__(self, weight=None, size_average=False, ignore_index=-100):
        super(InfillLossInduction, self).__init__(weight, size_average)
        self.ignore_index = ignore_index
        self.reduce = False
        self.size_average = size_average

    def forward(self, predict, true, input):
        _assert_no_grad(true)

        # Want three losses: non-dead, dead w/o charge, dead w/ charge
        L1loss=torch.nn.L1Loss(self.size_average)
        nondeadweight = 1.0
        deadnochargeweight = 500.0
        deadlowchargeweight = 100.0
        deadhighchargeweight = 100.0
        deadhighestchargeweight = 100.0

        # Identify dead channels
        goodch = (input != 0)
        for idx, col in enumerate(goodch[0, 0, :, :].T):
            if col.sum():
                goodch[0, 0, :, idx] = 1
        deadch = (~goodch).float()
        goodch = goodch.float()

        # Compute non-dead loss
        predictgood = goodch * predict
        adcgood = goodch * true
        totnondead = goodch.sum().float()
        if (totnondead == 0):
            totnondead = 1.0
        nondeadloss = (L1loss(predictgood, adcgood)*nondeadweight)/totnondead

        # Compute dead with highest true charge loss
        deadchhighestcharge = deadch * (true.abs() > 70).float()
        predictdeadhighestcharge = predict*deadchhighestcharge
        adcdeadhighestcharge = true*deadchhighestcharge
        totdeadhighestcharge = deadchhighestcharge.sum().float()
        if (totdeadhighestcharge == 0):
            totdeadhighestcharge = 1.0
        deadhighestchargeloss = (L1loss(predictdeadhighestcharge,adcdeadhighestcharge)*
            deadhighestchargeweight)/totdeadhighestcharge

        # Compute dead with high true charge loss
        deadchhighcharge = deadch * (true.abs() > 40).float() * (true.abs() < 70).float()
        predictdeadhighcharge = predict*deadchhighcharge
        adcdeadhighcharge = true*deadchhighcharge
        totdeadhighcharge = deadchhighcharge.sum().float()
        if (totdeadhighcharge == 0):
            totdeadhighcharge = 1.0
        deadhighchargeloss = (L1loss(predictdeadhighcharge,adcdeadhighcharge)*
            deadhighchargeweight)/totdeadhighcharge

        # Compute dead with low true chanrge and dead with no (< 10) true charge
        deadchlowcharge = deadch * (true.abs() > 10).float() * (true.abs() < 40).float()
        deadchnocharge = deadch * (true.abs() < 10).float()

        # < 10 charge regions where the bipolar signal crosses zero should be low true charge rather than dead.
        # Identify them and move them from deadnocharge to deadlowcharge
        dead_cols = torch.nonzero(deadch)
        for i, img in enumerate(true * deadch):
            dead_cols = torch.unique(torch.nonzero(img[0, :, :])[:, 1])
            for col in dead_cols:
                diffs = img[0, :, col].roll(1) - img[0, :, col].roll(-1)
                diffs[0] = 0
                diffs[-1] = 0
                deadchlowcharge_hidden = (diffs.abs() > 20).float() * deadchnocharge[i, 0, :, col].float()
                deadchnocharge[i, 0, :, col] -= deadchlowcharge_hidden
                deadchlowcharge[i, 0, :, col] += deadchlowcharge_hidden

        predictdeadlowcharge = predict*deadchlowcharge
        adcdeadlowcharge = true*deadchlowcharge
        totdeadlowcharge = deadchlowcharge.sum().float()
        if (totdeadlowcharge == 0):
            totdeadlowcharge = 1.0
        deadlowchargeloss = (L1loss(predictdeadlowcharge,adcdeadlowcharge)*deadlowchargeweight)/totdeadlowcharge

        predictdeadnocharge = predict*deadchnocharge
        adcdeadnocharge = true*deadchnocharge
        totdeadnocharge = deadchnocharge.sum().float()
        if (totdeadnocharge == 0):
            totdeadnocharge = 1.0
        deadnochargeloss = (L1loss(predictdeadnocharge,adcdeadnocharge)*deadnochargeweight)/totdeadnocharge

        totloss = nondeadloss + deadnochargeloss + deadlowchargeloss +deadhighchargeloss+deadhighestchargeloss
        return nondeadloss, deadnochargeloss, deadlowchargeloss, deadhighchargeloss,deadhighestchargeloss, totloss


class InfillLossCollection(nn.modules.loss._WeightedLoss):
    def __init__(self, weight=None, size_average=False, ignore_index=-100):
        super(InfillLossCollection, self).__init__(weight, size_average)
        self.ignore_index = ignore_index
        self.reduce = False
        self.size_average = size_average

    def forward(self, predict, true, input):
        _assert_no_grad(true)

        # Want three losses: non-dead, dead w/o charge, dead w/ charge
        L1loss=torch.nn.L1Loss(self.size_average)
        nondeadweight = 1.0
        deadnochargeweight = 500.0
        deadlowchargeweight = 100.0
        deadhighchargeweight = 100.0
        deadhighestchargeweight = 100.0

        # Identify dead channels
        goodch = (input != 0)
        for idx, col in enumerate(goodch[0, 0, :, :].T):
            if col.sum():
                goodch[0, 0, :, idx] = 1
        deadch = (~goodch).float()
        goodch = goodch.float()

        # Compute non-dead loss
        predictgood = goodch * predict
        adcgood = goodch * true
        totnondead = goodch.sum().float()
        if (totnondead == 0):
            totnondead = 1.0
        nondeadloss = (L1loss(predictgood, adcgood)*nondeadweight)/totnondead

        # Compute dead with highest true charge loss
        deadchhighestcharge = deadch * (true.abs() > 70).float()
        predictdeadhighestcharge = predict*deadchhighestcharge
        adcdeadhighestcharge = true*deadchhighestcharge
        totdeadhighestcharge = deadchhighestcharge.sum().float()
        if (totdeadhighestcharge == 0):
            totdeadhighestcharge = 1.0
        deadhighestchargeloss = (L1loss(predictdeadhighestcharge,adcdeadhighestcharge)*deadhighestchargeweight)/totdeadhighestcharge

        # Compute dead with high true charge loss
        deadchhighcharge = deadch * (true.abs() > 40).float()*(true.abs() < 70).float()
        predictdeadhighcharge = predict*deadchhighcharge
        adcdeadhighcharge = true*deadchhighcharge
        totdeadhighcharge = deadchhighcharge.sum().float()
        if (totdeadhighcharge == 0):
            totdeadhighcharge = 1.0
        deadhighchargeloss = (L1loss(predictdeadhighcharge,adcdeadhighcharge)*deadhighchargeweight)/totdeadhighcharge

        # Compute dead with low true charge loss
        deadchlowcharge = deadch * (true.abs() > 10).float() *(true.abs() < 40).float()
        deadchnocharge = deadch * (true.abs() < 10).float()
        predictdeadlowcharge = predict*deadchlowcharge
        adcdeadlowcharge = true*deadchlowcharge
        totdeadlowcharge = deadchlowcharge.sum().float()
        if (totdeadlowcharge == 0):
            totdeadlowcharge = 1.0
        deadlowchargeloss = (L1loss(predictdeadlowcharge,adcdeadlowcharge)*deadlowchargeweight)/totdeadlowcharge

        # Compute dead with no (< 10) true charge loss
        predictdeadnocharge = predict*deadchnocharge
        adcdeadnocharge = true*deadchnocharge
        totdeadnocharge = deadchnocharge.sum().float()
        if (totdeadnocharge == 0):
            totdeadnocharge = 1.0
        deadnochargeloss = (L1loss(predictdeadnocharge,adcdeadnocharge)*deadnochargeweight)/totdeadnocharge

        totloss = nondeadloss + deadnochargeloss + deadlowchargeloss +deadhighchargeloss+deadhighestchargeloss
        return nondeadloss, deadnochargeloss, deadlowchargeloss, deadhighchargeloss,deadhighestchargeloss, totloss