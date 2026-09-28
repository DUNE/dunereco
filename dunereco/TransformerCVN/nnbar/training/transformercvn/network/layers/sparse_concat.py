import torch
from MinkowskiCoordinateManager import CoordinateManager
from MinkowskiEngine import SparseTensor
from MinkowskiEngineBackend._C import CoordinateMapKey
from torch.autograd import Function
from torch.nn import Module


class MinkowskiConcatFunction(Function):
    @staticmethod
    def forward(
        ctx,
        in_coords_keys: list,
        out_coords_key: CoordinateMapKey,
        coordinate_manager: CoordinateManager,
        *in_feats,
    ):
        assert isinstance(in_feats, (list, tuple)), "Input must be a collection of Tensors"
        assert len(in_feats) > 1, "input must be a set with at least 2 Tensors"
        assert len(in_feats) == len(in_coords_keys), (
            "The input features and keys must have the same length"
        )

        union_maps = coordinate_manager.union_map(in_coords_keys, out_coords_key)
        feature_shapes = torch.tensor([0] + [feat.shape[1] for feat in in_feats], dtype=torch.long)
        feature_offsets = feature_shapes.cumsum(0)

        out_feat = torch.zeros(
            (coordinate_manager.size(out_coords_key), feature_shapes.sum()),
            dtype=in_feats[0].dtype,
            device=in_feats[0].device,
        )
        for i, (in_feat, union_map) in enumerate(zip(in_feats, union_maps)):
            out_feat[union_map[1], feature_offsets[i] : feature_offsets[i + 1]] += in_feat[
                union_map[0]
            ]

        ctx.keys = (in_coords_keys, coordinate_manager)
        ctx.save_for_backward(*union_maps)
        return out_feat

    @staticmethod
    def backward(ctx, grad_out_feat):
        if not grad_out_feat.is_contiguous():
            grad_out_feat = grad_out_feat.contiguous()

        union_maps = ctx.saved_tensors
        in_coords_keys, coordinate_manager = ctx.keys
        dtype, device = (
            grad_out_feat.dtype,
            grad_out_feat.device,
        )

        grad_in_feats = []
        current_offset = 0
        for in_coords_key, union_map, feats in zip(in_coords_keys, union_maps, grad_out_feat):
            num_ch = feats.shape[1]
            grad_in_feat = torch.zeros(
                (coordinate_manager.size(in_coords_key), num_ch),
                dtype=dtype,
                device=device,
            )
            grad_in_feat[union_map[0]] = grad_out_feat[
                union_map[1], current_offset : current_offset + num_ch
            ]
            grad_in_feats.append(grad_in_feat)
            current_offset += num_ch

        return (None, None, None, *grad_in_feats)


class MinkowskiConcat(Module):
    r"""Create a union of all sparse tensors and add overlapping features.
    Args:
        None
    .. warning::
       This function is experimental and the usage can be changed in the future updates.
    """

    def __init__(self):
        super(MinkowskiConcat, self).__init__()
        self.union = MinkowskiConcatFunction()

    def forward(self, *inputs):
        assert isinstance(inputs, (list, tuple)), "The input must be a list or tuple"
        for s in inputs:
            assert isinstance(s, SparseTensor), "Inputs must be sparse tensors."

        if len(inputs) == 1:
            return inputs[0]

        # Assert the same coordinate manager
        ref_coordinate_manager = inputs[0].coordinate_manager
        for s in inputs:
            assert ref_coordinate_manager == s.coordinate_manager, (
                "Invalid coordinate manager. All inputs must have the same coordinate manager."
            )

        in_coordinate_map_key = inputs[0].coordinate_map_key
        coordinate_manager = inputs[0].coordinate_manager
        out_coordinate_map_key = CoordinateMapKey(in_coordinate_map_key.get_coordinate_size())
        output = self.union.apply(
            [input.coordinate_map_key for input in inputs],
            out_coordinate_map_key,
            coordinate_manager,
            *[input.F for input in inputs],
        )
        return SparseTensor(
            output,
            coordinate_map_key=out_coordinate_map_key,
            coordinate_manager=coordinate_manager,
        )

    def __repr__(self):
        return self.__class__.__name__ + "()"
