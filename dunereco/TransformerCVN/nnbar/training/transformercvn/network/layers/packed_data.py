from typing import Tuple

import torch
from torch import Tensor, jit


def masked_pack(data: Tensor, mask: Tensor) -> Tensor:
    data_shape = data.shape

    while mask.ndim != data.ndim:
        mask = mask.unsqueeze(dim=-1)

    packed = torch.masked_select(data, mask)
    return packed.view(-1, *data_shape[2:])


def masked_pad(packed_data: Tensor, mask: Tensor) -> Tensor:
    output = torch.zeros(
        *mask.shape, *packed_data.shape[1:], dtype=packed_data.dtype, device=packed_data.device
    )

    while mask.ndim != output.ndim:
        mask = mask.unsqueeze(dim=-1)

    return torch.masked_scatter(output, mask, packed_data)


@jit.script
def masked_pack_1d(data: Tensor, mask: Tensor) -> Tensor:
    B, L, C = data.shape

    I1 = torch.arange(B).repeat_interleave(mask.sum(1))
    I2 = torch.masked_select(torch.arange(L).view(1, -1).repeat(B, 1), mask)
    return data[I1, I2]


@jit.script
def masked_pad_1d_v2(packed_data: Tensor, mask: Tensor) -> Tensor:
    B, L = mask.shape
    (
        _,
        C,
    ) = packed_data.shape

    output = torch.zeros(B, L, C, dtype=packed_data.dtype, device=packed_data.device)
    mask = mask.view(B, L, 1)

    return torch.masked_scatter(output, mask, packed_data)


@jit.script
def masked_pad_1d(packed_data: Tensor, mask: Tensor) -> Tensor:
    B, L = mask.shape
    (
        _,
        C,
    ) = packed_data.shape

    I1 = torch.arange(B).repeat_interleave(mask.sum(1))
    I2 = torch.masked_select(torch.arange(L).view(1, -1).repeat(B, 1), mask)
    output = torch.zeros(B, L, C, dtype=packed_data.dtype, device=packed_data.device)
    output[I1, I2] = packed_data

    return output


@jit.script
def masked_pack_1d_precomputed(data: Tensor, mask: Tensor) -> Tuple[Tensor, Tensor, Tensor]:
    B, L, C = data.shape
    device = data.device

    I1 = torch.arange(B, device=device).repeat_interleave(mask.sum(1))
    I2 = torch.masked_select(torch.arange(L, device=device).view(1, -1).repeat(B, 1), mask)
    return data[I1, I2], I1, I2


@jit.script
def masked_pad_1d_precomputed(
    packed_data: Tensor, I1: Tensor, I2: Tensor, batch_size: int, max_length: int
) -> Tensor:
    (
        _,
        C,
    ) = packed_data.shape

    output = torch.zeros(
        batch_size, max_length, C, dtype=packed_data.dtype, device=packed_data.device
    )
    output[I1, I2] = packed_data

    return output


@jit.script
def masked_pack_3d_precomputed(data: Tensor, mask: Tensor) -> Tuple[Tensor, Tensor, Tensor]:
    B, L, C, H, W = data.shape
    device = data.device

    I1 = torch.arange(B, device=device).repeat_interleave(mask.sum(1))
    I2 = torch.masked_select(torch.arange(L, device=device).view(1, -1).repeat(B, 1), mask)
    return data[I1, I2], I1, I2


@jit.script
def masked_pad_3d_precomputed(
    packed_data: Tensor, I1: Tensor, I2: Tensor, batch_size: int, max_length: int
) -> Tensor:
    _, C, H, W = packed_data.shape

    output = torch.zeros(
        batch_size, max_length, C, H, W, dtype=packed_data.dtype, device=packed_data.device
    )
    output[I1, I2] = packed_data
    return output


@jit.script
def masked_pack_3d(data: Tensor, mask: Tensor) -> Tensor:
    B, L, C, H, W = data.shape

    I1 = torch.arange(B).repeat_interleave(mask.sum(1))
    I2 = torch.masked_select(torch.arange(L).view(1, -1).repeat(B, 1), mask)
    return data[I1, I2]


@jit.script
def masked_pad_3d(packed_data: Tensor, mask: Tensor) -> Tensor:
    B, L = mask.shape
    _, C, H, W = packed_data.shape

    output = torch.zeros(B, L, C, H, W, dtype=packed_data.dtype, device=packed_data.device)
    mask = mask.view(B, L, 1, 1, 1)

    return torch.masked_scatter(output, mask, packed_data)
