"""
Unets for filling dead channels of collection view images and induction view images. Architectures
for the two models are the same but the padding is different as the views produce different sized images.

Architecture modelled on the submanifold sparse Unet (arXiv:1711.10275). Ideally this model would 
make use of sparse convolutions but current software restrictions means it does not.
"""

import torch
from torch import nn


class UnetInduction(nn.Module):
    def __init__(self, in_channels=1, out_channels=1):
        super(UnetInduction, self).__init__()

        self.conv_in = self.single_conv_in(1, 4, 3)
        self.convs1_L = self.conv_block(4, 4 ,3)
        self.down_conv1 = self.down_conv(4, 8, 3, 3)
        self.convs2_L = self.conv_block(8, 8, 3)
        self.down_conv2 = self.down_conv(8, 16, 3, 3)
        self.convs3_L = self.conv_block(16, 16, 3, padding=[(0,0), (0,0)])
        self.down_conv3 = self.down_conv(16, 32, 3, 3)
        self.convs4_L = self.conv_block(32, 32, 3, padding=[(0,0),(0,0)])

        self.down_conv_bottom = self.down_conv(32, 64, 3, 3)
        self.convs_bottom = self.conv_block(64, 64, 3, padding=[(1,1),(1,1)])
        self.up_conv_bottom = self.up_conv(64, 32, 3, 3, output_padding=(0,0))

        self.convs4_R = self.conv_block(32*2, 32, 3, padding=[(2,2),(2,2)])
        self.up_conv1 = self.up_conv(32, 16, 3, 3, output_padding=(2,0))
        self.convs3_R = self.conv_block(16*2, 16, 3, padding=[(2,2),(2,2)])
        self.up_conv2 = self.up_conv(16, 8, 3, 3, output_padding=(2,2))
        self.convs2_R = self.conv_block(8*2, 8, 3)
        self.up_conv3 = self.up_conv(8, 4, 3, 3, output_padding=(0,2))
        self.convs1_R = self.conv_block(4*2, 4, 3)
        self.conv_out = self.single_conv_out(4, 1, 3)

    def forward(self, conv1):   
        conv1 = self.conv_in(conv1)
        conv1 = self.convs1_L(conv1)
        conv2 = self.down_conv1(conv1)
        conv2 = self.convs2_L(conv2)
        conv3 = self.down_conv2(conv2)
        conv3 = self.convs3_L(conv3)
        conv4 = self.down_conv3(conv3)
        conv4 = self.convs4_L(conv4)

        conv_bottom = self.down_conv_bottom(conv4)
        conv_bottom = self.convs_bottom(conv_bottom)
        conv_bottom = self.up_conv_bottom(conv_bottom)

        conv4 = self.convs4_R(torch.cat([conv_bottom, conv4], 1))
        conv4 = self.up_conv1(conv4)
        conv3 = self.convs3_R(torch.cat([conv4, conv3], 1))
        conv3 = self.up_conv2(conv3)
        conv2 = self.convs2_R(torch.cat([conv3, conv2], 1))
        conv2 = self.up_conv3(conv2)
        conv1 = self.convs1_R(torch.cat([conv2, conv1], 1))
        conv1 = self.conv_out(conv1)

        return conv1

    def single_conv_in(self, in_channels, out_channels, kernel_size, padding=1):
        conv = nn.Sequential(
            nn.Conv2d(in_channels, out_channels, kernel_size, padding=padding))
        
        return conv

    def single_conv_out(self, in_channels, out_channels, kernel_size, padding=1):
        conv = nn.Sequential(
            nn.BatchNorm2d(in_channels),
            nn.ReLU(),
            nn.Conv2d(in_channels, out_channels, kernel_size, padding=padding))
        
        return conv

    def conv_block(self, in_channels, out_channels, kernel_size, padding=[(1,1),(1,1)]):
        conv_block = nn.Sequential(
            nn.BatchNorm2d(in_channels),
            nn.ReLU(),
            nn.Conv2d(in_channels, out_channels, kernel_size, padding=padding[0]),
            nn.BatchNorm2d(out_channels),
            nn.ReLU(),
            nn.Conv2d(out_channels, out_channels, kernel_size, padding=padding[1]))

        return conv_block

    def down_conv(self, in_channels, out_channels, kernel_size, stride):
        down_conv = nn.Sequential(
            nn.BatchNorm2d(in_channels),
            nn.ReLU(),
            nn.Conv2d(in_channels, out_channels, kernel_size, stride))

        return down_conv

    def up_conv(self, in_channels, out_channels, kernel_size, stride, output_padding=0):
        up_conv = nn.Sequential(
            nn.BatchNorm2d(in_channels),
            nn.ReLU(),
            nn.ConvTranspose2d(in_channels, out_channels, kernel_size, stride, output_padding=output_padding))
        
        return up_conv


class UnetCollection(nn.Module):
    def __init__(self, in_channels=1, out_channels=1):
        super(UnetCollection, self).__init__()

        self.conv_in = self.single_conv_in(1, 4, 3)
        self.convs1_L = self.conv_block(4, 4 ,3)
        self.down_conv1 = self.down_conv(4, 8, 3, 3)
        self.convs2_L = self.conv_block(8, 8, 3)
        self.down_conv2 = self.down_conv(8, 16, 3, 3)
        self.convs3_L = self.conv_block(16, 16, 3, padding=[(0,0), (0,0)])
        self.down_conv3 = self.down_conv(16, 32, 3, 3)
        self.convs4_L = self.conv_block(32, 32, 3, padding=[(0,0),(0,0)])

        self.down_conv_bottom = self.down_conv(32, 64, 3, 3)
        self.convs_bottom = self.conv_block(64, 64, 3, padding=[(1,1),(1,1)])
        self.up_conv_bottom = self.up_conv(64, 32, 3, 3, output_padding=(0,0))

        self.convs4_R = self.conv_block(32*2, 32, 3, padding=[(2,2),(2,2)])
        self.up_conv1 = self.up_conv(32, 16, 3, 3, output_padding=(2,1))
        self.convs3_R = self.conv_block(16*2, 16, 3, padding=[(2,2),(2,2)])
        self.up_conv2 = self.up_conv(16, 8, 3, 3, output_padding=(2,1))
        self.convs2_R = self.conv_block(8*2, 8, 3)
        self.up_conv3 = self.up_conv(8, 4, 3, 3, output_padding=(0,0))
        self.convs1_R = self.conv_block(4*2, 4, 3)
        self.conv_out = self.single_conv_out(4, 1, 3)

    def forward(self, conv1):        
        conv1 = self.conv_in(conv1)
        conv1 = self.convs1_L(conv1)
        conv2 = self.down_conv1(conv1)
        conv2 = self.convs2_L(conv2)
        conv3 = self.down_conv2(conv2)
        conv3 = self.convs3_L(conv3)
        conv4 = self.down_conv3(conv3)
        conv4 = self.convs4_L(conv4)

        conv_bottom = self.down_conv_bottom(conv4)
        conv_bottom = self.convs_bottom(conv_bottom)
        conv_bottom = self.up_conv_bottom(conv_bottom)

        conv4 = self.convs4_R(torch.cat([conv_bottom, conv4], 1))
        conv4 = self.up_conv1(conv4)
        conv3 = self.convs3_R(torch.cat([conv4, conv3], 1))
        conv3 = self.up_conv2(conv3)
        conv2 = self.convs2_R(torch.cat([conv3, conv2], 1))
        conv2 = self.up_conv3(conv2)
        conv1 = self.convs1_R(torch.cat([conv2, conv1], 1))
        conv1 = self.conv_out(conv1)

        return conv1

    def single_conv_in(self, in_channels, out_channels, kernel_size, padding=1):
        conv = nn.Sequential(
            nn.Conv2d(in_channels, out_channels, kernel_size, padding=padding))
        
        return conv

    def single_conv_out(self, in_channels, out_channels, kernel_size, padding=1):
        conv = nn.Sequential(
            nn.BatchNorm2d(in_channels),
            nn.ReLU(),
            nn.Conv2d(in_channels, out_channels, kernel_size, padding=padding))
        
        return conv

    def conv_block(self, in_channels, out_channels, kernel_size, padding=[(1,1),(1,1)]):
        conv_block = nn.Sequential(
            nn.BatchNorm2d(in_channels),
            nn.ReLU(),
            nn.Conv2d(in_channels, out_channels, kernel_size, padding=padding[0]),
            nn.BatchNorm2d(out_channels),
            nn.ReLU(),
            nn.Conv2d(out_channels, out_channels, kernel_size, padding=padding[1]))

        return conv_block

    def down_conv(self, in_channels, out_channels, kernel_size, stride):
        down_conv = nn.Sequential(
            nn.BatchNorm2d(in_channels),
            nn.ReLU(),
            nn.Conv2d(in_channels, out_channels, kernel_size, stride))

        return down_conv

    def up_conv(self, in_channels, out_channels, kernel_size, stride, output_padding=0):
        up_conv = nn.Sequential(
            nn.BatchNorm2d(in_channels),
            nn.ReLU(),
            nn.ConvTranspose2d(in_channels, out_channels, kernel_size, stride, output_padding=output_padding))
        
        return up_conv
        
