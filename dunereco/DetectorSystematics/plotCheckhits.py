import argparse

import ROOT
from matplotlib.lines import Line2D
from matplotlib import pyplot as plt
import numpy as np


def main(args):
    f = ROOT.TFile.Open(args.in_file)
    t = f.Get(args.ttree_loc)

    simchs = {}
    wires = {}
    hits = {}
    digits = {}
    ch_types = {}

    for i_event, event in enumerate(t):
        if i_event < args.n_skip:
            continue

        valid_chs = set()

        if args.wires:
            for i_resp, wire_resp in enumerate(event.wires):
                wire_resp = list(wire_resp)
                if sum(wire_resp) > args.wires_mag_thres:
                    wires[event.wires_chs[i_resp]] = wire_resp
                    ch_types[event.wires_chs[i_resp]] = event.ch_types[event.wires_chs[i_resp]]
                    valid_chs.add(event.wires_chs[i_resp])

        if args.simchannels:
            for i_resp, sc_resp in enumerate(event.simchannels):
                sc_resp = list(sc_resp)
                if sum(sc_resp) > args.simchannels_mag_thres:
                    #print("Got a simchannel above threshold with ", sum(sc_resp))
                    simchs[event.simchannels_chs[i_resp]] = sc_resp
                    ch_types[event.simchannels_chs[i_resp]] = (
                        event.ch_types[event.simchannels_chs[i_resp]]
                    )
                    valid_chs.add(event.simchannels_chs[i_resp])

        if args.simchannel_number:
            if args.simchannel_number in valid_chs:
                valid_chs = {args.simchannel_number}
            else:
                print("Channel ",args.simchannel_number," is not in the list of valid channels")
                valid_chs = set()

        if args.hits:
            for i_resp, hits_resp in enumerate(event.hits):
                hits_resp = list(hits_resp)
                if sum(hits_resp) > args.hits_mag_thres:
                    hits[event.hits_chs[i_resp]] = hits_resp
                    ch_types[event.hits_chs[i_resp]] = event.ch_types[event.hits_chs[i_resp]]
                    #valid_chs.add(event.hits_chs[i_resp])

        if args.digits:
            for i_resp, digit_resp in enumerate(event.digits):
                digit_resp = list(digit_resp)
                if event.digits_chs[i_resp] in valid_chs:
                    digits[event.digits_chs[i_resp]] = digit_resp

        # for i_ch, ch in enumerate(sorted(wires, key=lambda ch: sum(wires[ch]), reverse=True)):
        for i_ch, ch in enumerate(valid_chs):
            if args.n_plot != -1 and i_ch >= args.n_plot:
                break

            ticks = np.arange(1, 6001)

            _, ax = plt.subplots()

            ax.set_ylabel("signal")
            ax.set_xlabel("tick")

            if args.wires:
                ax.hist(
                    ticks, bins=len(ticks), weights=wires.get(ch, np.zeros(6000)), #weights=wires[ch],
                    histtype="step", label="recob::Wire", color='b'
                )

            if args.digits:
                ax.hist(
                    ticks, bins=len(ticks), weights=digits.get(ch, np.zeros(6000)),#weights=digits[ch],
                    histtype="step", label="raw::RawDigit", color='b'
                )

            if args.hits:
                ax.hist(
                    ticks, bins=len(ticks), weights=hits.get(ch, np.zeros(6000)), #weights=hits[ch],
                    histtype="step", label="recob::Hit", color='r'
                )

            ax.set_xlim(
                left=args.tick_low if args.tick_low != -1 else None,
                right=args.tick_high if args.tick_high != -1 else None
            )


            if args.simchannels:
                ax2 = ax.twinx()

                ax2.hist(
                    ticks, bins=len(ticks), weights=simchs.get(ch, np.zeros(6000)),
                    histtype="step", label="sim::SimChannels", color='g'
                )
                ax2.set_ylabel("energy")

                ax_ylims = ax.axes.get_ylim()
                ax_yratio = ax_ylims[0] / ax_ylims[1]
                ax2_ylims = ax2.axes.get_ylim()
                ax2_yratio = ax2_ylims[0] / ax2_ylims[1]
                if ax_yratio < ax2_yratio:
                    ax2.set_ylim(bottom = ax2_ylims[1]*ax_yratio)
                else:
                    ax.set_ylim(bottom = ax_ylims[1]*ax2_yratio)

                handles, labels = ax.get_legend_handles_labels()
                handles2, labels2 = ax2.get_legend_handles_labels()
                handles += handles2
                labels += labels2
                new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]
                plt.legend(handles=new_handles, labels=labels, prop={'size': 12})
            else:
                handles, labels = ax.get_legend_handles_labels()
                new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]
                plt.legend(handles=new_handles, labels=labels, prop={'size': 12})

            if args.tpc_idx != -1 and args.plane_idx != -1:
                plt.title(
                    "TPCIndex: {} PlaneIndex: {} ({}) Channel: {} ".format(
                        args.tpc_idx,
                        args.plane_idx,
                        "collection" if ch_types[ch] else "induction",
                        ch - args.ch_offset
                    )
                )
            else:
                plt.title(
                    "ch {} signal_type {}".format(
                        ch, "collection" if ch_types[ch] else "induction"
                    )
                )

            plt.tight_layout()
            plt.show()


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("in_file", type=str)
    parser.add_argument("ttree_loc", type=str)

    parser.add_argument("--wires", action="store_true")
    parser.add_argument("--hits", action="store_true")
    parser.add_argument("--digits", action="store_true")
    parser.add_argument("--simchannels", action="store_true")
    parser.add_argument("--n_plot", type=int, default=5)
    parser.add_argument("--n_skip", type=int, default=0)
    parser.add_argument("--wires_mag_thres", type=int, default=10)
    parser.add_argument("--simchannels_mag_thres", type=int, default=10)
    parser.add_argument("--hits_mag_thres", type=int, default=10)
    parser.add_argument("--tick_low", type=int, default=-1)
    parser.add_argument("--tick_high", type=int, default=-1)
    parser.add_argument("--ch_offset", type=int, default=0)
    parser.add_argument("--tpc_idx", type=int, default=-1)
    parser.add_argument("--plane_idx", type=int, default=-1)
    parser.add_argument("--simchannel_number", type=int, default=None)
    return parser.parse_args()


if __name__ == '__main__':
    main(parse_arguments())
