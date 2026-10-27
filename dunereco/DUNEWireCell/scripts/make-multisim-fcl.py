#!/usr/bin/env python3
"""Generate a "multi-sim" WC/LS FHiCL from a "single-sim" one.

The dune10kt-vd WC/LS job
(dune10kt-vd/wcls-sim-drift-simchannel-nf-sp.jsonnet) supports a
process_mode="multi-sim" in which the input depos are fanned out to one finger
per anode in structs.process_anode_indices, and EACH anode terminates in its own
trio of art-output nodes: a wclsDepoFluxWriter named "postdrift<id>"
(sim::SimChannel) and two wclsFrameSaver nodes named "simdigits<id>"
(raw::RawDigit) and "spsignals<id>" (recob::Wire).  For art to actually write
those products, the WireCellToolkit module's `outputers` list MUST name every
one of those nodes -- exactly the set implied by process_anode_indices.  Keeping
those lists in sync by hand is tedious and error-prone.

This script derives BOTH lists from a single user-supplied set of anode indices
and emits a small overlay FHiCL that #includes the original single-sim FHiCL and
overrides only the multi-sim keys.  Using FHiCL's native include+override keeps
the original file the single source of truth (no template to drift out of sync,
no fragile text surgery of the original).

Examples
--------
  # anodes 0,1,2,3 -> overlay on stdout
  make-multisim-fcl.py 0-3

  # a mixed spec, writing to a file
  make-multisim-fcl.py 0-3,10,20-22 -o wcls-multisim-0-3_10_20-22.fcl

  # a different base FHiCL / producer label / RawDigit tag pattern
  make-multisim-fcl.py --base my-single-sim.fcl --producer nfsp \\
                       --frame-tag-pattern 'daq%d' 0 5 17
"""

import argparse
import os
import re
import sys


def parse_anode_spec(tokens):
    """Parse anode-index tokens into a sorted, de-duplicated list of ints.

    Each token is a comma- and/or space-separated list of items, where an item
    is either a single index "N" or an inclusive range "A-B" (A<=B).
    """
    ids = set()
    for tok in tokens:
        for item in tok.replace(",", " ").split():
            m = re.fullmatch(r"(\d+)-(\d+)", item)
            if m:
                lo, hi = int(m.group(1)), int(m.group(2))
                if lo > hi:
                    raise ValueError(f"empty range '{item}' (lo > hi)")
                ids.update(range(lo, hi + 1))
                continue
            if re.fullmatch(r"\d+", item):
                ids.add(int(item))
                continue
            raise ValueError(f"bad anode-index item '{item}' "
                             "(want N or A-B, comma/space separated)")
    if not ids:
        raise ValueError("no anode indices given")
    return sorted(ids)


def detect_producer(base_path):
    """Return the FHiCL producer label whose block configures WireCellToolkit.

    Light line scan: track the most recent `<label> : {` seen inside a
    `producers:` table and return the label of the block that contains a
    `module_type ... WireCellToolkit` (or a `wcls_main`).  Returns None if the
    file can't be read or no such producer is found.
    """
    try:
        with open(base_path) as fp:
            text = fp.read()
    except OSError:
        return None
    # candidate producer labels: "<label> : {" or "<label>: {"
    label = None
    for line in text.splitlines():
        m = re.match(r"\s*([A-Za-z_]\w*)\s*:\s*\{", line)
        if m:
            label = m.group(1)
        if label and re.search(r"module_type\s*:\s*WireCellToolkit", line):
            return label
        if label and re.search(r"\bwcls_main\s*:", line):
            return label
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Generate a multi-sim WC/LS FHiCL overlay from a single-sim one.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    ap.add_argument("anodes", nargs="+", metavar="ANODESPEC",
                    help="anode indices: single indices and/or A-B ranges, "
                         "comma/space separated (e.g. '0-3,10,20-22')")
    ap.add_argument("--base", default="wcls-sim-drift-simchannel-nf-sp.fcl",
                    help="single-sim FHiCL to #include (a bare name resolved via "
                         "FHICL_FILE_PATH at run time, or a path to an existing "
                         "file -- a path is also read to auto-detect --producer). "
                         "[default: %(default)s]")
    ap.add_argument("--producer", default=None,
                    help="WireCellToolkit producer label to override "
                         "[default: auto-detect from --base, else 'tpcrawdecoder']")
    ap.add_argument("--frame-tag-pattern", default="daq%d",
                    help="per-anode raw::RawDigit product-instance pattern "
                         "(one %%d on the anode ident) [default: %(default)s]")
    ap.add_argument("--input-depos-tag", default=None,
                    help="override the input SimEnergyDeposit art tag "
                         "(e.g. 'IonAndScint' or 'plopper:bogus'); "
                         "omit to leave the base FHiCL's value")
    ap.add_argument("--process-name", default=None,
                    help="override process_name")
    ap.add_argument("--no-rawdigits", dest="save_rawdigits", action="store_false",
                    default=True,
                    help="do NOT save the per-anode raw::RawDigit (ADC) products: "
                         "sim feeds sigproc inside the WCT graph so the dense ADC "
                         "waveform never reaches the art::Event (big memory win "
                         "for many-anode jobs).  Omits the simdigits<id> outputers.")
    ap.add_argument("--max-anode", type=int, default=320,
                    help="warn if any index is >= this (FD-VD has 320 anodes) "
                         "[default: %(default)s]; 0 disables the check")
    ap.add_argument("-o", "--output", default=None,
                    help="output FHiCL file [default: stdout]")
    args = ap.parse_args(argv)

    try:
        ids = parse_anode_spec(args.anodes)
    except ValueError as e:
        ap.error(str(e))

    if args.max_anode and ids[-1] >= args.max_anode:
        print(f"WARNING: anode index {ids[-1]} >= --max-anode {args.max_anode}; "
              "may exceed the detector geometry.", file=sys.stderr)

    # Resolve the producer label and the #include name.
    include_name = os.path.basename(args.base)
    producer = args.producer
    if producer is None and os.path.isfile(args.base):
        producer = detect_producer(args.base)
    if producer is None:
        # try the default base sitting next to this script's dune10kt-vd sibling
        here = os.path.dirname(os.path.abspath(__file__))
        guess = os.path.join(here, os.pardir, "dune10kt-vd", include_name)
        if os.path.isfile(guess):
            producer = detect_producer(guess)
    if producer is None:
        producer = "tpcrawdecoder"

    indices = "[" + ", ".join(str(i) for i in ids) + "]"
    base_key = f"physics.producers.{producer}.wcls_main"

    lines = []
    lines.append("# Multi-sim WC/LS overlay -- GENERATED by make-multisim-fcl.py.")
    lines.append("# Do not edit by hand; re-run the script to change the anode set.")
    lines.append(f"#   anode indices: {indices}")
    lines.append(f"#   base:    {include_name}  (producer '{producer}')")
    lines.append("")
    lines.append(f'#include "{include_name}"')
    lines.append("")
    if args.process_name:
        lines.append(f"process_name: {args.process_name}")
        lines.append("")
    lines.append(f'{base_key}.params.process_mode: "multi-sim"')
    lines.append(f'{base_key}.params.frame_tag_pattern: "{args.frame_tag_pattern}"')
    if args.input_depos_tag is not None:
        lines.append(f'{base_key}.params.input_depos_tag: "{args.input_depos_tag}"')
    lines.append(f'{base_key}.structs.process_anode_indices: "{indices}"')
    lines.append(f'{base_key}.structs.save_rawdigits: '
                 f'{"true" if args.save_rawdigits else "false"}')
    lines.append("")
    if args.save_rawdigits:
        lines.append("# outputers: per anode, a SimChannel writer + a RawDigit saver +")
        lines.append("# a recob::Wire saver -- kept consistent with")
        lines.append("# process_anode_indices above.  In multi-sim each anode terminates")
        lines.append("# in its own trio named postdrift<id> (sim::SimChannel),")
        lines.append("# simdigits<id> (raw::RawDigit) and spsignals<id> (recob::Wire).")
    else:
        lines.append("# outputers: per anode, a SimChannel writer + a recob::Wire")
        lines.append("# saver (NO RawDigit -- save_rawdigits=false keeps ADC inside")
        lines.append("# the WCT graph): postdrift<id> (sim::SimChannel) and")
        lines.append("# spsignals<id> (recob::Wire).")
    lines.append(f"{base_key}.outputers: [")
    first = True
    for i in ids:
        lead = "      " if first else "    , "
        first = False
        trio = '"wclsDepoFluxWriter:postdrift%d"' % i
        if args.save_rawdigits:
            trio += ', "wclsFrameSaver:simdigits%d"' % i
        trio += ', "wclsFrameSaver:spsignals%d"' % i
        lines.append(lead + trio)
    lines.append("]")
    lines.append("")
    text = "\n".join(lines)

    if args.output:
        with open(args.output, "w") as fp:
            fp.write(text)
        print(f"wrote {args.output}  (anodes {indices}, producer '{producer}')",
              file=sys.stderr)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
