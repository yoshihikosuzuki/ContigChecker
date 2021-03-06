from collections import defaultdict
from logzero import logger
from .cc_map import run_mapping
from .cc_depth import calc_depth
from .cc_annotate import annotate_regions
from .cc_plot import Plotter
#from .cc_modify import modify_sequence

cfg_sections = ["Map", "Depth", "Annotate"]


def main():
    args = load_args()
    cfgs = load_cfgs(args.config_file)
    if "Map" in cfgs:
        run_mapping(**cfgs["Map"])
    if "Depth" in cfgs:
        calc_depth(**cfgs["Depth"])
    if "Annotate" in cfgs:
        #annotate_regions(**cfgs["Annotate"])
        annotate_regions()
    if "Plot" in cfgs:
        p = Plotter()
        p.plot(cfgs["Plot"]["contig_id"])
    if "Modify" in cfgs:
        #modify_sequence()
        pass


def load_cfgs(cfg_fname):
    import configparser
    parser = configparser.ConfigParser()
    parser.read(cfg_fname)

    cfgs = defaultdict(dict)
    for section in cfg_sections:
        if parser.has_section(section):
            cfgs[section] = dict(parser.items(section))
    return cfgs


def load_args():
    import argparse
    parser = argparse.ArgumentParser(description="Call methylation motifs for each genome bin")
    parser.add_argument('config_file',
                        type=str,
                        help="Config file")
    return parser.parse_args()


if __name__ == "__main__":
    main()
