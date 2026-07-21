#!/usr/bin/env python3
'''
tree_pipeline.py  —  Newick  ->  annotated tree SVG (and PNG/PDF), automatically.

This is a single-command wrapper around the existing tree-annotation toolkit.
It removes the manual GrapeTree-website step by generating a GrapeTree-compatible
JSON + SVG directly from a Newick file, then runs the established step-1
(clone labelling) and step-2 (annotation: #loci, multi-sample labels, node
colouring, concentric layers, batch marks, legend) logic without modifying it.

PIPELINE
  newick
    -> [newick_to_grapetree]  base.json  + base.svg
    -> [step 1: modify_new_samples_in_json]  labelled.json      (clone per node)
    -> [step 2: annotation functions]        annotated.svg
    -> [cairosvg]                            annotated.png / .pdf   (optional)

MINIMAL USAGE
  python tree_pipeline.py --newick tree.nwk --clone_map clones_map.tsv \
         --output annotated.svg

COMMON OPTIONS
  --color_map FILE        '#hex,clone' per line. A patterns SVG is auto-built.
  --cgmlst FILE           cgMLST tsv; adds a "#Loci: N" label.
  --batches FILE          'sample<TAB>batch' or 'sample,batch'; marks nodes.
  --legend_table FILE     tab-delimited legend layout (optional).
  --defaults FILE         override the bundled defaults.csv.
  --render png|pdf|both   also rasterise the final SVG.
  --width / --height      canvas size for the annotated SVG.
  --force                 overwrite existing outputs.

The script writes all intermediates next to --output (prefix = output stem)
so nothing is lost and every stage can be inspected or re-run manually.
'''

import argparse
import os
import sys
import csv
import json
from collections import Counter

import newick_to_grapetree as n2g
from grapetree_auxiliary_functions_v7 import *
import dan_auxiliary_functions as daf


# ----------------------------- small helpers ------------------------------- #

def sniff_delim(path):
    with open(path) as f:
        first = f.readline()
    return '\t' if first.count('\t') >= first.count(',') else ','


def to_csv_map(src, dst):
    '''Normalise a 2-column tsv/csv (sample,clone|batch) into a comma csv with
    header, as the downstream dictionary_from_csv_file helpers expect.'''
    delim = sniff_delim(src)
    rows = []
    with open(src) as f:
        r = csv.reader(f, delimiter=delim)
        rows = [row for row in r if row and len(row) >= 2]
    with open(dst, 'w', newline='') as f:
        w = csv.writer(f)
        for row in rows:
            w.writerow([row[0].strip(), row[1].strip()])
    return dst


def build_patterns_from_color_map(color_map_file, patterns_svg_out,
                                   pattern_clone_map_out):
    '''
    Convert a simple '#hex,clone' color map into:
      - an SVG <defs> patterns file (one solid-fill pattern per distinct colour), and
      - a pattern-clone map ('pattern,clone1,clone2,...').

    The colour map is the source of truth:
      * Clones with DIFFERENT hex values get different colours (own pattern each).
      * Clones sharing the SAME hex are grouped onto one pattern and are
        distinguished by concentric rings, in the order they appear in the file
        (first = plain node, second = one ring, third = two rings, ...).
    Clone names are never inspected, so cloneA1/cloneA2 only share a colour if
    you actually gave them the same hex.

    Returns (patterns_svg_out, pattern_clone_map_out).
    '''
    pairs = []
    with open(color_map_file) as f:
        for line in f:
            line = line.strip()
            if not line or ',' not in line:
                continue
            color, clone = [x.strip() for x in line.split(',', 1)]
            pairs.append((color, clone))

    # Group strictly by colour, preserving file order for both groups and members.
    groups = {}          # colour -> [clone, ...]  (ring order)
    order = []           # colours in first-seen order
    for color, clone in pairs:
        key = color.lower()
        if key not in groups:
            groups[key] = []
            order.append(key)
            
        groups[key].append(clone)

    # Pattern id: the single clone's name when a colour is unique, otherwise a
    # stable name derived from the group's first clone (keeps ids readable).
    pattern_id_of = {}
    for key in order:
        members = groups[key]
        pattern_id_of[key] = members[0] if len(members) == 1 else members[0] + '_grp'

    # patterns svg: one solid-fill pattern per distinct colour
    defs = ['<defs>']
    for key in order:
        defs.append(f'<pattern id="{pattern_id_of[key]}" width="10" height="10" '
                    f'patternUnits="userSpaceOnUse">')
        defs.append(f'<rect fill="{key}" height="100%" '
                    f'width="100%"/>')
        defs.append('</pattern>')
    defs.append('</defs>')
    with open(patterns_svg_out, 'w') as f:
        f.write('\n'.join(defs) + '\n')

    # pattern-clone map: pattern,clone1,clone2,...
    # Member order = file order, which sets the concentric-ring layer index
    # (first member plain, each subsequent member gains one more ring).
    with open(pattern_clone_map_out, 'w') as f:
        for key in order:
            f.write(','.join([pattern_id_of[key]] + groups[key]) + '\n')

    return patterns_svg_out, pattern_clone_map_out


def write_defaults_override(base_defaults, overrides, out_path):
    '''Read bundled defaults.csv, apply overrides, write a merged csv.'''
    d = daf.dictionary_from_csv_file(base_defaults)
    d.update({k: v for k, v in overrides.items() if v is not None})
    with open(out_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['argument', 'value'])
        for k, v in d.items():
            w.writerow([k, v])
    return out_path


def _fit_vis_to_canvas(svg_lines, target_w, target_h, legend_reserve_frac=0.35):
    '''
    The base tree is laid out at its own natural size (often much larger or
    smaller than the requested output canvas). Rather than just relabelling
    the outer <svg> width/height -- which clips or leaves nodes overlapping --
    rescale the <g id="vis" ...> content to fit inside the requested canvas,
    reserving a right-hand column (legend_reserve_frac of the width) so the
    legend has room and doesn't overlap the tree.
    '''
    src_w = src_h = None
    for line in svg_lines[:1]:
        if line.startswith('<svg'):
            import re
            mw = re.search(r'width="([\d.]+)"', line)
            mh = re.search(r'height="([\d.]+)"', line)
            if mw and mh:
                src_w, src_h = float(mw.group(1)), float(mh.group(1))
    if not src_w or not src_h:
        return svg_lines  # can't determine source size; leave as-is

    usable_w = target_w * (1 - legend_reserve_frac)
    usable_h = target_h
    scale = min(usable_w / src_w, usable_h / src_h)
    scale = max(scale, 1e-6)

    for i, line in enumerate(svg_lines):
        if line.startswith('<g id="vis"'):
            svg_lines[i] = f'<g id="vis" transform="translate(0,0) scale({scale})">'
            break
    return svg_lines

def run_step2(arguments, output_file):
    '''
    A faithful re-implementation of rearrange_tree_v15 step 2, driven by an
    already-assembled `arguments` dict (rather than defaults/differences csvs).
    Uses the unmodified grapetree_auxiliary_functions.
    '''
    json_data = get_json_data(arguments['json'])
    multi_sample_nodes = find_pie_charted_nodes(arguments['json'])
    svg_lines = split_svg_file(arguments['svg'])
    svg_lines = _fit_vis_to_canvas(svg_lines, arguments['svg_width'],
                                   arguments['svg_height'])
    clone_counter = None
    batches_marks = ['*', '%', '@', '$', '&', '|']

    if len(multi_sample_nodes) > 0:
        print(f"Adding multi-sample node labels ({len(multi_sample_nodes)} nodes)")
        svg_lines = add_multi_sample_nodes(svg_lines, multi_sample_nodes,
                                           arguments['spacing_factor'])
    else:
        print("No multi-sample nodes found.")

    if arguments['cgmlst'] != "NA":
        print("Adding #loci")
        svg_lines.insert(1, get_loci_string(arguments['cgmlst'],
                                            arguments['cgmlst_x'],
                                            arguments['cgmlst_y'],
                                            arguments['cgmlst_font_size']))
    else:
        print("No cgMLST file; skipping #loci.")

    node_idx = None
    if arguments['pattern_clone_legend_color_map'] != "NA":
        print("Colouring nodes by clone")
        sample_clone_map = _csv_map(arguments['clone_map'])
        clone_counter = Counter(sample_clone_map.values())
        clone_color_map = get_clone_color_map(
            arguments['pattern_clone_legend_color_map'])
        svg_lines = add_pattern_to_svg(svg_lines, arguments['patterns_file'])
        node_idx = get_mst_element_line_index(svg_lines)
        svg_lines = color_nodes_by_pattern(svg_lines, sample_clone_map,
                                           clone_color_map, node_idx)
        svg_lines = add_second_layer(svg_lines, sample_clone_map,
                                     clone_color_map, node_idx)
        node_idx = get_mst_element_line_index(svg_lines)
    else:
        print("No clone colour map; nodes left white.")

    batch_mark_map = None
    if arguments['sample_batch_map'] != "NA":
        print("Marking batches")
        sample_batch_map = _csv_map(arguments['sample_batch_map'])
        batch_mark_map = get_batch_mark_map(sample_batch_map, batches_marks)
        if node_idx is None:
            node_idx = get_mst_element_line_index(svg_lines)
        svg_lines = mark_batches(svg_lines, sample_batch_map,
                                 multi_sample_nodes, node_idx, batch_mark_map)
        node_idx = get_mst_element_line_index(svg_lines)
    else:
        print("No batches map; skipping batch marks.")

    if arguments['remove_samples_text'] is True:
        print("Removing sample text labels")
        svg_lines = remove_samples_text(svg_lines)

    have_legend = (arguments['legend_table'] != "NA"
                   or arguments['pattern_clone_legend_color_map'] != "NA"
                   or batch_mark_map is not None)
    if have_legend:
        legend_lines = get_legend_lines_with_table(arguments, clone_counter,
                                                   batch_mark_map)
    else:
        # Nothing to show; keep an empty legend group (avoids an upstream
        # add_batches() signature bug that only triggers on a fully-empty legend).
        legend_lines = ['<g class="legend" transform="translate('
                        + str(arguments['legend_x']) + ','
                        + str(arguments['legend_y']) + ')">', '</g>']
    legend_start = scale_start = -1
    for i, v in enumerate(svg_lines):
        if v.startswith('<g class="legend" transform'):
            legend_start = i
        if v.startswith('<g class="scale-bar"'):
            scale_start = i
    if legend_start > 0 and scale_start > 0:
        svg_lines = svg_lines[:legend_start] + legend_lines + svg_lines[scale_start:]
    svg_lines[0] = ('<svg xmlns="http://www.w3.org/2000/svg" width="'
                    + str(arguments['svg_width']) + '" height="'
                    + str(arguments['svg_height']) + '" id="mst-svg">')
    print_svg_lines(svg_lines, output_file)
    print("Created", output_file)


def _csv_map(csv_file):
    '''pandas-free reader for a header'd sample,value csv (mirrors the toolkit's
    dictionary_from_csv_file used inside step 2).'''
    import pandas as pd
    df = pd.read_csv(csv_file, header=0, index_col=None)
    c1, c2 = df.columns[0], df.columns[1]
    return {str(df.loc[i][c1]): df.loc[i][c2] for i in range(len(df))}


# --------------------------------- main ------------------------------------ #

def main():
    p = argparse.ArgumentParser(
        description="Newick -> annotated tree SVG, in one command.")
    p.add_argument('--newick', required=True, help="input Newick (.nwk/.tre) file")
    p.add_argument('--output', required=True, help="final annotated .svg path")
    p.add_argument('--clone_map', help="sample->clone map (tsv/csv)")
    p.add_argument('--color_map', help="'#hex,clone' colour map")
    p.add_argument('--patterns_file', help="pre-made SVG patterns (overrides color_map)")
    p.add_argument('--pattern_clone_map', help="pre-made pattern,clone,... map")
    p.add_argument('--cgmlst', help="cgMLST tsv for #loci label")
    p.add_argument('--batches', help="sample->batch map (tsv/csv)")
    p.add_argument('--legend_table', help="tab-delimited legend layout")
    default_defaults = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'defaults.csv')
    p.add_argument('--defaults', default=default_defaults, help="base defaults.csv")
    p.add_argument('--width', type=int, help="annotated svg width")
    p.add_argument('--height', type=int, help="annotated svg height")
    p.add_argument('--legend_x', type=int)
    p.add_argument('--legend_y', type=int)
    p.add_argument('--spacing_factor', type=int)
    p.add_argument('--remove_samples_text', action='store_true')
    p.add_argument('--legend_text', action='store_true')
    p.add_argument('--iterations', type=int, default=600,
                   help="layout iterations for the initial SVG")
    p.add_argument('--render', choices=['png', 'pdf', 'both'],
                   help="also rasterise the final svg")
    p.add_argument('--force', action='store_true')
    args = p.parse_args()

    if not args.output.endswith('.svg'):
        sys.exit("--output must end in .svg")
    if os.path.isfile(args.output) and not args.force:
        sys.exit("Output exists; use --force to overwrite.")

    outdir = os.path.dirname(os.path.abspath(args.output))
    stem = os.path.splitext(os.path.basename(args.output))[0]
    work = os.path.join(outdir, stem + "_work")
    os.makedirs(work, exist_ok=True)

    base_json = os.path.join(work, "base.json")
    base_svg = os.path.join(work, "base.svg")
    labelled_json = os.path.join(work, "labelled.json")

    # --- Stage A: Newick -> GrapeTree JSON + SVG ---
    print("== Stage A: Newick -> base tree ==")
    _, _, members = n2g.generate(args.newick, base_json, base_svg,
                                 layout_iterations=args.iterations)
    multi = {k: v for k, v in members.items() if len(v) > 1}
    print(f"  {len(members)} nodes ({len(multi)} multi-sample) laid out.")

    # --- Stage B: step 1 clone labelling (only if a clone map is given) ---
    json_for_step2 = base_json
    clone_csv = None
    if args.clone_map:
        print("== Stage B: labelling nodes with clones (step 1) ==")
        clone_csv = to_csv_map(args.clone_map, os.path.join(work, "clone_map.csv"))
        clone_dict = _csv_map(clone_csv)
        jd = get_json_data(base_json)
        jd = modify_new_samples_in_json(jd, clone_dict)
        output_json_data(jd, labelled_json)
        json_for_step2 = labelled_json
        print("  wrote", labelled_json)
    else:
        print("== Stage B: no clone map; skipping clone labelling ==")

    # --- prepare patterns/colour map for step 2 ---
    patterns_file = "NA"
    pattern_clone_map = "NA"
    if args.patterns_file and args.pattern_clone_map:
        patterns_file = args.patterns_file
        pattern_clone_map = args.pattern_clone_map
    elif args.color_map:
        print("  auto-building SVG patterns from colour map")
        patterns_file, pattern_clone_map = build_patterns_from_color_map(
            args.color_map,
            os.path.join(work, "patterns.svg"),
            os.path.join(work, "pattern_clone_map.txt"))

    batch_csv = "NA"
    if args.batches:
        batch_csv = to_csv_map(args.batches, os.path.join(work, "batch_map.csv"))

    # --- Stage C: annotation (step 2) ---
    print("== Stage C: annotating (step 2) ==")
    defaults = daf.dictionary_from_csv_file(args.defaults)

    def pick(key, cli, cast=str):
        if cli is not None:
            return cli
        return cast(defaults[key]) if key in defaults else None

    arguments = {
        'clone_map': clone_csv if clone_csv else "NA",
        'json': json_for_step2,
        'svg': base_svg,
        'spacing_factor': args.spacing_factor if args.spacing_factor is not None
                          else int(defaults.get('spacing_factor', 20)),
        'cgmlst': args.cgmlst if args.cgmlst else "NA",
        'cgmlst_x': int(defaults.get('cgmlst_x', 100)),
        'cgmlst_y': int(defaults.get('cgmlst_y', 50)),
        'cgmlst_font_size': int(defaults.get('cgmlst_font_size', 40)),
        'sample_batch_map': batch_csv,
        'remove_samples_text': args.remove_samples_text
                               or bool(int(defaults.get('remove_samples_text', 0))),
        'pattern_clone_legend_color_map': pattern_clone_map,
        'patterns_file': patterns_file,
        'legend_table': args.legend_table if args.legend_table else "NA",
        'legend_text': args.legend_text
                       or bool(int(defaults.get('legend_text', 0))),
        'svg_width': args.width if args.width else int(defaults.get('svg_width', 2500)),
        'svg_height': args.height if args.height else int(defaults.get('svg_height', 1500)),
        'legend_x': args.legend_x if args.legend_x is not None
                    else int(defaults.get('legend_x', 1100)),
        'legend_y': args.legend_y if args.legend_y is not None
                    else int(defaults.get('legend_y', 100)),
        'force': True,
    }
    run_step2(arguments, args.output)

    # --- Stage D: optional raster ---
    if args.render:
        try:
            import cairosvg
        except ImportError:
            print("cairosvg not installed; skipping raster.")
            return
        if args.render in ('png', 'both'):
            cairosvg.svg2png(url=args.output,
                             write_to=os.path.splitext(args.output)[0] + '.png',
                             output_width=1600)
            print("  wrote", os.path.splitext(args.output)[0] + '.png')
        if args.render in ('pdf', 'both'):
            cairosvg.svg2pdf(url=args.output,
                             write_to=os.path.splitext(args.output)[0] + '.pdf')
            print("  wrote", os.path.splitext(args.output)[0] + '.pdf')

    print("\nDone.  Final annotated tree:", args.output)
    print("Intermediates kept in:", work)


if __name__ == '__main__':
    main()
