'''
newick_to_grapetree.py
Replaces the manual GrapeTree-website step of the tree pipeline.

Given a Newick tree (as produced by chewBBACA / grapetree MSTree), this module:
  1. Parses the tree into samples (leaves) and internal (hypothetical) nodes.
  2. Collapses zero-length edges between real samples into single "pie" nodes,
     mimicking GrapeTree's node-collapsing behaviour. The retained node keeps
     one representative sample as its shown ID; the other members are recorded
     so downstream code can label them / draw pie slices.
  3. Computes a 2-D layout of the collapsed tree with a deterministic
     force-directed (spring) layout on the minimum-spanning structure.
  4. Emits:
       - a JSON file with a 'metadata' block in the exact shape the existing
         rearrange_tree step-1 expects: metadata[id] = {ID, __Node, id, __selected}
         plus 'newickTree', 'initial_category', 'layout_data'.
       - an SVG file in the exact structure the existing annotation code keys on:
         <g class="node mst-element fixed" id="..."> ... <path ... fill="white">,
         <g class="link mst-element"> ... <line> + distance <text>,
         <text class="node-group-number"> for the shown sample of each node,
         a placeholder <g class="legend"> and <g class="scale-bar">.

The output is intentionally faithful to GrapeTree's SVG conventions so that
grapetree_auxiliary_functions_v7 can operate on it without modification.
'''

import json
import math
import re
import argparse
from collections import defaultdict


# ----------------------------- Newick parsing ------------------------------ #

class TNode:
    __slots__ = ('name', 'length', 'children', 'parent', 'id')
    def __init__(self, name=None, length=0.0):
        self.name = name
        self.length = length
        self.children = []
        self.parent = None
        self.id = None


def parse_newick(text):
    '''A small robust Newick parser. Returns the root TNode.'''
    text = text.strip()
    if text.endswith(';'):
        text = text[:-1]
    pos = 0

    def parse_clade():
        nonlocal pos
        node = TNode()
        if text[pos] == '(':
            pos += 1  # consume '('
            while True:
                child = parse_clade()
                child.parent = node
                node.children.append(child)
                if text[pos] == ',':
                    pos += 1
                    continue
                if text[pos] == ')':
                    pos += 1
                    break
        # read label
        m = re.match(r"[^,():;]*", text[pos:])
        label = m.group(0)
        pos += len(label)
        length = 0.0
        if pos < len(text) and text[pos] == ':':
            pos += 1
            m2 = re.match(r"[-+0-9.eE]+", text[pos:])
            length = float(m2.group(0))
            pos += len(m2.group(0))
        if label != '':
            node.name = label
        node.length = length
        return node

    root = parse_clade()
    return root


def collect_nodes(root):
    nodes = []
    def walk(n):
        nodes.append(n)
        for c in n.children:
            walk(c)
    walk(root)
    return nodes


# --------------------- Collapse zero-length sample edges -------------------- #

def build_graph(root):
    '''
    Build an undirected weighted graph of *named* nodes (samples). GrapeTree's
    MSTree output has samples both as leaves and, effectively, at internal
    positions. We treat every named node as a real sample and connect it to its
    nearest named ancestor/descendant with the summed branch length between them.
    Hypothetical (unnamed) internal nodes are contracted away.
    '''
    named = [n for n in collect_nodes(root) if n.name]
    # For each named node, walk up to the nearest named ancestor, summing lengths.
    edges = []
    for n in named:
        if n.parent is None:
            continue
        dist = n.length
        anc = n.parent
        while anc is not None and not anc.name:
            dist += anc.length
            anc = anc.parent
        if anc is not None and anc.name:
            edges.append((n.name, anc.name, dist))
        else:
            # nearest named node is through the root among siblings' subtrees;
            # connect sibling-cluster representatives via the unnamed ancestor.
            pass
    # Also connect named siblings that share an unnamed parent (star topology)
    # so the graph stays connected when the nearest named node isn't a direct
    # ancestor. We add edges between each named node and the first named node
    # reachable in its parent's subtree list.
    # Build adjacency for connectivity check.
    return named, edges


def nearest_named_descendant(node):
    '''BFS for the first named node in this subtree (used for star centers).'''
    stack = list(node.children)
    while stack:
        c = stack.pop(0)
        if c.name:
            return c
        stack.extend(c.children)
    return None


def build_full_graph(root):
    '''
    Produce a graph over ALL tree nodes, keeping unnamed internal nodes as
    distinct vertices (named _hypo_N), exactly as GrapeTree does. Every parent
    -> child branch becomes one edge with the child's branch length. This
    preserves topology so that only genuinely adjacent zero-length pairs get
    merged later, never chained transitively through a hypothetical node.
    '''
    all_nodes = collect_nodes(root)
    hypo = 0
    for n in all_nodes:
        if not n.name:
            n.name = f"_hypo_{hypo}"
            hypo += 1
    edges = []
    for n in all_nodes:
        for c in n.children:
            edges.append((n.name, c.name, c.length))
    names = [n.name for n in all_nodes]
    return names, edges


# --------------------------- Zero-distance merging -------------------------- #

def merge_zero_distance(names, edges):
    '''
    Union-find merge of samples connected by distance-0 edges into pie nodes.
    Returns:
      shown_of[sample] -> representative sample id (the shown node)
      members[rep] -> list of all sample ids in that node (rep first)
      collapsed_edges -> list of (repA, repB, dist) with dist>0, deduped (min dist)
    '''
    parent = {s: s for s in names}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for a, b, d in edges:
        if d == 0:
            union(a, b)

    # Choose a real (non-hypothetical) sample as each cluster's representative
    # when one exists, matching GrapeTree's "shown sample" behaviour.
    groups = defaultdict(list)
    for s in names:
        groups[find(s)].append(s)
    rep_of_root = {}
    for rootrep, mem in groups.items():
        real = [m for m in mem if not m.startswith('_hypo_')]
        rep_of_root[rootrep] = (real[0] if real else mem[0])

    reps = {s: rep_of_root[find(s)] for s in names}
    members = defaultdict(list)
    for s in names:
        members[reps[s]].append(s)
    # Drop hypothetical members from the shown membership lists; they are not
    # real samples. A node that is purely hypothetical is kept as a branch point.
    for r in list(members.keys()):
        reals = [m for m in members[r] if not m.startswith('_hypo_')]
        if reals:
            members[r] = [r] + [m for m in reals if m != r]

    collapsed = {}
    for a, b, d in edges:
        ra, rb = reps[a], reps[b]
        if ra == rb:
            continue
        key = tuple(sorted((ra, rb)))
        if key not in collapsed or d < collapsed[key]:
            collapsed[key] = d
    collapsed_edges = [(k[0], k[1], v) for k, v in collapsed.items()]

    # Some hypothetical nodes have no zero-distance real neighbour and so
    # survive as their own cluster (rep == '_hypo_N'). GrapeTree never draws an
    # unlabeled node in the tree -- only real samples get circles -- so we
    # contract every such purely-hypothetical cluster out of the drawable graph:
    # its incident edges are spliced together (distances summed) to connect its
    # real-sample neighbours directly. If a hypothetical point has 3+ neighbours
    # (a true branch/Steiner point), each neighbour pair gets a spliced edge so
    # connectivity and total path length are preserved without adding a node.
    reps, members, collapsed_edges = _contract_hypothetical_clusters(
        reps, members, collapsed_edges)

    return reps, members, collapsed_edges


def _contract_hypothetical_clusters(reps, members, collapsed_edges):
    hypo_clusters = {r for r in members if r.startswith('_hypo_')}
    if not hypo_clusters:
        return reps, members, collapsed_edges

    adj = defaultdict(dict)  # node -> {neighbour: min_dist}
    for a, b, d in collapsed_edges:
        if b not in adj[a] or d < adj[a][b]:
            adj[a][b] = d
            adj[b][a] = d

    remaining_hypo = set(hypo_clusters)
    changed = True
    while changed and remaining_hypo:
        changed = False
        for h in list(remaining_hypo):
            neighbours = list(adj[h].items())
            # splice: connect every pair of real neighbours directly
            for i in range(len(neighbours)):
                for j in range(i + 1, len(neighbours)):
                    nb_a, d_a = neighbours[i]
                    nb_b, d_b = neighbours[j]
                    new_d = d_a + d_b
                    if nb_b not in adj[nb_a] or new_d < adj[nb_a][nb_b]:
                        adj[nb_a][nb_b] = new_d
                        adj[nb_b][nb_a] = new_d
            for nb, _ in neighbours:
                del adj[nb][h]
            del adj[h]
            remaining_hypo.discard(h)
            changed = True

    # Rebuild collapsed_edges / members / reps without the contracted nodes.
    new_members = {r: m for r, m in members.items() if r not in hypo_clusters}
    new_reps = {s: r for s, r in reps.items() if r not in hypo_clusters}
    seen = set()
    all_edges = []
    for a in adj:
        if a in hypo_clusters:
            continue
        for b, d in adj[a].items():
            if b in hypo_clusters:
                continue
            key = tuple(sorted((a, b)))
            if key in seen:
                continue
            seen.add(key)
            all_edges.append((d, key[0], key[1]))

    # Reduce to MST so spring_layout sees a sparse tree-topology graph
    # rather than the dense transitive-closure produced by contraction.
    all_edges.sort()
    uf = {n: n for n in new_members}
    def _find(x):
        while uf[x] != x:
            uf[x] = uf[uf[x]]
            x = uf[x]
        return x
    new_edges = []
    for d, a, b in all_edges:
        ra, rb = _find(a), _find(b)
        if ra != rb:
            uf[rb] = ra
            new_edges.append((a, b, d))
    return new_reps, new_members, new_edges


# ------------------------------- Layout ------------------------------------ #

def _tree_bfs_layout(rep_nodes, collapsed_edges, k):
    '''
    Radial BFS tree layout used as the initial positions for spring refinement.
    Each node's depth in a BFS traversal determines its radius; angular wedges
    are allocated proportionally to subtree size. This gives a topology-aware
    starting point so force-directed iterations refine rather than reshape.
    '''
    from collections import deque, defaultdict as _dd

    if len(rep_nodes) <= 1:
        return {r: [0.0, 0.0] for r in rep_nodes}

    adj = _dd(list)
    for a, b, _d in collapsed_edges:
        adj[a].append(b)
        adj[b].append(a)

    pos = {}
    unvisited = set(rep_nodes)
    x_off = 0.0  # horizontal offset between disconnected components

    while unvisited:
        root = next(iter(unvisited))
        parent_of = {}
        children_of = _dd(list)
        order = [root]
        comp = {root}
        q = deque([root])
        while q:
            v = q.popleft()
            for nb in adj[v]:
                if nb not in comp:
                    comp.add(nb)
                    parent_of[nb] = v
                    children_of[v].append(nb)
                    order.append(nb)
                    q.append(nb)
        unvisited -= comp

        # Subtree sizes for proportional wedge allocation
        sz = {v: 1 for v in comp}
        for v in reversed(order):
            for ch in children_of[v]:
                sz[v] += sz[ch]

        # Assign angular wedges top-down, place nodes radially
        wedge = {root: (0.0, 2 * math.pi)}
        depth = {root: 0}
        pos[root] = [x_off, 0.0]

        for v in order[1:]:
            par = parent_of[v]
            depth[v] = depth[par] + 1
            lo, hi = wedge[par]
            total = sum(sz[c] for c in children_of[par])
            cur = lo
            for ch in children_of[par]:
                span = (hi - lo) * sz[ch] / total
                if ch == v:
                    mid = cur + span / 2
                    wedge[v] = (cur, cur + span)
                    break
                cur += span
            r = k * depth[v]
            pos[v] = [x_off + r * math.cos(mid), r * math.sin(mid)]

        max_r = max(k * depth[v] for v in comp)
        x_off += max_r * 2 + k * 3

    return pos


def spring_layout(rep_nodes, collapsed_edges, iterations=600, seed=42):
    '''
    Fruchterman-Reingold layout seeded from a BFS radial tree layout so the
    tree topology is already reflected in the starting positions. Attraction
    forces are weighted by log(edge_distance+1) so longer branches produce
    greater separation. Iterations are auto-capped for large trees to keep
    O(n^2) repulsion time manageable.
    '''
    import random
    rng = random.Random(seed)
    n = len(rep_nodes)

    area = (n * 80.0) ** 2
    k = math.sqrt(area / max(1, n))

    # Topology-aware initial layout; much better than a circle for large trees
    pos = _tree_bfs_layout(rep_nodes, collapsed_edges, k)
    # Tiny jitter to break exact symmetry
    for r in rep_nodes:
        pos[r][0] += rng.uniform(-k * 0.005, k * 0.005)
        pos[r][1] += rng.uniform(-k * 0.005, k * 0.005)

    adj = defaultdict(list)
    for a, b, d in collapsed_edges:
        adj[a].append((b, d))
        adj[b].append((a, d))

    # Starting temperature: fraction of ideal separation so the BFS layout
    # is refined rather than scrambled in the first few iterations.
    temp = k * 0.2
    # Cap iterations for large n to keep O(n^2) repulsion time manageable;
    # the BFS start means fewer iterations are needed anyway.
    effective_iter = min(iterations, max(200, iterations * 60 // max(60, n)))

    for it in range(effective_iter):
        disp = {r: [0.0, 0.0] for r in rep_nodes}
        # repulsion (all pairs)
        for i in range(n):
            ri = rep_nodes[i]
            for j in range(i + 1, n):
                rj = rep_nodes[j]
                dx = pos[ri][0] - pos[rj][0]
                dy = pos[ri][1] - pos[rj][1]
                dist = math.hypot(dx, dy) or 0.01
                force = k * k / dist
                ux, uy = dx / dist, dy / dist
                disp[ri][0] += ux * force; disp[ri][1] += uy * force
                disp[rj][0] -= ux * force; disp[rj][1] -= uy * force
        # attraction along edges weighted by log(d+1)
        for a, b, d in collapsed_edges:
            dx = pos[a][0] - pos[b][0]
            dy = pos[a][1] - pos[b][1]
            dist = math.hypot(dx, dy) or 0.01
            ideal = k * (1 + math.log1p(d))
            force = dist * dist / ideal
            ux, uy = dx / dist, dy / dist
            disp[a][0] -= ux * force; disp[a][1] -= uy * force
            disp[b][0] += ux * force; disp[b][1] += uy * force
        # apply with temperature cap
        for r in rep_nodes:
            dl = math.hypot(*disp[r]) or 0.01
            pos[r][0] += disp[r][0] / dl * min(dl, temp)
            pos[r][1] += disp[r][1] / dl * min(dl, temp)
        temp *= 0.99

    return pos


def normalize_positions(pos, margin=80, target_span=1200):
    xs = [p[0] for p in pos.values()]
    ys = [p[1] for p in pos.values()]
    minx, miny = min(xs), min(ys)
    spanx = max(xs) - minx or 1.0
    spany = max(ys) - miny or 1.0
    # Scale so the largest dimension equals target_span, keeping aspect ratio.
    sc = target_span / max(spanx, spany)
    out = {}
    for r, p in pos.items():
        out[r] = [(p[0] - minx) * sc + margin,
                  (p[1] - miny) * sc + margin]
    xs2 = [p[0] for p in out.values()]
    ys2 = [p[1] for p in out.values()]
    w = max(xs2) - min(xs2) + 2 * margin
    h = max(ys2) - min(ys2) + 2 * margin
    return out, w, h


# ------------------------------ Emit JSON ---------------------------------- #

def build_json(root, reps, members, collapsed_edges, pos, newick_text):
    # metadata: one entry per sample. ID=sample, __Node=shown rep of its cluster.
    metadata = {}
    real_samples = [s for s in reps.keys() if not s.startswith('_hypo_')]
    for i, s in enumerate(real_samples):
        metadata[s] = {
            "ID": s,
            "__Node": reps[s],
            "__selected": False,
            "id": i,
        }
    rep_list = list(members.keys())
    node_index = {r: i for i, r in enumerate(rep_list)}
    nodes = rep_list
    links = []
    for a, b, d in collapsed_edges:
        links.append({"source": node_index[a], "target": node_index[b],
                      "distance": d})
    data = {
        "nodes": nodes,
        "links": links,
        "metadata": metadata,
        "initial_category": "ID",
        "category_num": 1,
        "newickTree": newick_text.strip(),
        "layout_data": {"display_category": "ID"},
        "metadata_options": {},
        "_positions": pos,       # helper for the SVG writer (not used by grapetree)
        "_members": members,
    }
    return data


# ------------------------------ Emit SVG ----------------------------------- #

def node_radius(n_members, base=28.4, per=6.0):
    return base + per * math.log2(n_members) if n_members > 1 else base


def circle_path(r):
    return f"M0,{r}A{r},{r} 0 1,1 0,-{r}A{r},{r} 0 1,1 0,{r}Z"


def build_svg(reps, members, collapsed_edges, pos, canvas_w, canvas_h):
    rep_list = list(members.keys())
    parts = []
    parts.append(f'<svg xmlns="http://www.w3.org/2000/svg" '
                 f'width="{int(canvas_w)}" height="{int(canvas_h)}" id="mst-svg">')
    parts.append(f'<rect pointer-events="all" width="{int(canvas_w)}" '
                 f'height="{int(canvas_h)}" style="fill: none;"></rect>')
    parts.append('<g id="vis" transform="translate(0,0) scale(1)">')

    # links first (drawn under nodes)
    for a, b, d in collapsed_edges:
        x1, y1 = pos[a]; x2, y2 = pos[b]
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        parts.append(f'<g id="{a}-{b}" class="link mst-element">')
        parts.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" '
                     f'stroke-dasharray="" stroke-width="3px" '
                     f'style="stroke: black; opacity: 1;"></line>')
        parts.append(f'<text class="distance-label" dy=".71em" '
                     f'text-anchor="middle" font-size="26px" '
                     f'font-family="sans-serif" x="{mx}" y="{my}" '
                     f'style="fill: rgb(102,102,102); stroke: white; '
                     f'stroke-width: 5px; opacity: 0.9;">{int(d)}</text>')  # self-contained
        parts.append('</g>')

    # nodes
    for r in rep_list:
        x, y = pos[r]
        rad = node_radius(len(members[r]))
        parts.append(f'<g class="node mst-element fixed" id="{r}" '
                     f'transform="translate({x},{y})">')
        parts.append(f'<path class="node-paths" d="{circle_path(rad)}" '
                     f'fill="white" style="stroke: black;"></path>')
        # shown sample label (grapetree shows the representative sample only)
        parts.append(f'<text class="node-group-number" dy=".71em" '
                     f'text-anchor="middle" font-size="12" '
                     f'font-family="sans-serif" '
                     f'transform="translate(0,-4)">{r}</text>')
        parts.append('</g>')

    # placeholder legend + scale-bar so the annotation code can splice them
    parts.append('<g class="legend" transform="translate(0,0)">')
    parts.append('</g>')
    parts.append('<g class="scale-bar" transform="translate(40,'
                 + str(int(canvas_h) - 40) + ')">')
    parts.append('</g>')
    parts.append('</g>')  # vis
    parts.append('</svg>')
    # GrapeTree emits the whole SVG on a single physical line; split_svg_file
    # (which joins physical lines with spaces then re-splits on '<') then yields
    # clean tokens with no trailing spaces. Join without newlines to match.
    return "".join(parts)


# ------------------------------- Driver ------------------------------------ #

def generate(newick_file, out_json, out_svg, layout_iterations=600):
    with open(newick_file) as f:
        newick_text = f.read()
    root = parse_newick(newick_text)
    names, edges = build_full_graph(root)
    reps, members, collapsed_edges = merge_zero_distance(names, edges)

    rep_list = list(members.keys())
    pos = spring_layout(rep_list, collapsed_edges, iterations=layout_iterations)
    pos, w, h = normalize_positions(pos)

    data = build_json(root, reps, members, collapsed_edges, pos, newick_text)
    # strip helper keys before writing json used by grapetree step-1
    json_out = {k: v for k, v in data.items()
                if k not in ('_positions', '_members')}
    with open(out_json, 'w') as f:
        json.dump(json_out, f)

    svg = build_svg(reps, members, collapsed_edges, pos, w, h)
    with open(out_svg, 'w') as f:
        f.write(svg)
    return out_json, out_svg, members


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--newick', required=True)
    p.add_argument('--out_json', required=True)
    p.add_argument('--out_svg', required=True)
    p.add_argument('--iterations', type=int, default=600)
    a = p.parse_args()
    generate(a.newick, a.out_json, a.out_svg, a.iterations)
    print("Wrote", a.out_json, "and", a.out_svg)
