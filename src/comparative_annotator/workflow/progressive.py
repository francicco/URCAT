from __future__ import annotations

from pathlib import Path
import json


def read_json(path):
    with open(path) as fh:
        return json.load(fh)


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True)


class TreeNode:
    def __init__(self, name=None, length=0.0):
        self.name = name
        self.length = float(length)
        self.children = []
        self.parent = None

    def add_child(self, child):
        child.parent = self
        self.children.append(child)


def collect_leaves(node: TreeNode) -> list[TreeNode]:
    if node.name and not node.children:
        return [node]
    out = []
    for child in node.children:
        out.extend(collect_leaves(child))
    return out


def parse_label(token: str):
    token = token.strip()
    if not token:
        return None, 0.0
    if ":" in token:
        name, length = token.split(":", 1)
        return name.strip() or None, float(length)
    return token, 0.0


def parse_newick_tree(newick: str) -> TreeNode:
    newick = newick.strip().rstrip(";")

    root = TreeNode()
    current = root
    stack = []
    token = ""

    def flush_token_as_child():
        nonlocal token, current
        tok = token.strip()
        token = ""
        if not tok:
            return
        name, length = parse_label(tok)
        child = TreeNode(name=name, length=length)
        current.add_child(child)

    i = 0
    while i < len(newick):
        ch = newick[i]

        if ch == "(":
            child = TreeNode()
            current.add_child(child)
            stack.append(current)
            current = child
            i += 1

        elif ch == ",":
            flush_token_as_child()
            i += 1

        elif ch == ")":
            flush_token_as_child()
            closed = current
            current = stack.pop()
            i += 1

            # parse optional internal node label / branch length
            label = []
            while i < len(newick) and newick[i] not in ",()":
                label.append(newick[i])
                i += 1
            label = "".join(label).strip()
            if label:
                name, length = parse_label(label)
                if name:
                    closed.name = name
                closed.length = float(length)

        else:
            token += ch
            i += 1

    if token.strip():
        flush_token_as_child()

    while root.parent is not None:
        root = root.parent
    return root


def get_leaf_nodes_by_name(root: TreeNode):
    leaves = collect_leaves(root)
    return {leaf.name: leaf for leaf in leaves if leaf.name is not None}


def path_to_root(node: TreeNode):
    path = []
    dist = 0.0
    while node is not None:
        path.append((node, dist))
        if node.parent is None:
            break
        dist += node.length
        node = node.parent
    return path


def pairwise_leaf_distance(a: TreeNode, b: TreeNode) -> float:
    path_a = path_to_root(a)
    path_b = path_to_root(b)

    dist_a = {node: d for node, d in path_a}
    dist_b = {node: d for node, d in path_b}

    common = set(dist_a).intersection(dist_b)
    lca_dist = min(dist_a[n] + dist_b[n] for n in common)
    return lca_dist


def compute_reference_order(seed_species: str, species_tree_newick: str, species_list: list[str]) -> list[str]:
    root = parse_newick_tree(species_tree_newick)
    leaf_map = get_leaf_nodes_by_name(root)

    if seed_species not in leaf_map:
        raise ValueError(f"Seed species '{seed_species}' not found in HAL tree")

    seed = leaf_map[seed_species]

    ranked = []
    for sp in species_list:
        if sp == seed_species:
            continue
        if sp not in leaf_map:
            raise ValueError(f"Species '{sp}' not found in HAL tree")
        d = pairwise_leaf_distance(seed, leaf_map[sp])
        ranked.append((d, sp))

    ranked.sort(key=lambda x: (x[0], x[1]))
    return [seed_species] + [sp for _, sp in ranked]


def pick_next_reference_from_order(reference_order, used_reference_species, pending_frontiers_by_species):
    used = set(used_reference_species)
    for sp in reference_order:
        if sp in used:
            continue
        if pending_frontiers_by_species.get(sp):
            return sp
    return None


def extract_missing_locus_payloads(round_merged_json: dict):
    missing_by_target_species = {}

    for target_block in round_merged_json["targets"]:
        target_species = target_block["target_species"]
        collected = []

        for r in target_block["results"]:
            missing = r.get("missing_annotations", {})
            if target_species not in missing:
                continue

            intervals = missing[target_species]
            if not intervals:
                continue

            collected.append(
                {
                    "source_species": r["source_species"],
                    "source_transcript": r["source_transcript"],
                    "target_species": target_species,
                    "intervals": intervals,
                }
            )

        if collected:
            missing_by_target_species[target_species] = collected

    return missing_by_target_species


def interval_overlap(a_start, a_end, b_start, b_end):
    return not (a_end < b_start or b_end < a_start)


def locus_overlaps_any_interval(locus, intervals):
    for iv in intervals:
        if locus.seqid != iv["seqid"]:
            continue
        if interval_overlap(locus.start, locus.end, iv["start"], iv["end"]):
            return True
    return False
