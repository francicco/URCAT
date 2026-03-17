from __future__ import annotations

from collections import defaultdict, deque
from pathlib import Path
import json
import re


def read_json(path):
    with open(path) as fh:
        return json.load(fh)


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True)


class TreeNode:
    def __init__(self, name=None):
        self.name = name
        self.children = []
        self.parent = None

    def add_child(self, child):
        child.parent = self
        self.children.append(child)


def collect_leaves(node: TreeNode) -> list[str]:
    if node.name and not node.children:
        return [node.name]
    out = []
    for child in node.children:
        out.extend(collect_leaves(child))
    return out


def parse_newick_tree(newick: str) -> TreeNode:
    newick = newick.strip().rstrip(";")
    stack = []
    current = TreeNode()
    root = current
    token = ""

    def flush_token():
        nonlocal token, current
        t = token.strip()
        token = ""
        if not t:
            return
        name = t.split(":")[0]
        node = TreeNode(name=name)
        current.add_child(node)

    for ch in newick:
        if ch == "(":
            node = TreeNode()
            current.add_child(node)
            stack.append(current)
            current = node
        elif ch == ",":
            flush_token()
        elif ch == ")":
            flush_token()
            current = stack.pop()
        else:
            token += ch
    flush_token()

    while root.parent is not None:
        root = root.parent
    return root


def build_leaf_graph(root: TreeNode):
    graph = defaultdict(set)

    def connect_subtree(node):
        if not node.children:
            return collect_leaves(node)

        child_leaf_sets = [connect_subtree(c) for c in node.children]
        for i in range(len(child_leaf_sets)):
            for j in range(i + 1, len(child_leaf_sets)):
                for a in child_leaf_sets[i]:
                    for b in child_leaf_sets[j]:
                        graph[a].add(b)
                        graph[b].add(a)

        merged = []
        for s in child_leaf_sets:
            merged.extend(s)
        return merged

    connect_subtree(root)
    return graph


def shortest_path_distances(graph, start):
    dist = {start: 0}
    q = deque([start])
    while q:
        u = q.popleft()
        for v in graph.get(u, []):
            if v not in dist:
                dist[v] = dist[u] + 1
                q.append(v)
    return dist


def compute_reference_order(seed_species: str, species_tree_newick: str, species_list: list[str]) -> list[str]:
    root = parse_newick_tree(species_tree_newick)
    graph = build_leaf_graph(root)
    dist = shortest_path_distances(graph, seed_species)

    ordered = sorted(
        [sp for sp in species_list if sp != seed_species],
        key=lambda sp: (dist.get(sp, 10**9), sp),
    )
    return [seed_species] + ordered


def pick_next_reference_from_order(reference_order, used_reference_species, pending_frontiers_by_species):
    used = set(used_reference_species)
    for sp in reference_order:
        if sp in used:
            continue
        if pending_frontiers_by_species.get(sp):
            return sp
    return None


def extract_missing_locus_payloads(round_merged_json: dict):
    missing_by_target_species = defaultdict(list)

    for target_block in round_merged_json["targets"]:
        target_species = target_block["target_species"]

        for r in target_block["results"]:
            missing = r.get("missing_annotations", {})
            if target_species not in missing:
                continue

            intervals = missing[target_species]
            if not intervals:
                continue

            missing_by_target_species[target_species].append(
                {
                    "source_species": r["source_species"],
                    "source_transcript": r["source_transcript"],
                    "target_species": target_species,
                    "intervals": intervals,
                }
            )

    return dict(missing_by_target_species)


def interval_overlap(a_start, a_end, b_start, b_end):
    return not (a_end < b_start or b_end < a_start)


def locus_overlaps_any_interval(locus, intervals):
    for iv in intervals:
        if locus.seqid != iv["seqid"]:
            continue
        if interval_overlap(locus.start, locus.end, iv["start"], iv["end"]):
            return True
    return False
