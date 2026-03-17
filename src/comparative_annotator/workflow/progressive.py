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


def parse_newick_leaf_names(newick: str) -> list[str]:
    tokens = re.findall(r"[A-Za-z0-9_.:-]+|[(),;]", newick)
    leaves = []
    prev = None
    for tok in tokens:
        if tok in {",", "(", ")", ";"}:
            prev = tok
            continue
        if prev in {"(", ",", None}:
            leaves.append(tok.split(":")[0])
        prev = tok
    return leaves


class TreeNode:
    def __init__(self, name=None):
        self.name = name
        self.children = []
        self.parent = None

    def add_child(self, child):
        child.parent = self
        self.children.append(child)


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


def build_undirected_graph(root: TreeNode):
    graph = defaultdict(set)

    def visit(node):
        for child in node.children:
            if node.name and child.name:
                graph[node.name].add(child.name)
                graph[child.name].add(node.name)
            visit(child)

            if node.name and not child.name:
                for leaf in collect_leaves(child):
                    graph[node.name].add(leaf)
                    graph[leaf].add(node.name)

    def connect_internal(node):
        child_leaves = [collect_leaves(c) for c in node.children]
        flat = [x for sub in child_leaves for x in sub]
        for i in range(len(flat)):
            for j in range(i + 1, len(flat)):
                if flat[j] not in graph[flat[i]]:
                    pass

    visit(root)

    leaves = collect_leaves(root)
    induced = defaultdict(set)

    def connect_subtree(node):
        subleaves = collect_leaves(node)
        for child in node.children:
            child_leaves = collect_leaves(child)
            other = [x for x in subleaves if x not in child_leaves]
            for a in child_leaves:
                for b in other:
                    induced[a].add(b)
                    induced[b].add(a)
            connect_subtree(child)

    connect_subtree(root)
    return induced


def collect_leaves(node: TreeNode) -> list[str]:
    if node.name and not node.children:
        return [node.name]
    out = []
    for child in node.children:
        out.extend(collect_leaves(child))
    return out


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


def choose_next_reference_species(
    current_reference: str,
    used_reference_species: set[str],
    candidate_species_with_new_loci: set[str],
    species_tree_newick: str,
):
    if not candidate_species_with_new_loci:
        return None

    root = parse_newick_tree(species_tree_newick)
    graph = build_undirected_graph(root)
    dist = shortest_path_distances(graph, current_reference)

    candidates = [
        sp for sp in candidate_species_with_new_loci
        if sp not in used_reference_species
    ]
    if not candidates:
        return None

    candidates.sort(key=lambda sp: (dist.get(sp, 10**9), sp))
    return candidates[0]


def extract_missing_locus_payloads(round_merged_json: dict):
    """
    Returns:
      missing_by_target_species[target_species] = [
          {
              "source_species": ...,
              "source_transcript": ...,
              "target_species": ...,
              "intervals": [...],
          }, ...
      ]
    """
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
