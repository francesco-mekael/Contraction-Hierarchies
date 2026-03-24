import heapq
import copy

# Maps shortcut edges to their middle nodes for path unpacking
SHORTCUT_MIDDLE_NODE = {}


def read_file(filename: str, directed: bool) -> dict[int, dict[int, int]]:
    """
    Reads a graph in DIMACS format.

    Constructs either a directed or undirected graph.
    Each edge line follows the format: 'a v1 v2 weight'

    Args:
        filename: Path to the DIMACS graph file.
        directed: If True, builds a directed graph; otherwise undirected.

    Returns:
        Adjacency list as dict[node -> {neighbor: weight}]
    """
    graph = {}

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if not line.startswith('a'):
                continue

            edge = line.split()
            u, v, weight = int(edge[1]), int(edge[2]), int(edge[3])

            graph.setdefault(u, {})
            graph.setdefault(v, {})
            graph[u][v] = weight

            if not directed:
                graph[v][u] = weight

    return graph


def visualize(g: dict[int, dict[int, int]], directed: bool = False, print_weights: bool = True):
    """
    Generates a DOT-format string for use with GraphViz.
    Intended for debugging and quick visualization of small graphs.
    """
    vis = ("graph" if not directed else "digraph") + " G {\n"
    for u in g:
        for v in g[u]:
            vis += f"\t{u} -- {v}" if not directed else f"\t{u} -> {v}"
            if print_weights:
                vis += f' [label="{g[u][v]}"]'
            vis += '\n'
    vis += "}\n"
    print(vis)


def dijkstra_s_t(g: dict[int, dict[int, int]], s: int, t: int) -> tuple[int, int, list[int]]:
    """
    Standard Dijkstra's algorithm from source s to target t.

    Used as a reference implementation to verify CH results.

    Args:
        g: Adjacency list graph.
        s: Source node.
        t: Target node.

    Returns:
        Tuple of (shortest distance, number of expanded nodes, path as list of nodes).
        Distance is -1 if t is unreachable from s.
    """
    priority_queue = [(0, s)]
    distances = {s: 0}
    predecessor = {s: None}
    expanded_nodes = 0
    final_distance = -1
    path = []

    while priority_queue:
        current_distance, current_node = heapq.heappop(priority_queue)

        if current_distance > distances.get(current_node, float('inf')):
            continue

        expanded_nodes += 1

        if current_node == t:
            final_distance = current_distance
            break

        for neighbor, weight in g.get(current_node, {}).items():
            new_distance = current_distance + weight
            if new_distance < distances.get(neighbor, float('inf')):
                distances[neighbor] = new_distance
                predecessor[neighbor] = current_node
                heapq.heappush(priority_queue, (new_distance, neighbor))

    if final_distance != -1:
        curr = t
        while curr is not None:
            path.append(curr)
            curr = predecessor[curr]
        path.reverse()

    return final_distance, expanded_nodes, path


def limited_dijkstra(g: dict[int, dict[int, int]], s: int, limit: int) -> tuple[dict[int, int], dict[int, int]]:
    """
    Distance-limited Dijkstra from source s. Stops once all remaining
    paths exceed the given distance limit.

    Used internally during node contraction to detect unnecessary shortcuts.

    Args:
        g: Adjacency list graph.
        s: Source node.
        limit: Maximum distance to explore.

    Returns:
        Tuple of (distance dict, predecessor dict).
    """
    priority_queue = [(0, s)]
    L = {node: float('inf') for node in g}
    L[s] = 0
    predecessor = {node: None for node in g}

    while priority_queue:
        current_distance, current_node = heapq.heappop(priority_queue)

        if current_distance > limit:
            break

        if current_distance > L.get(current_node, float('inf')):
            continue

        for neighbor, weight in g.get(current_node, {}).items():
            new_distance = current_distance + weight
            if new_distance > limit:
                continue
            if new_distance < L.get(neighbor, float('inf')):
                L[neighbor] = new_distance
                predecessor[neighbor] = current_node
                heapq.heappush(priority_queue, (new_distance, neighbor))

    return L, predecessor


def contract_node(g: dict[int, dict[int, int]], v: int) -> dict[tuple[int, int], int]:
    """
    Contracts a single node v from the graph.

    For each pair (u, w) where u -> v -> w, checks whether removing v
    creates a shortest-path gap. If so, a shortcut edge u -> w is added.

    Args:
        g: Current working graph.
        v: Node to contract.

    Returns:
        Dict of {(u, w): shortcut_weight} for all required shortcuts.
    """
    global SHORTCUT_MIDDLE_NODE

    contractions = {}
    incoming_nodes = [u for u in g if v in g[u]]
    outgoing_nodes = g.get(v, {})

    backup_outgoing = g[v]
    g[v] = {}  # Temporarily remove v's outgoing edges

    for u in incoming_nodes:
        dist_u_v = g[u][v]
        max_search_limit = 0
        potential_shortcuts = {}

        for w, dist_v_w in outgoing_nodes.items():
            if u == w:
                continue
            total_distance = dist_u_v + dist_v_w
            potential_shortcuts[w] = total_distance
            if total_distance > max_search_limit:
                max_search_limit = total_distance

        if not potential_shortcuts:
            continue

        distances, _ = limited_dijkstra(g, u, max_search_limit)

        for w, shortcut_distance in potential_shortcuts.items():
            if w in distances and distances[w] <= shortcut_distance:
                continue
            existing = contractions.get((u, w), float('inf'))
            if shortcut_distance < existing:
                contractions[(u, w)] = shortcut_distance
                SHORTCUT_MIDDLE_NODE[(u, w)] = v

    g[v] = backup_outgoing  # Restore v's outgoing edges
    return contractions


def contraction(g: dict[int, dict[int, int]]) -> tuple[dict[int, dict[int, int]], dict[int, int]]:
    """
    Full Contraction Hierarchies preprocessing.

    Contracts all nodes one by one in order of increasing degree.
    Ties in degree are broken by node ID (smaller ID first).
    Adds shortcut edges as needed to preserve shortest path distances.

    Args:
        g: Original directed graph as adjacency list.

    Returns:
        Tuple of:
            - g_ch: Augmented graph with original edges + shortcuts.
            - rank: Dict mapping each node to its contraction order (rank).
    """
    global SHORTCUT_MIDDLE_NODE
    SHORTCUT_MIDDLE_NODE.clear()

    g_ch = copy.deepcopy(g)
    working_g = copy.deepcopy(g)

    rank = {}
    rank_counter = 1

    # Build reverse graph for in-degree tracking
    g_rev_working = {}
    for u, neighbors in working_g.items():
        for v in neighbors:
            g_rev_working.setdefault(v, [])
            g_rev_working[v].append(u)

    while working_g:
        # Select node with minimum degree (in + out), break ties by node ID
        best_node = min(
            working_g,
            key=lambda v: (len(working_g[v]) + len(g_rev_working.get(v, [])), v)
        )

        rank[best_node] = rank_counter
        rank_counter += 1

        shortcuts = contract_node(working_g, best_node)

        # Add shortcuts to both graphs
        for (u, w), weight in shortcuts.items():
            g_ch.setdefault(u, {})
            if w not in g_ch[u] or weight < g_ch[u][w]:
                g_ch[u][w] = weight

            working_g.setdefault(u, {})
            if w not in working_g[u] or weight < working_g[u][w]:
                working_g[u][w] = weight

            g_rev_working.setdefault(w, [])
            if u not in g_rev_working[w]:
                g_rev_working[w].append(u)

        # Remove edges pointing to best_node
        if best_node in g_rev_working:
            for u in g_rev_working[best_node]:
                if u in working_g and best_node in working_g[u]:
                    del working_g[u][best_node]
            del g_rev_working[best_node]

        for neighbor in working_g[best_node]:
            if neighbor in g_rev_working and best_node in g_rev_working[neighbor]:
                g_rev_working[neighbor].remove(best_node)

        del working_g[best_node]

    return g_ch, rank


def bidirectional_dijkstra(g: dict[int, dict[int, int]], rank: dict[int, int], s: int, t: int) -> tuple[int, int, list[int]]:
    """
    CH-aware bidirectional Dijkstra for shortest path queries.

    Forward search expands nodes with higher rank than the current node.
    Backward search traverses the reverse graph upward in rank.

    Stopping criterion: min(L_up, L_down) > best_path_found

    Args:
        g: CH-augmented graph (original edges + shortcuts).
        rank: Node rank dict from contraction().
        s: Source node.
        t: Target node.

    Returns:
        Tuple of (shortest distance, expanded nodes, full unpacked path).
        Distance is -1 if unreachable.
    """
    # Build reverse graph
    g_rev = {}
    for u, neighbors in g.items():
        for v, w in neighbors.items():
            g_rev.setdefault(v, {})
            g_rev[v][u] = w

    pq_up = [(0, s)]
    pq_down = [(0, t)]
    dist_up = {s: 0}
    dist_down = {t: 0}
    parent_up = {s: None}
    parent_down = {t: None}
    settled_up = set()
    settled_down = set()

    best = float('inf')
    meeting_node = None
    expanded_nodes = 0

    while pq_up or pq_down:
        L_up = pq_up[0][0] if pq_up else float('inf')
        L_down = pq_down[0][0] if pq_down else float('inf')

        if min(L_up, L_down) > best:
            break

        if L_up <= L_down:
            d, u = heapq.heappop(pq_up)
            if u in settled_up:
                continue
            settled_up.add(u)
            expanded_nodes += 1

            if u in dist_down:
                total = d + dist_down[u]
                if total < best:
                    best = total
                    meeting_node = u

            for v, w in g.get(u, {}).items():
                if rank[v] > rank[u]:
                    new_d = d + w
                    if new_d < dist_up.get(v, float('inf')):
                        dist_up[v] = new_d
                        parent_up[v] = u
                        heapq.heappush(pq_up, (new_d, v))
        else:
            d, u = heapq.heappop(pq_down)
            if u in settled_down:
                continue
            settled_down.add(u)
            expanded_nodes += 1

            if u in dist_up:
                total = d + dist_up[u]
                if total < best:
                    best = total
                    meeting_node = u

            for v, w in g_rev.get(u, {}).items():
                if rank[v] > rank[u]:
                    new_d = d + w
                    if new_d < dist_down.get(v, float('inf')):
                        dist_down[v] = new_d
                        parent_down[v] = u
                        heapq.heappush(pq_down, (new_d, v))

    if meeting_node is None:
        return -1, expanded_nodes, []

    # Reconstruct compressed path
    path_ch = []
    curr = meeting_node
    while curr is not None:
        path_ch.append(curr)
        curr = parent_up.get(curr)
    path_ch.reverse()

    curr = parent_down.get(meeting_node)
    while curr is not None:
        path_ch.append(curr)
        curr = parent_down.get(curr)

    # Unpack shortcut edges into full path
    def unpack_edge(u, w):
        if (u, w) in SHORTCUT_MIDDLE_NODE:
            mid = SHORTCUT_MIDDLE_NODE[(u, w)]
            return unpack_edge(u, mid)[:-1] + unpack_edge(mid, w)
        return [u, w]

    full_path = [path_ch[0]]
    for i in range(len(path_ch) - 1):
        unpacked = unpack_edge(path_ch[i], path_ch[i + 1])
        full_path.extend(unpacked[1:])

    return best, expanded_nodes, full_path


if __name__ == '__main__':
    import time
    import random

    print("=" * 50)
    print("Contraction Hierarchies – Benchmark on rome99.gr")
    print("=" * 50)

    print("\n[1] Loading graph...")
    start = time.time()
    graph = read_file('rome99.gr', True)
    print(f"    Loaded in {time.time() - start:.4f}s | Nodes: {len(graph)}")

    print("\n[2] Preprocessing (contraction)...")
    start = time.time()
    g_ch, rank = contraction(graph)
    processing_time = time.time() - start
    print(f"    Processing time: {processing_time:.4f}s")

    print("\n[3] Shortcut statistics...")
    original_edges = sum(len(n) for n in graph.values())
    ch_edges = sum(len(n) for n in g_ch.values())
    print(f"    Original edges:  {original_edges}")
    print(f"    Edges with shortcuts: {ch_edges}")
    print(f"    Shortcuts inserted: {ch_edges - original_edges}")

    print("\n[4] Running 1000 random queries...")
    nodes = list(graph.keys())
    num_queries = 1000
    total_dijkstra, total_ch = 0, 0

    for _ in range(num_queries):
        s, t = random.choice(nodes), random.choice(nodes)
        _, exp_d, _ = dijkstra_s_t(graph, s, t)
        _, exp_ch, _ = bidirectional_dijkstra(g_ch, rank, s, t)
        total_dijkstra += exp_d
        total_ch += exp_ch

    print(f"\n    Results over {num_queries} queries:")
    print(f"    Avg expanded nodes (Dijkstra): {total_dijkstra / num_queries:.1f}")
    print(f"    Avg expanded nodes (CH Bidir): {total_ch / num_queries:.1f}")
    speedup = (total_dijkstra / num_queries) / max(total_ch / num_queries, 1)
    print(f"    Speedup factor: {speedup:.1f}x")
