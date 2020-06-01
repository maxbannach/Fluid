use crate::graph::{Graph, Decomposition, Bipgraph, GLOBAL_ABORT};
use std::collections::{HashMap, HashSet, BinaryHeap};
use rand::{Rng, seq::SliceRandom};
use std::cmp::Reverse;
use union_find::{UnionFind, QuickUnionUf, UnionBySize};
use std::sync::atomic::Ordering;
//use std::time::Instant;

// This code is ported from eppstein PADS.
// https://www.ics.uci.edu/~eppstein/PADS/
fn matching(g: &Bipgraph) -> HashMap<usize, usize> {
    // Find maximum cardinality matching of a bipartite graph (U,V,E).
    // The input format is a dictionary mapping members of U to lists
    // of their neighbors in V.  The output is a triple (M,A,B) where M is a
    // dictionary mapping members of V to their matches in U, A is the part
    // of the maximum independent set in U, and B is the part of the MIS in V.
    // The same object may occur in both U and V, and is treated as two
    // distinct vertices if this happens.

    // initialize greedy matching (redundant, but faster than full search)
    let mut matching: HashMap<usize, usize> = HashMap::new();
    for &u in &g.set_left {
        for &v in g.neighbors(u) {
            if !matching.contains_key(&v) {
                matching.insert(v, u);
                break;
            }
        }
    }

    loop {
        // structure residual graph into layers
        // pred[u] gives the neighbor in the previous layer for u in U
        // preds[v] gives a list of neighbors in the previous layer for v in V
        // unmatched gives a list of unmatched vertices in final layer of V,
        // and is also used as a flag value for pred[u] when u is in the first layer
        let mut preds: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut pred: HashMap<usize, usize> = HashMap::new();
        let mut unmatched: Vec<usize> = Vec::new();
        for &u in &g.set_left {
            pred.insert(u, u);
        }
        for v in matching.keys() {
            pred.remove(&matching[v]);
        }
        let mut layer: Vec<usize> = pred.keys().map(|v| *v).collect();

        // repeatedly extend layering structure by another pair of layers
        while layer.len() > 0 && unmatched.is_empty() {
            let mut new_layer: HashMap<usize, Vec<usize>> = HashMap::new();
            for u in layer {
                for v in g.neighbors(u) {
                    if !preds.contains_key(v) {
                        new_layer.entry(*v).or_insert(Vec::new()).push(u);
                    }
                }
            }
            layer = Vec::new();
            for &v in new_layer.keys() {
                preds.insert(v, new_layer[&v].clone());
                if let Some(&m) = matching.get(&v) {
                    layer.push(m);
                    pred.insert(m, v);
                } else {
                    unmatched.push(v);
                }
            }
        }

        // did we finish layering without finding any alternating paths?
        if unmatched.is_empty() {
            let mut unlayered = HashSet::new();
            for &u in &g.set_left {
                for v in g.neighbors(u) {
                    if !preds.contains_key(v) {
                        unlayered.insert(v);
                    }
                }
            }
            return matching;
        }

        // recursively search backward through layers to find alternating paths
        // recursion returns true if found path, false otherwise
        // Note: Rust does not allow to capture the context, as it was done in the python code.
        fn recurse(v: &usize, matching: &mut HashMap<usize, usize>, pred: &mut HashMap<usize, usize>, preds: &mut HashMap<usize, Vec<usize>>) -> bool {
            if let Some(v_list) = preds.remove(v) {
                for u in v_list.iter() {
                    if pred.contains_key(u) {
                        let pu = pred[u];
                        pred.remove(u);
                        if pu == *u || recurse(&pu, matching, pred, preds) {
                            matching.insert(*v, *u);
                            return true;
                        }
                    }
                }
            }
            return false;
        }
        
        for v in unmatched {
            recurse(&v, &mut matching, &mut pred, &mut preds);
        }
    }
}

// This is a Rust implementation of the algorithm described at
// <https://en.wikipedia.org/wiki/K%C3%B6nig%27s_theorem_%28graph_theory%29#Proof>.
fn to_vc(g: &Bipgraph, matching: &HashMap<usize, usize>) -> HashSet<usize> {
    let set_l = &g.set_left;
    let set_r = &g.set_right;
    let set_matching_l: HashSet<usize> = matching.values().cloned().collect();
    // Let U be the set of unmatched vertices in the left vertex set.
    let set_u: Vec<usize> = set_l.difference(&set_matching_l).cloned().collect();
    
    let _alternating_dfs = |matching_first: bool| {
        let mut visited = HashSet::new();
        let mut stack: Vec<(usize, bool)> = Vec::new();
        for &u in set_u.iter() {
            stack.push((u, matching_first));
            visited.insert(u);
        }
        while let Some((v, use_matching)) = stack.pop() {
            for &w in g.neighbors(v) {
                if visited.contains(&w) {
                    continue;
                }
                if !use_matching ^ (matching.get(&w) == Some(&v) || matching.get(&v) == Some(&w)) {
                    stack.push((w, !use_matching));
                    visited.insert(w);
                }
            }
        }
        return visited;
    };

    // Let Z be the set of vertices that are either in U or are connected to U
    // by alternating paths.
    let set_z = _alternating_dfs(true).union(&_alternating_dfs(false)).cloned().collect();
    // At this point, every edge either has a right endpoint in Z or a left
    // endpoint not in Z. This gives us the vertex cover.
    let diff_l_z: HashSet<usize> = set_l.difference(&set_z).cloned().collect();
    let inter_r_z: HashSet<usize> = set_r.intersection(&set_z).cloned().collect();
    return diff_l_z.union(&inter_r_z).cloned().collect();
}

// Given a graph and two communities, compute a minimum vertex separator.
// The two communities may separate a subgraph of g instead of the entire graph.
pub fn vertex_separator(g: &Graph, comm: &Vec<Vec<usize>>) -> HashSet<usize> {
    let bg = Bipgraph::new(&g, &comm[0], &comm[1]);
    let mat: HashMap<usize, usize> = matching(&bg);
    let vc = to_vc(&bg, &mat);
    vc
}

// Given a graph g, a set (here as vector) of vertices describing a subgraph, and the number of
// communities k, find k communities.
// Ported from networkx:
// https://networkx.github.io/documentation/stable/_modules/networkx/algorithms/community/asyn_fluid.html#asyn_fluidc
pub fn async_fluidc(g: &Graph, vert: &[usize], k: usize) -> Vec<Vec<usize>> {
    // Initialization
    let max_iter = 10000;
    let mut vertices = vert.to_vec(); // "Clone"
    let mut rng = rand::thread_rng();
    vertices.shuffle(&mut rng);
    let mut communities = HashMap::new();
    let mut com_to_numvertices = vec![0; k];
    let mut com_counter = HashMap::new();
    for v in vertices.iter() {
        com_counter.insert(v, vec![0; k]);
    }
    let mut touched = HashSet::new();
    for (i, &v) in vertices[0..k].iter().enumerate() {
        communities.insert(v, i);
        if let Some(vec) = com_counter.get_mut(&v) {
            vec[i] += 1;
        }
        touched.insert(v);
        for &w in g.neighbors[v].iter() {
            if let Some(vec) = com_counter.get_mut(&w) {
                vec[i] += 1;
            }
            touched.insert(w);
        }
        com_to_numvertices[i] = 1;
    }
    // Set up control variables and start iterating
    let mut iter_count = 0;
    let mut cont = true;
    let mut com_freq = vec![0.0; k];
    while cont {
        cont = false;
        iter_count += 1;
        // Loop over all vertices in graph in a random order
        let mut last_touched: Vec<usize> = touched.into_iter().collect();
        last_touched.shuffle(&mut rng);
        touched = HashSet::new();
        for vertex in last_touched {
            // Updating rule
            for com in 0..k {
                com_freq[com] = (com_counter[&vertex][com] as f64) / (com_to_numvertices[com] as f64);
            }
            // Check which is the community with highest density
            let max_freq = com_freq.iter().cloned().fold(-1./0. /* -inf */, f64::max);
            if max_freq > 0.0 {
                let mut best_communities = Vec::new();
                for com in 0..k {
                    if max_freq == com_freq[com] {
                        best_communities.push(com);
                    }
                }
                // If actual vertex com in best communities, it is preserved
                if let Some(com) = communities.get(&vertex) {
                    if best_communities.contains(&com) {
                        continue;
                    }
                }
                // If vertex community changes...
                // Set flag of non-convergence
                cont = true;
                // Randomly chose a new community from candidates
                let new_com = match best_communities.choose(&mut rng) {
                    Some(com) => *com,
                    None => continue // Should not happen.
                };

                // Update previous community status
                if let Some(&com) = communities.get(&vertex) {
                    com_to_numvertices[com] -= 1;
                    if let Some(vec) = com_counter.get_mut(&vertex) {
                        vec[com] -= 1;
                    }
                    for &v in g.neighbors[vertex].iter() {
                        if let Some(vec) = com_counter.get_mut(&v) {
                            vec[com] -= 1;
                        }
                    }
                }
                // Update new community status
                communities.insert(vertex, new_com);
                com_to_numvertices[new_com] += 1;
                if let Some(vec) = com_counter.get_mut(&vertex) {
                    vec[new_com] += 1;
                }
                touched.insert(vertex);
                for &v in g.neighbors[vertex].iter() {
                    touched.insert(v);
                    if let Some(vec) =  com_counter.get_mut(&v) {
                        vec[new_com] += 1;
                    }
                }
            }
        }
        // If maximum iterations reached --> output actual results
        if iter_count > max_iter {
            break;
        }
    }
    // Return results by grouping communities as list of vertices
    if iter_count > 9000 {
        eprintln!("{}", iter_count);
    }
    let mut res = vec![vec![]; k];
    for (&i, &com) in communities.iter() {
        res[com].push(i);
    }
    res
}

// Remove a separator set from g. Update decomposition tree.
// Returns the connected components in the subgraph given by "vertices".
fn remove_separator(g: &mut Graph,  vertices: &[usize], sep: &HashSet<usize>, dec: &mut Decomposition, parent: Option<usize>) -> (Vec<Vec<usize>>, Option<usize>) {
    let mut parent = parent;

    // Add separator as path and remove vertices from g.
    for &v in sep.iter() {
        if let Some(_) = dec.parent[v] {
            continue;
        }
        dec.parent[v] = parent;
        parent = Some(v);
        let (deg_zero, _deg_one) = g.remove_vertex(v);
        
        // Remove isolated vertices from g as well. Keep track of deg one vertices.
        for &w in deg_zero.iter() {
            dec.parent[w] = parent;
        }
        //remove_deg_one(&mut deg_one, g, dec, parent);
    }

    let remaining_vertices: Vec<usize> = vertices.iter().filter(|&v| g.neighbors[*v].len() > 0).map(|&v| v).collect();
    let cc = g.connected_components(&remaining_vertices);
    (cc, parent)
}

// Find a separator via the fluid algorithm or the greedy algorithm and remove it from g.
// Return the resulting connected components and decomposition parents.
pub fn greedy_fluid_separation(g: &mut Graph, vertices: &[usize], dec: &mut Decomposition, parent: Option<usize>, greedy_rep: usize, fluid_rep: usize) -> Result<(Vec<Vec<usize>>, Option<usize>), ()> {
    let mut sep = HashSet::new();
    for _ in 0..fluid_rep {
        if GLOBAL_ABORT.load(Ordering::Relaxed) {
            return Err(());
        }
        let cfcur = async_fluidc(&g, &vertices, 2);
        let sepcur = vertex_separator(&g, &cfcur);
        if sepcur.len() < sep.len() || sep.len() == 0 {
           sep = sepcur;
        }
    }

    let mut sep_greedy = HashSet::new();
    let mut best_score = -1.0;
    for _ in 0..greedy_rep {
        if GLOBAL_ABORT.load(Ordering::Relaxed) {
            return Err(());
        }
        let (cur_sep, cur_score) = greedy_separator(g, vertices);
        if best_score < cur_score {
            sep_greedy = cur_sep;
            best_score = cur_score;
        }
    }
    if sep_greedy.len() < sep.len() || sep.len() == 0 {
        sep = sep_greedy;
    }

    Ok(remove_separator(g, vertices, &sep, dec, parent))
}

fn remove_deg_one(deg_one_stack: &mut Vec<usize>, g: &mut Graph, dec: &mut Decomposition, parent: Option<usize>) {
    // Iteratively remove deg one vertices.
    while let Some(v) = deg_one_stack.pop() {
        if let Some(_) = dec.parent[v] {
            continue;
        }
        if g.neighbors[v].len() > 0 {
            dec.parent[v] = Some(g.neighbors[v][0]);
        } else {
            dec.parent[v] = parent;
        }
        let (deg_zero, mut deg_one) = g.remove_vertex(v);
        for &w in deg_zero.iter() {
            dec.parent[w] = parent;
        }
        deg_one_stack.append(&mut deg_one);
    }
}

pub fn solve(g_orig: &Graph, greedy_rep: usize, fluid_rep: usize) -> Result<Decomposition, ()> {
    let mut g = g_orig.myclone();
    let mut dec = Decomposition::new(&g);

    // Find all vertices with degree one.
    let mut deg_one_stack = vec![];
    for u in 0..g.n {
        if g.neighbors[u].len() == 1 {
            deg_one_stack.push(u);
        }
    }
    remove_deg_one(&mut deg_one_stack, &mut g, &mut dec, None);

    let vert: Vec<usize> = (0..g.n).filter(|v| g.neighbors[*v].len() > 0).collect();
    let mut cc_stack = vec![(vert, None)];
    
    while let Some((vert, mut parent)) = cc_stack.pop() {
        if GLOBAL_ABORT.load(Ordering::Relaxed) {
            return Err(());
        }

        if vert.len() <= 2 {
            for u in vert {
                dec.parent[u] = parent;
                parent = Some(u);
            }
            continue;
        }
        
        if let Ok((cc, parent)) = greedy_fluid_separation(&mut g, &vert, &mut dec, parent, greedy_rep, fluid_rep) {
            for comp in cc {
                cc_stack.push((comp, parent));
            }
        } else {
            return Err(());
        }
    }

    Ok(dec)
}

pub fn greedy_separator(g: &Graph, vertices: &[usize]) -> (HashSet<usize>, f64) {
    let mut set_a = HashSet::new();
    let mut set_b: HashSet<usize> = vertices.iter().cloned().collect();
    let mut set_c = HashSet::new(); // Separator

    let start_at = *vertices.choose(&mut rand::thread_rng()).unwrap();
    set_b.remove(&start_at);
    set_c.insert(start_at);
    let mut best = grow_separator(g, vertices, &mut set_a, &mut set_b, &mut set_c, true);
    
    // Move the separator to the other direction and check, if we find a better one.
    // This improves the results on many graphs.
    let best2 = grow_separator(g, vertices, &mut set_b, &mut set_a, &mut set_c, true);
    if best2.1 > best.1 {//best2.0.len() < best.0.len() { // && 
        //eprintln!("replace (n={}) before {:?}, after {:?}", vertices.len(), (best.0.len(), best.1), (best2.0.len(), best2.1));
        best = best2;
    } /*else {
        eprintln!("keep (n={}) before {:?}, after {:?}", vertices.len(), (best.0.len(), best.1), (best2.0.len(), best2.1));
    }*/
    best
}

pub fn grow_separator(g: &Graph, vertices: &[usize], set_a: &mut HashSet<usize>, set_b: &mut HashSet<usize>, set_c: &mut HashSet<usize>, allow_grow: bool) -> (HashSet<usize>, f64) {
    // For each vertex, compute the number of neighbors in set b.
    let mut neigh_in_b = HashMap::new();
    for &u in vertices.iter() {
        let deg = g.neighbors[u].len();
        neigh_in_b.insert(u, deg);
    }
    // Subtract vertices in A and C from neigh_in_b.
    for &u in set_a.union(&set_c) {
        for &v in g.neighbors[u].iter() {
            if let Some(val) = neigh_in_b.get_mut(&v) { *val -= 1 };
        }
    }
    // Keep a copy for later.
    let neigh_in_b_orig = neigh_in_b.clone();

    // All vertices in C go to the heap.
    let mut heap: BinaryHeap<(Reverse<usize>, usize, usize)> = BinaryHeap::new(); // Tuple (-deg, priority, vertex)
    for &u in set_c.iter() {
        heap.push((Reverse(neigh_in_b[&u]), rand::thread_rng().gen(), u));
    }
    
    let mut added_history = vec![];
    let mut removed_from_b = vec![];

    while (set_a.len() + set_c.len()) * 4 < vertices.len() * 3 {
        let mut added_to_a = vec![];
        let mut added_to_c = vec![];
        let (revdeg, _, u) = heap.pop().unwrap();
        if revdeg != Reverse(neigh_in_b[&u]) {
            continue;
        }
        if neigh_in_b[&u] > 1 && !allow_grow {
            break;
        }
        set_a.insert(u);
        added_to_a.push(u);
        set_c.remove(&u);
        for &v in g.neighbors[u].iter() {
            if !set_a.contains(&v) && !set_c.contains(&v) {
                // Vertex v has been in B and goes to C now.
                set_c.insert(v);
                added_to_c.push(v);
                set_b.remove(&v);
                removed_from_b.push(v);
                heap.push((Reverse(neigh_in_b[&v]), rand::thread_rng().gen(), v));
                // Neigh_in_b count need to be updated for the neighbors of v.
                for &w in g.neighbors[v].iter() {
                    if let Some(val) = neigh_in_b.get_mut(&w) { *val -= 1 };
                    // If vertex w is in c, the new degree value needs to be pushed.
                    if set_c.contains(&w) {
                        heap.push((Reverse(neigh_in_b[&w]), rand::thread_rng().gen(), w));
                    }
                }
            }
        }
        added_history.push((added_to_a, added_to_c));
    }

    // All vertices left in set_b go to the removal list as well.
    for &u in set_b.iter() {
        removed_from_b.push(u);
    }
    let (merge_tree, subtree_size, mut merge_roots) = component_tree(g, &removed_from_b);
    
    for (added_to_a, added_to_c) in added_history.iter().rev() {
        for &u in added_to_c.iter() {
            set_b.insert(u);
            set_c.remove(&u);
        }
        for &u in added_to_a.iter() {
            set_c.insert(u);
            set_a.remove(&u);
        }
    }

    // Restore the copy of the original hash map.
    neigh_in_b = neigh_in_b_orig;

    let mut best_score = -1.0;
    let mut best_in_valid_range = false;
    let mut best_c = set_c.clone();
    let mut set_c_queue_remove = vec![];
    let mut set_c_queue_add = vec![];

    for (added_to_a, added_to_c) in added_history.iter() {
        if (set_a.len() + set_c.len()) * 4 > vertices.len() * 3 {
            break;
        }
        for &u in added_to_a.iter() {
            if set_a.contains(&u) {
                continue;
            }
            set_a.insert(u);
            set_c.remove(&u);
            set_c_queue_remove.push(u);
        }
        
        let mut no_neigh_in_b = vec![];
        for &u in added_to_c.iter() {
            if set_a.contains(&u) {
                continue;
            }            

            set_c.insert(u);
            set_c_queue_add.push(u);
            set_b.remove(&u);
            
            // Update connected components list (via merge_roots).
            for &child in merge_tree[&u].iter() {
                merge_roots.insert(child);
            }
            merge_roots.remove(&u);
            
            // Neigh_in_b count need to be updated for the neighbors of u.
            for &w in g.neighbors[u].iter() {
                if let Some(val) = neigh_in_b.get_mut(&w) {
                    *val -= 1;
                    if *val == 0 {
                        no_neigh_in_b.push(w);
                    }
                }
            }
        }
        
        let mut components: Vec<(usize, usize)> = merge_roots.iter().map(|&u| (subtree_size[&u], u)).collect();
        components.sort();
        components.pop();
        for (_comp_size, comp_root) in components {
            // Move a component to set A.
            //eprintln!("Moving {} vertices", comp_size);
            merge_roots.remove(&comp_root);
            let mut stack = vec![comp_root];
            while let Some(node) = stack.pop() {
                set_a.insert(node);
                set_b.remove(&node);
                for &w in g.neighbors[node].iter() {
                    if let Some(val) = neigh_in_b.get_mut(&w) {
                        *val -= 1;
                        if *val == 0 {
                            no_neigh_in_b.push(w);
                        }
                    };
                }
                for &child in merge_tree[&node].iter() {
                    stack.push(child);
                }
            }
        }
        
        for &u in no_neigh_in_b.iter() {
            if (set_a.len() + set_c.len()) * 4 > vertices.len() * 3 {
                break;
            }
            set_a.insert(u);
            set_b.remove(&u);
            set_c.remove(&u);
            set_c_queue_remove.push(u);
            merge_roots.remove(&u);
        }

        let new_imbalance: f64 = if set_a.len() < set_b.len() {
            set_a.len() as f64 / set_b.len() as f64
        } else {
            set_b.len() as f64 / set_a.len() as f64
        };
        
        let new_score = new_imbalance.powf(0.8) / set_c.len() as f64;
        if set_a.len() * 4 > vertices.len() * 1 && new_score > best_score {
            //best_c = set_c.clone();
            while let Some(u) = set_c_queue_add.pop() {
                best_c.insert(u);
            }
            while let Some(u) = set_c_queue_remove.pop() {
                best_c.remove(&u);
            }
            best_score = new_score;
            best_in_valid_range = true;
        }
    }
    
    if !best_in_valid_range {
        let new_imbalance: f64 = if set_a.len() < set_b.len() {
            set_a.len() as f64 / set_b.len() as f64
        } else {
            set_b.len() as f64 / set_a.len() as f64
        };
        best_score = new_imbalance.powf(0.8) / set_c.len() as f64;
        best_c = set_c.clone();
    }
    
    (best_c, best_score)
}

fn component_tree(g: &Graph, removal_order: &[usize]) -> (HashMap<usize, Vec<usize>>, HashMap<usize, usize>, HashSet<usize>) {
    let mut merge_tree = HashMap::new(); // Key: vertex, value: vector with children
    let mut subtree_size = HashMap::new(); // Key: vertex, value: size of subtree in merge tree.
    let mut roots = HashSet::new(); // Roots in merge tree.

    let mut removal_order_ind = HashMap::new();
    for (i, &u) in removal_order.iter().enumerate() {
        removal_order_ind.insert(u, i);
    }
    let mut components: QuickUnionUf<UnionBySize> = QuickUnionUf::new(removal_order.len());
    let mut rep_to_root: HashMap<usize, usize> = HashMap::new(); // Key: vertex representant of a union find component, value: root vertex of that component in the merge tree.
    
    for &u in removal_order.iter().rev() {
        let mut children = vec![];
        let mut size_u = 1;
        for &v in g.neighbors[u].iter() {
            if subtree_size.contains_key(&v) { // Is v already added to the merge tree?
                let rep_v = components.find(removal_order_ind[&v]);
                if components.union(removal_order_ind[&u], removal_order_ind[&v]) {
                    let child = rep_to_root[&rep_v];
                    size_u += subtree_size[&child];
                    roots.remove(&child);
                    children.push(child);
                }
            }
        }
        merge_tree.insert(u, children);
        subtree_size.insert(u, size_u);
        roots.insert(u);
        rep_to_root.insert(components.find(removal_order_ind[&u]), u);
    }
    (merge_tree, subtree_size, roots)
}

