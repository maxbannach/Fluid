use std::collections::{VecDeque, HashMap};
use rand::Rng;
use crate::graph::{Graph, PickStrategy};
//use std::time::Instant;

pub fn maxdeg_elem(g: &Graph, pick: &PickStrategy) -> Vec<usize> {
    let mut degmap = HashMap::new();
    let mut elimination = Vec::with_capacity(g.n);

    // Compute degrees.
    let mut curdeg = vec![0; g.n];
    for u in 0..g.n {
        curdeg[u] = g.neighbors[u].len();
    }

    // Add the vertices to the degree list according to pick strategy.
    let search = match pick {
        PickStrategy::First => g.bfs_edges(0),
        _ => g.dfs_edges(0)
    };
    degmap.entry(curdeg[0]).or_insert(VecDeque::new()).push_back(0);
    for (_, v) in search {
        let d = curdeg[v];
        degmap.entry(d).or_insert(VecDeque::new()).push_back(v);
    }

    let mut md = match degmap.keys().max() {
        Some(mx) => *mx,
        None => 0
    };
    while md > 0 {
        match degmap.get_mut(&md) {
            None => {
                md -= 1;
                continue;
            },
            Some(vertices) => {
                if vertices.is_empty() {
                    md -= 1;
                    continue;
                }
                // Pick element from degree list according to pick strategy.
                let chosen = match pick {
                    PickStrategy::Random => {
                        let i = rand::thread_rng().gen_range(0, vertices.len());
                        vertices.swap_remove_back(i)
                    },
                    PickStrategy::First => vertices.pop_front(),
                    PickStrategy::Last => vertices.pop_back()
                };
                let u = match chosen {
                    None => {
                        md -= 1;
                        continue;
                    },
                    Some(uu) => uu
                };
                // Has the degree of the vertex changed, since it was added to the degmap list?
                if curdeg[u] != md {
                    continue;
                }
                elimination.push(u);
                // Update degrees
                curdeg[u] = 0;
                for &v in &g.neighbors[u] {
                    if curdeg[v] == 0 {
                        continue;
                    }
                    curdeg[v] -= 1;
                    let d = curdeg[v];
                    if d == 0 {
                        // Isolated.
                        elimination.push(v);
                    } else {
                        degmap.entry(d).or_insert(VecDeque::new()).push_back(v);
                    }
                }
            }
        }
    }
    elimination
}

