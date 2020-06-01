use std::collections::{HashSet, BinaryHeap};
use rand::Rng;
use crate::graph::{Graph, PickStrategy};

pub struct FillinSolver<'a> {
    pub g: &'a Graph,
    neighbors_as_set: Vec<HashSet<usize>>,
    pub fillin_values: Vec<usize>,
}

impl<'a> FillinSolver<'a> {
    pub fn new(g: &Graph) -> FillinSolver {
        let mut neighbors_as_set = vec![HashSet::new(); g.n];
        for u in 0..g.n {
            for &v in g.neighbors[u].iter() {
                neighbors_as_set[u].insert(v);
            }
        }

        let fillin_values = FillinSolver::fill_in_values(g, &neighbors_as_set);
        FillinSolver {
            g,
            fillin_values,
            neighbors_as_set
        }
    }

    fn fill_in_values(g: &Graph, nodeset: &Vec<HashSet<usize>>) -> Vec<usize> {
        let mut fiv = vec![0; g.n];
        for u in 0..g.n {
            let fill_in = FillinSolver::vertex_fill_in(g, u, nodeset);
            fiv[u] = fill_in;
        }
        fiv
    }

    fn vertex_fill_in(g: &Graph, u: usize, nodeset: &Vec<HashSet<usize>>) -> usize {
        let mut fill_in = 0;

        for &v in g.neighbors[u].iter() {
            fill_in += g.neighbors[u].len() - nodeset[u].intersection(&nodeset[v]).count() - 1;
        }
        fill_in /= 2;

        fill_in
    }

    fn fill_in_priority(pick: &PickStrategy, pushcnt: usize) -> isize {
        match pick {
            PickStrategy::First => -1 * pushcnt as isize,
            PickStrategy::Last => pushcnt as isize,
            PickStrategy::Random => rand::thread_rng().gen()
        }
    }

    pub fn fillin_elem(&self, pick: &PickStrategy) -> Vec<usize> {
        let mut nodeset = self.neighbors_as_set.clone();

        let mut heap = BinaryHeap::new(); // Contains tuples (fillin value, priority value, vertex number).
        let mut pushcnt = 0;
        let mut curval = self.fillin_values.clone();
        let mut elimination = Vec::with_capacity(self.g.n);

        // Add the vertices to the heap according to pick strategy.
        let search = match pick {
            PickStrategy::First => self.g.bfs_edges(0),
            _ => self.g.dfs_edges(0)
        };
        heap.push((curval[0], 0, 0));
        for (_, v) in search {
            let d = curval[v];
            heap.push((d, FillinSolver::fill_in_priority(&pick, pushcnt), v));
            pushcnt += 1;
        }

        let mut removed = vec![false; self.g.n];
        while !heap.is_empty() {
            let item = match heap.pop() {
                Some(item) => item,
                None => break
            };
            let val = item.0;
            let u = item.2;
            // Has the vertex been removed?
            if removed[u] {
                continue;
            }
            // Has the value of the vertex changed, since it was added to the heap?
            if curval[u] != val {
                continue;
            }
            elimination.push(u);
            // Update values
            for &v in self.g.neighbors[u].iter() {
                if removed[v] {
                    continue;
                }
                let v_diff = nodeset[v].len() - nodeset[v].intersection(&nodeset[u]).count() - 1;
                if v_diff > 0 {
                    curval[v] -= v_diff;
                    heap.push((curval[v], FillinSolver::fill_in_priority(&pick, pushcnt), v));
                    pushcnt += 1;
                }
                nodeset[v].remove(&u);
                if nodeset[v].len() == 0 {
                    // Isolated.
                    elimination.push(v);
                    removed[v] = true;
                }
            }
            //neighbors_as_set[u] = set() # The way it is, but we don't need to update this.
            removed[u] = true;
        }
        elimination
    }
}

