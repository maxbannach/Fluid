use std::collections::{HashSet, HashMap, VecDeque};
use std::mem;
use std::sync::atomic::AtomicBool;

pub static GLOBAL_ABORT: AtomicBool = AtomicBool::new(false);


#[derive(Debug)]
pub struct Graph {
    pub n:         usize,              // vertices are 0,1,...,n.1 (Some might get removed later on)
    pub neighbors: Vec<Vec<usize>>,    // neighbors of vertex i
    neighbors_pos: Vec<HashMap<usize, usize>>,// At which position of the neighbor[i] vector is vertex v?
}

impl Graph {

    /// Allocate memory for a graph with n vertices.
    ///
    pub fn new(n: usize) -> Graph {
        Graph{
            n: n,
            neighbors: vec![Vec::new(); n],
            neighbors_pos: vec![HashMap::new(); n],
        }
    }
    
    pub fn myclone(&self) -> Graph {
        Graph{
            n: self.n,
            neighbors: self.neighbors.clone(),
            neighbors_pos: self.neighbors_pos.clone(),
        }
    }

    /// Add the undirected edge {u,v} to the graph.
    ///
    /// *Assumes:*
    /// - u,v in {0,...,n-1}
    /// - {u,v} not already present in the graph
    ///
    pub fn add_edge(&mut self, u: usize, v: usize) {
        self.neighbors[u].push(v);
        self.neighbors[v].push(u);
        self.neighbors_pos[u].insert(v, self.neighbors[u].len() - 1);
        self.neighbors_pos[v].insert(u, self.neighbors[v].len() - 1);
    }

    // Remove vertex from graph.
    // Returns the neighbors, that have a degree of zero and of one after
    // this removal. Note: the information about deg zero nodes is crucial,
    // since there is no way to distinguish these isolated nodes from formerly
    // removed nodes.
    pub fn remove_vertex(&mut self, v: usize) -> (Vec<usize>, Vec<usize>) {
        let mut deg_zero = vec![];
        let mut deg_one = vec![];
        // To avoid mutability issues with the iterator, we move the neighbor vector of v
        // and replace it with an empty one (since v gets removed).
        let vneigh = mem::replace(&mut self.neighbors[v], Vec::new());
        // Remove v in the neighbor lists.
        for &w in vneigh.iter() {
            self.remove_neighbor(w, v);
            if self.neighbors[w].len() == 0 {
                deg_zero.push(w);
            } else if self.neighbors[w].len() == 1 {
                deg_one.push(w);
            }
        }
        (deg_zero, deg_one)
    }

    // Remove all edges that go to the given component.
    pub fn separate_vertex(&mut self, v: usize, component: &HashSet<usize>) -> (Vec<usize>, Vec<usize>) {
        let mut deg_zero = vec![];
        let mut deg_one = vec![];
        // To avoid mutability issues with the iterator, we move the neighbor vector of v
        // and replace it with an empty one (since v gets removed).
        let vneigh = mem::replace(&mut self.neighbors[v], Vec::new());
        //self.neighbors_pos[v] = HashMap::new();
        // Remove v in the neighbor lists.
        for &w in vneigh.iter() {
            if component.contains(&w) {
                self.remove_neighbor(w, v);
                if self.neighbors[w].len() == 0 {
                    deg_zero.push(w);
                } else if self.neighbors[w].len() == 1 {
                    deg_one.push(w);
                }
            } else {
                self.neighbors[v].push(w);
                self.neighbors_pos[v].insert(w, self.neighbors[v].len() - 1);
            }
        }
        if self.neighbors[v].len() == 1 {
            deg_one.push(v);
        }
        (deg_zero, deg_one)
    }
    
    // Remove neighbor v in the neighbor list of w.
    pub fn remove_neighbor(&mut self, w: usize, v: usize) {
        // Where is v in neighbors[w]? This position will be free.
        let v_pos = match self.neighbors_pos[w].remove(&v) {
            Some(v_pos) => v_pos,
            None => return
        };
        // Move the last element to v_pos...
        let last = match self.neighbors[w].pop() {
            Some(v) => v,
            None => return
        };
        // ... if it is not v itself.
        if last != v {
            self.neighbors[w][v_pos] = last;
            // Update position of "last" element.
            self.neighbors_pos[w].insert(last, v_pos);
        }
    }
    
    pub fn number_of_edges(&self, vertices: &[usize]) -> usize {
        let mut cnt = 0;
        for &v in vertices {
            cnt += self.neighbors[v].len();
        }
        cnt / 2
    }

    pub fn dfs_edges(&self, start_at: usize) -> Vec<(usize, usize)> {
        let mut edges: Vec<(usize,usize)> = Vec::new();
        let mut stack = Vec::new();
        let mut visited = HashSet::new();
        visited.insert(start_at);
        stack.push( (start_at, self.neighbors[start_at].iter()) );
        while let Some((u, mut neigh_iter)) = stack.pop() {
            // Find next child.
            let v = match neigh_iter.next() {
                Some(&v) => v,
                None => {
                    continue;
                }
            };
            // Put iterator back to the stack.
            stack.push((u, neigh_iter));
            
            if visited.contains(&v) {
                continue;
            }
            visited.insert(v);
            stack.push((v, self.neighbors[v].iter()));
            
            edges.push((u, v));
        }
        edges
    }

    pub fn bfs_edges(&self, start_at: usize) -> Vec<(usize, usize)> {
        let mut edges = Vec::new();
        let mut queue = VecDeque::with_capacity(self.n);
        let mut visited = HashSet::new();
        visited.insert(start_at);
        queue.push_back(start_at);
        while let Some(u) = queue.pop_front() {
            for &v in self.neighbors[u].iter() {
                if visited.contains(&v) {
                    continue;
                }
                visited.insert(v);
                queue.push_back(v);
                edges.push((u,v));
            }
        }
        edges
    }

    pub fn connected_components(&self, vertices: &[usize]) -> Vec<Vec<usize>> {
        let mut visited = HashSet::new();
        let mut comps = Vec::new();
        for &v in vertices {
            if visited.contains(&v) {
                continue;
            }
            let mut comp = vec![v];
            visited.insert(v);
            for (_v, w) in self.dfs_edges(v) {
                comp.push(w);
                visited.insert(w);
            }
            comps.push(comp);
        }
        comps
    }

    pub fn connected_subset_test(&self, vert_set: &HashSet<usize>, mut to_be_reached: HashSet<usize>) -> Vec<Vec<usize>> {
        if to_be_reached.len() <= 1 {
            return vec![];
        }
        let mut comps = Vec::new();
        let mut visited = HashSet::new();
        while let Some(&v) = to_be_reached.iter().next() {
            if visited.contains(&v) {
                continue;
            }
            let mut comp = vec![v];
            visited.insert(v);
            to_be_reached.remove(&v);

            let mut queue = VecDeque::new();
            queue.push_back(v);
            while let Some(u) = queue.pop_front() {
                for &w in self.neighbors[u].iter() {
                    if !vert_set.contains(&w) {
                        continue;
                    }
                    if visited.contains(&w) {
                        continue;
                    }
                    if to_be_reached.remove(&w) {
                        if to_be_reached.len() == 0 && comps.len() == 0 {
                            return vec![];
                        }
                    }
                    comp.push(w);
                    visited.insert(w);
                    queue.push_back(w);
                }
            }
            comps.push(comp);
        }
        comps
    }

    pub fn connected_components_subgraph(&self, vertices: &[usize], forbidden: HashSet<usize>) -> Vec<Vec<usize>> {
        let mut visited = forbidden;
        let mut comps = Vec::new();
        let mut queue = VecDeque::new();
        for &v in vertices {
            if visited.contains(&v) {
                continue;
            }
            let mut comp = vec![v];
            visited.insert(v);
            queue.push_back(v);
            while let Some(u) = queue.pop_front() {
                for &w in self.neighbors[u].iter() {
                    if visited.contains(&w) {
                        continue;
                    }
                    comp.push(w);
                    visited.insert(w);
                    queue.push_back(w);
                }
            }
            comps.push(comp);
        }
        comps
    }

    //fn articulation_points(vertices: &usize) -> Vec<usize> {
    //fn biconnected_components(vertices: &usize) -> Vec<Vec<usize>> {
    //    self._biconnected_dfs(vertices)
    //}
    
    pub fn biconnected_components(&self, vertices: &[usize]) -> (HashSet<usize>, Vec<HashSet<usize>>) {
        let (articulation, comp_edges) = self._biconnected_dfs(vertices, true);
        let art_set = articulation.iter().cloned().collect();
        let mut components = vec![];
        for edges in comp_edges.iter() {
            let mut comp = HashSet::new();
            for (u,v) in edges.iter() {
                comp.insert(*u);
                comp.insert(*v);
            }
            components.push(comp);
        }
        (art_set, components)
    }

    pub fn _biconnected_dfs(&self, vertices: &[usize], compute_components: bool) -> (Vec<usize>, Vec<Vec<(usize, usize)>>) {
        // depth-first search algorithm to generate articulation points
        // and biconnected components
        let mut components = vec![];
        let mut articulation = vec![];
        let mut visited = HashSet::new();
        for &start in vertices.iter() {
            if visited.contains(&start) {
                continue;
            }
            let mut discovery = HashMap::new();
            discovery.insert(start, 0); // "time" of first discovery of node during search
            let mut low = HashMap::new();
            low.insert(start, 0);
            let mut root_children = 0;
            visited.insert(start);
            let mut edge_stack = vec![];
            let mut stack = vec![(start, start, self.neighbors[start].iter())];
            while let Some((grandparent, parent, mut children)) = stack.pop() {
                // Find next child.
                let child = match children.next() {
                    Some(&v) => v,
                    None => {
                        if stack.len() > 1 {
                            if low[&parent] >= discovery[&grandparent] {
                                if compute_components {
                                    //eprintln!("suche {:?} in {:?}", (grandparent,parent), edge_stack);
                                    let ind = edge_stack.iter().position(|&elem| elem == (grandparent, parent)).unwrap();
                                    let suffix = edge_stack.split_off(ind);
                                    //eprintln!("split in {:?} und {:?}", edge_stack, suffix);
                                    components.push(suffix);
                                }
                                articulation.push(grandparent);
                            }
                            if low[&parent] < low[&grandparent] {
                                low.insert(grandparent, low[&parent]);
                            }
                        } else if stack.len() == 1 { // length 1 so grandparent is root
                            root_children += 1;
                            if compute_components {
                                let ind = edge_stack.iter().position(|&elem| elem == (grandparent, parent)).unwrap();
                                components.push(edge_stack[ind..].to_vec());
                            }
                        }
                        continue;
                    }
                };
                // Put iterator back to the stack.
                stack.push((grandparent, parent, children));
                if grandparent == child {
                    continue;
                }
                if visited.contains(&child) {
                    if discovery[&child] <= discovery[&parent] { // back edge
                        if discovery[&child] < low[&parent] {
                            low.insert(parent, discovery[&child]);
                        }
                        if compute_components {
                            //eprintln!("push {:?}", (parent, child));
                            edge_stack.push((parent, child));
                        }
                    }
                } else {
                    low.insert(child, discovery.len());
                    discovery.insert(child, discovery.len());
                    visited.insert(child);
                    stack.push((parent, child, self.neighbors[child].iter()));
                    if compute_components {
                        //eprintln!("push {:?}", (parent, child));
                        edge_stack.push((parent, child));
                    }
                }
            }
            if root_children > 1 {
                articulation.push(start);
            }
        }
        (articulation, components)
    }
}


pub struct Decomposition {
    pub parent:    Vec<Option<usize>>, // parent of each vertex in the treedepth decomposition (if set yet)
    pub orig_g:        Graph,
}

impl Decomposition {
    
    /// Allocate memory for a graph with n vertices.
    ///
    pub fn new(g: &Graph) -> Decomposition {
        let orig_g = g.myclone();
        Decomposition{
            parent:    vec![None; g.n],
            orig_g,
        }
    }

    /// Computes the depth of the tree stored in self.parent via path-compression.
    /// This operation takes time O(n) and does not modify the graph.
    ///
    pub fn get_depth(&self) -> usize {
        let n = self.parent.len();
        let mut s = Vec::with_capacity(n);                              // take clean copy of pointer structure
        for v in 0..n {                                                 // path-compression will modify this structure
            let p = match self.parent[v] {
                Some(p) => p,
                None    => v
            };
            s.insert(v,p);
        }

        let mut depth = vec![2; n];
        let mut stack = Vec::with_capacity(n);
        for mut v in 0..n {                                             // for every vertex v of the tree we do once:
            while v != s[v] {                                           // (i) crawl up the tree to find the root of v
                stack.push(v);
                v = s[v];
            }
            stack.pop();
            while let Some(w) = stack.pop() {                           // (ii) crawl down to compress the path and
                depth[w] = depth[s[w]] + 1;                             //      to compute the depth
                s[w] = v;
            }
        }
        
        // done
        return *depth.iter().max().unwrap();                            // maximum depth of any vertex is depth of the tree
    }
    
    pub fn subtrees(&self) -> Vec<(usize, Vec<usize>)> {
        let n = self.parent.len();
        
        // Compute children in decomposition tree T.
        let mut children = vec![Vec::new(); n];
        let mut child_cnt = vec![0; n];
        let mut root = 0;
        for v in 0..n {
            if let Some(w) = self.parent[v] {
                children[w].push(v);
                child_cnt[w] += 1;
            } else {
                root = v;
            }
        }

        // Build stack with children on top, and parents somewhere below.
        // Hence, for each node the children have been processed before the node
        // is popped.
        let mut children_first_stack = vec![root];
        {
            let mut helper_queue = VecDeque::new();
            helper_queue.push_back(root);
            while let Some(v) = helper_queue.pop_front() {
                for &w in children[v].iter() {
                    helper_queue.push_back(w);
                    children_first_stack.push(w);
                }
            }
        }
        // Compute, which subtrees of T induce a tree in G.
        let mut is_tree = vec![false; n];
        while let Some(v) = children_first_stack.pop() {
            if children[v].len() == 0 {
                is_tree[v] = true;
                continue;
            }
            let mut maybe_tree = true;
            for &w in children[v].iter() {
                if is_tree[w] == false {
                    maybe_tree = false;
                    break;
                }
            }
            if !maybe_tree { continue };
            // From now on, all children in this subtree of T are marked to be a tree.
            // All parents still are not marked as tree. Thus, is_tree can be used
            // as "visited" indicator.
            let mut cnt = 0;
            for &w in self.orig_g.neighbors[v].iter() {
                if is_tree[w] == true {
                    cnt += 1;
                }
            }
            if cnt == children[v].len() {
                is_tree[v] = true;
            }
        }

        // Find subtree roots in T, that induce a tree in G.
        let mut subtree_roots = vec![];
        for v in 0..n {
            if children[v].len() == 0 { continue; } // No leafs of T as subtree root.
            match self.parent[v] {
                None => if is_tree[v] { subtree_roots.push(v) },
                Some(p) => if is_tree[v] && !is_tree[p] { subtree_roots.push(v) },
            }
        }
        //eprintln!("{:?}", subtree_roots);

        let mut subtrees = vec![];
        for &st_root in subtree_roots.iter() {
            let mut subtree = vec![st_root];
            let mut stack = vec![st_root];
            while let Some(v) = stack.pop() {
                for &w in children[v].iter() {
                    stack.push(w);
                    subtree.push(w);
                }
            }
            //eprintln!("subtree with root {} is {:?}", st_root, subtree);
            if subtree.len() > 2 {
                subtrees.push((st_root, subtree));
            }
        }
        //eprintln!("{:?}", subtrees);

        subtrees
    }
}

pub struct Bipgraph {
    pub set_left: HashSet<usize>,
    pub set_right: HashSet<usize>,
    pub neighbors: HashMap<usize, Vec<usize>>,
}

impl Bipgraph {
    pub fn new(g: &Graph, left: &[usize], right: &[usize]) -> Bipgraph {
        let set_left: HashSet<usize> =  left.iter().cloned().collect();
        let set_right: HashSet<usize> =  right.iter().cloned().collect();
        let mut neighbors: HashMap<usize, Vec<usize>> = HashMap::new();
        for &u in left {
            neighbors.insert(u, g.neighbors[u].iter().filter(|&v| set_right.contains(v)).cloned().collect());
        }
        for &u in right {
            neighbors.insert(u, g.neighbors[u].iter().filter(|&v| set_left.contains(v)).cloned().collect());
        }
        Bipgraph {
            set_left,
            set_right,
            neighbors,
        }
    }
    
    //fn neighbors(&self, u: usize) -> impl Iterator<Item = &usize> {
    pub fn neighbors(&self, u: usize) -> &Vec<usize> {
        //self.g.neighbors[u].iter().filter(|&v| self.is_left[u] ^ self.is_left[*v])
        //self.g.neighbors[u].iter().filter(|&v| (self.is_left[u] && self.is_right[*v]) || (self.is_right[u] && self.is_left[*v])).collect()
        &self.neighbors[&u]
    }
}

#[derive(Debug)]
pub enum PickStrategy {
    First,
    Last,
    Random
}

