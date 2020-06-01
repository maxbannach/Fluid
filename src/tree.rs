use crate::graph::{Graph, Decomposition};
use std::collections::{HashSet, HashMap};

/// Simple macro to compute integer log2 – or, more precisely, max( floor(log2(x))+1, 0 ).
/// Only sensible defined for usize.
///
#[macro_export]
macro_rules! log2 {
    ($x:expr) => ( (0usize.leading_zeros() - ($x as usize).leading_zeros()) as usize )
}

pub fn solve_tree(g: &mut Graph, subgraph: &Vec<usize>, dec: &mut Decomposition, parent: Option<usize>) -> usize {
    let mut coloring = HashMap::new();
    let mut visited  = HashSet::new();
    if let Some(&start) = subgraph.iter().next() {
        rank_tree(g, &subgraph, start, &mut visited, &mut coloring);
        if let Some(&value) = coloring.values().max() {
            let value = value + 1;
            // extract coloring
            let mut stack = Vec::new();
            stack.push((subgraph.clone(), parent));
            while let Some((mut component, parent)) = stack.pop() {
                //eprintln!("{:?}", component);
                //eprintln!("{:?}", coloring);
                //for &v in component.iter() {
                //    eprintln!("neigh of {} is {:?}", v, g.neighbors[v]);
                //}
                let mut center = 0;
                let mut rank   = 0;
                for &v in component.iter() {
                    if coloring[&v] >= rank {
                        center = v;
                        rank   = coloring[&v];
                    }
                }
                //eprintln!("set parent of {} to {:?}", center, parent);
                dec.parent[center] = parent;
                g.remove_vertex(center);
                component.retain(|&x| x != center);
                for next in g.connected_components(&component) {
                    stack.push((next, Some(center)));
                }
            }
            return value;
        }
    }
    return 0;
}

/// Polynomial-Time tree ranking with the algorithm from
/// Ananth V Iyer, H Donald Ratliff, and Gopalakrishnan Vijayan: Optimal node ranking of trees.
///
fn rank_tree(g: &Graph, subgraph: &[usize], v: usize, visited: &mut HashSet<usize>, coloring: &mut HashMap<usize, usize>) -> Vec<usize> {
    visited.insert(v);
    let max_color    = log2!(subgraph.len()) + 1;
    let mut critical = vec![0; max_color];

    // compute critical list of children
    for &w in g.neighbors[v].iter() {
        if !subgraph.contains(&w) || visited.contains(&w) { continue; }
        let child_list = rank_tree(g, subgraph, w, visited, coloring);
        for &color in child_list.iter() {
            critical[color] = critical[color] + 1;
        }        
    }
    
    // extract own color
    let mut color = 0;
    for i in 0..max_color {
        if critical[i] > 1 { color = i + 1;}
    }
    while critical[color] > 0 { color = color + 1; } 
    coloring.insert(v, color);
    
    // create own critical list
    let mut result = Vec::new();
    let mut i = max_color - 1;
    while i > color {
        if critical[i] > 0 { result.push(i); }
        i = i - 1;
    }
    result.push(color);       
    return result;
}

