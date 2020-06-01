use std::collections::HashSet;
use crate::graph::Decomposition;
use crate::tree::solve_tree;

/// Tries to reduce the depth of the given treedepth decomposition by reordering
/// some subtrees.
///
pub fn optimize_dec(dec: &mut Decomposition) {
    //eprintln!("before {:?}", dec.get_depth());
    let mut gr = dec.orig_g.myclone();
    let subtrees = dec.subtrees();
    //eprintln!("subtree computed");
    let mut removed = HashSet::new();
    for (st_root, _subtree) in subtrees.iter() {
        let mut cur = *st_root;
        while let Some(parent) = dec.parent[cur] {
            if removed.contains(&parent) {
                break;
            }
            gr.remove_vertex(parent);
            removed.insert(parent);
            cur = parent;
        }
        //eprintln!("{:?} {:?}", subtree, parent);
    }
    //eprintln!("graph modified");
    //eprintln!("{:?}", subtrees);
    for (st_root, subtree) in subtrees.iter() {
        let parent = dec.parent[*st_root];
        let ccs = gr.connected_components(&subtree);
        if ccs.len() > 1 {
            // TODO should not happen in a well built tree. However, remove_separator does not check connectivity in between.
            //eprintln!("{:?}", ccs);
        }
        for cc in ccs {
            //eprintln!("{:?}, root {}", subtree, *st_root);
            solve_tree(&mut gr, &cc, dec, parent);
        }
    }
    //eprintln!("{:?}", gr);
    //eprintln!("after {:?}", dec.get_depth());
}
