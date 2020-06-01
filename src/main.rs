use signal_hook::{iterator::Signals, SIGTERM, SIGINT};
use std::{
    io::BufRead,
    error::Error,
    thread,
    sync::atomic::Ordering,
    cmp,
};
use union_find::{UnionFind, QuickUnionUf, UnionBySize};
use fluid::{
    graph::{Graph, Decomposition, PickStrategy, GLOBAL_ABORT},
    maxdeg::maxdeg_elem,
    fillin::FillinSolver,
    separation::solve,
    postprocessor::optimize_dec,
};
use std::time::Instant;

/// Read a graph in PACE2020 format from stdin.
///
fn read_graph() -> Result<Graph, Box<dyn Error>> {
    let mut g: Option<Graph> = None;
    let stdin = std::io::stdin(); // Compiler on optil.io complains, when we do not put this in a variable first.
    for line in stdin.lock().lines() {
        let line          = line?;
        let ll: Vec<&str> = line.split(" ").collect();
        match ll[0] {
            "c" => {},                                // skip comments
            "p" => {                                  // p-line
                let n = ll[2].parse::<usize>()?;      // get n and allocate memory
                g = Some(Graph::new(n));
            },
            _   => {                                  // edge line
                match g {
                    None => {                         // malformed input
                        return Err(From::from("c Found edge before p-line. Abort!"));
                    }
                    Some(ref mut g) => {              // add edge to graph
                        g.add_edge(
                            ll[0].parse::<usize>()? - 1,
                            ll[1].parse::<usize>()? - 1
                        );
                    }
                }
            }
        }
    }    
    match g {
        Some(g) => Ok(g),
        None    => Err(From::from("c Failed to parse input! Maybe it was empty?"))
    }
}

/// Print the elimination forest stored in g.parent in PACE2020 format.
/// Uses O(n) time to compute the depth of the forest and O(n) for printing.
///
/// *Assumes:*
/// - that the parents in g has been filled accordingly
///
fn print_elimination_forest(dec: &Decomposition) {
    println!("{}", dec.get_depth());
    for v in 0..dec.parent.len() {
        match dec.parent[v] {
            Some(p) => println!("{}", p+1),
            None    => println!("0")
        };
    }
}

fn get_decomposition(g: &Graph, elimination: &Vec<usize>) -> Decomposition {
    let mut dec = Decomposition::new(&g);
    
    // Initially, each vertex is in a single subtree.
    let mut subtrees: QuickUnionUf<UnionBySize> = QuickUnionUf::new(g.n);
    // The root of each subtree representant is the node itself.
    // When we merge two subtrees, we will update this array and
    // let the name-giving member of the new set point to the root
    // vertex.
    let mut root_from_rep: Vec<_> = (0..g.n).collect();
    // Has the vertex been addes to the partial forest? The parent might
    // still be None.
    let mut processed = vec![false; g.n];

    for &u in elimination.iter().rev() {
        processed[u] = true;
        let descendants: Vec<_> = g.neighbors[u].iter()
                                                .filter(|&v| processed[*v])
                                                .map(|x| *x)
                                                .collect();
        for &v in &descendants {
            // Set u as new root of the subtree containing v.
            let rep = subtrees.find(v);
            let root = root_from_rep[rep];
            dec.parent[root] = Some(u);
        }
        // Merge subtrees of u and all descendant vertices.
        for &v in &descendants {
            subtrees.union(u, v);
        }
        // Update pointer from rep to root of the new subtree.
        let new_rep = subtrees.find(u);
        root_from_rep[new_rep] = u;
    }
    
    dec
}

fn main() -> Result<(), Box<dyn Error>> {
    let now = Instant::now();
    let g = match read_graph() {
        Ok(g)  => g,
        Err(e) => {
            eprintln!("{}", e);
            std::process::exit(1);
        }
    };
    //eprintln!("read {}", now.elapsed().as_millis());

    let vertices: Vec<usize> = (0..g.n).collect();
    let g_number_of_edges = g.number_of_edges(&vertices);

    let mut best_depth = g.n + 1;
    let mut best_maxdeg = g.n + 1;
    let mut best_fillin = g.n + 1;
    let mut best_fluid = g.n + 1;
    let mut best_greedy = g.n + 1;
    let mut best_dec = Decomposition::new(&g); // Just avoid None stuff.

    let signals = Signals::new(&[SIGINT, SIGTERM])?; // TODO SIGINT might be removed here, but helps to test via CTRL-C.
    thread::spawn(move || {
        for sig in signals.forever() {
            GLOBAL_ABORT.store(true, Ordering::Relaxed);
            eprintln!("Received signal {:?}", sig);
        }
    });

    let strat_mix = vec![PickStrategy::First, PickStrategy::Last, PickStrategy::Random, PickStrategy::Random, PickStrategy::Random];

    for strat in strat_mix.iter() {
        if GLOBAL_ABORT.load(Ordering::Relaxed) {
            eprintln!("abort");
            break;
        }

        let md_elem = maxdeg_elem(&g, strat);
        let mut md_dec = get_decomposition(&g, &md_elem);
        optimize_dec(&mut md_dec);
        let md_depth = md_dec.get_depth();
        if md_depth < best_depth {
            best_depth = md_depth;
            best_dec = md_dec;
            eprintln!("Maxdeg {:?} {}", strat, best_depth);
        }
        if md_depth < best_maxdeg { best_maxdeg = md_depth }
    }

    let fi_solver = FillinSolver::new(&g);
    for strat in strat_mix.iter() {
        if GLOBAL_ABORT.load(Ordering::Relaxed) {
            eprintln!("abort");
            break;
        }

        let fi_elem = fi_solver.fillin_elem(strat);
        let mut fi_dec = get_decomposition(&g, &fi_elem);
        optimize_dec(&mut fi_dec);
        let fi_depth = fi_dec.get_depth();
        if fi_depth < best_depth {
            best_depth = fi_depth;
            best_dec = fi_dec;
            eprintln!("Fillin {:?} {}", strat, best_depth);
        }
        if fi_depth < best_fillin { best_fillin = fi_depth }
    }

    let mut edge_denominator = ((g_number_of_edges >> 10) as f64).powf(1.5) as usize / 50;
    if edge_denominator < 1 {
        edge_denominator = 1;
    } else if edge_denominator > 1_000 {
        edge_denominator = 1_000;
    }

    /*let rep_cnt = if g.n < 2000 { 100 }
                  else if g.n < 10000 { 30 }
                  else if g.n < 200000 { 3 }
                  else { 1 };*/
    let rep_cnt = cmp::min(1_000 / edge_denominator as usize, 30);

    /*if g.n < 1000 { 3000 }
                  else if g.n < 2000 { 500 }
                  else if g.n < 10000 { 200 }
                  else if g.n < 100000 { 30 }
                  else if g.n < 200000 { 10 }
                  else { 3 };*/

    let rep_cnt_greedy = 3_000 / edge_denominator as usize;

    eprintln!("Edge denominator {}, rep {}, rep_greedy {}", edge_denominator, rep_cnt, rep_cnt_greedy);

    for _ in 0..rep_cnt {
        if let Ok(mut fluid_dec) = solve(&g, 3, 3) {
            optimize_dec(&mut fluid_dec);
            let fluid_depth = fluid_dec.get_depth();
            if fluid_depth < best_depth {
                best_depth = fluid_depth;
                best_dec = fluid_dec;
                eprintln!("Fluid {}", best_depth);
            }
            if fluid_depth < best_fluid { best_fluid = fluid_depth }
        } else {
            break;
        }
    }

    for i in 0..rep_cnt_greedy {
        let sep_rep = if i > 5 && i < 100 { 20 } else { 3 };
        if let Ok(mut greedy_dec) = solve(&g, sep_rep, 0) {
            optimize_dec(&mut greedy_dec);
            let greedy_depth = greedy_dec.get_depth();
            if greedy_depth < best_depth {
                best_depth = greedy_depth;
                best_dec = greedy_dec;
                eprintln!("Greedy {}", best_depth);
            }
            //eprintln!("Greedy depth {}", greedy_depth);
            if greedy_depth < best_greedy { best_greedy = greedy_depth }
        } else {
            break;
        }
    }

    eprintln!("n={}, e={}, time={}, td={}, maxdeg +{}, fillin +{}, fluid +{}, greedy +{}", g.n, g.number_of_edges(&vertices), now.elapsed().as_secs(), best_depth,
                                                               best_maxdeg - best_depth,
                                                               best_fillin - best_depth,
                                                               best_fluid - best_depth,
                                                               best_greedy - best_depth);

    print_elimination_forest(&best_dec);

    Ok(())
}

