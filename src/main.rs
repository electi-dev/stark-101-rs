mod modules;

use crate::{
    modules::field::FieldElement,
    modules::polynomial::Polynomial,
    modules::channel::Channel,
    modules::merkle::{MerkleTree},
    modules::merkle::verify_decommitment, 
    modules::list_utils,
};

use rand::random;
use std::{collections::HashSet};

pub fn part_1() -> (
    Vec<FieldElement>, //a
    FieldElement, //g
    Vec<FieldElement>, //G 
    FieldElement, // h
    Vec<FieldElement>, // H
    Vec<FieldElement>, //eval_domain
    Polynomial, //f
    Vec<FieldElement>, // f_eval
    MerkleTree, //f_merkle
    Channel //channel
) {
    // Step 1: Construct a list of length 1023 whose first two elements are 1 and 3141592. The next 1021 elements will be the FibonacciSq sequence induced by these two elements
    // Modulus: 3221225473
    let mut a = vec![FieldElement::from(1), FieldElement::from(3141592)];
    while a.len() < 1023 {
        a.push(a[a.len() - 1].pow(2) + a[a.len() - 2].pow(2));
    }

    assert_eq!(a.len(), 1023);
    assert_eq!(a[0], FieldElement::one());
    assert_eq!(a[1022], FieldElement::from(2338775057));

    //Step 2: Find Group of Size 1024. Essentially this means to find an element s.t when raised to 3145728 (3 * 2**20) equals to 1 and has size 1024.
    let g = FieldElement::generator().pow(3145728);
    let mut G = Vec::<FieldElement>::new();
    for i in 0..1024 {
        G.push(g.pow(i));
    }

    assert!(g.is_order(1024));
    let mut b = FieldElement::one();
    for i in 0..1023 {
        assert!(b == G[i]);
        b = b * g;
        assert_ne!(b, FieldElement::one());
    }
    assert_eq!(b * g, FieldElement::one());

    // Step 3: Interpolate a polynomial where the x values are the vector `a` and y values the vector `G`
    let f = Polynomial::interpolate_poly(&G[..G.len() - 1].to_vec(), &a);
    assert_eq!(f.eval(&FieldElement::from(2)), FieldElement::from(1302089273));   
    for i in 0..1023 {
        assert_eq!(f.eval(&G[i]), a[i]);
    }

    // Step 4: Evaluate on a larger domain. We will create another subgroup H of order 8192. It exists because |H|(8192) divides |F|(3221225472). The generator is g**393216.
    // Note that H is a subgroup that (might) contain elements from G. The coset is built in such a way that is does *not* contain elements of G or H. The coset is not necessarily a subgroup
    // it is rather a shift of the subgroup
    let h = FieldElement::generator().pow(393216);
    let mut H = Vec::<FieldElement>::new();
    let mut eval_domain = Vec::<FieldElement>::new();
    for i in 0..8192 {
        H.push(h.pow(i));
        eval_domain.push(FieldElement::generator() * H[i as usize]);
    }

    assert!(h.is_order(8192));
    assert_eq!(HashSet::<i128>::from_iter(eval_domain.iter().map(|e| e.val)).len(), eval_domain.len());

    let w_inv = FieldElement::generator().inverse();
    for i in 0..8192 as usize{
        assert_eq!((w_inv * eval_domain[1]).pow(i as i128) * FieldElement::generator(), eval_domain[i]);
    }

    // Step 5: Evaluate the trace polynomial generated before over the coset generated in step 4.
    let mut f_eval = Vec::<FieldElement>::new();

    for i in 0..eval_domain.len() {
        f_eval.push(f.eval(&eval_domain[i]));
    }

    // Step 6: Create a Merkle Tree commitment
    let mut f_merkle = MerkleTree::new(f_eval.clone());
    f_merkle.build_tree();

    let leaf_id = random::<usize>() % f_eval.len();
    let leaf_data = f_eval[leaf_id];
    let decommitment = f_merkle.get_authentication_path(leaf_id);
    assert_eq!(verify_decommitment(leaf_id, &leaf_data, &decommitment, &f_merkle.root), true);

    // Step 7: Generate the interaction
    let mut channel = Channel::new();
    channel.send(&f_merkle.root);

    return(
        a,
        g,
        G, 
        h,
        H,
        eval_domain,
        f,
        f_eval,
        f_merkle,
        channel
    );

}

fn part_2() -> (
    Polynomial, 
    Vec<FieldElement>, 
    MerkleTree, 
    Channel, 
    Vec<FieldElement>,
    Vec<FieldElement>,
    MerkleTree) {
    let (_,
        g,
        G, 
        _,
        _,
        eval_domain,
        f,
        f_eval,
        f_merkle,
        mut channel) = part_1();

    // Step 1: Create the constraints

    //First Constraint
    let numer0 = f.clone() + Polynomial::from(-1);
    let denom0 = Polynomial::new(&vec![ - FieldElement::one(), FieldElement::one()], 'x');
    let p0 = numer0.clone() / denom0.clone();

    assert_eq!(Polynomial::new(&vec![], 'x'), numer0 % denom0);
    assert_eq!(p0.eval(&FieldElement::from(2718)), FieldElement::from(2509888982));

    //Second Constraint
    let numer1 = f.clone() + Polynomial::from(-2338775057);
    let denom1 = Polynomial::new(&vec![- G[1022], FieldElement::one()], 'x');
    let p1 = numer1.clone() / denom1.clone();

    assert_eq!(Polynomial::new(&vec![], 'x'), numer1 % denom1);
    assert_eq!(p1.clone().eval(&FieldElement::from(5772)), FieldElement::from(232961446));

    //Third Constraint
    let prod_vec = (0..1024)
    .into_iter()
    .map(|i| Polynomial::new(&vec![-g.pow(i), FieldElement::one()], 'x'))
    .collect::<Vec<Polynomial>>();
    
    let prod = list_utils::prod(&prod_vec);
    assert_eq!(prod.repr_latex(), "$3221225472 +x^{1024}$");

    let numer2 = f.compose(&(Polynomial::new(&vec![FieldElement::zero(), FieldElement::from(g.pow(2))], 'x'))) -
    f.compose(&(Polynomial::new(&vec![FieldElement::zero(), FieldElement::from(g)], 'x'))).pow(2) -
    f.pow(2);


    assert_eq!(numer2.eval(&g.pow(1020)), FieldElement::zero());
    assert_eq!(numer2.eval(&g.pow(1021)), FieldElement::from(230576507));

    let denom2 = prod / 
    (Polynomial::new(&vec![FieldElement::from(-g.pow(1021)), FieldElement::one()], 'x') *
    Polynomial::new(&vec![FieldElement::from(-g.pow(1022)), FieldElement::one()], 'x') *
    Polynomial::new(&vec![FieldElement::from(-g.pow(1023)), FieldElement::one()], 'x'));

    let p2 = numer2 / denom2;

    assert_eq!(p2.degree(), 1023);
    assert_eq!(p2.eval(&FieldElement::from(31415)), FieldElement::from(2090051528));

    //Step 2 Composition Polynomial
    //Compute 𝐶𝑃(𝑥)=𝛼0⋅𝑝0(𝑥)+𝛼1⋅𝑝1(𝑥)+𝛼2⋅𝑝2(𝑥) where  𝛼0,𝛼1,𝛼2 are random field elements obtained from the verifier, or in our case - from the channel.
    // This is to achieve succintness.
    let get_cp = |channel: &mut Channel| -> Polynomial {
        let alpha0 = channel.receive_random_field_element();
        let alpha1 = channel.receive_random_field_element();
        let alpha2 = channel.receive_random_field_element();
    
        return Polynomial::from(alpha0) * p0.clone() + Polynomial::from(alpha1) * p1.clone() + Polynomial::from(alpha2) * p2.clone();
    };

    // Step 3: Commit on the composition
    let cp_eval = |channel: &mut Channel| -> (Polynomial, Vec<FieldElement>) {
        let cp = get_cp(channel);

        let cp_eval = eval_domain
        .iter()
        .map(|d| cp.eval(d))
        .collect();

        return(cp, cp_eval);
    };

    let (cp, cp_eval) = cp_eval(&mut channel);

    assert_eq!(cp.degree(), 1023);

    let mut cp_merkle = MerkleTree::new(cp_eval.clone());
    cp_merkle.build_tree();
    channel.send(&cp_merkle.root);

    return (cp, cp_eval, cp_merkle, channel, eval_domain, f_eval, f_merkle)
        
}

fn part_3() -> (
    Vec<FieldElement>, 
    MerkleTree, 
    Vec<Polynomial>, 
    Vec<Vec<FieldElement>>, 
    Vec<Vec<FieldElement>>, 
    Vec<MerkleTree>,
    Channel) {
    let (
        cp, 
        cp_eval, 
        cp_merkle,
        mut channel,
        eval_domain,
        f_eval,
        f_merkle
        ) = part_2();

    let next_fri_domain = |eval_domain: &Vec<FieldElement>| -> Vec<FieldElement> {
        eval_domain[..(eval_domain.len() / 2)]
        .iter()
        .map(|e| e.pow(2))
        .collect()
    };

    let next_fri_polynomial = |poly: &Polynomial, beta: &FieldElement| -> Polynomial {
        let odd_coeffs = poly.poly.iter().enumerate().filter_map(|(index, value)| {
            if index % 2 == 1 {
                Some(*value) 
            } else {
                None
            }
        }).collect::<Vec<FieldElement>>();

        let even_coeffs = poly.poly.iter().enumerate().filter_map(|(index, value)| {
            if index % 2 == 0 {
                Some(*value) 
            } else {
                None
            }
        }).collect::<Vec<FieldElement>>();

        let odd = Polynomial::from(beta.clone()) * Polynomial::new(&odd_coeffs, 'x');
        let even = Polynomial::new(&even_coeffs, 'x');

        odd + even
    };
    
    let next_fri_layer = |poly: &Polynomial, domain: &Vec<FieldElement>, beta: &FieldElement| -> (Polynomial, Vec<FieldElement>, Vec<FieldElement>) {
        let next_poly = next_fri_polynomial(&poly, &beta);
        let next_domain = next_fri_domain(&domain);
        let next_layer = next_domain.iter().map(|e| next_poly.eval(e)).collect();

        (next_poly, next_domain, next_layer)
    };

    let test_poly = Polynomial::new(&vec![FieldElement::from(2), FieldElement::from(3), FieldElement::zero(), FieldElement::one()], 'x');
    let test_domain = vec![FieldElement::from(3), FieldElement::from(5)];
    let beta = FieldElement::from(7);
    let (next_p, next_d, next_l) = next_fri_layer(&test_poly, &test_domain, &beta);

    assert_eq!(next_p.poly, vec![FieldElement::from(23), FieldElement::from(7)]);
    assert_eq!(next_d, vec![FieldElement::from(9)]);
    assert_eq!(next_l, vec![FieldElement::from(86)]);

    let FriCommit = |cp: &Polynomial, domain: &Vec<FieldElement>, cp_eval: &Vec<FieldElement>, cp_merkle: &MerkleTree, channel: &mut Channel| 
    -> (Vec<Polynomial>, Vec<Vec<FieldElement>>, Vec<Vec<FieldElement>>, Vec<MerkleTree>){
        let mut fri_polys = vec![cp.clone()];
        let mut fri_domains = vec![domain.clone()];
        let mut fri_layers = vec![cp_eval.clone()];
        let mut fri_merkles = vec![cp_merkle.clone()];

        while fri_polys.last().unwrap().degree() > 0 {
            let beta = channel.receive_random_field_element();
            let (next_poly, next_domain, next_layer) = next_fri_layer(
                &fri_polys.last().unwrap(), 
                &fri_domains.last().unwrap(),
                 &beta);
            let mut next_merkle = MerkleTree::new(next_layer.clone());
            next_merkle.build_tree();
            fri_polys.push(next_poly);
            fri_domains.push(next_domain);
            fri_layers.push(next_layer);
            fri_merkles.push(next_merkle);
            channel.send(&fri_merkles.last().unwrap().root);
        }
        channel.send(&fri_polys.last().unwrap().poly[0].to_string());
        (fri_polys, fri_domains, fri_layers, fri_merkles)
    };

    let (fri_polys, fri_domains, fri_layers, fri_merkles) = FriCommit(&cp, &eval_domain, &cp_eval, &cp_merkle, &mut channel);
    assert_eq!(fri_layers.len(), 11);
    assert_eq!(fri_layers.last().unwrap().len(), 8);
    assert_eq!(fri_polys.last().unwrap().degree(), 0);

    return(f_eval, f_merkle, fri_polys, fri_domains, fri_layers, fri_merkles, channel)
    
}

fn part_4() {
    let (_, _, _, _, _, _, _, f_eval, f_merkle, _) = part_1();
    let (
        f_eval,
        f_merkle,
        fri_polys,
        fri_domains,
        fri_layers,
        fri_merkles,
        mut channel) = part_3();

    let decommit_on_fri_layers = |idx: usize, channel: &mut Channel| {
        let fri_layers_without_last = fri_layers.split_last().unwrap().1;
        let fri_merkles_without_last = fri_merkles.split_last().unwrap().1;
        for (layer, merkle) in fri_layers_without_last.iter().zip(fri_merkles_without_last.iter()) {
            let length = layer.len();
            let idx_mod = idx % length;
            let sib_idx = (idx + length / 2) % length;
            channel.send(&String::from(layer[idx_mod].val.to_string()));
            channel.send(&merkle.get_authentication_path(idx_mod).join(","));
            channel.send(&String::from(layer[sib_idx].to_string()));
            channel.send(&merkle.get_authentication_path(sib_idx).join(","));       
        }  
        channel.send(&fri_layers.last().unwrap()[0].val.to_string());
    };

    let decommit_on_query = |idx: usize, channel: &mut Channel| {
        assert!(idx + 16 < f_eval.len());
        channel.send(&f_eval[idx].val.to_string());
        channel.send(&f_merkle.get_authentication_path(idx).join(","));
        channel.send(&f_eval[idx + 8].val.to_string());
        channel.send(&f_merkle.get_authentication_path(idx + 8).join(","));
        channel.send(&f_eval[idx + 16].val.to_string());
        channel.send(&f_merkle.get_authentication_path(idx + 16).join(","));
        decommit_on_fri_layers(idx, channel);
    };

    let decommit_fri = |channel: &mut Channel| {
        for _ in 0..3 {
            let num = channel.receive_random_int(0, 8191-16, false).as_usize();
            decommit_on_query(num, channel);
        }
    };

    decommit_fri(&mut channel);

}

// The objective of this code is to prove that we know a field element X s.t the 1023rd element of the FibonacciSq sequence starting the 1, is 2338775057
// The FibonacciSq sequence is defined as a_{n+2} = a^2_{n+1} + a^2_n
fn main() {
    part_4();
}
