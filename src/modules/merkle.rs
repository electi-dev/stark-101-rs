use crate::FieldElement;
use sha256::digest;

/// # Fields
///
/// * `data` - A vector of FieldElement, representing the data leaves of the tree
/// * `facts` - A HashMap containing the facts about the tree (tree nodes and their children)
/// * `root` - The root hash of the MerkleTree
#[derive(Clone)]
pub struct MerkleTree {
    pub data: Vec<FieldElement>,
    pub facts: std::collections::HashMap<String, (String, String)>,
    pub root: String,
}

impl MerkleTree {
    /// Constructs a new MerkleTree with the given data.
    ///
    /// # Arguments
    ///
    /// * `data` - A vector of FieldElement
    ///
    /// # Returns
    ///
    /// A new MerkleTree with the given data
    pub fn new(data: Vec<FieldElement>) -> Self {
        assert!(data.len() > 0);

        let num_leaves = (data.len() as f64).log2().ceil().exp2() as usize;
        let mut padded_data = data.to_vec();
        padded_data.resize(num_leaves, FieldElement::zero());
        let facts = std::collections::HashMap::new();
        let root = "".to_string();

        MerkleTree { 
            data: padded_data,
            facts, 
            root 
        }
    }

    /// Returns the authentication path of the given leaf id.
    ///
    /// # Arguments
    ///
    /// * `leaf_id` - The index of the leaf in the tree
    ///
    /// # Returns
    ///
    /// A vector of String representing the authentication path
    pub fn get_authentication_path(&self, leaf_id: usize) -> Vec<String> {
        assert!(leaf_id < self.data.len());

        let node_id = leaf_id + self.data.len();
        let mut cur = self.root.clone();
        let mut decommitment = Vec::new();

        for bit in format!("{:b}", node_id).chars().skip(1) {
            let res = self.facts.get(&cur);
            cur = res.unwrap().0.clone();
            let mut auth = res.unwrap().1.clone();
            if bit == '1' {
                (auth, cur) = (cur, auth);
            }
            decommitment.push(auth.clone());
        }

        decommitment
    }

    /// Recursively builds the Merkle tree and populates the facts.
    ///
    /// # Arguments
    ///
    /// * `node_id` - The ID of the current node being processed
    ///
    /// # Returns
    ///
    /// A String representing the hash of the current node being processed
    pub fn recursive_build_tree(&mut self, node_id: usize) -> String  {
        if node_id >= self.data.len() {
            let id_in_data = node_id - self.data.len();
            let leaf_data = self.data[id_in_data].to_string();
            let h = digest(leaf_data.clone());
            self.facts.insert(h.clone(), (leaf_data, "".to_string()));
            return h;
        } else{
            let left = self.recursive_build_tree(node_id * 2);
            let right = self.recursive_build_tree(node_id * 2 + 1);
            let h = digest((format!("{}{}", left, right)));
            self.facts.insert(h.clone(), (left, right));
            return h;
        }
    }

    /// Builds the Merkle tree.
    ///
    /// # Arguments
    ///
    /// None
    ///
    /// # Returns
    ///
    /// None
    pub fn build_tree(&mut self) {
        self.root = self.recursive_build_tree(1).clone();
    }
}

/// Verifies if the given decommitment matches the provided leaf data and root.
///
/// # Arguments
///
/// * `leaf_id` - The index of the leaf in the tree
/// * `leaf_data` - A reference to a FieldElement representing the leaf data
/// * `decommitment` - A reference to a vector of String representing the decommitment
/// * `root` - A reference to a String representing the root hash
///
/// # Returns
///
/// A boolean, true if the decommitment is verified, false otherwise
pub fn verify_decommitment(leaf_id: usize, leaf_data: &FieldElement, decommitment: &Vec<String>, root: &String) -> bool {
    let leaf_num = 2_usize.pow(decommitment.len() as u32);
    let node_id = leaf_id + leaf_num;

    let mut cur = digest(leaf_data.clone().to_string());
    let bits = format!("{:b}", node_id).chars().rev().collect::<Vec<_>>();

    let decommitment_iter = decommitment.iter().rev(); 
    let mut h = "".to_string();

    for (bit, auth) in bits.into_iter().zip(decommitment_iter) {
        if bit == '0' {
            h = cur.clone() + auth;
        } else {
            h = auth.clone() + &cur;
        }
        cur = digest(h.clone());
    }

    cur == *root
}

#[cfg(test)]
mod merkle {
    use crate::{
        FieldElement,
        MerkleTree,
        *
    };

    #[test]
    fn test_get_authentication_path() {
        let data = vec![
            FieldElement::zero(), 
            FieldElement::one(),
            FieldElement::from(2),
            FieldElement::from(3),
            FieldElement::from(4),
            FieldElement::from(5),
            FieldElement::from(6),
            FieldElement::from(7)];

        let mut merkle_tree = MerkleTree::new(data.clone());
        merkle_tree.build_tree();
        
        let leaf_id = 3;
        let decommitment = merkle_tree.get_authentication_path(leaf_id);
        let content = data[leaf_id];
        assert!(verify_decommitment(leaf_id, &content, &decommitment, &merkle_tree.root));
        

    }
}