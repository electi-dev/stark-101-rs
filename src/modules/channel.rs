use std::{
    any::type_name, 
    fmt::Debug
};
use crate::FieldElement;
use sha256::digest;
use serde::Serialize;
use primitive_types::U256;

pub fn serialize<T: Serialize + Debug>(obj: &T) -> String {
    serde_json::to_string(obj).unwrap()
}
#[derive(Clone)]
pub struct Channel{
    pub state: String,
    pub proof: Vec<String>
}

impl Channel {
    pub fn new() -> Self {
        Channel {
            state: "0".to_string(),
            proof: Vec::<String>::new()
        }
    }

    pub fn send(&mut self, s: &String) {
        self.state = digest(self.state.clone() + s);
        self.proof.push(format!("{}:{}", type_name::<Self>(), s));
    }

    pub fn receive_random_int(&mut self, min: usize, max: usize, show_in_proof: bool) -> U256 {
        let state = U256::from_str_radix(&self.state, 16).unwrap();
        let num = U256::from(min) + state % U256::from(max - min + 1);
        self.state = digest(self.state.clone());
        if show_in_proof {
            self.proof.push(format!("{}:{}", type_name::<Self>(), num));
        }
        num
    }

    pub fn receive_random_field_element(&mut self) -> FieldElement {
        let num = self.receive_random_int(0, (FieldElement::modulus() - 1) as usize, false);
        self.proof.push(format!("{}:{}", type_name::<Self>(), num));
        FieldElement::from(num)
    }
}