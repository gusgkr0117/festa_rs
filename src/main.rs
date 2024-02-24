use festa_rs::{
    discrete_log::xgcd_big,
    festa::{FESTACrypto, TrapdoorInput},
    fields::FpFESTAExt::{D1, D2, L_POWER},
};
use num_bigint::{BigUint, RandBigInt};
use rand::thread_rng;
use std::time::Instant;

fn gen_random_trapdoor_input() -> TrapdoorInput {
    // Pick random constants
    let d1 = BigUint::from_slice(&D1);
    let d2 = BigUint::from_slice(&D2);
    let mut rng = thread_rng();
    let s1 = rng.gen_biguint_range(&BigUint::from(0u32), &d1);
    let s2 = rng.gen_biguint_range(&BigUint::from(0u32), &d2);
    let mut b11 = rng.gen_biguint((L_POWER - 1) as u64);
    b11.set_bit(0, true);
    let l_power = BigUint::from(1u32) << L_POWER;
    let b22 = xgcd_big(&b11, &l_power) % &l_power;
    let trapdoor_input = TrapdoorInput::new(&s1, &s2, &(b11, b22));
    trapdoor_input
}

fn festa_trapdoor_bench() {
    let festa = FESTACrypto::new();

    const TOTAL_CASES: usize = 10;
    for case_num in 0..TOTAL_CASES {
        println!("=========TEST #{case_num}=======");
        let timer = Instant::now();
        let (pubkey, seckey) = festa.keygen();
        println!("keygen elapsed : {} ms", timer.elapsed().as_millis());

        let trapdoor_input = gen_random_trapdoor_input();

        let timer = Instant::now();
        let trapdoor_output = festa.trapdoor_eval(&pubkey, &trapdoor_input);
        println!("eval elapsed : {} ms", timer.elapsed().as_millis());

        let timer = Instant::now();
        let inverted = festa.trapdoor_inverse(&trapdoor_output, &seckey);
        println!("inverse elapsed : {} ms", timer.elapsed().as_millis());

        assert_eq!(trapdoor_input, inverted);
    }
}
fn main() {
    festa_trapdoor_bench();
}
