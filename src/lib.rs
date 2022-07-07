//! The Poseidon algebraic hash function lightweight implementation

#![no_std]

use core::iter;

mod constants;
use constants::{MDS, ROUND_CONSTANTS};
use pasta_curves::arithmetic::FieldExt;
use pasta_curves::group::ff::Field;
use pasta_curves::pallas::Base as Fp;
#[cfg(test)]
mod test_vectors;

/// The type used to hold permutation state.
pub(crate) type State = [Fp; 3];

#[inline]
fn sbox(x: Fp) -> Fp {
    x.pow_vartime(&[5])
}

#[inline]
fn full_round_sbox_layer(state: &mut State) {
    state[0] = sbox(state[0]);
    state[1] = sbox(state[1]);
    state[2] = sbox(state[2]);
}

#[inline]
fn half_round_sbox_layer(state: &mut State) {
    state[0] = sbox(state[0]);
}

#[inline]
fn add_round_constants(state: &mut State, rcs: &State) {
    state[0] += rcs[0];
    state[1] += rcs[1];
    state[2] += rcs[2];
}

#[inline]
fn apply_mds(state: &mut State) {
    let mut new_state = [Fp::zero(); 3];
    for i in 0..3 {
        for j in 0..3 {
            new_state[i] += MDS[i][j] * state[j];
        }
    }
    *state = new_state;
}

/// Runs the Poseidon permutation on the given state.
pub(crate) fn permute(state: &mut State) {
    let rounds = iter::empty()
        .chain(iter::repeat(true).take(4))
        .chain(iter::repeat(false).take(56))
        .chain(iter::repeat(true).take(4));

    for (is_full_round, rcs) in rounds.zip(ROUND_CONSTANTS.iter()) {
        add_round_constants(state, rcs);
        if is_full_round {
            full_round_sbox_layer(state);
        } else {
            half_round_sbox_layer(state);
        }
        apply_mds(state);
    }
}

/// Poseidon Hash
pub fn hash(x: Fp, y: Fp) -> Fp {
    let mut state = [x, y, Fp::from_u128((2 as u128) << 64)];
    permute(&mut state);
    state[0]
}

#[cfg(test)]
mod tests {
    use super::*;
    use ff::PrimeField;
    #[test]
    fn permute_test_vectors() {
        for tv in test_vectors::PERMUTE_TEST_VECTORS.iter() {
            let mut state = [
                Fp::from_repr(tv.initial_state[0]).unwrap(),
                Fp::from_repr(tv.initial_state[1]).unwrap(),
                Fp::from_repr(tv.initial_state[2]).unwrap(),
            ];

            permute(&mut state);

            for (expected, actual) in tv.final_state.iter().zip(state.iter()) {
                assert_eq!(&actual.to_repr(), expected);
            }
        }
    }

    #[test]
    fn hash_test_vectors() {
        for tv in test_vectors::HASH_TEST_VECTORS.iter() {
            let x = Fp::from_repr(tv.input[0]).unwrap();
            let y = Fp::from_repr(tv.input[1]).unwrap();

            let result = hash(x, y);

            assert_eq!(result.to_repr(), tv.output);
        }
    }
}
