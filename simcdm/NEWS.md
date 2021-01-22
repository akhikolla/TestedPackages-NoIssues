# simcdm 0.1.1

## Changes

- Addressed a signed vs. unsigned integer comparison.

## Documentation

- Fixed spacing issues brought on from an organization of _C++_ code.

## Deployment

- Changed unit tests to use R 3.5.0's RNG setup in anticipation for the RNG
  change in R 3.6.0.

# simcdm 0.1.0

## Changes

- Renamed `sim_attribute_classes()` to `attribute_classes()`.
- Addressed ambiguous calls to `std::pow(<int>, <int>)`

## Documentation

- Added a `CITATION` file for the package
- Improved Vignette Examples
- Improved README contents

## Deployment

- Added Unit Tests to verify simulation routines and attribute generations.
- Added testing on Travis-CI for the previous release of _R_, e.g. the `oldrel`.

# simcdm 0.0.5

## Features

- Added _C++_ and _R_ functions for simulation of:
    - Deterministic Input, Noisy "And" Gate (DINA)
        - Item Response: `sim_dina_items()`
        - $\eta$ Response: `sim_dina_attributes()`.
    - reduced Reparameterized Unified Model (rRUM):
        - Item Reponse: `sim_rrum_items()`
    - Matrices:
        - Random Q Matrix: `sim_q_matrix()`
        - ETA Matrix: `sim_eta_matrix()`
        - Latent Attribute Profiles for Subjects: `sim_subject_attributes()`
    - Attributes:
        - Latent Attribute Profile Matrix: `sim_attribute_classes()`
        - Attribute Bijection: `attribute_bijection()`
        - Attribute Inverse Bijection: `attribute_inv_bijection()`