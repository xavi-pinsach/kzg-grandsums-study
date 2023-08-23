# Grand Sums Study
![Tests](https://github.com/xavi-pinsach/kzg-grandsums-study/actions/workflows/tests.yml/badge.svg)

This repository is a study about the grand-sum usage in multiset equality arguments proposed by [Ulrich Haböck](https://eprint.iacr.org/2022/1530.pdf) in the multivariate setting and adapted by [Héctor Masip](https://hecmas.github.io/) / Polygon zkEvm to the univariate setting. It explores the concept and implementation of grand-sums instead grand-product.

To be able to compare the performance of grand-sums and grand-product, we have implemented both. The grand-product protocol is described [here](https://hackmd.io/-5URcycYTlOgJTt3rCMrpw). The grand-sum protocol is described [here](https://hackmd.io/D3-fws5dQkiUQK9kSHABbA?view).

⚠️ **Warning: This Repository is for Study Purposes Only!** ⚠️

Please note that this repository is intended for study and experimental purposes. It is not meant to be used in a production environment. The code provided here may not have undergone thorough testing, and it may contain bugs or security vulnerabilities.

## Overview

The purpose of this study is to analyze and understand the grand-sums concept proposed by [Héctor Masip](https://hecmas.github.io/). It provides a implementation of grand sums in JavaScript, allowing developers to explore and experiment with the concept.

## Contents

The repository contains the following files:

- `grand-sums.js`: The JavaScript file containing the implementation of the grand-sums concept.
- `example.js`: An example usage of the grand-sums implementation.
- `LICENSE`: The license file specifying the terms of use for this repository.

## Usage

Before running the tests, you need to create a `tmp` folder in the root of the repository. Then, you need to download the `powersOfTau28_hez_final_11.ptau` file from [this link](https://github.com/iden3/snarkjs#7-prepare-phase-2) and place it in the `tmp` folder.

## TODO List
- [ ] Benchmarking: Compare grand-product vs grand-sum (either single ones and multiple ones).
- [x] Implement NTT-based polynomial arithmetic to be able to operate over the coefficients. Also $O(n\log(n))$ vs $O(n^2)$.
- [ ] Fix the edge case where `nBits = 1`.
- [x] Generalize the API/protocol to run in the case of a vector multiset equality.
- [ ] Generalize the API/protocol to run in the case of selected vector multiset equality.
- [ ] Generalize the API/protocol to run in the case of multiple multiset equality checks.

## License

This study is released under the [MIT License](LICENSE). Please review the license file for more details.

