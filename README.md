# BIE-periodic

Code implementing the Boundary Integral Equation (BIE) methods with the doubly-periodic Green's function from *"Layer Potential Methods for Doubly-Periodic Harmonic Functions"* by Bohyun Kim and Braxton Osting (preprint, [arXiv:2505.03074](https://arxiv.org/abs/2505.03074). Please contact hyunbrianna@gmail.com for any inquiries.

## Acknowledgements

This software utilizes several classes and functions from `mpspack`, authored by Alex Barnett and Timo Betcke (see [https://github.com/ahbarnett/mpspack](https://github.com/ahbarnett/mpspack)). I acknowledge partial support from the National Science Foundation under grant NSF DMS-2136198.

## Requirements

- Matlab 2008a or newer (No additional toolbox required).

## Installation

1. **Install Git** (if not already installed).
2. Clone the repository with the following command:

    ```bash
    $ git clone https://github.com/BohyunKim92/BIE-periodic
    ```

    Alternatively, click the `<> Code` tab in the upper-right corner of the repository page and select **Download ZIP**.

3. In MATLAB, open the `BIE-periodic` directory and add the necessary path by running:

    ```matlab
    addpath('mps_helpers')
    ```

4. **Example Usage**: 
   - You can run the `dirichlet_single_hole` example to observe how this software works. This example corresponds to **EXAMPLE 1** in the manuscript.

## Examples

This repository includes MATLAB implementations for the following problems:

- **Dirichlet Boundary Value Problem (BVP)**
- **Neumann Boundary Value Problem (BVP)**
- **Steklov Eigenvalue Problem (EVP)**
  
These examples feature both single and multiple holes.
