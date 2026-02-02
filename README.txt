# Gravity and Gravity Gradient Forward Modeling Using Trapezoidal Cells

## 1. Purpose of this repository

This repository provides MATLAB implementations for 2D forward modeling of gravity and gravity gradient tensor fields produced by a density contrast model discretized into trapezoidal cells.

The code computes the following physical quantities:

* Gravity components:

  * **Gx**: horizontal component
  * **Gz**: vertical component

* Gravity gradient tensor components:

  * **Gxz**, **Gzz**

In addition to the field responses, the corresponding kernel matrices (**Fx, Fz, Fxz, Fzz**) are explicitly constructed, allowing direct integration into gravity inversion or joint gravity–seismic inversion workflows.

The formulation supports **general trapezoidal discretization**, with rectangular cells treated as a special case.

---

## 2. Software requirements and setup

* MATLAB (standard distribution)
* No external toolboxes are required

The repository is self-contained. Users only need to download or clone the repository and set the working directory accordingly in MATLAB.

---

## 3. Input data

A test model is provided in `model.mat` to facilitate verification and performance evaluation. The file contains:

* `mod`:
  Density contrast matrix defining the subsurface model, where rows and columns correspond to the discretized model grid.

* `nodex`, `nodez`:
  Node coordinates (x and z) defining the geometry of each trapezoidal cell.

* `gobs`:
  Observation point coordinates:

  * column 1: x-coordinates
  * column 2: z-coordinates

Users may replace the test model with their own datasets, provided that the variable structure and dimensional consistency are preserved.

---

## 4. Usage

1. Launch MATLAB and set the repository folder as the current working directory.

2. Ensure that `model.mat` (or a user-defined input file with the same structure) is present.

3. Run the main program:

   ```matlab
   main.m
   ```

4. The program performs forward gravity modeling and outputs gravity and gravity gradient responses computed using different implementation strategies.

---

## 5. Numerical implementations

To evaluate computational efficiency, the forward modeling is implemented using four progressively optimized strategies:

（1）. **Method 1 – Baseline implementation**
   Direct calculation without optimization.

（2）. **Method 2 – Edge-based simplification**
   Simplified analytical expressions for vertical and horizontal edges.

（3）. **Method 3 – Edge contribution precomputation**
   Further optimization based on Method 2, using precomputation and storage of edge-related terms.

（4）. **Method 4 – Node contribution precomputation**
   Additional optimization based on Method 3, incorporating precomputation and reuse of node-related terms.

These implementations are provided for methodological comparison.
For practical applications, **Method 4 is recommended**, as it offers the highest computational efficiency while maintaining numerical accuracy.

---

## 6. Output

Depending on the selected implementation, the program outputs:

* Gravity components: **Gx**, **Gz**
* Gravity gradient tensor components: **Gxz**, **Gzz**
* Kernel matrices: **Fx**, **Fz**, **Fxz**, **Fzz**

The outputs can be directly used for further analysis, validation studies, or inversion and joint inversion applications.

---

## 7. Reproducibility

The inclusion of a complete test model and explicit kernel construction ensures that all numerical results produced by this code are fully reproducible.
Users are encouraged to report the MATLAB version used when publishing results based on this software.

---

## 8. Contact

For questions or comments regarding the code, please contact:

**[Xiaogang Zhu]**
[Ocean University of China]
[zhuxiaogang@88.com]

## 9. License

This software is released under the **GNU General Public License (GPL), version 3 or later**.

The program is distributed in the hope that it will be useful, but **WITHOUT ANY WARRANTY**; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

A copy of the license should be included with this repository.
If not, see [https://www.gnu.org/licenses/](https://www.gnu.org/licenses/).
