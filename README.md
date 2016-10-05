# effective_mass_tensor

## Example : A band minimum at L point in reciprocal space of FCC material

## Basic idea
  i. Calculate energy on a mesh around L point using Quantum Espresso. <br />
  ii. Interplate the energy value at +dx, -dx, +dy, -dy, +dz, -dz, *etc*. <br />
  iii. Use finite method to calculate the tensor of effective mass. <br />
