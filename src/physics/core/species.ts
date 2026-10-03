/**
 * Canonical species set shared by every physics module.
 *
 * A "composition vector" is a Float64Array of length NS indexed by these
 * indices. Whether it holds moles, mole fractions or mass fractions is always
 * stated by the variable name (N / X / Y).
 *
 * Index layout (do not reorder — other modules rely on it):
 *   [0, N_EQ)       equilibrium product species (burned-gas set, Olikara-Borman + Ar)
 *   [N_EQ, NS)      fuel species (present only in unburned gas; never in equilibrium)
 */
export const SPECIES = [
  'N2', 'O2', 'AR', 'CO2', 'H2O', 'CO', 'H2', 'OH', 'H', 'O', 'NO', 'N',
  'IC8H18', 'NC7H16', 'CH4', 'C3H8', 'C2H5OH',
] as const;

export type SpeciesName = (typeof SPECIES)[number];

export const NS = SPECIES.length;
/** Number of equilibrium (burned-gas) species; they occupy indices [0, N_EQ). */
export const N_EQ = 12;

export const SP: { readonly [K in SpeciesName]: number } = Object.freeze(
  Object.fromEntries(SPECIES.map((s, i) => [s, i])) as { [K in SpeciesName]: number },
);

export const FUEL_SPECIES = ['IC8H18', 'NC7H16', 'CH4', 'C3H8', 'C2H5OH'] as const;
export type FuelSpeciesName = (typeof FUEL_SPECIES)[number];

/** Element order for element-abundance vectors (length NE). */
export const ELEMENTS = ['C', 'H', 'O', 'N', 'AR'] as const;
export type ElementName = (typeof ELEMENTS)[number];
export const NE = ELEMENTS.length;
export const EL: { readonly [K in ElementName]: number } = Object.freeze(
  Object.fromEntries(ELEMENTS.map((e, i) => [e, i])) as { [K in ElementName]: number },
);

/**
 * Atom counts: ELEMENT_COUNTS[k][e] = atoms of element e in species k.
 * Order of e follows ELEMENTS.
 */
export const ELEMENT_COUNTS: readonly (readonly number[])[] = [
  //C  H  O  N  AR
  [0, 0, 0, 2, 0], // N2
  [0, 0, 2, 0, 0], // O2
  [0, 0, 0, 0, 1], // AR
  [1, 0, 2, 0, 0], // CO2
  [0, 2, 1, 0, 0], // H2O
  [1, 0, 1, 0, 0], // CO
  [0, 2, 0, 0, 0], // H2
  [0, 1, 1, 0, 0], // OH
  [0, 1, 0, 0, 0], // H
  [0, 0, 1, 0, 0], // O
  [0, 0, 1, 1, 0], // NO
  [0, 0, 0, 1, 0], // N
  [8, 18, 0, 0, 0], // IC8H18 (2,2,4-trimethylpentane)
  [7, 16, 0, 0, 0], // NC7H16 (n-heptane)
  [1, 4, 0, 0, 0], // CH4
  [3, 8, 0, 0, 0], // C3H8
  [2, 6, 1, 0, 0], // C2H5OH
];

/** Allocate a zeroed composition vector. */
export const newComposition = (): Float64Array => new Float64Array(NS);
