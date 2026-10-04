/**
 * Combustion sub-models: laminar burning velocity, flame geometry (sphere ∩ disc chamber), the chamber
 * interface (flat disc / side-valve L-head: volume, front area, wetted walls), 0-D K–k turbulence, and the
 * Keck / Tabaczynski entrainment–burn-up model.
 */
export * from './flame-geometry';
export * from './chamber';
export * from './turbulence';
export * from './entrainment';
export * from './laminar-flame-speed';
