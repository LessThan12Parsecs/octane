/**
 * Mean of a per-snapshot signal over a gated window (e.g. intake manifold pressure while
 * the intake valve is open), in SIMULATED time so it is independent of playback speed.
 *
 * Why gated: on a single-cylinder engine with a small intake plenum the manifold refills
 * almost to ambient while the intake valve is shut and is drawn down during the intake
 * stroke (≈ 0.4–1.0 bar within one cycle at 20 % throttle on the CFR). The time-averaged
 * MAP therefore hides most of the throttling; the pressure the cylinder actually breathes
 * is the mean over the valve-open period.
 *
 * The value is published when the window closes and held until the next window closes.
 * Stepping backwards abandons the window in progress.
 */
export class GatedMean {
  private sum = 0;
  private dur = 0;
  private tLast = NaN;
  private xLast = NaN;
  private open = false;
  private held = NaN;

  /** Feed sample x at simulated time t (s) with gate state `gate`; returns the held mean (NaN until the first window closes). */
  update(x: number, t: number, gate: boolean): number {
    const dt = t - this.tLast;
    if (!(dt >= 0)) {
      // first sample or stepped backwards: restart the window
      this.sum = 0;
      this.dur = 0;
      this.open = gate;
    } else if (this.open) {
      this.sum += 0.5 * (x + this.xLast) * dt;
      this.dur += dt;
      if (!gate) {
        if (this.dur > 0) this.held = this.sum / this.dur;
        this.open = false;
      }
    } else if (gate) {
      this.open = true;
      this.sum = 0;
      this.dur = 0;
    }
    this.tLast = t;
    this.xLast = x;
    return this.held;
  }

  get value(): number {
    return this.held;
  }

  reset(): void {
    this.sum = 0;
    this.dur = 0;
    this.tLast = NaN;
    this.xLast = NaN;
    this.open = false;
    this.held = NaN;
  }
}
