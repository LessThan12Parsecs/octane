/**
 * The App's animation loop, free of DOM and WebGL so its error policy can be unit-tested with a fake
 * frame scheduler:
 *
 *  - start() is idempotent (one pending frame at most) and does nothing after dispose();
 *  - every frame schedules the next one first, then runs the frame callback;
 *  - a frame that throws HALTS the loop (its pending frame is cancelled) and reports the error with
 *    the number of consecutive failed frames — a deterministic throw never spins;
 *  - resume() restarts a halted loop, and only a halted one that its owner started: the App calls it
 *    after a successful engine swap and on Reset, the recovery paths its failure card offers. If the
 *    first frame after a resume throws again the loop halts again (consecutive count 2, 3, …); a frame
 *    that completes resets the count.
 */

/** requestAnimationFrame / cancelAnimationFrame (injectable for tests). */
export interface FrameScheduler {
  request(cb: (nowMs: number) => void): number;
  cancel(handle: number): void;
}

const browserScheduler: FrameScheduler = {
  request: (cb) => requestAnimationFrame(cb),
  cancel: (h) => cancelAnimationFrame(h),
};

export class AnimationLoop {
  private handle = 0;
  /** The owner asked the loop to run (start() and no stop() since). */
  private wanted = false;
  /** A frame threw: the loop stays stopped until resume() or start(). */
  private halted = false;
  private failures = 0;
  private disposed = false;

  /**
   * @param frame   one animation frame (rAF timestamp, ms)
   * @param onError a frame threw; `consecutive` = failed frames since the last one that completed (≥ 1)
   */
  constructor(
    private readonly frame: (nowMs: number) => void,
    private readonly onError: (err: unknown, consecutive: number) => void,
    private readonly scheduler: FrameScheduler = browserScheduler,
  ) {}

  /** A frame is scheduled. */
  get running(): boolean {
    return this.handle !== 0;
  }

  /** The loop was started and is stopped by a frame error (resume() restarts it). */
  get haltedOnError(): boolean {
    return this.halted;
  }

  /** Run the loop (idempotent; clears a halt). */
  start(): void {
    if (this.disposed) return;
    this.wanted = true;
    this.halted = false;
    if (!this.handle) this.handle = this.scheduler.request(this.tick);
  }

  /** Restart the loop if a frame error halted it; no-op otherwise (running, stopped by stop(), never started). */
  resume(): void {
    if (this.halted && this.wanted) this.start();
  }

  /** Stop the loop (until start()). */
  stop(): void {
    this.wanted = false;
    this.halted = false;
    this.cancel();
  }

  dispose(): void {
    this.stop();
    this.disposed = true;
  }

  private cancel(): void {
    if (this.handle) this.scheduler.cancel(this.handle);
    this.handle = 0;
  }

  private readonly tick = (nowMs: number): void => {
    if (this.disposed) return;
    this.handle = this.scheduler.request(this.tick);
    try {
      this.frame(nowMs);
    } catch (err) {
      this.cancel();
      this.halted = true;
      this.onError(err, ++this.failures);
      return;
    }
    this.failures = 0;
  };
}
