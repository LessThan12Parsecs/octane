/**
 * AnimationLoop: the App's frame loop and its render-error policy (fake scheduler, no DOM).
 * Regression for the review finding "the animation loop is never restarted after a render-loop error":
 * after a frame error the failure card offers "Back to <engine>" and Reset, and both must bring the
 * loop back (App.afterSwap / the Reset hook call resume()), while a deterministic error must not spin.
 */
import { describe, expect, it } from 'vitest';
import { AnimationLoop, type FrameScheduler } from './loop';

/** Manual rAF: frames run only when the test flushes them. */
class FakeScheduler implements FrameScheduler {
  private next = 1;
  readonly pending = new Map<number, (now: number) => void>();
  now = 0;
  request(cb: (now: number) => void): number {
    const h = this.next++;
    this.pending.set(h, cb);
    return h;
  }
  cancel(h: number): void {
    this.pending.delete(h);
  }
  /** Run the frames pending now (not the ones they schedule). */
  flush(): void {
    const due = [...this.pending.entries()];
    this.pending.clear();
    this.now += 16;
    for (const [, cb] of due) cb(this.now);
  }
}

function setup(failOn: (frame: number) => boolean = () => false) {
  const sched = new FakeScheduler();
  const frames: number[] = [];
  const errors: [string, number][] = [];
  let n = 0;
  const loop = new AnimationLoop(
    (now) => {
      n++;
      if (failOn(n)) throw new Error(`frame ${n}`);
      frames.push(now);
    },
    (err, consecutive) => errors.push([(err as Error).message, consecutive]),
    sched,
  );
  return { sched, frames, errors, loop };
}

describe('AnimationLoop', () => {
  it('start is idempotent: one pending frame, every frame schedules the next', () => {
    const { sched, frames, loop } = setup();
    loop.start();
    loop.start();
    expect(sched.pending.size).toBe(1);
    sched.flush();
    sched.flush();
    expect(frames).toEqual([16, 32]);
    expect(sched.pending.size).toBe(1);
    expect(loop.running).toBe(true);
  });

  it('a frame error halts the loop (no pending frame) and is reported once', () => {
    const { sched, errors, loop } = setup((k) => k === 2);
    loop.start();
    sched.flush();
    sched.flush(); // throws
    expect(errors).toEqual([['frame 2', 1]]);
    expect(loop.running).toBe(false);
    expect(loop.haltedOnError).toBe(true);
    expect(sched.pending.size).toBe(0);
    sched.flush();
    expect(errors).toHaveLength(1); // nothing runs any more
  });

  it('resume() restarts a halted loop (session swap / Reset), and is a no-op otherwise', () => {
    const { sched, frames, errors, loop } = setup((k) => k === 1);
    loop.resume(); // never started (the App constructor's afterSwap): stays stopped
    expect(sched.pending.size).toBe(0);
    loop.start();
    sched.flush(); // frame 1 throws → halted
    expect(errors).toEqual([['frame 1', 1]]);
    loop.resume();
    expect(loop.running).toBe(true);
    expect(loop.haltedOnError).toBe(false);
    loop.resume(); // already running: still one pending frame
    expect(sched.pending.size).toBe(1);
    sched.flush();
    sched.flush();
    expect(frames).toHaveLength(2);
  });

  it('a second consecutive error halts again (count 2) instead of spinning; a completed frame resets the count', () => {
    let failing = true;
    const { sched, errors, loop } = setup(() => failing);
    loop.start();
    sched.flush();
    loop.resume();
    sched.flush(); // throws again right after the restart
    expect(errors).toEqual([['frame 1', 1], ['frame 2', 2]]);
    expect(sched.pending.size).toBe(0);
    for (let k = 0; k < 5; k++) sched.flush();
    expect(errors).toHaveLength(2);
    failing = false;
    loop.resume();
    sched.flush(); // completes
    failing = true;
    sched.flush();
    expect(errors[2]).toEqual(['frame 4', 1]);
  });

  it('stop() cancels and is not undone by resume(); dispose() ends it for good', () => {
    const { sched, loop } = setup((k) => k === 1);
    loop.start();
    sched.flush(); // halted
    loop.stop();
    loop.resume();
    expect(sched.pending.size).toBe(0);
    loop.start();
    expect(sched.pending.size).toBe(1);
    loop.dispose();
    expect(sched.pending.size).toBe(0);
    loop.start();
    loop.resume();
    expect(sched.pending.size).toBe(0);
  });
});
