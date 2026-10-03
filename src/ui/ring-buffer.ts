/** Small allocation-free containers used by the UI. */

/** Fixed-capacity FIFO ring buffer. Index 0 is the oldest element. */
export class RingBuffer<T> {
  private readonly buf: (T | undefined)[];
  private head = 0; // index of the oldest element
  private len = 0;

  constructor(readonly capacity: number) {
    if (!(capacity >= 1) || !Number.isInteger(capacity)) throw new RangeError(`RingBuffer capacity must be a positive integer, got ${capacity}`);
    this.buf = new Array<T | undefined>(capacity);
  }

  get length(): number {
    return this.len;
  }

  /** Append; returns the evicted (oldest) element when full, else undefined. */
  push(v: T): T | undefined {
    if (this.len < this.capacity) {
      this.buf[(this.head + this.len) % this.capacity] = v;
      this.len++;
      return undefined;
    }
    const evicted = this.buf[this.head];
    this.buf[this.head] = v;
    this.head = (this.head + 1) % this.capacity;
    return evicted;
  }

  /** Element i (0 = oldest). */
  get(i: number): T | undefined {
    if (i < 0 || i >= this.len) return undefined;
    return this.buf[(this.head + i) % this.capacity];
  }

  /** k-th newest element (0 = newest). */
  latest(k = 0): T | undefined {
    return this.get(this.len - 1 - k);
  }

  clear(): void {
    this.buf.fill(undefined);
    this.head = 0;
    this.len = 0;
  }

  /** Oldest → newest. */
  toArray(): T[] {
    const out: T[] = [];
    for (let i = 0; i < this.len; i++) out.push(this.buf[(this.head + i) % this.capacity] as T);
    return out;
  }

  *[Symbol.iterator](): IterableIterator<T> {
    for (let i = 0; i < this.len; i++) yield this.buf[(this.head + i) % this.capacity] as T;
  }
}

/** Growable Float64Array-backed column (amortised O(1) push, no per-push allocation). */
export class Float64Column {
  data: Float64Array;
  length = 0;

  constructor(initialCapacity = 1024) {
    this.data = new Float64Array(Math.max(1, initialCapacity));
  }

  push(v: number): void {
    if (this.length === this.data.length) this.grow(this.data.length * 2);
    this.data[this.length++] = v;
  }

  grow(capacity: number): void {
    if (capacity <= this.data.length) return;
    const next = new Float64Array(capacity);
    next.set(this.data.subarray(0, this.length));
    this.data = next;
  }

  clear(): void {
    this.length = 0;
  }

  /** View of the filled part (shares memory; invalidated by the next grow). */
  view(n: number = this.length): Float64Array {
    return this.data.subarray(0, Math.min(n, this.length));
  }
}

/** Index of the first element of sorted `arr[0..n)` that is > x (n if none). */
export function upperBound(arr: ArrayLike<number>, n: number, x: number): number {
  let lo = 0;
  let hi = n;
  while (lo < hi) {
    const mid = (lo + hi) >>> 1;
    if (arr[mid] <= x) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

/** Index of the first element of sorted `arr[0..n)` that is ≥ x (n if none). */
export function lowerBound(arr: ArrayLike<number>, n: number, x: number): number {
  let lo = 0;
  let hi = n;
  while (lo < hi) {
    const mid = (lo + hi) >>> 1;
    if (arr[mid] < x) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}
