/**
 * Per-viewer UI preferences (panel open, active tab, collapsed cards) in
 * localStorage. Every access is guarded: storage may be unavailable (private
 * mode, sandboxed iframe, node tests) and the UI must work without it.
 */
const PREFIX = 'octane.ui.';

export function loadPref<T>(key: string, fallback: T): T {
  try {
    const raw = globalThis.localStorage?.getItem(PREFIX + key);
    if (raw == null) return fallback;
    return JSON.parse(raw) as T;
  } catch {
    return fallback;
  }
}

export function savePref(key: string, value: unknown): void {
  try {
    globalThis.localStorage?.setItem(PREFIX + key, JSON.stringify(value));
  } catch {
    /* storage unavailable: preferences are best-effort */
  }
}
