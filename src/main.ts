/**
 * Octane entry point: mounts the app into index.html's viewport and UI hosts. The engine comes from
 * ?engine=<id>, the remembered choice or the default; the app switches engines in place (App.setEngine).
 */
import { startApp, type App } from './app/index';

const viewport = document.getElementById('viewport');
const uiContainer = document.getElementById('ui');
if (!viewport || !uiContainer) throw new Error('index.html must provide #viewport and #ui');

const app: App | null = startApp({ viewport, uiContainer, search: location.search });

if (import.meta.env.DEV && app) {
  // Debug handle for the browser console (dev server only).
  (globalThis as { octane?: App }).octane = app;
}

if (import.meta.hot) {
  import.meta.hot.dispose(() => app?.dispose());
}
