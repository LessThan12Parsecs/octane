import { defineConfig } from 'vite';

export default defineConfig({
  // Relative asset URLs so the build works from a subpath (GitHub Pages: /octane/).
  base: './',
  worker: { format: 'es' },
  build: { target: 'es2022', sourcemap: true },
});
